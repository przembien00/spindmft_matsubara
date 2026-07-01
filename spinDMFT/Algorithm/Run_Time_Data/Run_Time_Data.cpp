#include"Run_Time_Data.h"
#include<fstream>
#include<iomanip>
#include<mpi.h>
#include<Globals/MPI_Types.h>

#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;

#include"RTD_Error_Handling.h"
namespace error = spinDMFT::Run_Time_Data::Error_Handling;

namespace spinDMFT::Run_Time_Data
{
// ============================================================================
// =========================== RUN TIME DATA CLASS ============================
// ============================================================================
// constructor
RunTimeData::RunTimeData( const ps::ParameterSpace& pspace, const int my_rank ):
    num_Iterations(                         0 ),
    my_rank(                                my_rank ),
    num_PrintDigits(                        pspace.num_PrintDigits ),
    critical_eigenvalue_ratio(              pspace.critical_eigenvalue_ratio ),
    seed(                                   pspace.seed ),
    num_Samples(                            pspace.num_Samples ),
    num_SamplesPerCore(                     pspace.num_SamplesPerCore ),
    num_Cores(                              pspace.get_num_Cores() ),
    self_consistency(                       pspace.self_consistency ),
    absolute_iteration_error_threshold(     pspace.absolute_iteration_error_threshold ),
    Iteration_Limit(                        pspace.Iteration_Limit )
{
    sample_sqsum_Re = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    sample_sqsum_Im = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    sample_stds_Re  = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    sample_stds_Im  = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    tau_int_Re      = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    tau_int_Im      = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    spin_expval_sqsum = FieldVector{0.,0.,0.};
    spin_expval_stds  = FieldVector{0.,0.,0.};

    error_method = pspace.error_method;
    num_blocks   = pspace.mh_num_blocks;
    cur_block_Re = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    cur_block_Im = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };

    // set the truncation scheme for negative eigenvalues:
    if( pspace.truncation_scheme_negative_eigenvalues == "set_zero" )
    {
        truncate_if_negative = []( RealType& eig ) -> void
        { 
            if( eig < RealType{0.} )
            {
                eig = RealType{0.};
            }
        };
    }
    else
    {
        error::EIGENVALUE_TRUNCATION( pspace.truncation_scheme_negative_eigenvalues, __PRETTY_FUNCTION__ );
    }
}

// return the seed string
std::string RunTimeData::get_seed_str()
{
    return seed;
}

// search for negative eigenvalues and replace them
void RunTimeData::process_and_check_eigenvalues( EigenValuesList& eigenvalues )
{
    // 1.) compute positive and negative eigenvalues sums blockwise and truncate negative eigenvalues 
    std::vector<RealType> sum_pos{}, sum_neg{};
    std::for_each( eigenvalues.begin(), eigenvalues.end(),
    [this, &sum_pos, &sum_neg]( auto& eig_block )
    {
        // A.) compute the sum of positive and negative eigenvalues of the block
        sum_pos.emplace_back( std::accumulate( eig_block.cbegin(), eig_block.cend(), RealType{0.}, []( RealType sum, const RealType& eig )
        {  
            RealType add = ( eig > 0 ) ? eig : RealType{0.};
            return std::move( sum ) + add;
        } ) );
        sum_neg.emplace_back( std::accumulate( eig_block.cbegin(), eig_block.cend(), RealType{0.}, []( RealType sum, const RealType& eig )
        {  
            RealType add = ( eig < 0 ) ? eig : RealType{0.};
            return std::move( sum ) + add;
        } ) );
        
        // B.) set the negative eigenvalues of the block to zero
        std::for_each( eig_block.begin(), eig_block.end(), truncate_if_negative );
    } );

    // 2.) compute ratio of negative and positive eigenvalues sum
    RealType psum = std::accumulate( sum_pos.cbegin(), sum_pos.cend(), RealType{0.} );
    RealType nsum = std::accumulate( sum_neg.cbegin(), sum_neg.cend(), RealType{0.} );
    RealType ratio = std::abs( nsum / psum ); 

    // 3.) terminate if negative eigenvalue ratio is too large
    if( ratio > critical_eigenvalue_ratio )
    {
        error::CRITICAL_EIGENVALUE_RATIO( ratio, critical_eigenvalue_ratio, __PRETTY_FUNCTION__ );
    }

    // 4.) save negative eigenvalue ratio
    negative_eigenvalue_ratios.emplace_back( ratio );
}

// return number of samples
size_t RunTimeData::get_num_SamplesPerCore() const
{
    return num_SamplesPerCore;
}
size_t RunTimeData::get_num_Samples() const
{
    return this->get_num_SamplesPerCore() * num_Cores;
}

// Estimate the lag-1 autocorrelation rho1 of the pCN chain from the latest acceptance rate
// and the step size, and return the AR(1) variance-inflation factor sqrt( (1+rho1)/(1-rho1) ).
// A rejected pCN step copies the state, an accepted one retains a fraction sqrt(1-beta^2) of
// it, so rho1 ~ 1 - a*(1 - sqrt(1-beta^2)). This is a crude single-exponential estimate: it
// assumes one relaxation timescale and underestimates the true autocorrelation when the chain
// has a slow / trapped mode. rho1 is clamped below 1 so the factor stays finite when the chain
// barely moves.
RealType RunTimeData::pcn_autocorrelation_factor( RealType pcn_step_size ) const
{
    if( acceptance_rates.empty() ) return RealType{1.};
    const RealType a = acceptance_rates.back();
    const RealType retained = std::sqrt( std::max( RealType{0.}, RealType{1.} - pcn_step_size * pcn_step_size ) );
    RealType rho1 = RealType{1.} - a * ( RealType{1.} - retained );
    rho1 = std::min( std::max( rho1, RealType{0.} ), RealType{1.} - RealType{1e-12} );
    return std::sqrt( ( RealType{1.} + rho1 ) / ( RealType{1.} - rho1 ) );
}

// Dispatch to the configured standard-error estimator. "blocking" (batch means, default) is
// robust to autocorrelation; "ar1" is the legacy acceptance-based single-mode factor that
// underestimates the error for slow/trapped chains. `sample_mean_*` holds <o> (already divided
// by N in main.cpp); for the AR(1) path `sample_sqsum_*` holds the MPI-reduced sum_i o_i^2.
void RunTimeData::compute_sample_stds( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im,
                                       const FieldVector& spin_mean, RealType N, RealType pcn_step_size )
{
    if( error_method == "ar1" )
        compute_sample_stds_ar1( sample_mean_Re, sample_mean_Im, spin_mean, N, pcn_step_size );
    else
        compute_sample_stds_blocking( sample_mean_Re, sample_mean_Im, N );
}

// Legacy estimator: stderr( <o> ) = sqrt( ( <o^2> - <o>^2 ) / N ), scaled by the AR(1)
// autocorrelation factor from pcn_autocorrelation_factor.
void RunTimeData::compute_sample_stds_ar1( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im,
                                           const FieldVector& spin_mean, RealType N, RealType pcn_step_size )
{
    const RealType autocorr_factor = pcn_autocorrelation_factor( pcn_step_size );
    auto per_point_stderr = [N, autocorr_factor]( const RealType& sample_mean, const RealType& sample_sqsum ) -> RealType
    {
        const RealType var = std::abs( sample_sqsum / N - sample_mean * sample_mean );
        return autocorr_factor * std::sqrt( var / N );
    };

    auto stderr_correlation = [&]( const CorrTen& mean, const CorrTen& sqsum, CorrTen& out )
    {
        auto mean_dir = mean.cbegin();
        auto sqsum_dir = sqsum.cbegin();
        auto out_dir = out.begin();
        for( ; mean_dir != mean.cend(); ++mean_dir, ++sqsum_dir, ++out_dir )
        {
            Corr buf( mean_dir->size() );
            auto m_it = mean_dir->cbegin();
            auto q_it = sqsum_dir->cbegin();
            auto b_it = buf.begin();
            for( ; m_it != mean_dir->cend(); ++m_it, ++q_it, ++b_it )
                *b_it = per_point_stderr( *m_it, *q_it );
            *out_dir = std::move( buf );
        }
    };

    stderr_correlation( sample_mean_Re, sample_sqsum_Re, sample_stds_Re );
    stderr_correlation( sample_mean_Im, sample_sqsum_Im, sample_stds_Im );

    for( size_t i = 0; i < 3; ++i )
        spin_expval_stds[i] = per_point_stderr( spin_mean[i], spin_expval_sqsum[i] );

    // reset per-iteration accumulators
    sample_sqsum_Re = CorrTen{ sample_mean_Re.get_symmetry(), sample_mean_Re[0].size() };
    sample_sqsum_Im = CorrTen{ sample_mean_Re.get_symmetry(), sample_mean_Re[0].size() };
    spin_expval_sqsum = FieldVector{ 0., 0., 0. };
}

// Allocate the batch-mean blocks for one production sweep and zero all per-iteration
// accumulators. block_length = num_SamplesPerCore / num_blocks; any remainder samples
// (< num_blocks) still enter the global mean but not the variance estimate.
void RunTimeData::init_blocks( const CorrTen& template_corr, size_t num_SamplesPerCore )
{
    const char sym  = template_corr.get_symmetry();
    const size_t ntp = template_corr[0].size();
    if( num_blocks == 0 ) num_blocks = 1;
    block_length = num_SamplesPerCore / num_blocks;
    if( block_length == 0 ) { block_length = 1; num_blocks = num_SamplesPerCore; }
    blocks_filled = 0;

    cur_block_Re   = CorrTen{ sym, ntp };
    cur_block_Im   = CorrTen{ sym, ntp };
    cur_block_spin = FieldVector{ 0., 0., 0. };
    block_means_Re.assign(   num_blocks, CorrTen{ sym, ntp } );
    block_means_Im.assign(   num_blocks, CorrTen{ sym, ntp } );
    block_means_spin.assign( num_blocks, FieldVector{ 0., 0., 0. } );

    // also reset the AR(1) accumulators (kept available for errmethod=ar1)
    sample_sqsum_Re = CorrTen{ sym, ntp };
    sample_sqsum_Im = CorrTen{ sym, ntp };
    spin_expval_sqsum = FieldVector{ 0., 0., 0. };
}

// Turn the current block's running sums into block means, store them, and reset the running
// sums for the next block. Excess calls (beyond num_blocks) are ignored.
void RunTimeData::close_block()
{
    if( blocks_filled >= num_blocks ) return;
    const RealType inv = RealType{1.} / static_cast<RealType>( block_length );
    auto finalize = []( CorrTen& cur, CorrTen& dst, const RealType inv )
    {
        auto c = cur.begin(); auto d = dst.begin();
        for( ; c != cur.end(); ++c, ++d )
        {
            auto ci = c->begin(); auto di = d->begin();
            for( ; ci != c->end(); ++ci, ++di ) { *di = *ci * inv; *ci = RealType{0.}; }
        }
    };
    finalize( cur_block_Re, block_means_Re[blocks_filled], inv );
    finalize( cur_block_Im, block_means_Im[blocks_filled], inv );
    for( size_t i = 0; i < 3; ++i ) { block_means_spin[blocks_filled][i] = cur_block_spin[i] * inv; cur_block_spin[i] = RealType{0.}; }
    ++blocks_filled;
}

// Batch-means standard error of the global mean. Each core contributes num_blocks block means;
// these are pooled across cores via MPI_Allgather, and the error of the global mean is
// std(block means)/sqrt(num_blocks*num_cores) per correlation point. Valid when
// block_length >> autocorrelation time (see the blocking-curve diagnostic).
void RunTimeData::compute_sample_stds_blocking( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im, RealType N )
{
    int num_cores = 1;
    MPI_Comm_size( MPI_COMM_WORLD, &num_cores );

    // number of correlation points (sum of per-direction lengths)
    size_t P = 0;
    for( const auto& dir : sample_stds_Re ) P += dir.size();

    // flatten this core's block means: layout [block][point]
    std::vector<RealType> local_Re( num_blocks * P ), local_Im( num_blocks * P ), local_spin( num_blocks * 3 );
    for( size_t b = 0; b < num_blocks; ++b )
    {
        size_t p = 0;
        for( const auto& dir : block_means_Re[b] ) for( size_t t = 0; t < dir.size(); ++t ) local_Re[b * P + (p++)] = dir[t];
        p = 0;
        for( const auto& dir : block_means_Im[b] ) for( size_t t = 0; t < dir.size(); ++t ) local_Im[b * P + (p++)] = dir[t];
        for( size_t i = 0; i < 3; ++i ) local_spin[b * 3 + i] = block_means_spin[b][i];
    }

    // pool all blocks from all cores
    const size_t n_tot = num_blocks * static_cast<size_t>( num_cores );
    std::vector<RealType> all_Re( n_tot * P ), all_Im( n_tot * P ), all_spin( n_tot * 3 );
    MPI_Allgather( local_Re.data(),   num_blocks * P, MPI_REALTYPE, all_Re.data(),   num_blocks * P, MPI_REALTYPE, MPI_COMM_WORLD );
    MPI_Allgather( local_Im.data(),   num_blocks * P, MPI_REALTYPE, all_Im.data(),   num_blocks * P, MPI_REALTYPE, MPI_COMM_WORLD );
    MPI_Allgather( local_spin.data(), num_blocks * 3, MPI_REALTYPE, all_spin.data(), num_blocks * 3, MPI_REALTYPE, MPI_COMM_WORLD );

    (void) n_tot;
    const size_t nc = static_cast<size_t>( num_cores );

    // Batch-means standard error of the global mean at block-merge factor g, for one point:
    // merge g consecutive blocks WITHIN each core, then take the spread of the (num_blocks/g)*nc
    // merged block means divided by sqrt of their count (Welford, single pass, no allocation).
    auto sigma_at = [&]( const std::vector<RealType>& pool, size_t stride, size_t idx, size_t g ) -> RealType
    {
        const size_t bpc  = num_blocks / g;
        const size_t ntot = bpc * nc;
        if( ntot < 2 ) return RealType{0.};
        RealType mean = RealType{0.}, M2 = RealType{0.};
        size_t count = 0;
        for( size_t c = 0; c < nc; ++c )
            for( size_t k = 0; k < bpc; ++k )
            {
                RealType s = RealType{0.};
                for( size_t e = 0; e < g; ++e ) s += pool[ ( c * num_blocks + k * g + e ) * stride + idx ];
                const RealType val = s / static_cast<RealType>( g );
                ++count;
                const RealType d = val - mean; mean += d / static_cast<RealType>( count );
                M2 += d * ( val - mean );
            }
        const RealType var = M2 / static_cast<RealType>( count - 1 );
        return std::sqrt( var / static_cast<RealType>( count ) );
    };

    // The reported error is the blocking-curve plateau: the largest batch-means estimate over
    // the merge factors g whose pooled block count stays >= MIN_BLOCKS (g=1 is always allowed).
    // Reporting the finest binning alone underestimates the error when the block length is not
    // yet >> the autocorrelation time; taking the plateau (the maximum) avoids that bias.
    constexpr size_t MIN_BLOCKS = 8;
    auto plateau = [&]( const std::vector<RealType>& pool, size_t stride, size_t idx ) -> RealType
    {
        RealType best = RealType{0.};
        for( size_t g = 1; num_blocks / g >= 2; g *= 2 )
        {
            if( g > 1 && ( num_blocks / g ) * nc < MIN_BLOCKS ) break;
            best = std::max( best, sigma_at( pool, stride, idx, g ) );
        }
        return best;
    };

    auto fill = [&]( const std::vector<RealType>& pool, CorrTen& out )
    {
        size_t p = 0;
        for( auto& dir : out ) for( auto it = dir.begin(); it != dir.end(); ++it, ++p ) *it = plateau( pool, P, p );
    };
    fill( all_Re, sample_stds_Re );
    fill( all_Im, sample_stds_Im );
    for( size_t i = 0; i < 3; ++i ) spin_expval_stds[i] = plateau( all_spin, 3, i );

    // aggregate blocking curve over the Re correlation points (diagnostic, stored to HDF5)
    blocking_curve_len.clear();
    blocking_curve_sigma.clear();
    for( size_t g = 1; num_blocks / g >= 2; g *= 2 )
    {
        RealType acc = RealType{0.};
        size_t counted = 0;
        for( size_t p = 0; p < P; ++p )
        {
            const RealType sg = sigma_at( all_Re, P, p, g );
            if( sg == RealType{0.} ) continue;   // skip points with no signal (not-yet-symmetrized upper half)
            acc += sg;
            ++counted;
        }
        blocking_curve_len.push_back( static_cast<RealType>( block_length * g ) );
        blocking_curve_sigma.push_back( counted ? acc / static_cast<RealType>( counted ) : RealType{0.} );
    }

    // Integrated autocorrelation time (in pCN steps), recovered per point from the variance
    // inflation between the blocking error and the i.i.d. error: with sigma_iid^2 = s^2 / N
    // (s^2 the marginal sample variance from sample_sqsum) and sigma_blocking^2 the reported
    // error, sigma_blocking^2 / sigma_iid^2 = 2 tau_int, so tau_int = 0.5 sigma_blocking^2 N / s^2.
    auto fill_tau = [N]( const CorrTen& mean, const CorrTen& sqsum, const CorrTen& stds, CorrTen& out )
    {
        auto m = mean.cbegin(); auto q = sqsum.cbegin(); auto sd = stds.cbegin(); auto o = out.begin();
        for( ; m != mean.cend(); ++m, ++q, ++sd, ++o )
        {
            auto mi = m->cbegin(); auto qi = q->cbegin(); auto si = sd->cbegin(); auto oi = o->begin();
            for( ; mi != m->cend(); ++mi, ++qi, ++si, ++oi )
            {
                const RealType s2 = std::abs( *qi / N - ( *mi ) * ( *mi ) ); // marginal sample variance
                *oi = ( s2 > RealType{0.} ) ? RealType{0.5} * ( *si ) * ( *si ) * N / s2 : RealType{0.};
            }
        }
    };
    fill_tau( sample_mean_Re, sample_sqsum_Re, sample_stds_Re, tau_int_Re );
    fill_tau( sample_mean_Im, sample_sqsum_Im, sample_stds_Im, tau_int_Im );
}


// compute the absolute iteration error Delta I_abs = max_ab{ timeav(Delta I^ab) }
void RunTimeData::compute_iteration_error( const CorrTen& CT, const CorrTen& new_CT )
{
    std::vector<RealType> time_average{};
    std::transform( CT.cbegin(), CT.cend(), new_CT.cbegin(), std::back_inserter(time_average), []( const Corr& C, const Corr& new_C )
    {
        Corr diff = C - new_C; // difference vector
        RealType abs_sum = std::accumulate( diff.cbegin(), diff.cend(), RealType{0.}, []( RealType sum, const RealType& d )
        { 
            return sum + std::abs(d);
        } ); // adding up the absolute difference
        return abs_sum/static_cast<RealType>(diff.size()); // divide by size and store
    } );
    absolute_iteration_errors.emplace_back( *std::max_element( time_average.cbegin(), time_average.cend() ) );
}

// finalize iteration step
void RunTimeData::finalize_iteration_step()
{
    size_t pc_space = 30;
    if( my_rank == 0 )
    {
        std::cout << "Iteration step finished:\n"
        << print::quantity_to_output_line( pc_space, "current absolute iteration error",    print::round_value_to_string( absolute_iteration_errors.back(), num_PrintDigits ) )
        << print::quantity_to_output_line( pc_space, "absolute iteration error threshold",  print::round_value_to_string( absolute_iteration_error_threshold, num_PrintDigits ) );
        if( !acceptance_rates.empty() )
        {
            std::cout << print::quantity_to_output_line( pc_space, "MH acceptance rate",    print::round_value_to_string( acceptance_rates.back(), num_PrintDigits ) );
        }
        if( error_method != "ar1" && !blocking_curve_sigma.empty() )
        {
            // plateau ratio = std at the largest block length / std at the smallest. ~1 means
            // the blocks are long enough and the error bar is trustworthy; >>1 means the run is
            // under-sampled for its autocorrelation and the reported error is a lower bound.
            const RealType plateau = ( blocking_curve_sigma.front() > RealType{0.} )
                ? blocking_curve_sigma.back() / blocking_curve_sigma.front() : RealType{1.};
            std::cout << print::quantity_to_output_line( pc_space, "blocking plateau ratio", print::round_value_to_string( plateau, num_PrintDigits ) );
            // worst-mixed Re correlation point: max tau_int (in pCN steps)
            RealType tau_max = RealType{0.};
            for( const auto& dir : tau_int_Re ) for( size_t t = 0; t < dir.size(); ++t ) tau_max = std::max( tau_max, dir[t] );
            std::cout << print::quantity_to_output_line( pc_space, "max tau_int (steps)", print::round_value_to_string( tau_max, num_PrintDigits ) );
        }
    }
    num_Iterations++;
}

// MPI-reduce the per-core MH acceptance counters and append the global rate.
void RunTimeData::record_mh_acceptance()
{
    long long local[2] = { static_cast<long long>(mh_accepted_count),
                           static_cast<long long>(mh_proposed_count) };
    long long global[2] = { 0, 0 };
    MPI_Allreduce( local, global, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
    RealType rate = ( global[1] > 0 ) ? static_cast<RealType>(global[0]) / static_cast<RealType>(global[1])
                                      : RealType{0.};
    acceptance_rates.emplace_back( rate );
    mh_accepted_count = 0;
    mh_proposed_count = 0;
}

// determine whether the iteration should be terminated
bool RunTimeData::terminate()
{
    if( !self_consistency )
    {
        termination = "no self consistency";
        return true;
    }
    else if( absolute_iteration_errors.back() < absolute_iteration_error_threshold )
    {
        print::print_R0( my_rank, "\033[1;32mTerminating regularly due to converged self-consistency.\033[0m\nThe iteration has been stopped and the current data will be stored regularly.\n" );
        termination = "by convergence";
        return true;
    }
    else if( num_Iterations+1 >= Iteration_Limit )
    {
        print::print_R0( my_rank, "\033[1;31mTerminating because I have reached the preset iteration limit of " + std::to_string(Iteration_Limit) + ".\033[0m\nThe iteration has been stopped and the current data will be stored regularly.\n" );
        termination = "by iteration limit";
        return true;
    }
    return false;
}

};


