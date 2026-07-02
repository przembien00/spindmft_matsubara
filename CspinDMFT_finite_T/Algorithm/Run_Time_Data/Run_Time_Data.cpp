#include"Run_Time_Data.h"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<filesystem>
#include<Globals/MPI_Types.h>

#include<Standard_Algorithms/Standard_Algorithms.h>
namespace stda = Standard_Algorithms;

#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;

#include"Error_Handling.h"
namespace error = Run_Time_Data::Error_Handling;

namespace Run_Time_Data
{
// ============================================================================
// =========================== RUN TIME DATA CLASS ============================
// ============================================================================
// constructor
RunTimeData::RunTimeData(const ps::ParameterSpace& pspace, const int my_rank ):
    num_Iterations(                     0 ),
    my_rank(                            my_rank ),
    Run_ID(                             pspace.Run_ID ),
    num_PrintDigits(                    pspace.num_PrintDigits ),
    eigenvalue_ratio_tolerance(         pspace.eigenvalue_ratio_tolerance ),
    seed(                               pspace.seed ),
    num_Samples(                        pspace.num_Samples ),
    num_SamplesPerCore(                 pspace.num_SamplesPerCore ),
    num_SamplesPerSet(                  pspace.num_SamplesPerSet ),
    num_Cores(                          pspace.get_num_Cores() ),
    adaptive_sample_size(               pspace.adaptive_sample_size ),
    statistical_error_tolerance(        pspace.statistical_error_tolerance ),
    Iteration_Limit(                    pspace.Iteration_Limit ),
    iteration_error_mode(               pspace.iteration_error_mode ),
    absolute_iteration_error_tolerance( pspace.absolute_iteration_error_tolerance ),
    relative_iteration_error_tolerance( pspace.relative_iteration_error_tolerance )
{
    // INITIALIZATION
    // ...concerning the eigenvalues:
    if( pspace.truncation_scheme_negative_eigenvalues == "set_zero" )
    {
        truncate_if_negative = []( RealType& eig ) -> void
        { 
            if( eig < RealType{0.} ){ eig = RealType{0.}; }
        };
    }
    else if( pspace.truncation_scheme_negative_eigenvalues == "abs" )
    {
        truncate_if_negative = []( RealType& eig ) -> void
        { 
            if( eig < RealType{0.} ){ eig = std::abs( eig ); }
        };
    }
    else
    {
        error::EIGENVALUE_TRUNCATION( pspace.truncation_scheme_negative_eigenvalues, __PRETTY_FUNCTION__ );
    }

    // ...concerning the statistics
    sample_sqsum_Re = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    sample_sqsum_Im = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    sample_stds_Re  = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    sample_stds_Im  = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    tau_int_Re      = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    tau_int_Im      = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    spin_expval_sqsum = std::vector<FieldVector>( pspace.num_Spins, FieldVector{0.,0.,0.} );
    spin_expval_stds  = std::vector<FieldVector>( pspace.num_Spins, FieldVector{0.,0.,0.} );
    if( pspace.adaptive_sample_size ){ adaptive_num_SamplesPerCore.emplace_back( num_SamplesPerCore ); }

    // ...concerning the batch-means (blocking) error estimation
    error_method = pspace.error_method;
    num_blocks   = pspace.mh_num_blocks;
}

// process the eigenvalues, search for negative ones and compare the ratio to a predefined threshold
void RunTimeData::process_and_check_eigenvalues( mvgb::EigenValuesBlocks& EVBs )
{
    // 1.) compute the largest positive and the largest negative eigenvalues
    Vec mins{};
    std::transform( EVBs.cbegin(), EVBs.cend(), std::back_inserter(mins), 
    []( const auto & ev_block ) -> RealType
    {
        return *std::min_element( ev_block.cbegin(), ev_block.cend() );
    } );
    largest_negative_eigen_value_list.emplace_back( std::min(RealType{0.}, *std::min_element(mins.cbegin(), mins.cend())) );

    Vec maxs{};
    std::transform( EVBs.cbegin(), EVBs.cend(), std::back_inserter(maxs), 
    []( const auto & ev_block ) -> RealType
    {
        return *std::max_element( ev_block.cbegin(), ev_block.cend() );
    } );
    largest_positive_eigen_value_list.emplace_back( std::max(RealType{0.}, *std::max_element(maxs.cbegin(), maxs.cend())) );

    // 2.) compute the sum of positive and the sum of negative eigenvalues and their ratio
    RealType positive_eigen_values_sum{};
    RealType negative_eigen_values_sum{};
    for( const auto & ev_block : EVBs )
    {
        for( const auto & ev : ev_block )
        {
            if( ev > RealType{0.} )
            {
                positive_eigen_values_sum += ev;
            }
            else 
            {
                negative_eigen_values_sum += ev;
            }
        }
    }
    positive_eigen_values_sum_list.emplace_back( positive_eigen_values_sum );
    negative_eigen_values_sum_list.emplace_back( negative_eigen_values_sum );
    negative_eigenvalue_ratio_list.emplace_back( negative_eigen_values_sum / positive_eigen_values_sum );

    // 3.) truncate negative eigevalues
    std::for_each( EVBs.begin(), EVBs.end(), [this]( auto & EVB )
    {
        std::for_each( EVB.begin(), EVB.end(), truncate_if_negative );
    } );

    // 4.) compare negative eigenvalue ratio to threshold
    eigenvalue_threshold_violated_list.emplace_back( std::abs(negative_eigenvalue_ratio_list.back()) > eigenvalue_ratio_tolerance );
    if( eigenvalue_threshold_violated_list.back() && my_rank == 0 ) // then print warning to terminal
    {
        std::cout << "\033[1;31mWarning: Negative eigenvalue ratio violates the tolerance according to " << std::abs(negative_eigenvalue_ratio_list.back()) << " > " << eigenvalue_ratio_tolerance << "\033[0m\n";
        std::cout << "The simulation continues regularly.\n";
    }
}

std::string RunTimeData::get_seed_str()
{
    return seed;
}
size_t RunTimeData::get_num_SamplesPerCore() const
{
    if( adaptive_sample_size ) 
    {
        return adaptive_num_SamplesPerCore.back(); // current sample size
    }
    return num_SamplesPerCore; // initially inserted, constant sample size
}
size_t RunTimeData::get_num_Samples() const
{
    return get_num_SamplesPerCore() * num_Cores;
}
size_t RunTimeData::get_num_SetsPerCore() const
{
    return (size_t) std::ceil( static_cast<double>(get_num_SamplesPerCore()) / static_cast<double>(num_SamplesPerSet) ); 
}

// Estimate the lag-1 autocorrelation rho1 of the pCN chain from the latest acceptance rate
// and the step size, and return the AR(1) variance-inflation factor sqrt( (1+rho1)/(1-rho1) ).
// A rejected pCN step copies the state, an accepted one retains a fraction sqrt(1-beta^2) of
// it, so rho1 ~ 1 - a*(1 - sqrt(1-beta^2)). This is a crude single-exponential estimate: it
// assumes one relaxation timescale and underestimates the true autocorrelation when the chain
// has a slow / trapped mode. rho1 is clamped below 1 so the factor stays finite when the chain
// barely moves.
RealType RunTimeData::pcn_autocorrelation_factor( const RealType pcn_step_size ) const
{
    if( acceptance_rates.empty() ) return RealType{1.};
    const RealType a = acceptance_rates.back();
    const RealType retained = std::sqrt( std::max( RealType{0.}, RealType{1.} - pcn_step_size * pcn_step_size ) );
    RealType rho1 = RealType{1.} - a * ( RealType{1.} - retained );
    rho1 = std::min( std::max( rho1, RealType{0.} ), RealType{1.} - RealType{1e-12} );
    return std::sqrt( ( RealType{1.} + rho1 ) / ( RealType{1.} - rho1 ) );
}

// pCN estimators are plain chain means; the i.i.d. standard error of the mean is
// stderr( <o> ) = sqrt( ( <o^2> - <o>^2 ) / N ). Because the samples are autocorrelated this
// is scaled by the AR(1) factor from pcn_autocorrelation_factor. sample_mean_* hold the
// already normalized means and sample_sqsum_* the MPI-reduced raw sum_i o_i^2.
void RunTimeData::compute_and_process_sample_stds( const CluCorrTen& sample_mean_Re, const CluCorrTen& sample_mean_Im, const RealType N, const RealType pcn_step_size )
{
    const RealType autocorr_factor = pcn_autocorrelation_factor( pcn_step_size );
    auto per_point_stderr = [N, autocorr_factor]( const RealType sample_mean, const RealType sample_sqsum ) -> RealType
    {
        const RealType var = std::abs( sample_sqsum / N - sample_mean * sample_mean );
        return autocorr_factor * std::sqrt( var / N );
    };

    auto fill_stderr = [&]( const CluCorrTen& mean_CCT, CluCorrTen& sqsum_CCT, CluCorrTen& std_CCT )
    {
        std_CCT.iterate2( sqsum_CCT, [&]( auto& std_CT, auto& sqsum_CT, const auto& ij )
        {
            const auto& mean_CT = mean_CCT( ij[0], ij[1] );
            std_CT.iterate2( sqsum_CT, [&]( auto& std_C, auto& sqsum_C, const auto& alphabeta )
            {
                const auto mean_C = mean_CT( alphabeta[0], alphabeta[1] );
                for( size_t tau_index = 0; tau_index < std_C.size(); ++tau_index )
                {
                    std_C[tau_index] = per_point_stderr( mean_C[tau_index], sqsum_C[tau_index] );
                }
            } );
        } );
    };

    fill_stderr( sample_mean_Re, sample_sqsum_Re, sample_stds_Re );
    fill_stderr( sample_mean_Im, sample_sqsum_Im, sample_stds_Im );

    // reset the per-iteration accumulators
    sample_sqsum_Re = CluCorrTen{ sample_mean_Re.get_site_pairs(), sample_mean_Re[0].get_symmetry(), sample_mean_Re[0][0].size() };
    sample_sqsum_Im = CluCorrTen{ sample_mean_Im.get_site_pairs(), sample_mean_Im[0].get_symmetry(), sample_mean_Im[0][0].size() };

    // determine a new sample size for the next iteration step in case of adaptive settings
    if( adaptive_sample_size )
    {
        std::vector<RealType> max_std_per_CCT{};
        std::transform( sample_stds_Re.cbegin(), sample_stds_Re.cend(), std::back_inserter(max_std_per_CCT), []( const CorrTen& std_CT )
        {
            std::vector<RealType> max_std_per_CT{};
            std::transform( std_CT.cbegin(), std_CT.cend(), std::back_inserter(max_std_per_CT), []( const Corr& std_C )
            {
                return *std::max_element( std_C.cbegin(), std_C.cend() );
            } );
            return *std::max_element( max_std_per_CT.cbegin(), max_std_per_CT.cend() );
        } );
        RealType max_std_of_sum = *std::max_element( max_std_per_CCT.cbegin(), max_std_per_CCT.cend() );

        size_t new_num_SamplesPerCore{};
        if( max_std_of_sum < statistical_error_tolerance )
        {
            sample_size_updated = false;
            new_num_SamplesPerCore = this->get_num_SamplesPerCore();
        }
        else
        {
            sample_size_updated = true;
            RealType ratio = max_std_of_sum / statistical_error_tolerance;
            new_num_SamplesPerCore = (size_t) ( static_cast<RealType>( this->get_num_SamplesPerCore() ) * pow(ratio,2) + 1.0 );
        }
        adaptive_num_SamplesPerCore.emplace_back( new_num_SamplesPerCore );
    }
}

void RunTimeData::compute_and_process_spin_expval_stds( const std::vector<FieldVector>& spin_mean, const RealType N, const RealType pcn_step_size )
{
    const RealType autocorr_factor = pcn_autocorrelation_factor( pcn_step_size );
    for( size_t site = 0; site < spin_mean.size(); ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            const RealType mean = spin_mean[site][alpha];
            const RealType sqsum = spin_expval_sqsum[site][alpha];
            const RealType var = std::abs( sqsum / N - mean * mean );
            spin_expval_stds[site][alpha] = autocorr_factor * std::sqrt( var / N );
        }
    }

    spin_expval_sqsum = std::vector<FieldVector>( spin_mean.size(), FieldVector{0.,0.,0.} );
}

// Dispatch to the configured standard-error estimator. "blocking" (batch means, default) is
// robust to autocorrelation; "ar1" is the legacy acceptance-based single-mode factor that
// underestimates the error for slow/trapped chains. sample_mean_* hold <o> (already divided by
// N in main.cpp); for the AR(1) path sample_sqsum_* hold the MPI-reduced sum_i o_i^2.
void RunTimeData::compute_sample_stds( const CluCorrTen& sample_mean_Re, const CluCorrTen& sample_mean_Im, const std::vector<FieldVector>& spin_mean, const RealType N, const RealType pcn_step_size )
{
    if( error_method == "ar1" )
    {
        compute_and_process_spin_expval_stds( spin_mean, N, pcn_step_size );
        compute_and_process_sample_stds( sample_mean_Re, sample_mean_Im, N, pcn_step_size );
    }
    else
    {
        compute_sample_stds_blocking( sample_mean_Re, sample_mean_Im, N );
    }
}

// Allocate the batch-mean blocks for one production sweep and zero all per-iteration
// accumulators. block_length = num_SamplesPerCore / num_blocks; any remainder samples
// (< num_blocks) still enter the global mean but not the variance estimate.
void RunTimeData::init_blocks( const CluCorrTen& template_corr, size_t num_SamplesPerCore )
{
    const auto& site_pairs = template_corr.get_site_pairs();
    const char  sym        = template_corr[0].get_symmetry();
    const size_t ntp       = template_corr[0][0].size();
    const size_t num_Spins = spin_expval_sqsum.size();
    if( num_blocks == 0 ) num_blocks = 1;
    block_length = num_SamplesPerCore / num_blocks;
    if( block_length == 0 ) { block_length = 1; num_blocks = num_SamplesPerCore; }
    blocks_filled = 0;

    const CluCorrTen zero_corr{ site_pairs, sym, ntp };
    const std::vector<FieldVector> zero_spin( num_Spins, FieldVector{0.,0.,0.} );
    cur_block_Re   = zero_corr;
    cur_block_Im   = zero_corr;
    cur_block_spin = zero_spin;
    block_means_Re.assign(   num_blocks, zero_corr );
    block_means_Im.assign(   num_blocks, zero_corr );
    block_means_spin.assign( num_blocks, zero_spin );

    // also reset the AR(1) accumulators (kept available for errmethod=ar1 and the tau_int recovery)
    sample_sqsum_Re   = zero_corr;
    sample_sqsum_Im   = zero_corr;
    spin_expval_sqsum = zero_spin;
}

// Turn the current block's running sums into block means, store them, and reset the running
// sums for the next block. Excess calls (beyond num_blocks) are ignored.
void RunTimeData::close_block()
{
    if( blocks_filled >= num_blocks ) return;
    const RealType inv = RealType{1.} / static_cast<RealType>( block_length );
    for( size_t p = 0; p < cur_block_Re.size(); ++p )
    {
        for( size_t d = 0; d < cur_block_Re[p].size(); ++d )
        {
            corr::CorrelationVector& cur_re = cur_block_Re[p][d];
            corr::CorrelationVector& cur_im = cur_block_Im[p][d];
            corr::CorrelationVector& dst_re = block_means_Re[blocks_filled][p][d];
            corr::CorrelationVector& dst_im = block_means_Im[blocks_filled][p][d];
            for( size_t t = 0; t < cur_re.size(); ++t )
            {
                dst_re[t] = cur_re[t] * inv; cur_re[t] = RealType{0.};
                dst_im[t] = cur_im[t] * inv; cur_im[t] = RealType{0.};
            }
        }
    }
    for( size_t site = 0; site < cur_block_spin.size(); ++site )
    {
        for( size_t alpha = 0; alpha < 3; ++alpha )
        {
            block_means_spin[blocks_filled][site][alpha] = cur_block_spin[site][alpha] * inv;
            cur_block_spin[site][alpha] = RealType{0.};
        }
    }
    ++blocks_filled;
}

// Batch-means standard error of the global mean. Each core contributes num_blocks block means;
// these are pooled across cores via MPI_Allgather, and the error of the global mean is
// std(block means)/sqrt(num_blocks*num_cores) per correlation point. Valid when
// block_length >> autocorrelation time (see the blocking-curve diagnostic).
void RunTimeData::compute_sample_stds_blocking( const CluCorrTen& sample_mean_Re, const CluCorrTen& sample_mean_Im, RealType N )
{
    int num_cores = 1;
    MPI_Comm_size( MPI_COMM_WORLD, &num_cores );
    const size_t nc = static_cast<size_t>( num_cores );

    // number of correlation points (sum over site pairs of sum over direction pairs of tau-length)
    size_t P = 0;
    for( auto& CT : sample_stds_Re ) for( auto& C : CT ) P += C.size();
    const size_t num_Spins = spin_expval_stds.size();

    // flatten this core's block means in a fixed (site-pair, direction-pair, tau) order
    std::vector<RealType> local_Re( num_blocks * P ), local_Im( num_blocks * P ), local_spin( num_blocks * num_Spins * 3 );
    for( size_t b = 0; b < num_blocks; ++b )
    {
        size_t p = 0;
        for( auto& CT : block_means_Re[b] ) for( auto& C : CT ) for( size_t t = 0; t < C.size(); ++t ) local_Re[b * P + (p++)] = C[t];
        p = 0;
        for( auto& CT : block_means_Im[b] ) for( auto& C : CT ) for( size_t t = 0; t < C.size(); ++t ) local_Im[b * P + (p++)] = C[t];
        for( size_t site = 0; site < num_Spins; ++site ) for( size_t a = 0; a < 3; ++a ) local_spin[b * num_Spins * 3 + site * 3 + a] = block_means_spin[b][site][a];
    }

    // pool all blocks from all cores
    std::vector<RealType> all_Re( nc * num_blocks * P ), all_Im( nc * num_blocks * P ), all_spin( nc * num_blocks * num_Spins * 3 );
    MPI_Allgather( local_Re.data(),   num_blocks * P,             MPI_REALTYPE, all_Re.data(),   num_blocks * P,             MPI_REALTYPE, MPI_COMM_WORLD );
    MPI_Allgather( local_Im.data(),   num_blocks * P,             MPI_REALTYPE, all_Im.data(),   num_blocks * P,             MPI_REALTYPE, MPI_COMM_WORLD );
    MPI_Allgather( local_spin.data(), num_blocks * num_Spins * 3, MPI_REALTYPE, all_spin.data(), num_blocks * num_Spins * 3, MPI_REALTYPE, MPI_COMM_WORLD );

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

    auto fill = [&]( const std::vector<RealType>& pool, CluCorrTen& out )
    {
        size_t p = 0;
        for( auto& CT : out ) for( auto& C : CT ) for( auto it = C.begin(); it != C.end(); ++it, ++p ) *it = plateau( pool, P, p );
    };
    fill( all_Re, sample_stds_Re );
    fill( all_Im, sample_stds_Im );
    for( size_t site = 0; site < num_Spins; ++site ) for( size_t a = 0; a < 3; ++a ) spin_expval_stds[site][a] = plateau( all_spin, num_Spins * 3, site * 3 + a );

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
            if( sg == RealType{0.} ) continue;   // skip points with no signal
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
    auto fill_tau = [N]( const CluCorrTen& mean, const CluCorrTen& sqsum, const CluCorrTen& stds, CluCorrTen& out )
    {
        for( size_t p = 0; p < out.size(); ++p )
        {
            for( size_t d = 0; d < out[p].size(); ++d )
            {
                const corr::CorrelationVector& m  = mean[p][d];
                const corr::CorrelationVector& q  = sqsum[p][d];
                const corr::CorrelationVector& sd = stds[p][d];
                corr::CorrelationVector& o        = out[p][d];
                for( size_t t = 0; t < o.size(); ++t )
                {
                    const RealType s2 = std::abs( q[t] / N - m[t] * m[t] ); // marginal sample variance
                    o[t] = ( s2 > RealType{0.} ) ? RealType{0.5} * sd[t] * sd[t] * N / s2 : RealType{0.};
                }
            }
        }
    };
    fill_tau( sample_mean_Re, sample_sqsum_Re, sample_stds_Re, tau_int_Re );
    fill_tau( sample_mean_Im, sample_sqsum_Im, sample_stds_Im, tau_int_Im );

    // determine a new sample size for the next iteration step in case of adaptive settings
    if( adaptive_sample_size )
    {
        RealType max_std_of_sum = RealType{0.};
        for( auto& CT : sample_stds_Re ) for( auto& C : CT ) for( auto& v : C ) max_std_of_sum = std::max( max_std_of_sum, v );

        size_t new_num_SamplesPerCore{};
        if( max_std_of_sum < statistical_error_tolerance )
        {
            sample_size_updated = false;
            new_num_SamplesPerCore = this->get_num_SamplesPerCore();
        }
        else
        {
            sample_size_updated = true;
            RealType ratio = max_std_of_sum / statistical_error_tolerance;
            new_num_SamplesPerCore = (size_t) ( static_cast<RealType>( this->get_num_SamplesPerCore() ) * pow(ratio,2) + 1.0 );
        }
        adaptive_num_SamplesPerCore.emplace_back( new_num_SamplesPerCore );
    }
}

// compute the absolute time average of a correlation
RealType average_of_absolute( const Corr& C )
{
    return std::accumulate( C.cbegin(), C.cend(), RealType{0.}, [](RealType sum, const RealType& value) 
    {
        return sum + std::abs(value);
    } ) / static_cast<RealType>(C.size());
}

// compute the iteration error, i.e., the deviation between current and previous iteration step
void RunTimeData::compute_iteration_error( const CluCorrTen& new_CCT, const CluCorrTen& CCT )
{
    // compute the deviation CCT 
    CluCorrTen deviation_CCT{CCT.get_site_pairs()};
    std::transform( new_CCT.cbegin(), new_CCT.cend(), CCT.cbegin(), deviation_CCT.begin(), []( const CorrTen& new_CT, const CorrTen& CT )
    {
        CorrTen deviation_CT{CT.get_symmetry(), CT[0].get_num_TimePoints()};
        std::transform( new_CT.cbegin(), new_CT.cend(), CT.cbegin(), deviation_CT.begin(), []( const Corr& new_C, const Corr& C )
        {
            return new_C - C;
        } );
        return deviation_CT;
    } );

    // compute the relative iteration error Delta I_rel = max_ij{ max_ab{ timeav(Delta I^ab_ij) / timeav(sigma^ab_ij) } }
    // The sample standard deviations already carry the 1/sqrt(M) normalization.
    Vec max_rel_avdevs{};
    std::transform( sample_stds_Re.cbegin(), sample_stds_Re.cend(), deviation_CCT.cbegin(), std::back_inserter( max_rel_avdevs ), []( const CorrTen& std_CT, const CorrTen& dev_CT )
    {
        Vec rel_avdevs{};
        std::transform( std_CT.cbegin(), std_CT.cend(), dev_CT.cbegin(), std::back_inserter( rel_avdevs ), []( const Corr& std_C, const Corr& dev_C )
        {   
            return average_of_absolute(dev_C) / average_of_absolute(std_C);
        } );
        return *max_element( rel_avdevs.cbegin(), rel_avdevs.cend() );
    } );
    relative_iteration_error_list.emplace_back( *max_element( max_rel_avdevs.cbegin(), max_rel_avdevs.cend() ) );

    // compute the absolute iteration error Delta I_abs = max_ij{ max_ab{ timeav(Delta I^ab_ij) } }
    Vec max_abs_devs{};
    std::transform( deviation_CCT.cbegin(), deviation_CCT.cend(), std::back_inserter( max_abs_devs ), []( const CorrTen& dev_CT )
    {
        Vec abs_avdevs{};
        std::transform( dev_CT.cbegin(), dev_CT.cend(), std::back_inserter( abs_avdevs ), []( const Corr& dev_C )
        {   
            return average_of_absolute(dev_C);
        } );
        return *max_element( abs_avdevs.cbegin(), abs_avdevs.cend() );
    } );
    absolute_iteration_error_list.emplace_back( *max_element( max_abs_devs.cbegin(), max_abs_devs.cend() ) );
}

// produce some output and increment the iteration counter
void RunTimeData::finalize_iteration_step()
{
    size_t pc_space = 30;
    if( my_rank == 0 )
    {
        std::cout << "Iteration step finished:\n";
        if( adaptive_sample_size )
        {
            if( sample_size_updated )
            {
                std::cout << "total sample size is updated to " << std::to_string(get_num_Samples()) << "\n";
            }
            else
            {
                std::cout << "total sample size remains constant at " << std::to_string(get_num_Samples()) << "\n";
            }
        }
        std::cout << print::quantity_to_output_line( pc_space, "current relative iteration error", print::round_value_to_string( relative_iteration_error_list.back(), num_PrintDigits ) )
        << print::quantity_to_output_line( pc_space, "relative iteration error tolerance" + regarded("relative"), print::round_value_to_string( relative_iteration_error_tolerance, num_PrintDigits ) )
        << print::quantity_to_output_line( pc_space, "current absolute iteration error", print::round_value_to_string( absolute_iteration_error_list.back(), num_PrintDigits ) )
        << print::quantity_to_output_line( pc_space, "absolute iteration error tolerance" + regarded("absolute"), print::round_value_to_string( absolute_iteration_error_tolerance, num_PrintDigits ) );
        if( !acceptance_rates.empty() )
        {
            std::cout << print::quantity_to_output_line( pc_space, "pCN acceptance rate", print::round_value_to_string( acceptance_rates.back(), num_PrintDigits ) );
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
            for( auto& CT : tau_int_Re ) for( auto& C : CT ) for( auto& v : C ) tau_max = std::max( tau_max, v );
            std::cout << print::quantity_to_output_line( pc_space, "max tau_int (steps)", print::round_value_to_string( tau_max, num_PrintDigits ) );
        }
    }

    num_Iterations++;
}

// MPI-reduce the per-core pCN acceptance counters and append the global acceptance rate.
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

// check whether the results are converged
bool RunTimeData::is_converged() const
{
    if( iteration_error_mode == "relative" )
    {
        return relative_iteration_error_list.back() < relative_iteration_error_tolerance;
    }
    else if( iteration_error_mode == "absolute" )
    {
        return absolute_iteration_error_list.back() < absolute_iteration_error_tolerance;
    }
    else if( iteration_error_mode == "either" )
    {
        return relative_iteration_error_list.back() < relative_iteration_error_tolerance
            || absolute_iteration_error_list.back() < absolute_iteration_error_tolerance;
    }
    else
    {
        error::ITERATION_ERROR_MODE( iteration_error_mode, __PRETTY_FUNCTION__ ); 
        return true;
    }
}

// the iteration can be terminated manually by creating a signal file named terminate_<Run_ID> in the directory of the executable
// this function searches for the signal file at the end of each iteration step
bool RunTimeData::signal_file_exists() const
{
    bool exists = false;
    if( my_rank == 0 ) // core 0: search for signal file
    {
        exists = std::filesystem::exists("terminate_" + std::to_string(Run_ID));
        if( exists )
        {
            std::filesystem::remove("terminate_" + std::to_string(Run_ID)); // remove termination file
        }
    }
    int manual_termination_int = static_cast<int>(exists);
    MPI_Bcast( &manual_termination_int, 1, MPI_INT, 0, MPI_COMM_WORLD ); // broadcast manual_termination_int to all the other ranks
    exists = static_cast<bool>(manual_termination_int);
    return exists;
}

// check all possible termination conditions: convergence, manual termination, iteration limit
bool RunTimeData::terminate()
{
    if( is_converged() )
    {
        print::print_R0( my_rank, "\033[1;32mTerminating regularly due to converged self-consistency.\033[0m\nThe iteration has been stopped and the current data will be stored regularly.\n" );
        termination = "by convergence";
        return true;
    }
    else if( signal_file_exists() ) 
    {
        print::print_R0( my_rank, "\033[1;31mReceived termination signal at rank 0 and broadcasted it to all other cores.\033[0m\nThe iteration has been stopped and the current data will be stored regularly.\n" );
        termination = "by hand";
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

}
