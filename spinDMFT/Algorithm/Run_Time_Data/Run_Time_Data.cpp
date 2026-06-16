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
    spin_expval_sqsum = FieldVector{0.,0.,0.};
    spin_expval_stds  = FieldVector{0.,0.,0.};

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

// Standard error of the mean: stderr( <o> ) = sqrt( ( <o^2> - <o>^2 ) / N ), scaled by the
// AR(1) autocorrelation factor from pcn_autocorrelation_factor since the pCN samples are
// correlated. `sample_mean_*` holds <o> (already divided by N in main.cpp via normalize()),
// `sample_sqsum_*` holds sum_i o_i^2 (MPI-reduced, raw).
void RunTimeData::compute_sample_stds( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im,
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


