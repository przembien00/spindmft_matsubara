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
    sample_sqsum = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    sample_stds  = CluCorrTen{ pspace.correlation_categories, pspace.symmetry_type, pspace.num_TimePoints };
    if( pspace.adaptive_sample_size ){ adaptive_num_SamplesPerCore.emplace_back( num_SamplesPerCore ); }
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

// compute the standard deviation of a single correlation Monte-Carlo sample, adapt the sample size if desired
void RunTimeData::compute_and_process_sample_stds( const CluCorrTen& sample_sum )
{
    // compute sample stds:
    RealType M = static_cast<RealType>( this->get_num_Samples() );
    stda::for_3each( sample_stds.begin(), sample_stds.end(), sample_sum.cbegin(), sample_sqsum.cbegin(), 
    [&M]( auto& std_CT, const auto& av_CT, const auto& sqav_CT )
    {
        stda::for_3each( std_CT.begin(), std_CT.end(), av_CT.cbegin(), sqav_CT.cbegin(), 
        [&M]( auto& std_C, const auto& av_C, const auto& sqav_C )
        {
            stda::for_3each( std_C.begin()+1, std_C.end(), av_C.cbegin()+1, sqav_C.cbegin()+1, // at time zero the std is zero
            [&M]( RealType & std, const RealType & av, const RealType & sqav ) 
            {
                std = std::sqrt( std::abs( sqav / M - std::pow( av / M, 2 ) ) );
            } );
        } );
    } );

    // reset the square sum for the next iteration step:
    sample_sqsum = CluCorrTen{ sample_sum.get_site_pairs(), sample_sum[0].get_symmetry(), sample_sum[0][0].size() }; 

    // determine a new sample size for the next iteration step in case of adaptive settings
    if( adaptive_sample_size ) 
    {
        // determine the maximum standard deviation sigma_max = max_{ij} max_{ab} max_{t} sigma^ab_ij(t,0)
        std::vector<RealType> max_std_per_CCT{};
        std::transform( sample_stds.cbegin(), sample_stds.cend(), std::back_inserter(max_std_per_CCT), []( const CorrTen& std_CT )
        {
            std::vector<RealType> max_std_per_CT{};
            std::transform( std_CT.cbegin(), std_CT.cend(), std::back_inserter(max_std_per_CT), []( const Corr& std_C )
            {
                return *std::max_element( std_C.cbegin(), std_C.cend() );
            } );
            return *std::max_element( max_std_per_CT.cbegin(), max_std_per_CT.cend() );
        } );
        RealType max_std = *std::max_element( max_std_per_CCT.cbegin(), max_std_per_CCT.cend() );
        RealType max_std_of_sum = max_std / std::sqrt( static_cast<RealType>( this->get_num_Samples() ) ); // Sigma_max = sigma_max/sqrt(M)

        // check if the previous sample size was large enough and update the sample size
        size_t new_num_SamplesPerCore{};
        if( max_std_of_sum < statistical_error_tolerance ) // then leave the number of samples constant
        {
            sample_size_updated = false;
            new_num_SamplesPerCore = this->get_num_SamplesPerCore();
        }
        else // otherwise set the new sample size so that the threshold will be approximately fulfilled in the next iteration
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

    // compute the relative iteration error Delta I_rel = sqrt(M) * max_ij{ max_ab{ timeav(Delta I^ab_ij) / timeav(sigma^ab_ij) } }
    Vec max_rel_avdevs{};
    std::transform( sample_stds.cbegin(), sample_stds.cend(), deviation_CCT.cbegin(), std::back_inserter( max_rel_avdevs ), []( const CorrTen& std_CT, const CorrTen& dev_CT )
    {
        Vec rel_avdevs{};
        std::transform( std_CT.cbegin(), std_CT.cend(), dev_CT.cbegin(), std::back_inserter( rel_avdevs ), []( const Corr& std_C, const Corr& dev_C )
        {   
            return average_of_absolute(dev_C) / average_of_absolute(std_C);
        } );
        return *max_element( rel_avdevs.cbegin(), rel_avdevs.cend() );
    } );
    relative_iteration_error_list.emplace_back( *max_element( max_rel_avdevs.cbegin(), max_rel_avdevs.cend() ) * std::sqrt(get_num_Samples()) );

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
    }

    num_Iterations++; 
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