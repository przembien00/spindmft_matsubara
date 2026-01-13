#include"Run_Time_Data.h"

#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<filesystem>
#include<Globals/MPI_Types.h>

#include<Standard_Algorithms/Standard_Algorithms.h>
namespace stda = Standard_Algorithms;

#include"Error_Handling.h"
namespace error = Run_Time_Data::Error_Handling;

namespace Run_Time_Data
{
// ============================================================================
// =========================== RUN TIME DATA CLASS ============================
// ============================================================================
// constructor
RunTimeData::RunTimeData(const ps::ParameterSpace& pspace, const int my_rank ):
    my_rank(                            my_rank ),
    eigenvalue_ratio_tolerance(         pspace.eigenvalue_ratio_tolerance ),
    seed(                               pspace.seed ),
    num_Samples(                        pspace.num_Samples ),
    num_SamplesPerCore(                 pspace.num_SamplesPerCore ),
    num_Cores(                          pspace.get_num_Cores() )
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

    // ...concerning the statistics:
    sample_sqsum = CluCorrTen{ pspace.compute_only_categories, pspace.symmetry_type, pspace.num_TimePoints };
    sample_stds  = CluCorrTen{ pspace.compute_only_categories, pspace.symmetry_type, pspace.num_TimePoints };
}

// process the eigenvalues, search for negative ones and compare the ratio to a predefined threshold
void RunTimeData::process_and_check_eigenvalues( mvgb::EigenValuesBlocks& EVBs )
{
    // 1.) compute the largest positive and the largest negative eigenvalues
    std::vector<RealType> mins{};
    std::transform( EVBs.cbegin(), EVBs.cend(), std::back_inserter(mins), 
    []( const auto & ev_block ) -> RealType
    {
        return *std::min_element( ev_block.cbegin(), ev_block.cend() );
    } );
    largest_negative_eigen_value = std::min( RealType{0.}, *std::min_element( mins.cbegin(), mins.cend() ) );

    std::vector<RealType> maxs{};
    std::transform( EVBs.cbegin(), EVBs.cend(), std::back_inserter(maxs), 
    []( const auto & ev_block ) -> RealType
    {
        return *std::max_element( ev_block.cbegin(), ev_block.cend() );
    } );
    largest_positive_eigen_value = std::max( RealType{0.}, *std::max_element( maxs.cbegin(), maxs.cend() ) );

    // 2.) compute the sum of positive and the sum of negative eigenvalues and their ratio
    positive_eigen_values_sum = double{0.};
    negative_eigen_values_sum = double{0.};
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
    negative_eigenvalue_ratio = negative_eigen_values_sum / positive_eigen_values_sum;

    // 3.) truncate negative eigevalues
    std::for_each( EVBs.begin(), EVBs.end(), [this]( auto & EVB )
    {
        std::for_each( EVB.begin(), EVB.end(), truncate_if_negative );
    } );

    // 4.) compare negative eigenvalue ratio to threshold
    eigenvalue_threshold_violated = std::abs(negative_eigenvalue_ratio) > eigenvalue_ratio_tolerance;
    if( eigenvalue_threshold_violated && my_rank == 0 ) // then print warning to terminal
    {
        std::cout << "\033[1;31mWarning: Negative eigenvalue ratio violates the tolerance according to " << std::abs(negative_eigenvalue_ratio) << " > " << eigenvalue_ratio_tolerance << "\033[0m\n";
        std::cout << "The simulation continues regularly.\n";
    }
}

std::string RunTimeData::get_seed_str()
{
    return seed;
}
size_t RunTimeData::get_num_SamplesPerCore() const
{
    return num_SamplesPerCore;
}
size_t RunTimeData::get_num_Samples() const
{
    return this->get_num_SamplesPerCore() * num_Cores;
}

// compute the standard deviation of a single correlation Monte-Carlo sample
void RunTimeData::compute_sample_stds( const CluCorrTen& sample_sum )
{
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
}

}