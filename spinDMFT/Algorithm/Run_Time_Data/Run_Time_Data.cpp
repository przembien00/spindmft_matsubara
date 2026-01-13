#include"Run_Time_Data.h"
#include<fstream>
#include<iomanip>

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
    sample_cov_Re   = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    sample_cov_Im   = CorrTen{ pspace.correlation_symmetry_type, pspace.num_TimePoints };
    Z_sqsum = RealType{0.};

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

template <class InputIterator1, class InputIterator2, class InputIterator3,
          class OutputIterator, class TrenaryOperation>
  OutputIterator transform3(InputIterator1 first1, InputIterator1 last1,
                            InputIterator2 first2, InputIterator3 first3, OutputIterator result,
                            TrenaryOperation trenary_op)
{
  while (first1 != last1) {
    *result = trenary_op(*first1, *first2, *first3);
    ++result; ++first1; ++first2; ++first3;
  }
  return result;
}

// compute the single sample standard deviations
void RunTimeData::compute_sample_stds( const CorrTen& sample_sum_Re, const CorrTen& sample_sum_Im, const RealType& Z )
{
    // RealType M = static_cast<RealType>( this->get_num_Samples() );

    transform3( sample_sum_Re.cbegin(), sample_sum_Re.cend(), sample_sqsum_Re.cbegin(), sample_cov_Re.cbegin(), sample_stds_Re.begin(), [&]( const Corr& sample_sum_C, const Corr& sample_sqsum_C, const Corr& sample_cov_CZ )
    {
        Corr std_C( sample_sum_C.size() );
        transform3( sample_sum_C.cbegin(), sample_sum_C.cend(), sample_sqsum_C.cbegin(), sample_cov_CZ.cbegin(), std_C.begin(), [&]( const auto& sample_sum, const auto& sample_sqsum, const auto& sample_cov )
        {
            return std::sqrt(std::abs( sample_sqsum + std::pow(sample_sum,2)*Z_sqsum - 2 * sample_sum * sample_cov ) / std::pow(Z,2));
        } );
        return std_C;
    } );

    transform3( sample_sum_Im.cbegin(), sample_sum_Im.cend(), sample_sqsum_Im.cbegin(), sample_cov_Im.cbegin(), sample_stds_Im.begin(), [&]( const Corr& sample_sum_C, const Corr& sample_sqsum_C, const Corr& sample_cov_CZ )
    {
        Corr std_C( sample_sum_C.size() );
        transform3( sample_sum_C.cbegin(), sample_sum_C.cend(), sample_sqsum_C.cbegin(), sample_cov_CZ.cbegin(), std_C.begin(), [&]( const auto& sample_sum, const auto& sample_sqsum, const auto& sample_cov )
        {
            return std::sqrt(std::abs( sample_sqsum + std::pow(sample_sum,2)*Z_sqsum - 2 * sample_sum * sample_cov ) / std::pow(Z,2));
        } );
        return std_C;
    } );

    sample_sqsum_Re = CorrTen{ sample_sum_Re.get_symmetry(), sample_sum_Re[0].size() }; // reset the square sum for the next iteration
    sample_sqsum_Im = CorrTen{ sample_sum_Re.get_symmetry(), sample_sum_Re[0].size() }; // reset the square sum for the next iteration
    sample_cov_Re = CorrTen{ sample_sum_Re.get_symmetry(), sample_sum_Re[0].size() }; // reset the covariance for the next iteration
    sample_cov_Im = CorrTen{ sample_sum_Re.get_symmetry(), sample_sum_Re[0].size() }; // reset the covariance for the next iteration
    Z_sqsum = RealType{0.}; // reset the Z square sum for the next iteration
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
    }
    num_Iterations++; 
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


