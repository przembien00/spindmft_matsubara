#pragma once

#include<vector>
#include<string>
#include<functional>
#include<iostream>
#include<blaze/Math.h>
#include<Globals/Types.h>
#include<Observables/Tensors.h>
#include<Observables/Correlations.h>
#include"../Parameter_Space/Parameter_Space.h"

namespace spinDMFT::Run_Time_Data
{

namespace ps = spinDMFT::Parameter_Space;
namespace ten = Observables::Tensors;
namespace corr = Observables::Correlations;
  
using EigenValuesList = std::vector<blaze::StaticVector<RealType,3UL,blaze::columnVector>>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;

// ============================================================================
// =========================== RUN TIME DATA CLASS ============================
// ============================================================================
class RunTimeData
{ 
 public:
   // CONSTRUCTORS
   RunTimeData() = default;
   RunTimeData( const ps::ParameterSpace& pspace, const int my_rank );

   // PUBLIC MEMBERS
   // ...concerning the eigenvalues
   std::vector<RealType> negative_eigenvalue_ratios{};     // over the iterations : ratio relative to the positive_eigenvalues_sum

   // ...concerning the statistics
   size_t generated_seed{};
   CorrTen sample_sqsum_Re{};    // M*<x^2> (<x> is not saved in runtimedata)
   CorrTen sample_sqsum_Im{};    // M*<x^2> (<x> is not saved in runtimedata)
   RealType Z_sqsum{};          // M*<Z^2> (<Z> is not saved in runtimedata)
   CorrTen sample_cov_Re{};    // covariance between real part and Z
   CorrTen sample_cov_Im{};    // covariance between imag part and Z
   CorrTen sample_stds_Re{};     // sqrt(<x^2> - <x>^2)
   CorrTen sample_stds_Im{};     // sqrt(<x^2> - <x>^2)

   // ...concerning the iteration
   size_t num_Iterations{};
   std::vector<RealType> absolute_iteration_errors{};
   std::string termination{};

   // PUBLIC METHODS
   // ...concerning the eigenvalues
   void process_and_check_eigenvalues( EigenValuesList& eigenvalues );

   // ...concerning the statistics
   std::string get_seed_str();
   void compute_sample_stds( const CorrTen& sample_sum_Re, const CorrTen& sample_sum_Im, const RealType& Z );
   size_t get_num_SamplesPerCore() const;
   size_t get_num_Samples() const;

   // ...concerning the iteration
   void compute_iteration_error( const CorrTen& CT, const CorrTen& new_CT );
   void finalize_iteration_step();
   bool terminate();

 private:
   const int my_rank{};

   // MEMBERS IMPORTED FROM A PARAMETER SPACE (AND THEREFORE CONST)
   const size_t num_PrintDigits{};

   // ...concerning the eigenvalues
   const RealType critical_eigenvalue_ratio{ 0.01 }; // should be set by the statistics
   std::function<void(RealType&)> truncate_if_negative;

   // ...concerning the statistics
   const std::string seed{};
   const size_t num_Samples{};           // (initial) number of samples
   const size_t num_SamplesPerCore{};    // (initial) number of samples per core
   const size_t num_Cores{};

   // ...concerning the iteration
   const bool self_consistency{};
   const RealType absolute_iteration_error_threshold{};
   const size_t Iteration_Limit{};
};

}
