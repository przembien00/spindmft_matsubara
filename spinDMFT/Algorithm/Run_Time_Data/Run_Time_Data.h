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
   CorrTen sample_sqsum_Re{};       // sum_i o_i^2 (raw, MPI-reduced before compute_sample_stds)
   CorrTen sample_sqsum_Im{};
   CorrTen sample_stds_Re{};        // standard error of the mean, sqrt( (<o^2> - <o>^2) / N )
   CorrTen sample_stds_Im{};
   // ...spin expectation value statistics (S_x, S_y, S_z)
   FieldVector spin_expval_sqsum{}; // sum_i s_i^2 (raw, MPI-reduced before compute_sample_stds)
   FieldVector spin_expval_stds{};  // standard error of the mean for <S^alpha>

   // ...concerning the iteration
   size_t num_Iterations{};
   std::vector<RealType> absolute_iteration_errors{};
   std::string termination{};

   // ...concerning the Metropolis-Hastings chain diagnostics
   // mh_accepted_count and mh_proposed_count are accumulated per MPI core during one
   // iteration; record_mh_acceptance() reduces them across MPI ranks and appends the
   // resulting global acceptance ratio to acceptance_rates.
   size_t mh_accepted_count{};
   size_t mh_proposed_count{};
   std::vector<RealType> acceptance_rates{};
   void record_mh_acceptance();

   // PUBLIC METHODS
   // ...concerning the eigenvalues
   void process_and_check_eigenvalues( EigenValuesList& eigenvalues );

   // ...concerning the statistics
   std::string get_seed_str();
   // Given the normalized means (sample_mean_*) and the global sample count N, compute the
   // standard error of the mean for each correlation and each spin component. The pCN samples
   // are autocorrelated, so the i.i.d. stderr is scaled by the AR(1) variance-inflation factor
   // sqrt( (1+rho1)/(1-rho1) ) (see pcn_autocorrelation_factor); pcn_step_size is the pCN beta
   // used to estimate rho1.
   void compute_sample_stds( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im,
                             const FieldVector& spin_mean, RealType N, RealType pcn_step_size );

   // AR(1) autocorrelation correction factor sqrt( (1+rho1)/(1-rho1) ) for the pCN stds.
   // rho1 is estimated from the latest acceptance rate a and the pCN step size beta via the
   // pCN move: a rejected step copies the state, an accepted one retains a fraction
   // sqrt(1-beta^2) of it, giving rho1 ~ 1 - a*(1 - sqrt(1-beta^2)). This is a crude single-
   // exponential estimate: it assumes one relaxation timescale and UNDERESTIMATES the true
   // autocorrelation when the chain has a slow / trapped mode.
   RealType pcn_autocorrelation_factor( RealType pcn_step_size ) const;
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
