#pragma once

#include<vector>
#include<string>
#include<functional>
#include<iostream>
#include<blaze/Math.h>
#include<Globals/Types.h>
#include<Observables/Tensors.h>
#include<Observables/Correlations.h>
#include<Observables/Magnetization.h>
#include"../Parameter_Space/Parameter_Space.h"

namespace spinDMFT::Run_Time_Data
{

namespace ps = spinDMFT::Parameter_Space;
namespace ten = Observables::Tensors;
namespace corr = Observables::Correlations;
namespace mag = Observables::Magnetization;

using EigenValuesList = std::vector<blaze::StaticVector<RealType,3UL,blaze::columnVector>>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using MagVec = mag::MagnetizationVector;

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
   CorrTen sample_sqsum_Re{};       // sum_i o_i^2 (raw, MPI-reduced; used by the legacy AR(1) estimator)
   CorrTen sample_sqsum_Im{};
   CorrTen sample_stds_Re{};        // standard error of the mean of the correlations
   CorrTen sample_stds_Im{};
   // Integrated autocorrelation time tau_int (in pCN steps) per correlation point, recovered
   // from the variance inflation tau_int = 0.5 * sigma_blocking^2 / sigma_iid^2. Only filled
   // for errmethod=blocking; the blocking error already captures the autocorrelation, this just
   // exposes it as a number (burn-in / chain-length rule of thumb: a few x tau_int).
   CorrTen tau_int_Re{}, tau_int_Im{};
   // ...spin expectation value statistics (symmetry-permitted components only)
   MagVec spin_expval_sqsum{}; // sum_i s_i^2 (raw, MPI-reduced; used by the legacy AR(1) estimator)
   MagVec spin_expval_stds{};  // standard error of the mean for <S^alpha>

   // ...concerning batch-means (blocking) error estimation
   // The production samples of one core are split into num_blocks consecutive blocks of
   // block_length samples each. compute_spin_correlations sums the observable into the
   // current block (cur_block_*); close_block() turns a finished block sum into a block
   // mean and stores it. The error of the global mean is then the spread of all
   // (num_blocks * num_cores) pooled block means, which captures the full autocorrelation
   // (unlike the lag-1 AR(1) factor) as long as block_length >> autocorrelation time.
   std::string error_method{ "blocking" };  // "blocking" or "ar1"
   size_t num_blocks{};
   size_t block_length{};
   size_t blocks_filled{};
   CorrTen cur_block_Re{}, cur_block_Im{};       // running per-block sums of the observable
   MagVec cur_block_spin{};
   std::vector<CorrTen> block_means_Re{}, block_means_Im{};
   std::vector<MagVec> block_means_spin{};
   // Flyvbjerg-Petersen blocking-curve diagnostic: Re-correlation std (averaged over the
   // computed correlation points) vs block length. A plateau means block_length >>
   // autocorrelation time and the error bar is trustworthy; a curve still rising at the
   // largest block length means the run is too short / too autocorrelated for its length.
   std::vector<RealType> blocking_curve_len{};
   std::vector<RealType> blocking_curve_sigma{};

   void init_blocks( const CorrTen& template_corr, size_t num_SamplesPerCore );
   void close_block();

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
                             const MagVec& spin_mean, RealType N, RealType pcn_step_size );

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
   // The iteration error covers both self-consistency channels: the correlations (time-averaged
   // absolute deviation per direction pair) and the first moments (absolute deviation per
   // symmetry-permitted component); the reported error is the maximum over both channels.
   void compute_iteration_error( const CorrTen& CT, const CorrTen& new_CT, const MagVec& mag_old, const MagVec& mag_new );
   void finalize_iteration_step();
   bool terminate();

 private:
   const int my_rank{};

   // ...standard-error estimators (selected by error_method in compute_sample_stds)
   void compute_sample_stds_ar1( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im,
                                 const MagVec& spin_mean, RealType N, RealType pcn_step_size );
   void compute_sample_stds_blocking( const CorrTen& sample_mean_Re, const CorrTen& sample_mean_Im, RealType N );

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
