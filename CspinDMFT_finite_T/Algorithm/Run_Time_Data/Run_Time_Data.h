#pragma once

#include<vector>
#include<iostream>
#include<Globals/Types.h>
#include"../Parameter_Space/Parameter_Space.h"
#include<Multivariate_Gaussian/Multivariate_Gaussian_Blocks.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Clusters.h>

namespace Run_Time_Data
{
    
namespace ps = DMFT_parameter_space;
namespace mvgb = Multivariate_Gaussian::Blocks;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace clu = Observables::Clusters;

using Vec = std::vector<RealType>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using CluCorrTen = clu::CorrelationCluster<CorrTen>;

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
    Vec positive_eigen_values_sum_list{};
    Vec negative_eigen_values_sum_list{}; 
    Vec largest_positive_eigen_value_list{};
    Vec largest_negative_eigen_value_list{};
    Vec negative_eigenvalue_ratio_list{}; // in the latest iteration step; sum of negative_eigenvalue_ratios
    std::vector<bool> eigenvalue_threshold_violated_list{}; // is the threshold violated?
  
    // ...concerning the statistics
    size_t generated_seed{};
    CluCorrTen sample_sqsum_Re{}; // sum of squared real-part correlation samples
    CluCorrTen sample_sqsum_Im{}; // sum of squared imaginary-part correlation samples
    CluCorrTen sample_stds_Re{};  // std of the single sample real-part correlations
    CluCorrTen sample_stds_Im{};  // std of the single sample imaginary-part correlations
    std::vector<FieldVector> spin_expval_sqsum{}; // sum of squared single sample spin expectation values
    std::vector<FieldVector> spin_expval_stds{};   // std of the single sample spin expectation values
    std::vector<size_t> adaptive_num_SamplesPerCore{}; // number of samples per core in case of adaptive sample size
    bool sample_size_updated = {false}; // was the sample size updated in the last iteration step?

    // ...concerning the iterations
    size_t num_Iterations{};
    std::string termination{};
    Vec relative_iteration_error_list{};
    Vec absolute_iteration_error_list{};

    // ...concerning the preconditioned Crank-Nicolson chain diagnostics
    // mh_accepted_count and mh_proposed_count are accumulated per MPI core during one
    // iteration; record_mh_acceptance() reduces them across MPI ranks and appends the
    // resulting global acceptance ratio to acceptance_rates.
    size_t mh_accepted_count{};
    size_t mh_proposed_count{};
    std::vector<RealType> acceptance_rates{};
    void record_mh_acceptance();

    // PUBLIC METHODS
    // ...concerning the eigenvalues
    void process_and_check_eigenvalues( mvgb::EigenValuesBlocks& EVB ); // how to make the EigenValues const?

    // ...concerning the statistics
    std::string get_seed_str();
    size_t get_num_SamplesPerCore() const;
    size_t get_num_Samples() const;
    size_t get_num_SetsPerCore() const;
    // pCN estimators are plain chain means, so the i.i.d. standard error of the mean is
    //   stderr(<o>) = sqrt( ( <o^2> - <o>^2 ) / N ).
    // The pCN samples are autocorrelated, so this underestimates the true error. We multiply
    // it by sqrt( (1+rho1)/(1-rho1) ), the AR(1) variance-inflation factor, where rho1 is the
    // lag-1 autocorrelation (see pcn_autocorrelation_factor). sample_mean_* hold the already-
    // normalized means (sum_i o_i / N); sample_sqsum_* hold the MPI-reduced raw sum_i o_i^2;
    // N is the global sample count; pcn_step_size is the pCN beta used to estimate rho1.
    void compute_and_process_sample_stds( const CluCorrTen& sample_mean_Re, const CluCorrTen& sample_mean_Im, const RealType N, const RealType pcn_step_size );
    void compute_and_process_spin_expval_stds( const std::vector<FieldVector>& spin_mean, const RealType N, const RealType pcn_step_size );

    // AR(1) autocorrelation correction factor sqrt( (1+rho1)/(1-rho1) ) for the pCN stds.
    // rho1 is estimated from the latest acceptance rate a and the pCN step size beta via the
    // pCN move: a rejected step copies the state, an accepted one retains a fraction
    // sqrt(1-beta^2) of it, giving rho1 ~ 1 - a*(1 - sqrt(1-beta^2)). This is a crude single-
    // exponential estimate: it assumes one relaxation timescale and UNDERESTIMATES the true
    // autocorrelation when the chain has a slow / trapped mode.
    RealType pcn_autocorrelation_factor( const RealType pcn_step_size ) const;

    // ...concerning the iterations
    std::string regarded(const std::string& mode){return (mode==iteration_error_mode) ? " (regarded)" : "";}
    void compute_iteration_error( const CluCorrTen& new_CCT, const CluCorrTen& CCT );
    void finalize_iteration_step();
    bool is_converged() const;
    bool signal_file_exists() const;
    bool terminate();

private:
    const int my_rank{};

    // PRIVATE MEMBERS IMPORTED FROM A PARAMETER SPACE (AND THEREFORE CONST)
    const int Run_ID{};
    const size_t num_PrintDigits{};

    // ...concerning the eigenvalues
    const RealType eigenvalue_ratio_tolerance{};
    std::function<void(RealType&)> truncate_if_negative;

    // ...concerning the statistics
    const std::string seed{};
    const size_t num_Samples{};           // (initial) number of samples
    const size_t num_SamplesPerCore{};    // (initial) number of samples per core
    const size_t num_SamplesPerSet{};     // number of samples per set
    const size_t num_Cores{};
    const bool adaptive_sample_size{};
    const RealType statistical_error_tolerance{};
    
    // ...concerning the iterations
    const size_t Iteration_Limit{};
    const std::string iteration_error_mode{};
    const RealType absolute_iteration_error_tolerance{};
    const RealType relative_iteration_error_tolerance{};
};

}
