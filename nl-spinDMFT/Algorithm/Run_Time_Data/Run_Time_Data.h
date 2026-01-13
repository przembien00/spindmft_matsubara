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

using CluCorrTen = clu::CorrelationCluster<ten::CorrelationTensor<corr::CorrelationVector>>;

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
    RealType positive_eigen_values_sum{};
    RealType negative_eigen_values_sum{}; 
    RealType largest_positive_eigen_value{};
    RealType largest_negative_eigen_value{};
    RealType negative_eigenvalue_ratio{}; // in the latest iteration step; sum of negative_eigenvalue_ratios
    bool eigenvalue_threshold_violated{false};  // is the soft threshold violated?

    // ...concerning the statistics
    size_t generated_seed{};
    CluCorrTen sample_sqsum{};    // sum of the squared single sample correlations
    CluCorrTen sample_stds{};     // std of the single sample correlations

    // PUBLIC METHODS
    // ...concerning the eigenvalues
    void process_and_check_eigenvalues( mvgb::EigenValuesBlocks& EVB ); // how to make the EigenValues const?

    // ...concerning the statistics
    std::string get_seed_str();
    void compute_sample_stds( const CluCorrTen& sample_sum );
    size_t get_num_SamplesPerCore() const;
    size_t get_num_Samples() const;

private:
    const int my_rank{};

    // PRIVATE MEMBERS IMPORTED FROM A PARAMETER SPACE (AND THEREFORE CONST)
    // ...concerning the eigenvalues
    const RealType eigenvalue_ratio_tolerance{};
    std::function<void(RealType&)> truncate_if_negative;

    // ...concerning the statistics
    const std::string seed{};
    const size_t num_Samples{};           // (initial) number of samples
    const size_t num_SamplesPerCore{};    // (initial) number of samples per core
    const size_t num_Cores{};
};

}
