#pragma once

#include<complex>
#include<string>
#include<vector>

#include<Globals/Types.h>
#include<Observables/Correlations.h>
#include<Observables/Magnetization.h>
#include<Observables/Tensors.h>

#include"../Contour/Contour_Types.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace spinDMFT::Run_Time_Data
{
namespace ps = spinDMFT::Parameter_Space;
namespace ten = Observables::Tensors;
namespace corr = Observables::Correlations;
namespace mag = Observables::Magnetization;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using MagVec = mag::MagnetizationVector;
using IndexPair = ten::IndexPair;
using CorrelationSet = spinDMFT::Contour::CorrelationSet;

class RunTimeData
{
 public:
    RunTimeData() = default;
    RunTimeData( const ps::ParameterSpace& pspace, int my_rank );

    size_t generated_seed{};
    CorrelationSet contour_sample_stds{};
    CorrelationSet contour_tau_int{};
    std::vector<FieldVector> magnetization_time_Re{};
    std::vector<FieldVector> magnetization_time_Im{};
    std::vector<FieldVector> magnetization_time_Re_stds{};
    std::vector<FieldVector> magnetization_time_Im_stds{};
    std::vector<RealType> closed_contour_ratio_Re{};
    std::vector<RealType> closed_contour_ratio_Im{};
    std::vector<RealType> closed_contour_residual_Re_sample_stds{};
    std::vector<RealType> closed_contour_residual_Im_sample_stds{};
    std::vector<RealType> closed_contour_residual_abs_sample_stds{};

    std::vector<RealType> covariance_symmetry_errors{};
    std::vector<RealType> branch_identity_errors{};
    std::vector<RealType> gaussian_factor_reconstruction_errors{};
    std::vector<RealType> gaussian_covariance_approximation_errors{};
    std::vector<unsigned int> gaussian_factor_latent_dimensions{};
    std::vector<unsigned int> gaussian_largest_factorization_dimensions{};
    std::vector<RealType> average_phase_magnitudes{};
    std::vector<RealType> effective_sample_sizes{};
    std::vector<RealType> partition_sum_Re{};
    std::vector<RealType> partition_sum_Im{};
    std::vector<RealType> partition_abs_sums{};
    std::vector<RealType> denominator_constancy_residuals{};
    std::vector<RealType> mh_acceptance_rates{};
    std::vector<RealType> mh_nonpositive_rejection_rates{};
    std::vector<RealType> maximum_relative_imaginary_partitions{};
    std::vector<RealType> blocking_curve_block_lengths{};
    std::vector<RealType> blocking_curve_mean_errors{};
    std::vector<RealType> blocking_curve_max_errors{};
    std::vector<RealType> maximum_tau_int{};
    std::vector<RealType> block_length_to_tau_ratios{};
    std::vector<RealType> absolute_iteration_errors{};
    std::vector<RealType> standardized_iteration_errors{};
    size_t num_Iterations{};
    std::string termination{};

    std::string get_seed_str() const { return m_seed; }
    void begin_iteration_accumulation();
    void begin_sample( ComplexType Z, RealType observable_normalization=RealType{1.} );
    void accumulate_closed_contour_trace( size_t real_time, ComplexType value );
    void accumulate_magnetization( size_t real_time, size_t component,
                                   ComplexType numerator );
    void accumulate_edge_correlation( size_t real_time, size_t component,
                                      size_t imaginary_edge, ComplexType numerator );
    void end_sample();

    size_t num_magnetization_components() const { return m_mag_directions.size(); }
    size_t magnetization_direction( size_t c ) const { return m_mag_directions[c]; }
    size_t num_correlation_components() const { return m_corr_directions.size(); }
    IndexPair correlation_direction( size_t p ) const { return m_corr_directions[p]; }
    size_t num_real_time_points() const { return m_num_real_points; }
    size_t num_imaginary_edge_points() const { return m_num_imaginary_edge_points; }

    void record_complex_field_diagnostics( RealType covariance_symmetry_error,
                                           RealType branch_identity_error,
                                           RealType factor_reconstruction_error,
                                           RealType covariance_approximation_error,
                                           size_t latent_dimension,
                                           size_t largest_factorization_dimension );
    void record_pcn_diagnostics( size_t accepted, size_t proposed,
                                 size_t rejected_nonpositive,
                                 RealType maximum_relative_imaginary_partition );
    void mpi_reduce_and_finalize( CorrelationSet& correlations,
                                  CorrelationSet& standard_errors );
    void record_iteration_error( RealType absolute_error,
                                 RealType standardized_error );
    void finalize_iteration_step();
    bool terminate();
    size_t num_blocks() const { return m_num_blocks; }
    size_t samples_per_block() const { return m_samples_per_block; }

 private:
    struct ComplexBlockSums
    {
        ComplexBlockSums() = default;
        ComplexBlockSums( char symmetry, size_t num_imaginary_edge_points,
                          size_t num_real_points );
        ComplexType partition{};
        RealType partition_abs{};
        RealType partition_abs_sq{};
        CorrelationSet correlations{};
        std::vector<MagVec> mag_Re{};
        std::vector<MagVec> mag_Im{};
        std::vector<ComplexType> closure{};
    };

    int m_my_rank{};
    size_t m_num_print_digits{};
    char m_symmetry{};
    size_t m_num_imaginary_edge_points{};
    size_t m_num_real_points{};
    size_t m_num_samples_per_core{};
    size_t m_num_cores{};
    size_t m_num_blocks{};
    size_t m_samples_per_block{};
    size_t m_samples_seen{};
    size_t m_current_block{};
    std::string m_seed{};
    bool m_self_consistency{};
    bool m_harmonic_bath{};
    bool m_pcn{};
    bool m_antithetic_pairs{};
    RealType m_iteration_error_sigma_threshold{};
    size_t m_iteration_limit{};
    RealType m_covariance_tolerance{};
    RealType m_branch_identity_tolerance{};
    RealType m_takagi_tolerance{};
    RealType m_minimum_phase_magnitude{};
    RealType m_denominator_constancy_tolerance{};
    RealType m_imaginary_magnetization_sigma{};
    RealType m_partition_imaginary_tolerance{};
    std::vector<IndexPair> m_corr_directions{};
    std::vector<size_t> m_mag_directions{};
    ComplexBlockSums m_total{};
    std::vector<ComplexBlockSums> m_blocks{};
    CorrelationSet m_sample_squares{};
    CorrelationSet m_current_sample{};
    RealType m_current_observable_normalization{RealType{1.}};

    static ComplexType edge_value( const CorrelationSet& values,
                                   size_t t, size_t p, size_t tau );
    bool diagnostics_pass() const;
};
}
