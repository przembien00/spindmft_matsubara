#pragma once

#include<complex>
#include<cmath>
#include<array>
#include<stdexcept>
#include<string>
#include<vector>

#include<Globals/Types.h>
#include<Observables/Correlations.h>
#include<Observables/Magnetization.h>
#include<Observables/Tensors.h>

#include"../Contour/Contour_Types.h"
#include"../Parameter_Space/Parameter_Space.h"

namespace spinDMFT::Functions { struct MeasuredSample; }

namespace spinDMFT::Run_Time_Data
{
namespace ps = spinDMFT::Parameter_Space;
namespace ten = Observables::Tensors;
namespace corr = Observables::Correlations;
namespace mag = Observables::Magnetization;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using MagVec = mag::MagnetizationVector;
using MagTen = mag::MagnetizationTensor<MagVec>;
using IndexPair = ten::IndexPair;
using CorrelationSet = spinDMFT::Contour::CorrelationSet;
using Magnetization = spinDMFT::Contour::Magnetization;

// Results retained independently for each self-consistency loop.
struct LoopStatistics
{
    CorrelationSet contour_sample_stds{};
    CorrelationSet contour_tau_int{};
    Magnetization magnetization_time_stds{};
    std::vector<RealType> closed_contour_ratio_Re{};
    std::vector<RealType> closed_contour_ratio_Im{};
    std::vector<RealType> closed_contour_residual_Re_sample_stds{};
    std::vector<RealType> closed_contour_residual_Im_sample_stds{};
    std::vector<RealType> closed_contour_residual_abs_sample_stds{};

    std::vector<RealType> covariance_symmetry_errors{};
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
    std::vector<RealType> blocking_curve_block_lengths{};
    std::vector<RealType> blocking_curve_mean_errors{};
    std::vector<RealType> blocking_curve_max_errors{};
    std::vector<RealType> maximum_tau_int{};
    std::vector<RealType> block_length_to_tau_ratios{};
    std::vector<RealType> absolute_iteration_errors{};
    std::vector<RealType> standardized_iteration_errors{};
    size_t num_Iterations{};
    std::string termination{};
    size_t num_blocks{}, samples_per_block{};
};

class RunTimeData
{
 public:
    RunTimeData() = default;
    RunTimeData( const ps::ParameterSpace& pspace, int my_rank );

    size_t generated_seed{};
    LoopStatistics preparation{}, real_time{};
    // Real-time feedback residuals against the frozen Matsubara covariance.
    std::vector<RealType> frozen_errors{};

    std::string get_seed_str() const { return m_seed; }
    // Zero the caller-owned observable sums and accumulate directly into them.
    // They must remain alive until mpi_reduce_and_normalize() completes.
    void begin_iteration_accumulation(CorrelationSet& correlations,
        Magnetization& magnetization,bool matsubara_preparation=false);
    void begin_sample( ComplexType Z, RealType observable_normalization=RealType{1.} );
    void accumulate_closed_contour_trace( size_t real_time, ComplexType value );
    void accumulate_magnetization( size_t real_time, size_t component,
                                   ComplexType numerator );
    void accumulate_edge_correlation( size_t real_time, size_t component,
                                      size_t imaginary_edge, ComplexType numerator );
    void end_sample();
    void accumulate_sample( const spinDMFT::Functions::MeasuredSample& sample, RealType normalization=RealType{1.} );

    // Independent measurement writes directly into total and block sums.
    // Keep the contraction inlined into this loop, with no numerator buffer.
    void begin_independent_sample(ComplexType Z);
    template<typename Evaluate>
    void accumulate_independent_correlation_row(size_t t,size_t c,const Evaluate& evaluate)
    {
        if(t>=m_num_real_points||c>=m_corr_directions.size())
            throw std::invalid_argument("independent correlation row exceeds accumulator");
        RealType* total_re=m_correlations->Re[t][c].data();
        RealType* total_im=m_correlations->Im[t][c].data();
        RealType* block_re=m_blocks[m_current_block].correlations.Re[t][c].data();
        RealType* block_im=m_blocks[m_current_block].correlations.Im[t][c].data();
        for(size_t tau=0;tau<m_num_imaginary_edge_points;++tau)
        {
            const ComplexType value=evaluate(tau);
            const RealType re=std::real(value),im=std::imag(value);
            if(!std::isfinite(re)||!std::isfinite(im))
                throw std::runtime_error("non-finite edge contour-correlation numerator");
            total_re[tau]+=re; total_im[tau]+=im;
            block_re[tau]+=re; block_im[tau]+=im;
        }
    }

    size_t num_magnetization_components() const { return m_mag_directions.size(); }
    size_t magnetization_direction( size_t c ) const { return m_mag_directions[c]; }
    size_t num_correlation_components() const { return m_corr_directions.size(); }
    IndexPair correlation_direction( size_t p ) const { return m_corr_directions[p]; }
    size_t num_real_time_points() const { return m_num_real_points; }
    size_t num_imaginary_edge_points() const { return m_num_imaginary_edge_points; }

    void record_complex_field_diagnostics( RealType covariance_symmetry_error,
                                           RealType factor_reconstruction_error,
                                           RealType covariance_approximation_error,
                                           size_t latent_dimension,
                                           size_t largest_factorization_dimension );
    void record_pcn_diagnostics( size_t accepted, size_t proposed,
                                 size_t rejected_zero_weight );
    // Reduce the bound sums, calculate errors, then normalize those same variables.
    void mpi_reduce_and_normalize();
    void record_iteration_error( RealType absolute_error,
                                 RealType standardized_error );
    void record_matsubara_iteration_error( const CorrelationSet& old_correlations,
        const CorrelationSet& raw_correlations, const MagTen& old_magnetization_Re,
        const MagTen& raw_magnetization_Re, RealType standardized_error );
    RealType matsubara_iteration_error_threshold() const
    { return m_matsubara_iteration_error_threshold; }
    void finalize_iteration_step();
    bool terminate( const MagTen& magnetization_Im );

 private:
    struct ComplexBlockSums
    {
        ComplexBlockSums() = default;
        ComplexBlockSums( char symmetry, size_t num_imaginary_edge_points,
                          size_t num_real_points );
        ComplexType partition{};
        RealType partition_abs{};
        CorrelationSet correlations{};
        Magnetization magnetization{};
        std::vector<ComplexType> closure{};
    };

    LoopStatistics& loop() { return m_matsubara_preparation?preparation:real_time; }
    const LoopStatistics& loop() const { return m_matsubara_preparation?preparation:real_time; }
    size_t m_real_time_points{}, m_real_samples_per_core{}, m_real_iteration_limit{};
    size_t m_preparation_samples_per_core{}, m_preparation_iteration_limit{};
    size_t m_requested_num_blocks{};

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
    bool m_matsubara_preparation{};
    RealType m_matsubara_iteration_error_threshold{};
    RealType m_iteration_error_sigma_threshold{};
    RealType m_real_time_residual_scale{};
    size_t m_iteration_limit{};
    RealType m_covariance_tolerance{};
    RealType m_takagi_tolerance{};
    RealType m_minimum_phase_magnitude{};
    std::vector<IndexPair> m_corr_directions{};
    std::vector<size_t> m_mag_directions{};
    struct TotalSampleSums
    {
        ComplexType partition{};
        RealType partition_abs{},partition_abs_sq{};
        std::vector<ComplexType> closure{};
    };
    TotalSampleSums m_total{};
    // The caller owns the observable sums; runtime keeps only references.
    CorrelationSet* m_correlations{};
    Magnetization* m_magnetization{};
    std::vector<ComplexBlockSums> m_blocks{};
    CorrelationSet m_sample_squares{};
    CorrelationSet m_current_sample{};
    // Cross moments of the reweighted numerator n and denominator z:
    // n_Re*n_Im, n_Re*z_Re, n_Re*z_Im, n_Im*z_Re, n_Im*z_Im.
    // These give the iid variance of (n-g*z)/mean(z), used to infer tau_int
    // from ratio batch errors in modulus/phase-reweighted staged sampling.
    std::vector<std::array<RealType,5>> m_ratio_cross_moments;
    std::array<RealType,3> m_partition_second_moments{};
    ComplexType m_current_partition{};
    RealType m_current_observable_normalization{RealType{1.}};

    void configure_iteration(bool matsubara_preparation);
    void normalize_observables();
    void accumulate_ratio_cross_moments(size_t flat,RealType re,RealType im);

    static ComplexType edge_value( const CorrelationSet& values,
                                   size_t t, size_t p, size_t tau );
    bool diagnostics_pass( const MagTen& magnetization_Im ) const;
};
}
