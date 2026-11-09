#pragma once

#include"Complex_Gaussian.h"

namespace spinDMFT::Functions
{
// Exact stationary lag tables and unrefined mixed primitive. Component
// blocks are in row-major (a,b) order. All transposes are unconjugated.
class StructuredContourCovariance
{
 public:
    using Block=std::array<ComplexType,9>;
    using FieldVector=JointComplexGaussianSampler::FieldVector;
    StructuredContourCovariance(size_t imaginary_points,size_t measured_intervals,
        size_t substeps,std::vector<Block> matsubara_lags,
        std::vector<Block> eta_lags,std::vector<Block> response_lags,
        std::vector<Block> mixed_primitive);
    ~StructuredContourCovariance();
    FieldVector apply(const FieldVector& x) const;
    size_t size() const;
    // Certified absolute row-sum upper bound on ||Gamma||_2, including
    // absolute interpolation weights. Also identifies exact zero rows.
    RealType spectral_bound() const;
    const std::vector<bool>& active_rows() const;
    size_t stored_complex_values() const;
 private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

// Fixed polynomial approximation to sqrt(max(A,0)). It never constructs
// Gamma, the real lift, eigenvectors, or a dense sampling factor.
class ChebyshevComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    ChebyshevComplexGaussianSampler(std::shared_ptr<StructuredContourCovariance> covariance,
        RealType relative_moment_tolerance=RealType{1e-3},size_t max_degree=16384);
    LatentVector draw_latent(std::mt19937& engine) override;
    FieldVector field_from_latent(const LatentVector& latent) override;
    FieldVector draw(std::mt19937& engine) override;
    // No factorization is performed. Polynomial covariance error is reported
    // separately as an absolute operator-norm bound for EACH complex moment.
    RealType reconstruction_error() const override { return RealType{}; }
    RealType covariance_approximation_error() const override { return m_moment_bound; }
    size_t latent_dimension() const override { return 2*size(); }
    size_t largest_factorization_dimension() const override { return 0; }
    size_t size() const override { return m_covariance->size(); }
    size_t polynomial_degree() const override { return m_coefficients.empty()?0:m_coefficients.size()-1; }
    RealType spectral_bound() const override { return m_covariance->spectral_bound(); }
 private:
    LatentVector apply_scaled_lift(const LatentVector& x) const;
    std::shared_ptr<StructuredContourCovariance> m_covariance;
    std::vector<RealType> m_coefficients;
    RealType m_moment_bound{};
    std::normal_distribution<RealType> m_normal{RealType{},RealType{1}};
};
}
