#pragma once

#include<array>
#include<complex>
#include<cstddef>
#include<memory>
#include<vector>
#include<random>
#include<string>
#include<stdexcept>

#include<blaze/Math.h>
#include<Globals/Types.h>

namespace spinDMFT::Functions
{

using ComplexFieldVector = blaze::StaticVector<ComplexType,3UL,blaze::columnVector>;
using ComplexDynamicMatrix = blaze::DynamicMatrix<ComplexType,blaze::rowMajor>;

struct TakagiFactor
{
    // Gamma = L L^T. Columns belonging to numerically zero singular values
    // are omitted, so the factor can be rectangular.
    ComplexDynamicMatrix L{};
    RealType reconstruction_error{};
    size_t numerical_rank{};
    std::vector<RealType> singular_values{};
};

struct GaussianBlockFactors;

// Autonne--Takagi factorization of a complex symmetric pseudo-covariance.
// A nonsymmetric input is rejected instead of being projected. The factor is
// constructed from the equivalent real-symmetric eigenproblem. Nonfinite or
// excessive reconstruction residuals are rejected before returning the factor.
TakagiFactor autonne_takagi(const ComplexDynamicMatrix& Gamma,
    RealType takagi_tolerance=RealType{1e-8});

class JointComplexGaussianSampler
{
 public:
    using FieldVector = blaze::DynamicVector<ComplexType,blaze::columnVector>;
    using LatentVector = blaze::DynamicVector<RealType,blaze::columnVector>;

    struct ContourFieldSample
    {
        // With staggered Strang this stores [V_M, eta_midpoint, nu_boundary]
        // in StaggeredLayout order, from one joint draw.
        FieldVector edge_field{};
    };

    virtual ~JointComplexGaussianSampler() = default;
    // Every factorization exposes the same independent N(0,1) latent state.
    // pCN must act here: applying it directly to the complex physical field
    // would not, in general, preserve the Gaussian pseudo-covariance prior.
    virtual LatentVector draw_latent( std::mt19937& engine ) = 0;
    virtual FieldVector field_from_latent( const LatentVector& latent ) = 0;
    virtual FieldVector draw( std::mt19937& engine ) = 0;
    virtual ContourFieldSample contour_field_from_latent(
        const LatentVector& latent );
    ContourFieldSample draw_contour_field(
        std::mt19937& engine );
    virtual std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine, size_t count );
    virtual RealType reconstruction_error() const = 0;
    virtual size_t latent_dimension() const = 0;
    virtual size_t largest_factorization_dimension() const { return size(); }
    virtual size_t size() const = 0;
};

// Dense joint-contour sampler. Independent real standard normals r are mapped
// through L to z=Lr, giving E[z z^T]=Gamma.
class DenseComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    using LatentVector = JointComplexGaussianSampler::LatentVector;
    using FieldVector = JointComplexGaussianSampler::FieldVector;

    // Reject nonfinite or excessive physical reconstruction residuals at construction.
    explicit DenseComplexGaussianSampler(const ComplexDynamicMatrix& covariance,
        RealType takagi_tolerance=RealType{1e-8});
    // Preserve a diagonal sampling objective in another field basis. The
    // returned draws reconstruct the original, unweighted pseudo-covariance.
    DenseComplexGaussianSampler(const ComplexDynamicMatrix& covariance,
        const std::vector<RealType>& weights, RealType takagi_tolerance=RealType{1e-8});

    LatentVector draw_latent( std::mt19937& engine ) override;
    FieldVector field_from_latent( const LatentVector& latent ) override;
    FieldVector draw( std::mt19937& engine ) override;

    RealType reconstruction_error() const override;
    size_t numerical_rank() const { return latent_dimension(); }
    size_t latent_dimension() const override;
    size_t size() const override;
    size_t largest_factorization_dimension() const override;
    // Selected leading physical rows, with the actual block assembly and
    // inverse weights used by field_from_latent (no new factorization).
    ComplexDynamicMatrix leading_factor_rows(size_t rows) const;
    std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine,size_t count ) override;

 private:
    void validate_reconstruction(RealType tolerance) const;
    std::shared_ptr<GaussianBlockFactors> m_factors;
    std::vector<RealType> m_weight_roots{};
    RealType m_unweighted_reconstruction_error{};
    std::normal_distribution<RealType> m_standard_normal{ RealType{0.}, RealType{1.} };
};

}
