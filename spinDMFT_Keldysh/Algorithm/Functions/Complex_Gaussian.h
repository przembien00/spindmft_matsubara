#pragma once

#include<array>
#include<complex>
#include<cstddef>
#include<memory>
#include<functional>
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

// Read-only physical covariance view used during covariance construction.
struct CovarianceSource
{
    size_t dimension{};
    std::function<ComplexType(size_t,size_t)> value;
    CovarianceSource() = default;
    CovarianceSource(size_t n, std::function<ComplexType(size_t,size_t)> entry)
        : dimension(n),value(std::move(entry)) {}
    CovarianceSource(const ComplexDynamicMatrix& matrix)
        : dimension(matrix.rows()),value([&matrix](size_t i,size_t j){return matrix(i,j);})
    {
        if(matrix.rows()!=matrix.columns()) throw std::invalid_argument("covariance must be square");
    }
    size_t rows() const { return dimension; }
    size_t columns() const { return dimension; }
    ComplexType operator()(size_t i,size_t j) const { return value(i,j); }
};
struct GaussianBlockFactors;

// Autonne--Takagi factorization of a complex symmetric pseudo-covariance.
// A nonsymmetric input is rejected instead of being projected. The factor is
// constructed from the equivalent real-symmetric eigenproblem.
TakagiFactor autonne_takagi( const ComplexDynamicMatrix& Gamma );

class JointComplexGaussianSampler
{
 public:
    using FieldVector = blaze::DynamicVector<ComplexType,blaze::columnVector>;
    using LatentVector = blaze::DynamicVector<RealType,blaze::columnVector>;

    struct ContourFieldSample
    {
        // With staggered Strang this vector instead stores the joint
        // [V_M,eta_midpoint,nu_boundary] draw in StaggeredLayout order.
        FieldVector edge_field{};
        // For CF4, node 0/1 is the earlier/later Gauss--Legendre node in
        // every real-time propagation subinterval. Each vector is laid out as
        // [interval][forward x,y,z; backward x,y,z].
        std::array<FieldVector,2> real_gauss_fields{};

        bool has_real_gauss_fields() const
        {
            return real_gauss_fields[0].size()!=0
                && real_gauss_fields[1].size()!=0;
        }
    };

    virtual ~JointComplexGaussianSampler() = default;
    // Every factorization exposes the same independent N(0,1) latent state.
    // pCN must act here: applying it directly to the complex physical field
    // would not, in general, preserve the Gaussian pseudo-covariance prior.
    virtual LatentVector draw_latent( std::mt19937& engine ) = 0;
    virtual FieldVector field_from_latent( const LatentVector& latent ) = 0;
    virtual FieldVector draw( std::mt19937& engine ) = 0;
    virtual ContourFieldSample contour_field_from_latent(
        const LatentVector& latent, bool include_real_gauss_fields );
    ContourFieldSample draw_contour_field(
        std::mt19937& engine, bool include_real_gauss_fields );
    virtual std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine, size_t count, bool include_real_gauss_fields );
    virtual RealType reconstruction_error() const = 0;
    virtual RealType covariance_approximation_error() const { return RealType{}; }
    virtual size_t latent_dimension() const = 0;
    virtual size_t largest_factorization_dimension() const { return size(); }
    virtual size_t polynomial_degree() const { return 0; }
    virtual RealType spectral_bound() const { return RealType{}; }
    virtual size_t size() const = 0;
};

// Dense joint-contour sampler. Independent real standard normals r are mapped
// through L to z=Lr, giving E[z z^T]=Gamma.
class DenseComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    using LatentVector = JointComplexGaussianSampler::LatentVector;
    using FieldVector = JointComplexGaussianSampler::FieldVector;

    explicit DenseComplexGaussianSampler( const ComplexDynamicMatrix& covariance );
    // Diagonal noise objective in an already transformed field basis.
    DenseComplexGaussianSampler(const ComplexDynamicMatrix& covariance,
        const std::vector<RealType>& weights);

    LatentVector draw_latent( std::mt19937& engine ) override;
    FieldVector field_from_latent( const LatentVector& latent ) override;
    FieldVector draw( std::mt19937& engine ) override;

    RealType reconstruction_error() const override;
    size_t numerical_rank() const { return latent_dimension(); }
    size_t latent_dimension() const override;
    size_t size() const override;
    size_t largest_factorization_dimension() const override;
    std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine,size_t count,bool include_real_gauss_fields ) override;

 private:
    std::shared_ptr<GaussianBlockFactors> m_factors;
    std::vector<RealType> m_weight_roots{};
    RealType m_unweighted_reconstruction_error{};
    std::normal_distribution<RealType> m_standard_normal{ RealType{0.}, RealType{1.} };
};

// Minimize tr(W E[z z^dagger]) at fixed E[v v^T], where
// z=(V_M,eta=(V_++V_-)/2,kappa=(V_+-V_-)/2) and W has positive
// sector weights (M,eta,kappa). Samples are returned in the physical v basis.
// (1,2,2) recovers the canonical dense ensemble; (1,2,8) uses the
// Schmitz--Stockburger relative real-time weights, since nu=2*kappa.
class WeightedDenseComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    WeightedDenseComplexGaussianSampler(
        const ComplexDynamicMatrix& covariance, size_t num_matsubara_intervals,
        size_t num_real_points,
        const std::array<RealType,3>& weights={RealType{1.},RealType{2.},RealType{8.}} );
    LatentVector draw_latent( std::mt19937& engine ) override;
    FieldVector field_from_latent( const LatentVector& latent ) override;
    FieldVector draw( std::mt19937& engine ) override;
    std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine,size_t count,bool include_real_gauss_fields ) override;
    RealType reconstruction_error() const override { return m_reconstruction_error; }
    size_t latent_dimension() const override;
    size_t size() const override;
    size_t largest_factorization_dimension() const override;

 private:
    void to_physical_field( FieldVector& field ) const;
    std::shared_ptr<GaussianBlockFactors> m_factors;
    size_t m_matsubara_size{},m_real_size{};
    std::array<RealType,3> m_roots{};
    RealType m_reconstruction_error{};
    std::normal_distribution<RealType> m_standard_normal{RealType{0.},RealType{1.}};
};

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const std::string& algorithm,
    const ComplexDynamicMatrix& covariance,
    size_t num_matsubara_intervals,
    size_t num_real_points,
    const std::array<RealType,3>& weights={RealType{1.},RealType{2.},RealType{8.}} );

}
