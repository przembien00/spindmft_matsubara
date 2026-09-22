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

// Reflection tau -> beta-tau on an edge grid with distinct 0+ and beta-
// endpoints. num_intervals is the index of the beta endpoint.
size_t one_sided_edge_reflection_index( size_t point, size_t num_intervals );

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

// Read-only physical covariance view. The callable is used only during
// factor construction, permitting FFT setup without a dense physical matrix.
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

// Canonical Takagi factorization obtained from a full complex SVD.  Degenerate
// singular-value subspaces receive a small Takagi phase correction.  The
// returned factor has both the same pseudo-covariance and the same Hermitian
// covariance as autonne_takagi():
//
//   L L^T      = Gamma,
//   L L^dagger = sqrt(Gamma Gamma^dagger).
TakagiFactor svd_takagi( const ComplexDynamicMatrix& Gamma );

class JointComplexGaussianSampler
{
 public:
    using FieldVector = blaze::DynamicVector<ComplexType,blaze::columnVector>;
    using LatentVector = blaze::DynamicVector<RealType,blaze::columnVector>;

    struct ContourFieldSample
    {
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
    std::normal_distribution<RealType> m_standard_normal{ RealType{0.}, RealType{1.} };
};

// Alternative canonical sampler using the complex-SVD Takagi factor.  It
// samples exactly the same complete Gaussian ensemble as the dense real-lift
// algorithm; only factor construction differs.
class SVDComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    using LatentVector = JointComplexGaussianSampler::LatentVector;
    using FieldVector = JointComplexGaussianSampler::FieldVector;

    explicit SVDComplexGaussianSampler( const ComplexDynamicMatrix& covariance );

    LatentVector draw_latent( std::mt19937& engine ) override;
    FieldVector field_from_latent( const LatentVector& latent ) override;
    FieldVector draw( std::mt19937& engine ) override;
    RealType reconstruction_error() const override;
    size_t latent_dimension() const override;
    size_t size() const override;
    size_t largest_factorization_dimension() const override;
    std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine,size_t count,bool include_real_gauss_fields ) override;

 private:
    std::shared_ptr<GaussianBlockFactors> m_factors;
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

// Frequency-space sampler. The first N_tau Matsubara values (0+ through
// beta-dtau) and the doubled-real axis are transformed with unitary FFTs. The
// distinct beta- endpoint remains untransformed in the dense Matsubara block.
// After each draw the inverse FFT is applied and only the physical contour is
// returned.
class FFTDenseComplexGaussianSampler : public JointComplexGaussianSampler
{
 public:
    using LatentVector = JointComplexGaussianSampler::LatentVector;
    using FieldVector = JointComplexGaussianSampler::FieldVector;

    FFTDenseComplexGaussianSampler( const CovarianceSource& covariance,
                                    size_t num_matsubara_intervals,
                                    size_t num_real_points,
                                    RealType delta_real_time,
                                    RealType cross_frequency_cutoff, size_t real_time_substeps=1 );
    ~FFTDenseComplexGaussianSampler() override;

    LatentVector draw_latent( std::mt19937& engine ) override;
    FieldVector field_from_latent( const LatentVector& latent ) override;
    FieldVector draw( std::mt19937& engine ) override;
    ContourFieldSample contour_field_from_latent(
        const LatentVector& latent, bool include_real_gauss_fields ) override;
    std::vector<ContourFieldSample> draw_contour_batch(
        std::mt19937& engine,size_t count,bool include_real_gauss_fields ) override;
    RealType reconstruction_error() const override { return m_reconstruction_error; }
    RealType covariance_approximation_error() const override
    { return m_covariance_approximation_error; }
    size_t latent_dimension() const override { return m_latent_dimension; }
    size_t largest_factorization_dimension() const override
    { return m_largest_factorization_dimension; }
    size_t size() const override { return m_physical_size; }

 private:
    struct FFTPlans;
    struct FrequencyFactors;
    ContourFieldSample contour_field_from_frequency(
        size_t sample,bool include_real_gauss_fields );

    size_t m_num_matsubara_intervals{};
    size_t m_num_matsubara_points{};
    size_t m_num_real_points{};
    size_t m_embedded_real_points{};
    size_t m_real_time_substeps{1};
    size_t m_physical_size{};
    size_t m_latent_dimension{};
    RealType m_reconstruction_error{};
    RealType m_covariance_approximation_error{};
    size_t m_largest_factorization_dimension{};
    std::unique_ptr<FrequencyFactors> m_frequency_factors{};
    std::unique_ptr<FFTPlans> m_fft{};
    std::normal_distribution<RealType> m_standard_normal{ RealType{0.}, RealType{1.} };
};

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const CovarianceSource& covariance, size_t num_matsubara_intervals,
    size_t num_real_points, RealType delta_real_time, RealType cross_frequency_cutoff, size_t real_time_substeps=1 );

std::unique_ptr<JointComplexGaussianSampler> make_complex_gaussian_sampler(
    const std::string& algorithm,
    const ComplexDynamicMatrix& covariance,
    size_t num_matsubara_intervals,
    size_t num_real_points,
    RealType delta_real_time=RealType{1.},
    RealType cross_frequency_cutoff=RealType{-1.},
    const std::array<RealType,3>& weights={RealType{1.},RealType{2.},RealType{8.}},
    size_t real_time_substeps=1 );

}
