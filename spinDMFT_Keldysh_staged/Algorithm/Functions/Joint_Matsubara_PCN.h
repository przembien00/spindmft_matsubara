#pragma once
#include "Complex_Gaussian.h"

namespace spinDMFT::Functions
{
// Imaginary-only view of the CURRENT joint Takagi law. U spans the real
// row space of [Re L_M; Im L_M]. M=L_M U a; fresh real normals are projected
// onto its orthogonal complement for real-time production draws.
class JointMatsubaraSampler : public JointComplexGaussianSampler
{
 public:
    JointMatsubaraSampler(const DenseComplexGaussianSampler& joint,size_t matsubara_size);
    LatentVector draw_latent(std::mt19937& engine) override;
    FieldVector field_from_latent(const LatentVector& a) override;
    FieldVector draw(std::mt19937& engine) override;
    size_t size() const override { return m_factor.rows(); }
    size_t latent_dimension() const override { return m_basis.columns(); }
    RealType reconstruction_error() const override { return m_error; }
    size_t residual_dimension() const { return m_basis.rows()-m_basis.columns(); }
    // Deterministic map exposed for covariance-level verification.
    LatentVector joint_latent(const LatentVector& a,const LatentVector& fresh) const;
    LatentVector draw_joint_latent(const LatentVector& a,std::mt19937& engine);
 private:
    blaze::DynamicMatrix<RealType,blaze::rowMajor> m_basis;
    ComplexDynamicMatrix m_factor;
    RealType m_error{};
    std::normal_distribution<RealType> m_normal{RealType{},RealType{1.}};
};
}
