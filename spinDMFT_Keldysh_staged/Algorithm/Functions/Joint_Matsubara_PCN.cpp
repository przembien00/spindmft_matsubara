#include "Joint_Matsubara_PCN.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace spinDMFT::Functions
{
JointMatsubaraSampler::JointMatsubaraSampler(
    const DenseComplexGaussianSampler& joint,size_t matsubara_size)
{
    const auto M=joint.leading_factor_rows(matsubara_size);
    const size_t q=M.columns();
    m_basis.resize(q,0);m_factor.resize(matsubara_size,0);
    if(q==0||matsubara_size==0)return;
    using Matrix=blaze::DynamicMatrix<RealType,blaze::rowMajor>;
    Matrix real(2*matsubara_size,q),left,right;
    RealType scale{};
    for(size_t i=0;i<matsubara_size;++i)for(size_t j=0;j<q;++j)
    {
        real(i,j)=std::real(M(i,j));real(matsubara_size+i,j)=std::imag(M(i,j));
        scale=std::max(scale,std::abs(M(i,j)));
    }
    if(scale==RealType{})return;
    // Normalize for scale-independent rank decisions, including tiny baths.
    real/=scale;
    blaze::DynamicVector<RealType> singular;
    blaze::svd(real,left,singular,right);
    const RealType eps=std::numeric_limits<RealType>::epsilon();
    const RealType cutoff=eps*static_cast<RealType>(std::max(real.rows(),q))*singular[0];
    size_t rank{};while(rank<singular.size()&&singular[rank]>cutoff)++rank;
    m_basis.resize(q,rank);
    // Blaze returns V^T: right singular vectors are its rows.
    for(size_t j=0;j<rank;++j)for(size_t i=0;i<q;++i)m_basis(i,j)=right(j,i);
    m_factor=M*m_basis;
    const ComplexDynamicMatrix restored=m_factor*blaze::trans(m_basis);
    RealType error{},norm{};
    for(size_t i=0;i<M.rows();++i)for(size_t j=0;j<q;++j)
    {
        error+=std::norm((M(i,j)-restored(i,j))/scale);
        norm+=std::norm(M(i,j)/scale);
    }
    m_error=std::sqrt(error/norm);
    if(!std::isfinite(m_error)||m_error>RealType{100.}*eps*static_cast<RealType>(std::max(real.rows(),q)))
        throw std::runtime_error("joint Matsubara row-space reconstruction failed");
}

JointMatsubaraSampler::LatentVector JointMatsubaraSampler::draw_latent(std::mt19937& engine)
{
    LatentVector a(latent_dimension());for(auto& x:a)x=m_normal(engine);return a;
}
JointMatsubaraSampler::FieldVector JointMatsubaraSampler::field_from_latent(const LatentVector& a)
{
    if(a.size()!=latent_dimension())throw std::invalid_argument("wrong joint Matsubara latent dimension");
    return m_factor*a;
}
JointMatsubaraSampler::FieldVector JointMatsubaraSampler::draw(std::mt19937& engine)
{return field_from_latent(draw_latent(engine));}

JointMatsubaraSampler::LatentVector JointMatsubaraSampler::joint_latent(
    const LatentVector& a,const LatentVector& fresh) const
{
    if(a.size()!=latent_dimension()||fresh.size()!=m_basis.rows())
        throw std::invalid_argument("wrong joint Matsubara or fresh latent dimension");
    // No full q-by-q null-space basis is needed. This is Ua+(I-UU^T)xi.
    const LatentVector correction=a-blaze::trans(m_basis)*fresh;
    return fresh+m_basis*correction;
}
JointMatsubaraSampler::LatentVector JointMatsubaraSampler::draw_joint_latent(
    const LatentVector& a,std::mt19937& engine)
{
    LatentVector fresh(m_basis.rows());for(auto& x:fresh)x=m_normal(engine);
    return joint_latent(a,fresh);
}
}
