#include "../Functions/Joint_Matsubara_PCN.h"
#include <iostream>
#include <cmath>

namespace f=spinDMFT::Functions;
using Matrix=f::ComplexDynamicMatrix;
using Vector=f::JointComplexGaussianSampler::LatentVector;
void require(bool ok,const char* message)
{if(!ok)throw std::runtime_error(message);}
RealType distance(const Matrix& a,const Matrix& b)
{
    RealType result{};
    for(size_t i=0;i<a.rows();++i)for(size_t j=0;j<a.columns();++j)
        result=std::max(result,std::abs(a(i,j)-b(i,j)));
    return result;
}
void check(const Matrix& gamma,size_t nm,const std::vector<RealType>& weights={})
{
    auto joint=weights.empty()?f::DenseComplexGaussianSampler(gamma)
                              :f::DenseComplexGaussianSampler(gamma,weights);
    f::JointMatsubaraSampler split(joint,nm);
    const size_t q=joint.latent_dimension(),r=split.latent_dimension();
    require(r<=2*nm&&r<=q,"M-active rank is bounded by real M dimension and joint rank");
    Matrix original(joint.size(),q),rotated(joint.size(),r+q);
    for(size_t j=0;j<q;++j)
    {
        Vector e(q,0.);e[j]=1.;const auto z=joint.field_from_latent(e);
        for(size_t i=0;i<joint.size();++i)original(i,j)=z[i];
    }
    require(distance(original,joint.leading_factor_rows(joint.size()))<1e-14,
        "exported factor includes physical inverse weights and disconnected block assembly");
    for(size_t j=0;j<r+q;++j)
    {
        Vector a(r,0.),xi(q,0.);if(j<r)a[j]=1.;else xi[j-r]=1.;
        const auto z=joint.field_from_latent(split.joint_latent(a,xi));
        const auto M=split.field_from_latent(a);
        for(size_t i=0;i<joint.size();++i)rotated(i,j)=z[i];
        for(size_t i=0;i<nm;++i)require(std::abs(z[i]-M[i])<1e-11,
            "fresh projected coordinates do not change M, including its imaginary part");
    }
    require(distance(rotated*blaze::trans(rotated),original*blaze::trans(original))<1e-10,
        "every physical pseudo-covariance block is preserved");
    require(distance(rotated*blaze::ctrans(rotated),original*blaze::ctrans(original))<1e-10,
        "every Hermitian sampling covariance block is preserved");
    require(split.reconstruction_error()<1e-11,"M row-space reconstruction");
    std::mt19937 rng(734);
    const auto a=split.draw_latent(rng);
    const auto M=split.field_from_latent(a);
    for(size_t k=0;k<8;++k)
    {
        const auto z=joint.field_from_latent(split.draw_joint_latent(a,rng));
        for(size_t i=0;i<nm;++i)require(std::abs(z[i]-M[i])<1e-11,"fixed M across fresh production draws");
    }
}
int main()
{
    try
    {
        Matrix gamma(7,7,ComplexType{});
        for(size_t i=0;i<7;++i)for(size_t j=i;j<7;++j)
            gamma(i,j)=gamma(j,i)=ComplexType{std::cos(1.+i+3*j),std::sin(2.+3*i+j)};
        check(gamma,2);check(gamma,2,{1.,1.,2.,2.,0.5,0.5,0.5});
        check(gamma,7); // No null space: the entire field is M.
        gamma=ComplexType{};check(gamma,2); // Zero bath and zero latent rank.
        gamma(5,5)=1.;check(gamma,2); // All real-time noise is freshly drawn.
        gamma(0,5)=gamma(5,0)=0.4;check(gamma,2); // A=0, B!=0 is supported.
        gamma(0,0)=1.;gamma(1,1)=1e-10;check(gamma,2); // Weak/singular A.
        Matrix L(7,2,ComplexType{});
        for(size_t i=0;i<7;++i){L(i,0)={0.1*(i+1),0.2};L(i,1)={0.3,0.07*i};}
        check(L*blaze::trans(L),2); // Low-rank full contour.
        std::cout<<"Joint Matsubara split preserves both Gaussian covariances and fresh residual support\n";
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
