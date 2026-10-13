#include"../Functions/Complex_Gaussian.h"
#include<cmath>
#include<iostream>
#include<stdexcept>

namespace f=spinDMFT::Functions;
using Matrix=f::ComplexDynamicMatrix;
void require(bool condition,const char* message)
{if(!condition)throw std::runtime_error(message);}
RealType distance(const Matrix& a,const Matrix& b)
{
    RealType sum{};
    for(size_t i=0;i<a.rows();++i)for(size_t j=0;j<a.columns();++j)
        sum+=std::norm(a(i,j)-b(i,j));
    return std::sqrt(sum);
}

// Recover the actual sampling map, including disconnected-block assembly,
// from its action on unit latent vectors, rather than testing a separate factor.
Matrix sampling_factor(f::DenseComplexGaussianSampler& sampler)
{
    Matrix L(sampler.size(),sampler.latent_dimension(),ComplexType{});
    for(size_t j=0;j<L.columns();++j)
    {
        f::JointComplexGaussianSampler::LatentVector e(L.columns(),RealType{});
        e[j]=RealType{1.};const auto v=sampler.field_from_latent(e);
        for(size_t i=0;i<L.rows();++i)L(i,j)=v[i];
    }
    return L;
}
Matrix check(const Matrix& covariance)
{
    f::DenseComplexGaussianSampler sampler(covariance);
    const auto L=sampling_factor(sampler);
    require(sampler.reconstruction_error()<1e-10,"reported joint reconstruction");
    require(distance(L*blaze::trans(L),covariance)<1e-10,"actual draw reconstructs every joint block");
    require(L.columns()<=L.rows(),"no conditional completion columns");
    return L;
}
Matrix joint(const Matrix& A,RealType scale)
{
    const size_t nm=A.rows(),nr=4;Matrix g(nm+nr,nm+nr,ComplexType{});
    for(size_t i=0;i<nm;++i)for(size_t j=0;j<nm;++j)g(i,j)=A(i,j);
    // Two real branches, with mixed couplings also to weak/null M directions.
    for(size_t i=0;i<nm;++i)for(size_t j=0;j<nr;++j)
        g(i,nm+j)=g(nm+j,i)=scale*ComplexType{0.1*(i+1),0.07*(j+1)};
    for(size_t i=0;i<nr;++i)for(size_t j=0;j<nr;++j)
        g(nm+i,nm+j)={i==j?0.7:0.1,0.03*(i+j)};
    return g;
}
int main()
{
    try
    {
        Matrix A(2,2,ComplexType{});
        A(0,0)={1.2,0.3};A(1,1)={-0.8,0.2};A(0,1)=A(1,0)={0.1,-0.4};
        check(joint(A,RealType{0.4}));check(joint(A,RealType{1.2}));
        A=ComplexType{};A(0,0)=RealType{1.};A(1,1)=RealType{1e-10};check(joint(A,RealType{1.}));
        A(1,1)=ComplexType{};check(joint(A,RealType{1.}));
        // Even A=0, B!=0 is representable as a complex Gaussian: its M
        // pseudo-covariance is zero, but M is not deterministically zero.
        A=ComplexType{};check(joint(A,RealType{1.}));
        Matrix column(3,1);column(0,0)={1.,0.2};column(1,0)={0.3,-0.7};column(2,0)={0.2,0.5};
        check(joint(column*blaze::trans(column),RealType{1.}));

        // Gamma^2=(1+b^2) I: the exact canonical Hermitian covariance
        // equals sqrt(1+b^2) I despite the unchanged Gamma_MM=1.
        for(const RealType b:{RealType{},RealType{0.4},RealType{1.2}})
        {
            Matrix g(2,2,ComplexType{});g(0,0)=RealType{1.};g(1,1)=RealType{-1.};g(0,1)=g(1,0)=b;
            const auto L=check(g);
            const Matrix H=L*blaze::ctrans(L);
            Matrix expected(2,2,ComplexType{});expected(0,0)=expected(1,1)=std::sqrt(RealType{1.}+b*b);
            require(distance(H,expected)<1e-11,"fixed pseudo-covariance permits a changing Hermitian marginal");
            const RealType conditional_trace=RealType{2.}+RealType{2.}*b*b;
            require(std::real(H(0,0)+H(1,1))<=conditional_trace+1e-11,"canonical field norm bound");
        }
        Matrix zero(7,7,ComplexType{});f::DenseComplexGaussianSampler empty(zero);
        require(empty.latent_dimension()==0,"fully zero bath has no artificial fluctuations");
        std::mt19937 rng(17);for(const auto value:empty.draw(rng))require(value==ComplexType{},"zero bath draw");
        std::cout<<"Joint Takagi blocks, singular preparation, and changing auxiliary marginal passed\n";
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
