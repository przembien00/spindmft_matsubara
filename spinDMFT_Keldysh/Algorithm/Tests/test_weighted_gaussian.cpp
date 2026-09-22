#include "../Functions/Complex_Gaussian.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace f=spinDMFT::Functions;
using Matrix=f::ComplexDynamicMatrix;
using Weights=std::array<RealType,3>;

namespace
{
int require(bool condition,const char* message)
{
    if(condition)return 0;
    std::cerr<<"FAILED: "<<message<<'\n';return 1;
}
RealType residual(const Matrix& a,const Matrix& b)
{
    RealType difference{},scale{};
    for(size_t i=0;i<a.rows();++i)for(size_t j=0;j<a.columns();++j)
    { difference+=std::norm(a(i,j)-b(i,j));scale+=std::norm(b(i,j)); }
    return std::sqrt(difference/(scale>RealType{}?scale:RealType{1.}));
}
Matrix factor(f::JointComplexGaussianSampler& sampler)
{
    Matrix F(sampler.size(),sampler.latent_dimension());
    f::JointComplexGaussianSampler::LatentVector unit(sampler.latent_dimension(),RealType{});
    for(size_t j=0;j<unit.size();++j)
    {
        unit[j]=RealType{1.};const auto value=sampler.field_from_latent(unit);unit[j]=RealType{};
        for(size_t i=0;i<value.size();++i)F(i,j)=value[i];
    }
    return F;
}
// Independent explicit change-of-basis matrices, with two distinct M edges,
// three real-time points and all spin components.
constexpr size_t M=6,R=9,N=M+2*R;
Matrix rotation(bool inverse)
{
    Matrix T(N,N,ComplexType{});
    for(size_t i=0;i<M;++i)T(i,i)=RealType{1.};
    for(size_t i=0;i<R;++i)
    {
        const RealType a=inverse?RealType{1.}:RealType{0.5};
        T(M+i,M+i)=T(M+i,M+R+i)=T(M+R+i,M+i)=a;
        T(M+R+i,M+R+i)=-a;
    }
    return T;
}
Matrix kernel()
{
    Matrix G(N,N,ComplexType{});
    for(size_t i=0;i<M;++i)for(size_t j=i;j<M;++j)
        G(i,j)=G(j,i)=i==j?ComplexType{RealType{0.8}+RealType{0.1}*i,RealType{0.02}*i}
                            :ComplexType{RealType{0.03}/(1+j-i),RealType{0.01}};
    for(size_t i=0;i<R;++i)for(size_t j=0;j<R;++j)
    {
        G(M+i,M+j)=RealType{0.4}*std::exp(-RealType{0.2}*std::abs(RealType(i)-RealType(j)));
        if(i/3>j/3)
            G(M+i,M+R+j)=G(M+R+j,M+i)=ComplexType{0.,RealType{0.08}/(1+i-j)};
    }
    for(size_t i=0;i<M;++i)for(size_t j=0;j<R;++j)
        G(i,M+j)=G(M+j,i)=ComplexType{RealType{0.02}*(1+i)/(1+j),RealType{0.015}*(1+j)/(1+i)};
    const Matrix inverse=rotation(true);
    return inverse*G*blaze::trans(inverse);
}
RealType cost(const Matrix& F,const Weights& w)
{
    const Matrix Z=rotation(false)*F;
    RealType value{};
    for(size_t i=0;i<N;++i)for(size_t j=0;j<Z.columns();++j)
        value+=w[i<M?0:i<M+R?1:2]*std::norm(Z(i,j));
    return value;
}
RealType lower_bound(const Matrix& covariance,const Weights& w)
{
    Matrix B=rotation(false);
    for(size_t i=0;i<N;++i)for(size_t j=0;j<N;++j)
        B(i,j)*=std::sqrt(w[i<M?0:i<M+R?1:2]);
    Matrix P=B*covariance*blaze::trans(B),U,V;
    blaze::DynamicVector<RealType> singular;
    blaze::svd(P,U,singular,V);
    RealType sum{};for(auto value:singular)sum+=value;
    return sum;
}
template<class Action> bool rejects(Action action)
{
    try{action();}catch(const std::exception&){return true;}return false;
}
}

int main()
{
    int failures{};
    const Matrix G=kernel();
    f::DenseComplexGaussianSampler dense(G);
    const Matrix F0=factor(dense),C0=F0*blaze::ctrans(F0);
    for(const Weights w:{Weights{1.,2.,2.},Weights{1.,2.,8.},Weights{1.,8.,2.},Weights{0.01,2.,2.}})
    {
        auto sampler=f::make_complex_gaussian_sampler("weighted-dense",G,1,3,0.1,-1.,w);
        const Matrix F=factor(*sampler),actual=F*blaze::trans(F),C=F*blaze::ctrans(F);
        const RealType error=residual(actual,G),bound=lower_bound(G,w);
        failures+=require(error<1e-11,"weighted factor preserves every physical covariance block");
        failures+=require(std::abs(error-sampler->reconstruction_error())<1e-12,
                          "reported reconstruction residual is in the physical basis");
        failures+=require(std::abs(cost(F,w)-bound)/bound<1e-11,"weighted cost attains independent SVD nuclear-norm bound");
        failures+=require(cost(F,w)<=cost(F0,w)+1e-10,"weighted cost improves or matches canonical cost");
        if(w[2]==2.&&w[1]==2.&&w[0]==1.)
            failures+=require(residual(C,C0)<1e-11,"(1,2,2) recovers complete canonical Gaussian distribution");
        else
            failures+=require(residual(C,C0)>1e-3,"noncanonical weights change the conjugated covariance");
        const Matrix T=rotation(false),rotated=T*actual*blaze::trans(T);
        for(size_t i=0;i<R;++i)for(size_t j=0;j<R;++j)
        {
            failures+=require(std::abs(rotated(M+R+i,M+R+j))<1e-11,"zero kappa pseudo-covariance survives");
            if(i/3<=j/3)failures+=require(std::abs(rotated(M+i,M+R+j))<1e-11,"causal response survives");
        }
        for(size_t i=0;i<M;++i)for(size_t j=0;j<R;++j)
            failures+=require(std::abs(rotated(i,M+R+j))<1e-11,"zero Matsubara-kappa block survives");

        Weights scaled=w;for(auto& value:scaled)value*=1e-80;
        f::WeightedDenseComplexGaussianSampler rescaled(G,1,3,scaled);
        const Matrix Fs=factor(rescaled);
        failures+=require(residual(Fs*blaze::ctrans(Fs),C)<1e-11,"common weight scale leaves the ensemble unchanged");
        failures+=require(rescaled.latent_dimension()==sampler->latent_dimension(),"common weight scale leaves numerical rank unchanged");
    }
    f::WeightedDenseComplexGaussianSampler single(G,1,3),batched(G,1,3);
    std::mt19937 e1{752},e2{752};
    for(size_t count:{size_t{0},size_t{1},size_t{7},size_t{32},size_t{3}})
    {
        const auto fields=batched.draw_contour_batch(e2,count,true);
        for(const auto& sample:fields)
        {
            const auto expected=single.draw(e1);
            failures+=require(!sample.has_real_gauss_fields(),"weighted CF4 uses edge interpolation");
            for(size_t i=0;i<expected.size();++i)
                failures+=require(std::abs(expected[i]-sample.edge_field[i])<1e-12,"batch and scalar draws agree");
        }
    }
    failures+=require(e1==e2,"batch and scalar draws consume identical RNG streams");
    auto r=single.draw_latent(e1);const auto plus=single.field_from_latent(r);
    for(auto& value:r)value=-value;
    const auto minus=single.field_from_latent(r);
    for(size_t i=0;i<N;++i)failures+=require(plus[i]==-minus[i],"latent sign reversal negates all branches jointly");
    failures+=require(rejects([&]{single.field_from_latent(f::JointComplexGaussianSampler::LatentVector(1));}),
                      "wrong latent dimension is rejected");
    f::WeightedDenseComplexGaussianSampler zero(Matrix(N,N,ComplexType{}),1,3);
    failures+=require(zero.latent_dimension()==0&&zero.reconstruction_error()==0.,"zero covariance has zero rank and residual");
    for(const auto& sample:zero.draw_contour_batch(e1,4,false))for(auto value:sample.edge_field)
        failures+=require(value==ComplexType{},"zero-rank batch returns zero fields");
    for(const RealType invalid:{RealType{0.},RealType{-1.},std::numeric_limits<RealType>::infinity(),
                               std::numeric_limits<RealType>::quiet_NaN()})
        for(size_t i=0;i<3;++i)
        {
            Weights w{1.,2.,8.};w[i]=invalid;
            failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(G,1,3,w);}),"invalid weights rejected in every sector");
        }
    failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(G,2,3);}),"incompatible contour dimensions rejected");
    failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(G,1,0);}),"empty real branch rejected");
    Matrix bad=G;bad(0,1)+=1.;
    failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(bad,1,3);}),"nonsymmetric covariance rejected");
    bad=G;bad(0,0)=std::numeric_limits<RealType>::quiet_NaN();
    failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(bad,1,3);}),"nonfinite covariance rejected");
    failures+=require(rejects([&]{f::WeightedDenseComplexGaussianSampler s(G,1,3,Weights{1e-100,1.,1.});}),
                      "weight contrast that loses physical covariance is rejected");
    return failures==0?0:1;
}
