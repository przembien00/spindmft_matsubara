#include"../Functions/Functions.h"
#include<cmath>
#include<iostream>

namespace f=spinDMFT::Functions;
namespace c=spinDMFT::Contour;
namespace ps=spinDMFT::Parameter_Space;
namespace
{
void require(bool ok,const char* message)
{if(!ok)throw std::runtime_error(message);}
RealType norm(const f::ComplexDynamicMatrix& a)
{RealType s{};for(size_t i=0;i<a.rows();++i)for(size_t j=0;j<a.columns();++j)s+=std::norm(a(i,j));return std::sqrt(s);}
RealType spectral_norm(f::ComplexDynamicMatrix a)
{
    f::ComplexDynamicMatrix u,v;blaze::DynamicVector<RealType> singular;
    blaze::svd(a,u,singular,v);
    RealType result{};for(auto s:singular)result=std::max(result,s);return result;
}
f::ComplexDynamicMatrix factor(f::JointComplexGaussianSampler& sampler)
{
    f::ComplexDynamicMatrix L(sampler.size(),sampler.latent_dimension());
    for(size_t j=0;j<L.columns();++j)
    {
        f::JointComplexGaussianSampler::LatentVector e(L.columns(),RealType{});e[j]=1;
        const auto z=sampler.field_from_latent(e);
        for(size_t i=0;i<L.rows();++i)L(i,j)=z[i];
    }
    return L;
}
ps::ParameterSpace parameters(size_t m,size_t N,size_t q)
{
    ps::ParameterSpace p;
    p.spin_float=.5;p.num_HilbertSpaceDimension=2;
    p.spin_model=Physics::SpinModel{"DRF"};p.correlation_symmetry_type='D';
    p.B=Physics::MagneticField{"z",.3,0.,0.};
    p.noise=Physics::Noise{"z",.12};p.JQ=.7;p.JL=-.4;
    p.beta=.8;p.num_TimeSteps=m-1;p.num_TimePoints=m;p.delta_t=p.beta/(m-1);
    p.Tmax=.5;p.num_RealTimeSteps=N;p.num_RealTimePoints=N+1;p.delta_real_t=p.Tmax/N;
    p.real_time_substeps=q;p.gaussian_factorization="chebyshev";
    return p;
}
void check(size_t m,size_t N,size_t q,bool moments)
{
    auto p=parameters(m,N,q);
    c::CorrelationSet primitive{'D',m,N+1};
    // Intentionally noisy non-Hermitian component data test the ORIGINAL
    // canonical-triangle prescription, not additional physical symmetrizing.
    for(size_t t=0;t<=N;++t)for(size_t ab=0;ab<primitive.Re[t].size();++ab)
    {
        const auto pair=primitive.Re[t].get_direction_pair(ab);
        for(size_t tau=0;tau<m;++tau)
        {
            const RealType x=RealType(t)/N,y=RealType(tau)/(m-1);
            primitive.Re[t][ab][tau]=.12+.01*(pair[0]+2*pair[1])+.03*x*x+.015*y;
            primitive.Im[t][ab][tau]=.004*(RealType(pair[0])-RealType(pair[1]))+.02*x*(y-.5);
        }
    }
    f::MagTen re{'D',N+1},im{'D',N+1};
    const auto reference=f::self_consistent_equations(p,primitive,re,im,true);
    const auto field=f::self_consistent_equations(p,primitive,re,im,false);
    require(field.covariance.rows()==0,"compressed path allocated a dense covariance");
    require(field.structured_covariance!=nullptr,"compressed operator missing");
    const auto& op=*field.structured_covariance;
    std::mt19937 engine{891};std::normal_distribution<RealType> normal;
    for(size_t trial=0;trial<4;++trial)
    {
        f::JointComplexGaussianSampler::FieldVector x(op.size());
        for(auto& z:x)z={normal(engine),normal(engine)};
        const auto actual=op.apply(x);
        const f::JointComplexGaussianSampler::FieldVector expected=reference.covariance*x;
        RealType error{},scale{};
        for(size_t i=0;i<x.size();++i)
        {error=std::max(error,std::abs(actual[i]-expected[i]));scale=std::max(scale,std::abs(expected[i]));}
        require(error<2e-12*std::max(RealType{1},scale),"compressed action differs from canonical dense matrix");
    }
    RealType row_bound{};
    for(size_t i=0;i<op.size();++i)
    {RealType sum{};for(size_t j=0;j<op.size();++j)sum+=std::abs(reference.covariance(i,j));row_bound=std::max(row_bound,sum);}
    require(op.spectral_bound()+1e-13>=row_bound,"spectral bound underestimates a row sum");
    require(op.stored_complex_values()<op.size()*op.size(),"compressed operator unexpectedly uses dense storage");
    if(!moments)return;
    f::DenseComplexGaussianSampler dense(reference.covariance);
    const auto L=factor(dense);
    const f::ComplexDynamicMatrix targetH=L*ctrans(L);
    RealType previous_error=std::numeric_limits<RealType>::infinity();
    for(const RealType tolerance:{RealType{.04},RealType{.01},RealType{.001}})
    {
        f::ChebyshevComplexGaussianSampler sampler(field.structured_covariance,tolerance,16384);
        const auto K=factor(sampler);
        const f::ComplexDynamicMatrix deltaGamma=K*trans(K)-reference.covariance;
        const f::ComplexDynamicMatrix deltaH=K*ctrans(K)-targetH;
        RealType largest{};
        for(size_t i=0;i<K.rows();++i)for(size_t j=0;j<K.rows();++j)
            largest=std::max({largest,std::abs(deltaGamma(i,j)),std::abs(deltaH(i,j))});
        require(largest<sampler.covariance_approximation_error()+1e-11,"actual covariance error exceeds reported bound");
        require(std::max(spectral_norm(deltaGamma),spectral_norm(deltaH))
            <sampler.covariance_approximation_error()+1e-11,"covariance operator-norm error exceeds reported bound");
        const RealType error=std::max(norm(deltaGamma),norm(deltaH));
        require(error<previous_error,"tightening tolerance did not improve BOTH-moment comparison");previous_error=error;
        const c::StaggeredLayout layout{m,N*q};
        for(size_t a=0;a<3;++a)for(size_t j=0;j<K.columns();++j)
            require(K(layout.flat(c::StaggeredSector::Nu,N*q,a),j)==ComplexType{},"null response endpoint acquired noise");
        std::cout<<"degree="<<sampler.polynomial_degree()<<" two-moment Frobenius error="<<error
                 <<" maximum entry error="<<largest<<" bound="<<sampler.covariance_approximation_error()<<'\n';
        // Check the map stays linear (including batch use) for fixed setup.
        auto a=sampler.draw_latent(engine),b=sampler.draw_latent(engine);
        const f::JointComplexGaussianSampler::LatentVector sum=a+b;
        const auto combined=sampler.field_from_latent(sum),za=sampler.field_from_latent(a),zb=sampler.field_from_latent(b);
        for(size_t i=0;i<za.size();++i)require(std::abs(combined[i]-za[i]-zb[i])<1e-10,"sample map is not fixed and linear");
    }
    bool rejected=false;
    try{f::ChebyshevComplexGaussianSampler impossible(field.structured_covariance,RealType{1e-10},16);}
    catch(const std::runtime_error&){rejected=true;}
    require(rejected,"insufficient maximum degree was silently accepted");
}
}
int main()
{
    try
    {
        check(5,4,1,true);check(5,2,3,false);check(20,17,2,false);check(4,1,1,false);
        auto p=parameters(4,1,1);p.JQ=0;p.noise=Physics::Noise{"none",0.};
        c::CorrelationSet primitive{'D',4,2};f::MagTen re{'D',2},im{'D',2};
        const auto zero=f::self_consistent_equations(p,primitive,re,im,false);
        f::ChebyshevComplexGaussianSampler sampler(zero.structured_covariance);
        std::mt19937 rng{1};const auto z=sampler.draw(rng);
        for(auto v:z)require(v==ComplexType{},"zero covariance produced a nonzero sample");
        std::cout<<"Structured Gaussian action, moments, linearity, zero modes and failure checks passed\n";
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
