#include"../Functions/Functions.h"
#include<cmath>
#include<iostream>

namespace func=spinDMFT::Functions;
namespace ps=spinDMFT::Parameter_Space;
namespace contour=spinDMFT::Contour;
namespace {
int failures{};
void require(bool ok,const char* message) {
    if(!ok) { ++failures; std::cerr<<"FAILED: "<<message<<'\n'; }
}
RealType frobenius(const Operator& x) {
    RealType sum{};
    for(size_t i=0;i<x.rows();++i)for(size_t j=0;j<x.columns();++j)sum+=std::norm(x(i,j));
    return std::sqrt(sum);
}

void check_spectral_sampling() {
    constexpr size_t nm=3,nr=4,L=2*nr,q=3;
    constexpr size_t dimension=3*(nm+1)+6*nr;
    func::ComplexDynamicMatrix gamma(dimension,dimension,ComplexType{});
    for(size_t i=0;i<dimension;++i)gamma(i,i)=1.;
    func::FFTDenseComplexGaussianSampler native(gamma,nm,nr,0.2,-1.,0);
    func::FFTDenseComplexGaussianSampler refined(gamma,nm,nr,0.2,-1.,q);
    require(native.latent_dimension()==refined.latent_dimension(),"substeps preserve latent rank");
    require(native.largest_factorization_dimension()==refined.largest_factorization_dimension(),
            "substeps preserve factorization size");
    std::vector<RealType> times;
    const RealType c=RealType{0.5}-std::sqrt(RealType{3.})/6;
    for(size_t i=0;i<nr;++i)times.push_back(RealType(i));
    for(size_t node=0;node<2;++node)for(size_t i=0;i<(nr-1)*q;++i)
        times.push_back((RealType(i)+(node?1-c:c))/q);
    func::ComplexDynamicMatrix B(times.size(),refined.latent_dimension());
    // Recover the exact linear sampling map, avoiding Monte-Carlo tolerances.
    for(size_t k=0;k<refined.latent_dimension();++k) {
        func::JointComplexGaussianSampler::LatentVector r(refined.latent_dimension(),0.);
        r[k]=1.;
        const auto edges=refined.contour_field_from_latent(r,false);
        const auto gauss=refined.contour_field_from_latent(r,true);
        const auto original=native.contour_field_from_latent(r,false);
        for(size_t i=0;i<dimension;++i)
            require(std::abs(edges.edge_field[i]-original.edge_field[i])<1e-13,
                    "q=3 retains every native field including Matsubara endpoints");
        require(!original.has_real_gauss_fields(),"q=0 samples only native endpoints");
        size_t row=0;
        for(size_t i=0;i<nr;++i)B(row++,k)=edges.edge_field[3*(nm+1)+3*i];
        for(size_t node=0;node<2;++node)for(size_t i=0;i<(nr-1)*q;++i)
            B(row++,k)=gauss.real_gauss_fields[node][6*i];
    }
    const func::ComplexDynamicMatrix actual=B*trans(B);
    const func::ComplexDynamicMatrix actual_hermitian=B*ctrans(B);
    // Independent direct Fourier interpolation of white noise on the full
    // embedded grid, including negative frequencies and the positive Nyquist.
    func::ComplexDynamicMatrix P(times.size(),L,ComplexType{});
    const RealType pi=std::acos(RealType{-1.});
    for(size_t i=0;i<times.size();++i)for(size_t j=0;j<L;++j)for(size_t k=0;k<L;++k) {
        const auto signed_k=k<=L/2?int(k):int(k)-int(L);
        P(i,j)+=std::exp(ComplexType{0.,2*pi*signed_k*(times[i]-RealType(j))/L})/RealType(L);
    }
    const func::ComplexDynamicMatrix expected=P*trans(P),expected_hermitian=P*ctrans(P);
    RealType error{},herror{};
    for(size_t i=0;i<times.size();++i)for(size_t j=0;j<times.size();++j) {
        error=std::max(error,std::abs(actual(i,j)-expected(i,j)));
        herror=std::max(herror,std::abs(actual_hermitian(i,j)-expected_hermitian(i,j)));
    }
    require(error<1e-11,"substep pseudo-covariance equals direct signed-frequency interpolation");
    require(herror<1e-11,"substep Hermitian covariance equals direct signed-frequency interpolation");
}

void check_dense_sampling() {
    constexpr size_t nm=3,nr=4,dimension=3*(nm+1)+6*nr;
    func::ComplexDynamicMatrix gamma(dimension,dimension,ComplexType{});
    for(size_t i=0;i<dimension;++i)gamma(i,i)=ComplexType{1.,0.01*RealType(i)};
    auto native=func::make_complex_gaussian_sampler("dense",gamma,nm,nr,0.2,-1.,{1.,1.,1.},0);
    for(size_t q:{1,2,3,4}) {
        auto refined=func::make_complex_gaussian_sampler("dense",gamma,nm,nr,0.2,-1.,{1.,1.,1.},q);
        require(native->latent_dimension()==refined->latent_dimension(),"dense substeps preserve latent rank");
        require(native->largest_factorization_dimension()==refined->largest_factorization_dimension(),
                "dense substeps preserve factorization size");
        std::mt19937 original_engine(734),refined_engine(734);
        const auto original=native->draw_contour_field(original_engine,false);
        const auto sample=refined->draw_contour_field(refined_engine,true);
        require(original_engine==refined_engine,"dense substeps preserve RNG consumption");
        require(!sample.has_real_gauss_fields(),"dense substeps use edge interpolation");
        for(size_t i=0;i<dimension;++i)
            require(original.edge_field[i]==sample.edge_field[i],"dense substeps preserve sampled edge fields exactly");
    }
}

func::ComplexFieldVector field(RealType t,bool backward,bool equal,bool polynomial=false) {
    const RealType sign=backward&&!equal?RealType{-0.7}:RealType{1.};
    // Cubic interpolation is exact for this noncommuting complex field,
    // including the one-sided stencils at the first and last intervals.
    if(polynomial)
        return {ComplexType{sign*(1.+0.3*t*t*t),0.12*t},
                ComplexType{0.4*t*t,-sign*0.13},ComplexType{0.3,sign*0.2*(1.-t*t)}};
    return {ComplexType{sign*std::cos(t),0.12*std::sin(t)},
            ComplexType{0.4*std::sin(2*t),-sign*0.13},ComplexType{0.3,sign*0.2*std::cos(t)}};
}
func::ContourTrajectory trajectory(ps::ParameterSpace p,size_t q,bool equal,
                                   bool polynomial=false,bool gauss_nodes=true) {
    p.real_time_substeps=q;
    func::JointComplexGaussianSampler::ContourFieldSample sample;
    const contour::ContourLayout layout{p.num_TimePoints,p.num_RealTimePoints};
    sample.edge_field.resize(layout.dimension()); reset(sample.edge_field);
    for(size_t t=0;t<p.num_RealTimePoints;++t)for(size_t b=0;b<2;++b)for(size_t c=0;c<3;++c)
        sample.edge_field[layout.flat(b?contour::Branch::Backward:contour::Branch::Forward,t,c)]
            =field(t*p.delta_real_t,b,equal,polynomial)[c];
    const size_t n=p.num_RealTimeSteps*p.real_time_steps_per_interval();
    const RealType h=p.delta_real_t/p.real_time_steps_per_interval();
    if(p.uses_cf4()&&gauss_nodes)for(auto& g:sample.real_gauss_fields)g.resize(6*n);
    if(p.uses_cf4()&&gauss_nodes)
    for(size_t t=0;t<n;++t)for(size_t node=0;node<2;++node)for(size_t b=0;b<2;++b)for(size_t c=0;c<3;++c) {
        const RealType fraction=0.5+(node?1.:-1.)*std::sqrt(RealType{3.})/6;
        sample.real_gauss_fields[node][6*t+3*b+c]=field((t+fraction)*h,b,equal,polynomial)[c];
    }
    func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    for(size_t t=0;t<mean.size();++t)for(size_t c=0;c<3;++c)
        mean[t][c]=0.07*(c+1)*(1+t*p.delta_real_t
            +(polynomial?0.2*std::pow(t*p.delta_real_t,3):0.));
    return func::compute_contour_trajectory(p,sample,mean);
}
std::pair<Operator,Operator> total(const func::ContourTrajectory& tr) {
    Operator f=IDENTITY,b=IDENTITY;
    for(size_t t=1;t<tr.forward_steps.size();++t){f=tr.forward_steps[t]*f;b=b*tr.backward_steps[t];}
    return {f,b};
}
void check_propagation() {
    ps::ParameterSpace p;
    p.spin_float=0.5;p.num_HilbertSpaceDimension=2;
    p.spin_model=Physics::SpinModel{"ISO"};
    p.B=Physics::MagneticField{"z",0.2,0.,0.};
    p.noise=Physics::Noise{"none",0.};p.extra_interaction=Physics::ExtraInteraction{"none",0.};
    p.beta=0.8;p.num_TimeSteps=4;p.num_TimePoints=5;p.delta_t=0.2;
    p.Tmax=1.2;p.num_RealTimeSteps=4;p.num_RealTimePoints=5;p.delta_real_t=0.3;
    func::initialize_matrices(p);
    const auto endpoint=total(trajectory(p,0,true));
    require(frobenius(endpoint.second*endpoint.first-IDENTITY)<1e-12,
            "q=0 endpoint propagation closes equal complex branch fields");
    for(bool dense:{false,true}) {
        p.gaussian_factorization=dense?"dense":"fft";
        const auto coarse=trajectory(p,3,false,dense,!dense);
        auto fine=p;fine.num_RealTimeSteps*=3;fine.num_RealTimePoints=fine.num_RealTimeSteps+1;fine.delta_real_t/=3;
        // Supply analytically evaluated Gauss nodes independently of interpolation.
        fine.gaussian_factorization="fft";
        const auto reference=trajectory(fine,1,false,dense);
        require(coarse.forward_steps.size()==p.num_RealTimePoints,"measurements retain native grid");
        require(frobenius(coarse.imaginary_density_operator-reference.imaginary_density_operator)<1e-13,
                "Matsubara preparation is unchanged");
        for(size_t t=1;t<p.num_RealTimePoints;++t) {
            Operator f=IDENTITY,b=IDENTITY;
            for(size_t j=1;j<=3;++j) {
                f=reference.forward_steps[3*(t-1)+j]*f;
                b=b*reference.backward_steps[3*(t-1)+j];
            }
            require(frobenius(f-coarse.forward_steps[t])<1e-12,"forward microstep multiplication order");
            require(frobenius(b-coarse.backward_steps[t])<1e-12,"backward microstep multiplication order");
        }
        const auto same=total(trajectory(p,3,true,dense,!dense));
        require(frobenius(same.second*same.first-IDENTITY)<1e-12,"equal complex branch fields close algebraically");
        const auto ref=total(trajectory(p,128,false,dense,!dense));
        RealType previous{};
        for(size_t q:{1,2,4}) {
            const auto result=total(trajectory(p,q,false,dense,!dense));
            const RealType error=frobenius(result.first-ref.first)+frobenius(result.second-ref.second);
            if(q>1)require(previous/error>12.,"substep convergence matches propagator order");
            previous=error;
        }
    }
}
}
int main(){check_spectral_sampling();check_dense_sampling();check_propagation();return failures?1:0;}
