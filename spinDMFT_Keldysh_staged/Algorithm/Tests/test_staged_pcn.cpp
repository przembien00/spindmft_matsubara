#include"../Functions/Staged_Solver.h"
#include"../Functions/Joint_Matsubara_PCN.h"
#include<cmath>
#include<iostream>
#include<stdexcept>
#include<mpi.h>

namespace f=spinDMFT::Functions;
namespace ps=spinDMFT::Parameter_Space;
namespace {
void require(bool ok,const char* message)
{if(!ok)throw std::runtime_error(message);}

ps::ParameterSpace parameters()
{
    ps::ParameterSpace p;
    p.spin_float=0.5;p.num_HilbertSpaceDimension=2;
    p.spin_model=Physics::SpinModel{"ISO"};p.correlation_symmetry_type='D';
    p.B=Physics::MagneticField{"z",0.,0.,0.};
    p.noise=Physics::Noise{"none",0.};p.extra_interaction=Physics::ExtraInteraction{"none",0.};
    p.beta=1.;p.num_TimeSteps=4;p.num_TimePoints=5;p.delta_t=0.25;
    p.Tmax=0.2;p.num_RealTimeSteps=2;p.num_RealTimePoints=3;p.delta_real_t=0.1;
    p.sampling_strategy="pcn";p.mh_step_size=0.65;
    return p;
}

// A commuting complex field gives Z_M(x)=2 cosh(k x) exactly.
// This test includes both negative Re Z_M and a nontrivial complex phase.
class ScalarSampler : public f::JointComplexGaussianSampler
{
 public:
    ScalarSampler(size_t points,ComplexType k):m_points(points),m_k(k){}
    LatentVector draw_latent(std::mt19937& engine) override
    {LatentVector x(1);x[0]=m_normal(engine);return x;}
    FieldVector field_from_latent(const LatentVector& x) override
    {
        FieldVector z(size(),ComplexType{});
        for(size_t t=0;t<m_points;++t)z[3*t+2]=RealType{2.}*m_k*x[0];
        return z;
    }
    FieldVector draw(std::mt19937& engine) override{return field_from_latent(draw_latent(engine));}
    RealType reconstruction_error() const override{return 0.;}
    size_t latent_dimension() const override{return 1;}
    size_t size() const override{return 3*m_points;}
 private:
    size_t m_points;
    ComplexType m_k;
    std::normal_distribution<RealType> m_normal{0.,1.};
};

void check_chain()
{
    auto p=parameters();f::initialize_matrices(p);
    const ComplexType k{0.4,1.2};
    ScalarSampler sampler(p.num_TimePoints,k);
    f::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    std::mt19937 rng(2718);
    f::MatsubaraPCNChain chain(p,sampler,mean,rng);
    for(size_t i=0;i<2048;++i)chain.step();
    constexpr size_t blocks=32,per_block=2048;
    std::vector<ComplexType> ratios;
    ComplexType total_n{},total_d{};
    size_t negative{},rejected{};
    for(size_t b=0;b<blocks;++b)
    {
        ComplexType n{},d{};
        for(size_t i=0;i<per_block;++i)
        {
            const auto previous=chain.latent();
            if(!chain.step())
            {
                ++rejected;
                require(chain.latent()==previous,"rejection retains Matsubara state");
            }
            const auto z=RealType{2.}*std::cosh(k*chain.latent()[0]);
            require(std::abs(z-chain.partition())<1e-10,"acceptance uses only Matsubara trace");
            require(std::abs(std::abs(z)-chain.weight())<1e-10,"pCN weight is |Z_M|");
            if(std::real(z)<0.)++negative;
            const auto phase=z/chain.weight();
            d+=phase;n+=phase*chain.latent()[0]*chain.latent()[0];
        }
        ratios.push_back(n/d);total_n+=n;total_d+=d;
    }
    const auto estimate=total_n/total_d,exact=RealType{1.}+k*k;
    RealType variance{};
    for(auto value:ratios)variance+=std::norm(value-estimate);
    const RealType stderr=std::sqrt(variance/(blocks*(blocks-1)));
    require(std::abs(estimate-exact)<6*stderr+0.01,"phase-reweighted Gaussian moment matches analytic result");
    require(negative>100,"negative real partition functions remain in sampling support");
    require(rejected>100,"test exercises rejected pCN steps");
    require(chain.rejected_zero_weight()==0,"negative real weights are not zero-modulus rejections");
    std::cout<<"Reweighted x^2: "<<estimate<<" exact "<<exact<<" batch error "<<stderr<<'\n';

    // The same M chain is independent of real-time grids and mean values
    // beyond the origin, so no real propagation belongs in acceptance.
    auto other=p;other.real_time_evolution="staggered-strang";
    other.num_RealTimeSteps=8;other.num_RealTimePoints=9;
    f::MeanFieldTrajectory other_mean(other.num_RealTimePoints,FieldVector{9.,8.,7.});
    other_mean.front()=mean.front();
    std::mt19937 a_rng(19),b_rng(19);
    ScalarSampler a_sampler(p.num_TimePoints,k),b_sampler(p.num_TimePoints,k);
    f::MatsubaraPCNChain a(p,a_sampler,mean,a_rng),b(other,b_sampler,other_mean,b_rng);
    for(size_t i=0;i<32;++i)
    {
        require(a.step()==b.step(),"real evolution cannot change acceptance decisions");
        require(a.latent()==b.latent(),"real grid cannot change M states");
    }
}

void check_residual_refresh()
{
    f::ComplexDynamicMatrix gamma(6,6,ComplexType{});
    gamma(0,0)=1.;gamma(0,3)=gamma(3,0)=0.25;
    for(size_t i=3;i<6;++i)gamma(i,i)=0.7;
    f::DenseComplexGaussianSampler sampler(gamma);
    f::JointMatsubaraSampler m_sampler(sampler,3);
    std::mt19937 rng(17);
    const auto x=m_sampler.draw_latent(rng);
    const auto a=m_sampler.draw_joint_latent(x,rng);
    const auto b=m_sampler.draw_joint_latent(x,rng);
    require(a!=b,"fresh residuals when the Matsubara state is unchanged");
    const auto M=m_sampler.field_from_latent(x),va=sampler.field_from_latent(a),vb=sampler.field_from_latent(b);
    for(size_t i=0;i<3;++i)require(std::abs(va[i]-M[i])<1e-14&&std::abs(vb[i]-M[i])<1e-14,"fresh residuals leave M unchanged");
    require(std::abs(va[3]-vb[3])>1e-6,"real fields refresh at fixed M");
    require(sampler.reconstruction_error()<1e-12,"mixed and real pseudo-covariance reconstruction");
    bool rejected=false;
    try{m_sampler.draw_joint_latent(f::JointComplexGaussianSampler::LatentVector{},rng);}
    catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"shared latent dimension is checked");
}

void check_rejected_joint_refresh()
{
    auto p=parameters();f::initialize_matrices(p);
    const size_t nm=3*p.num_TimePoints;
    f::ComplexDynamicMatrix L(nm+3,2,ComplexType{});
    for(size_t t=0;t<p.num_TimePoints;++t)L(3*t+2,0)=ComplexType{0.8,2.4};
    L(nm,1)=1.;
    f::DenseComplexGaussianSampler joint(L*blaze::trans(L));
    f::JointMatsubaraSampler marginal(joint,nm);
    f::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    std::mt19937 rng(734);
    f::MatsubaraPCNChain chain(p,marginal,mean,rng);
    size_t rejected{};
    auto previous=joint.field_from_latent(marginal.draw_joint_latent(chain.latent(),rng));
    for(size_t i=0;i<512;++i)
    {
        const bool accepted=chain.step();
        const auto current=joint.field_from_latent(marginal.draw_joint_latent(chain.latent(),rng));
        if(!accepted)
        {
            ++rejected;
            for(size_t j=0;j<nm;++j)require(std::abs(current[j]-previous[j])<1e-12,
                "rejected joint pCN proposal retains every Matsubara field");
            require(std::abs(current[nm]-previous[nm])>1e-12,
                "rejected joint pCN proposal still refreshes real-only noise");
        }
        previous=current;
    }
    require(rejected>10,"joint test exercises actual pCN rejections");
}

void check_ratio_statistics(int world,int rank)
{
    char name[]="test_staged_pcn";char* argv[]{name};
    ps::ParameterSpace p(1,argv,world,rank);
    p.sampling_strategy="pcn";p.correlation_normalization="partition-function";
    p.correlation_symmetry_type='C';p.num_TimeSteps=0;p.num_TimePoints=1;
    p.num_RealTimeSteps=0;p.num_RealTimePoints=1;
    p.num_SamplesPerCore=512;p.num_Samples=512*world;p.num_blocks=32;
    f::rtd::RunTimeData cached(p,rank),scaled(p,rank),streamed(p,rank);
    std::vector<ComplexType> numerators,phases;
    for(size_t i=0;i<p.num_SamplesPerCore;++i)
    {
        const RealType s=static_cast<RealType>(i);
        const auto z=std::exp(ComplexType{0.,0.6*std::sin(0.17*s)+0.35});
        const auto n=ComplexType{0.2,-0.1}*z+ComplexType{0.3*std::sin(0.31*s),0.2*std::cos(0.27*s)};
        const RealType w=2.+(i%3);
        phases.push_back(z);numerators.push_back(n);
        f::rtd::MeasuredSample sample;
        sample.partition=w*z;sample.closure={w*z};
        sample.correlations.assign(cached.num_correlation_components(),w*n);
        sample.magnetization.assign(cached.num_magnetization_components(),w*n);
        cached.accumulate_sample(sample,1./w);
        // Overall likelihood normalization cannot change ratios, errors or tau.
        scaled.accumulate_sample(sample,7./w);
        streamed.begin_sample(sample.partition,1./w);
        streamed.accumulate_closed_contour_trace(0,w*z);
        for(size_t c=0;c<streamed.num_magnetization_components();++c)
            streamed.accumulate_magnetization(0,c,w*n);
        for(size_t c=0;c<streamed.num_correlation_components();++c)
            streamed.accumulate_edge_correlation(0,c,0,w*n);
        streamed.end_sample();
    }
    f::CorrelationSet g,g_scaled,g_streamed;f::MagTen re,im;
    cached.mpi_reduce_and_finalize(g,re,im);
    scaled.mpi_reduce_and_finalize(g_scaled,re,im);
    streamed.mpi_reduce_and_finalize(g_streamed,re,im);
    ComplexType mean_z{};for(auto z:phases)mean_z+=z;mean_z/=RealType(phases.size());
    const ComplexType estimate{g.Re[0][0][0],g.Im[0][0][0]};
    RealType var_re{},var_im{};
    for(size_t i=0;i<phases.size();++i)
    {
        const auto influence=(numerators[i]-estimate*phases[i])/mean_z;
        var_re+=std::real(influence)*std::real(influence);
        var_im+=std::imag(influence)*std::imag(influence);
    }
    var_re/=phases.size();var_im/=phases.size();
    const RealType count=RealType(p.num_Samples);
    const RealType error_re=cached.contour_sample_stds.Re[0][0][0];
    const RealType error_im=cached.contour_sample_stds.Im[0][0][0];
    require(std::abs(cached.contour_tau_int.Re[0][0][0]-0.5*error_re*error_re*count/var_re)<1e-10,
        "Re tau uses complex ratio influence variance");
    require(std::abs(cached.contour_tau_int.Im[0][0][0]-0.5*error_im*error_im*count/var_im)<1e-10,
        "Im tau uses complex ratio influence variance");
    require(std::abs(cached.maximum_tau_int.back()-scaled.maximum_tau_int.back())<1e-10,
        "ratio tau is invariant to overall likelihood normalization");
    require(std::abs(cached.maximum_tau_int.back()-streamed.maximum_tau_int.back())<1e-10,
        "cached and streamed samples have identical phase cross moments");
}
}
int main()
{
    MPI_Init(nullptr,nullptr);int world{},rank{},status{};
    MPI_Comm_size(MPI_COMM_WORLD,&world);MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    try{check_chain();check_residual_refresh();check_rejected_joint_refresh();check_ratio_statistics(world,rank);}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';status=1;}
    MPI_Finalize();return status;
}
