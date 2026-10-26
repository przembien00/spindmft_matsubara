#include"../Functions/Functions.h"
#include"../Functions/Staged_Solver.h"
#include<cmath>
#include<iostream>

namespace f=spinDMFT::Functions;
namespace c=spinDMFT::Contour;
namespace ps=spinDMFT::Parameter_Space;
using Sector=c::StaggeredSector;
namespace {
int failures{};
void require(bool ok,const char* message)
{ if(!ok){++failures;std::cerr<<"FAILED: "<<message<<'\n';} }
RealType frobenius(const Operator& x)
{ RealType sum{};for(size_t i=0;i<x.rows();++i)for(size_t j=0;j<x.columns();++j)sum+=std::norm(x(i,j));return std::sqrt(sum); }
ps::ParameterSpace parameters()
{
    ps::ParameterSpace p;
    p.spin_float=0.5;p.num_HilbertSpaceDimension=2;
    p.spin_model=Physics::SpinModel{"ISO"};p.correlation_symmetry_type='D';
    p.B=Physics::MagneticField{"z",0.2,0.,0.};
    p.noise=Physics::Noise{"none",0.};p.extra_interaction=Physics::ExtraInteraction{"none",0.};
    p.JQ=1.3;p.JL=-0.4;
    p.beta=0.8;p.num_TimeSteps=4;p.num_TimePoints=5;p.delta_t=0.2;
    p.Tmax=1.2;p.num_RealTimeSteps=4;p.num_RealTimePoints=5;p.delta_real_t=0.3;
    p.real_time_evolution="staggered-strang";
    return p;
}
ComplexType primitive_value(RealType t,size_t a,size_t b,size_t tau)
{
    return {0.2+0.01*(a+b)+0.03*t*t+0.02*t*t*t+0.004*tau,
            0.003*(RealType(a)-RealType(b))+0.002*t*RealType(tau)};
}
void check_covariance()
{
    auto p=parameters();p.real_time_substeps=3;
    c::CorrelationSet primitive{'D',p.num_TimePoints,p.num_RealTimePoints};
    for(size_t t=0;t<p.num_RealTimePoints;++t)for(size_t i=0;i<primitive.Re[t].size();++i)
    {
        const auto ab=primitive.Re[t].get_direction_pair(i);
        for(size_t tau=0;tau<p.num_TimePoints;++tau)
        {
            auto z=primitive_value(RealType(t),ab[0],ab[1],tau);
            primitive.Re[t][i][tau]=std::real(z);primitive.Im[t][i][tau]=std::imag(z);
        }
    }
    f::MagTen re{'D',p.num_RealTimePoints},im{'D',p.num_RealTimePoints};
    for(size_t t=0;t<re.size();++t)
    {
        re[t]=Observables::Magnetization::MagnetizationVector{'D',FieldVector{0.1,0.02,0.03}};
        im[t]=Observables::Magnetization::MagnetizationVector{'D',FieldVector{0.01,0.,-0.01}};
    }
    const auto connected=f::connected_contour_primitive(primitive, spinDMFT::Contour::Magnetization{re,im});
    const auto field=f::self_consistent_equations(p, primitive, spinDMFT::Contour::Magnetization{re,im});
    c::CorrelationSet preparation{'D',p.num_TimePoints,1};
    preparation.Re.front()=primitive.Re.front();
    preparation.Im.front()=primitive.Im.front();
    c::Magnetization preparation_mean{'D',1};
    preparation_mean.Re.front()=re.front();
    preparation_mean.Im.front()=im.front();
    const auto prepared=f::self_consistent_equations(p,preparation,preparation_mean);
    const size_t nm=3*p.num_TimePoints;
    require(prepared.covariance.rows()==nm&&prepared.covariance.columns()==nm,
        "preparation constructs only the Matsubara covariance");
    require(prepared.mean_time.size()==1,"preparation mean has only the origin");
    for(size_t i=0;i<nm;++i)for(size_t j=0;j<nm;++j)
        require(prepared.covariance(i,j)==field.covariance(i,j),
            "direct preparation covariance equals the joint Matsubara block");
    const c::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps*p.real_time_substeps};
    require(field.covariance.rows()==layout.dimension(),"staggered joint dimension");
    RealType error{};
    for(size_t i=0;i<layout.dimension();++i)for(size_t j=0;j<layout.dimension();++j)
    {
        const auto a=layout.decode(i),b=layout.decode(j);
        const auto actual=field.covariance(i,j);
        const auto expected=p.JQ*p.JQ*c::staggered_correlation(connected,a,b,p.real_time_substeps);
        error=std::max(error,std::abs(actual-expected));
        require(actual==field.covariance(j,i),"joint covariance uses transpose symmetry");
        if((a.sector==Sector::Nu&&b.sector!=Sector::Eta)
            ||(a.sector==Sector::Eta&&b.sector==Sector::Nu&&b.point>a.point))
            require(actual==ComplexType{},"forbidden response contractions are exactly zero");
    }
    require(error<1e-13,"connected primitive and JQ rotation retained");
    const ComplexType ma{0.1,0.01},mb{0.02,0.};
    const size_t eta=layout.flat(Sector::Eta,2,0),m=layout.flat(Sector::Matsubara,3,1);
    require(std::abs(field.covariance(eta,m)-p.JQ*p.JQ*(primitive_value(2.5/3,0,1,3)-ma*mb))<1e-13,
        "mixed kernel is interpolated directly at the midpoint including complex mean subtraction");
    const size_t nu=layout.flat(Sector::Nu,2,1);
    require(std::abs(field.covariance(eta,nu)-p.JQ*p.JQ*
        (primitive_value(0.5/3,0,1,4)-primitive_value(0.5/3,0,1,0)))<1e-13,
        "earlier response contact uses the smooth half-integer lag");
    {
        const auto make_sampler=[&]() -> std::unique_ptr<f::JointComplexGaussianSampler>
        {
            return f::make_staged_real_sampler(p,field.covariance);
        };
        auto sampler=make_sampler();
        f::ComplexDynamicMatrix L(layout.dimension(),sampler->latent_dimension());
        for(size_t k=0;k<L.columns();++k)
        {
            f::JointComplexGaussianSampler::LatentVector r(L.columns(),0.);r[k]=1.;
            const auto z=sampler->field_from_latent(r);
            for(size_t i=0;i<L.rows();++i)L(i,k)=z[i];
        }
        const f::ComplexDynamicMatrix reconstructed=L*trans(L);
        RealType residual{};
        for(size_t i=0;i<L.rows();++i)for(size_t j=0;j<L.rows();++j)
            residual=std::max(residual,std::abs(reconstructed(i,j)-field.covariance(i,j)));
        require(residual<1e-11,"latent basis reconstructs full joint staggered covariance");
        {
            // Unconjugated reconstruction alone cannot detect a wrong
            // sampling gauge. Check the conjugated covariance as well,
            // against an independent Takagi factor in the branch norm.
            const auto weights=f::staggered_sampling_weights(p);
            f::ComplexDynamicMatrix gamma=field.covariance;
            for(size_t i=0;i<gamma.rows();++i)for(size_t j=0;j<gamma.columns();++j)
                gamma(i,j)*=std::sqrt(weights[i]*weights[j]);
            auto expected=f::autonne_takagi(gamma).L;
            for(size_t i=0;i<expected.rows();++i)for(size_t j=0;j<expected.columns();++j)
                expected(i,j)/=std::sqrt(weights[i]);
            const f::ComplexDynamicMatrix actual_H=L*ctrans(L),expected_H=expected*ctrans(expected);
            RealType gauge_error{};
            for(size_t i=0;i<L.rows();++i)for(size_t j=0;j<L.rows();++j)
                gauge_error=std::max(gauge_error,std::abs(actual_H(i,j)-expected_H(i,j)));
            require(gauge_error<1e-10,"production staggered sampler preserves the branch-equivalent conjugated covariance");
        }
        require(sampler->reconstruction_error()<1e-11,"original-basis reconstruction diagnostic");
        auto reference_sampler=make_sampler();
        std::mt19937 batch_engine{321},single_engine{321};
        const auto batch=sampler->draw_contour_batch(batch_engine,3);
        for(const auto& sample:batch)
        {
            const auto single=reference_sampler->draw(single_engine);
            for(size_t i=0;i<single.size();++i)
                require(std::abs(single[i]-sample.edge_field[i])<1e-12,"batch and single latent draws agree");
        }
    }
    // Coupling rotation and classical noise do not create response noise.
    p.spin_model=Physics::SpinModel{"DRF"};p.noise=Physics::Noise{"z",0.1};
    const auto rotated=f::self_consistent_equations(p, primitive, spinDMFT::Contour::Magnetization{re,im});
    const auto& D=p.spin_model.coupling_matrix;
    ComplexType expected{};
    for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)
        expected+=p.JQ*p.JQ*D(0,a)*D(1,b)*c::staggered_correlation(connected,
            {Sector::Eta,2,a},{Sector::Nu,2,b},3);
    require(std::abs(rotated.covariance(eta,nu)-expected)<1e-13,"spin coupling rotation retains causal response");
}
f::ComplexFieldVector eta(RealType t)
{return {{0.3+0.2*t,0.13*t},{0.17*t*t,-0.11},{0.2,0.07*t}};}
f::ComplexFieldVector nu(RealType t)
{return {{0.12,0.05*t},{0.13*t,-0.07},{-0.09*t*t,0.08}};}
Operator hamiltonian(const f::ComplexFieldVector& field)
{return field[0]*S_X+field[1]*S_Y+field[2]*S_Z;}
void check_exponential_pair()
{
    auto p=parameters();f::initialize_matrices(p);
    // Includes a nonzero nilpotent field with q^2=0, the small-q series,
    // non-Hermitian fields, and complex scalar shifts.
    const std::array<f::ComplexFieldVector,4> fields{
        f::ComplexFieldVector{},
        f::ComplexFieldVector{ComplexType{1.,0.},ComplexType{0.,1.},ComplexType{}},
        f::ComplexFieldVector{ComplexType{1e-7,2e-8},ComplexType{3e-8,-1e-7},ComplexType{2e-7,0.}},
        eta(0.7)};
    for(const auto& field:fields)
        for(ComplexType scalar:{ComplexType{},ComplexType{0.3,0.2}})
            for(ComplexType step:{ComplexType{},ComplexType{-0.2,0.},ComplexType{0.,-0.3},ComplexType{0.2,-0.4}})
            {
                const auto pair=f::spin_half_field_exponential_pair(field,scalar,step);
                const Operator H=hamiltonian(field)+scalar*IDENTITY;
                require(frobenius(pair.first-f::general_matrix_exponential(step*H))<1e-12,
                    "paired forward exponential agrees with independent matrix exponential");
                require(frobenius(pair.second-f::general_matrix_exponential(-step*H))<1e-12,
                    "paired inverse includes the complex scalar and is not an adjoint");
                require(frobenius(pair.second*pair.first-IDENTITY)<1e-12,
                    "paired complex exponentials are mutual inverses");
            }
}

void check_skipped_final_trace()
{
    for(const std::string evolution:{"staggered-strang"})
    {
        auto p=parameters();p.real_time_evolution=evolution;
        p.num_SamplesPerCore=p.num_Samples=p.num_blocks=1;
        f::initialize_matrices(p);
        f::JointComplexGaussianSampler::ContourFieldSample field;
        const size_t dimension=p.uses_staggered_strang()
            ?c::StaggeredLayout{p.num_TimePoints,p.num_RealTimeSteps}.dimension()
            :c::ContourLayout{p.num_TimePoints,p.num_RealTimePoints}.dimension();
        field.edge_field.resize(dimension);
        for(size_t i=0;i<dimension;++i)
            field.edge_field[i]=ComplexType{0.1*std::sin(RealType(i)),0.05*std::cos(RealType(i))};
        const f::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
        spinDMFT::Run_Time_Data::RunTimeData data(p,0);
        spinDMFT::Contour::CorrelationSet data_sums;
        spinDMFT::Contour::Magnetization data_magnetization;
        data.begin_iteration_accumulation(data_sums,data_magnetization);
        spinDMFT::Functions::MeasuredSample full,skipped;
        f::MeasurementWorkspace measurement;
        f::TrajectoryWorkspace workspace;
        f::ContourTrajectory tr;
        f::build_contour_trajectory(p,field,mean,tr,workspace);
        const ComplexType final_trace=tr.final_closed_contour_trace;
        f::measure_contour_observables(data,tr,full,measurement);
        // Reuse the same trajectory to detect an accidentally stale final trace.
        f::build_contour_trajectory(p,field,mean,tr,workspace,false,f::FinalTraceMode::Skip);
        require(std::isnan(std::real(tr.final_closed_contour_trace)),"skipped trace is marked unavailable");
        f::measure_contour_observables(data,tr,skipped,measurement);
        require(full.partition==skipped.partition&&full.correlations==skipped.correlations
            &&full.magnetization==skipped.magnetization&&full.closure==skipped.closure,
            "skipping the redundant trace preserves every independent observable");
        require(std::abs(skipped.closure.back()-final_trace)<1e-12,
            "measurement still supplies the independently computed final closure");
        f::complete_contour_trajectory(p,field,mean,tr,workspace);
        require(std::abs(tr.final_closed_contour_trace-final_trace)<1e-12,
            "default completion continues to populate the final trace");
    }
}

void check_static_trajectory()
{
    for(const std::string evolution:{"staggered-strang"})
        for(size_t substeps:{size_t{1},size_t{3}})
            for(bool preparation:{false,true})
            {
                auto p=parameters();p.real_time_evolution=evolution;p.real_time_substeps=substeps;
                p.num_SamplesPerCore=p.num_Samples=p.num_blocks=1;
                if(preparation){p.num_RealTimeSteps=0;p.num_RealTimePoints=1;p.Tmax=0.;}
                f::initialize_matrices(p);
                f::JointComplexGaussianSampler::ContourFieldSample field;
                const size_t dimension=p.uses_staggered_strang()
                    ?c::StaggeredLayout{p.num_TimePoints,p.num_RealTimeSteps*substeps}.dimension()
                    :c::ContourLayout{p.num_TimePoints,p.num_RealTimePoints}.dimension();
                field.edge_field.resize(dimension);
                for(size_t i=0;i<dimension;++i)
                    field.edge_field[i]=ComplexType{0.2*std::sin(RealType(i)),0.08*std::cos(RealType(i))};
                f::MeanFieldTrajectory mean(p.num_RealTimePoints);
                for(size_t t=0;t<mean.size();++t)mean[t]=FieldVector{0.01*RealType(t),-0.02,0.03};
                f::ContourTrajectory dynamic;
                f::TrajectoryWorkspace dynamic_workspace;
                f::SpinHalfContourTrajectory fixed;
                f::SpinHalfTrajectoryWorkspace fixed_workspace;
                f::build_contour_trajectory(p,field,mean,dynamic,dynamic_workspace);
                // Exercise split build/completion as well as the stored static prefixes.
                f::build_contour_trajectory(p,field,mean,fixed,fixed_workspace,true);
                f::complete_contour_trajectory(p,field,mean,fixed,fixed_workspace);
                const auto same_matrix=[&](const auto& a,const auto& b)
                    {require(frobenius(Operator(a)-Operator(b))<2e-12,"fixed and dynamic trajectory matrices agree");};
                same_matrix(fixed.imaginary_density_operator,dynamic.imaginary_density_operator);
                for(size_t t=0;t<p.num_RealTimePoints;++t)
                {
                    same_matrix(fixed.forward_steps[t],dynamic.forward_steps[t]);
                    same_matrix(fixed.backward_steps[t],dynamic.backward_steps[t]);
                }
                for(size_t a=0;a<3;++a)for(size_t tau=0;tau<p.num_TimePoints;++tau)
                    same_matrix(fixed.imaginary_edge_insertions[a][tau],dynamic.imaginary_edge_insertions[a][tau]);
                require(std::abs(fixed.partition_function-dynamic.partition_function)<2e-12
                    &&std::abs(fixed.final_closed_contour_trace-dynamic.final_closed_contour_trace)<2e-12,
                    "fixed-size partition and final trace agree");
                spinDMFT::Run_Time_Data::RunTimeData data(p,0);
                spinDMFT::Contour::CorrelationSet data_sums;
                spinDMFT::Contour::Magnetization data_magnetization;
                data.begin_iteration_accumulation(data_sums,data_magnetization);
                spinDMFT::Functions::MeasuredSample a,b;
                f::MeasurementWorkspace measurement;
                {
                    f::measure_contour_observables(data,dynamic,a,measurement);
                    f::measure_contour_observables(data,fixed,b,measurement);
                    const auto same_values=[&](const auto& x,const auto& y)
                    {
                        require(x.size()==y.size(),"fixed-size observable shape");
                        for(size_t i=0;i<x.size();++i)
                            require(std::abs(x[i]-y[i])<2e-12,"fixed-size observables agree with dynamic reference");
                    };
                    same_values(a.correlations,b.correlations);
                    same_values(a.magnetization,b.magnetization);
                    same_values(a.closure,b.closure);
                }
            }
}

void check_prepared_mean()
{
    const auto polynomial=[](RealType t)
        {return FieldVector{0.01*(t*t*t+2*t),0.02*t*t-0.03,0.03*t+0.01};};
    for(const std::string evolution:{"staggered-strang"})
        for(size_t substeps:{size_t{1},size_t{3}})
            for(size_t points:{size_t{1},size_t{2},size_t{3},size_t{5}})
            {
                auto p=parameters();p.real_time_evolution=evolution;p.real_time_substeps=substeps;
                p.num_RealTimePoints=points;p.num_RealTimeSteps=points-1;p.Tmax=(points-1)*p.delta_real_t;
                f::MeanFieldTrajectory mean(points);
                for(size_t t=0;t<points;++t)mean[t]=polynomial(RealType(t));

                f::initialize_matrices(p);
                const f::PreparedMeanField prepared(p,mean);
                const size_t q=p.real_time_steps_per_interval();
                for(size_t t=0;t<p.num_RealTimeSteps;++t)for(size_t j=0;j<q;++j)
                {
                    if(p.uses_staggered_strang())
                    {
                        const RealType u=(RealType(j)+0.5)/RealType(q);
                        const FieldVector expected=points>=4?polynomial(RealType(t)+u)
                            :FieldVector((1-u)*mean[t]+u*mean[t+1]);
                        require(blaze::length(prepared.midpoints()[t*q+j]-expected)<1e-13,
                            "cached midpoints reproduce cubic or short-grid linear interpolation");
                    }

                }
                const auto original=mean;
                mean.front()[0]+=0.2; // Must not alter the immutable prepared snapshot.
                const f::PreparedMeanField refreshed(p,mean);
                const size_t dimension=p.uses_staggered_strang()
                    ?c::StaggeredLayout{p.num_TimePoints,p.num_RealTimeSteps*q}.dimension()
                    :c::ContourLayout{p.num_TimePoints,p.num_RealTimePoints}.dimension();
                f::JointComplexGaussianSampler::ContourFieldSample field;
                field.edge_field.resize(dimension);
                f::SpinHalfContourTrajectory reference,cached;
                f::SpinHalfTrajectoryWorkspace reference_workspace,cached_workspace;
                const auto same=[&]()
                {
                    require(frobenius(Operator(reference.imaginary_density_operator-cached.imaginary_density_operator))<1e-13,
                        "prepared mean preserves the thermal density");
                    require(std::abs(reference.final_closed_contour_trace-cached.final_closed_contour_trace)<1e-13,
                        "prepared mean preserves final trace");
                    for(size_t t=0;t<points;++t)
                    {
                        require(frobenius(Operator(reference.forward_steps[t]-cached.forward_steps[t]))<1e-13,
                            "cached mean preserves forward propagation");
                        require(frobenius(Operator(reference.backward_steps[t]-cached.backward_steps[t]))<1e-13,
                            "cached mean preserves backward propagation");
                    }
                };
                for(size_t sample=0;sample<2;++sample)
                {
                    for(size_t i=0;i<dimension;++i)
                        field.edge_field[i]=ComplexType{0.1*std::sin(RealType(i+sample)),0.03*std::cos(RealType(i))};
                    f::build_contour_trajectory(p,field,original,reference,reference_workspace);
                    f::build_contour_trajectory(p,field,prepared,cached,cached_workspace,true);
                    f::complete_contour_trajectory(p,field,prepared,cached,cached_workspace);
                    same();
                    f::build_contour_trajectory(p,field,mean,reference,reference_workspace);
                    f::build_contour_trajectory(p,field,refreshed,cached,cached_workspace);
                    same();
                }
                auto changed=p;++changed.real_time_substeps;
                bool rejected=false;
                try{prepared.validate_grid(changed);}catch(const std::invalid_argument&){rejected=true;}
                require(rejected,"prepared mean rejects a changed substep grid");
            }
}
f::ContourTrajectory trajectory(ps::ParameterSpace p,bool zero_response=false,bool constant=false)
{
    const size_t n=p.num_RealTimeSteps*p.real_time_substeps;
    const RealType h=p.Tmax/RealType(n);
    const c::StaggeredLayout layout{p.num_TimePoints,n};
    f::JointComplexGaussianSampler::ContourFieldSample sample;
    sample.edge_field.resize(layout.dimension());reset(sample.edge_field);
    for(size_t i=0;i<n;++i)for(size_t a=0;a<3;++a)
        sample.edge_field[layout.flat(Sector::Eta,i,a)]=eta(constant?0.4:(RealType(i)+0.5)*h)[a];
    for(size_t i=0;i<=n;++i)for(size_t a=0;a<3;++a)
        sample.edge_field[layout.flat(Sector::Nu,i,a)]=zero_response?ComplexType{}:nu(constant?0.4:RealType(i)*h)[a];
    return f::compute_contour_trajectory(p,sample,f::MeanFieldTrajectory(p.num_RealTimePoints,FieldVector{}));
}
std::pair<Operator,Operator> total(const f::ContourTrajectory& tr)
{
    Operator U=IDENTITY,B=IDENTITY;
    for(size_t i=1;i<tr.forward_steps.size();++i){U=tr.forward_steps[i]*U;B=B*tr.backward_steps[i];}
    return {U,B};
}
// Independent very fine unsplit midpoint reference for H+=Heta+Hnu/2
// and H-=Heta-Hnu/2, including the anti-time ordering of B.
std::pair<Operator,Operator> reference(const ps::ParameterSpace& p)
{
    constexpr size_t n=8192;const RealType h=p.Tmax/RealType(n);
    Operator U=IDENTITY,B=IDENTITY;
    for(size_t i=0;i<n;++i)
    {
        const RealType t=(RealType(i)+0.5)*h;
        auto e=eta(t);for(size_t a=0;a<3;++a)e[a]+=p.B.m_h[a];
        U=f::general_matrix_exponential(ComplexType{0.,-h}*hamiltonian(f::ComplexFieldVector(e+0.5*nu(t))))*U;
        B=B*f::general_matrix_exponential(ComplexType{0.,h}*hamiltonian(f::ComplexFieldVector(e-0.5*nu(t))));
    }
    return {U,B};
}
void check_propagation()
{
    auto p=parameters();f::initialize_matrices(p);
    const auto ref=reference(p);RealType previous{};
    for(size_t q:{1,2,4,8})
    {
        p.real_time_substeps=q;
        const auto result=total(trajectory(p));
        const RealType error=frobenius(result.first-ref.first)+frobenius(result.second-ref.second);
        if(q>1)require(previous/error>3.8&&previous/error<4.2,"second-order prescribed-field convergence");
        previous=error;
        const auto equal=total(trajectory(p,true));
        require(frobenius(equal.second*equal.first-IDENTITY)<1e-12,"zero nu closes complex precession using actual inverse");
    }
    p.real_time_substeps=3;
    const auto coarse=trajectory(p);
    auto fine=p;fine.num_RealTimeSteps*=3;fine.num_RealTimePoints=fine.num_RealTimeSteps+1;
    fine.delta_real_t/=3;fine.real_time_substeps=1;
    const auto refined=trajectory(fine);
    require(coarse.forward_steps.size()==p.num_RealTimePoints,"native measurement grid retained");
    for(size_t i=1;i<p.num_RealTimePoints;++i)
    {
        Operator U=IDENTITY,B=IDENTITY;
        for(size_t j=1;j<=3;++j){U=refined.forward_steps[3*(i-1)+j]*U;B=B*refined.backward_steps[3*(i-1)+j];}
        require(frobenius(U-coarse.forward_steps[i])+frobenius(B-coarse.backward_steps[i])<1e-12,
            "microsteps compose with shared response boundaries");
    }
    // Freeze noncommuting complex fields: explicitly traverse half/full/half
    // kicks and independently check that adjacent half kicks combine.
    const auto tr=trajectory(p,false,true);
    auto e=eta(0.4);for(size_t a=0;a<3;++a)e[a]+=p.B.m_h[a];
    const RealType h=p.delta_real_t/3;
    const Operator Q=f::general_matrix_exponential(ComplexType{0.,-h/4}*hamiltonian(nu(0.4)));
    const Operator P=f::general_matrix_exponential(ComplexType{0.,-h}*hamiltonian(e));
    Operator U=Q,B=Q;
    for(size_t i=0;i<p.num_RealTimeSteps*3;++i)
    {U=(i+1==p.num_RealTimeSteps*3?Q:Operator(Q*Q))*P*U;B=B*blaze::inv(P)*(i+1==p.num_RealTimeSteps*3?Q:Operator(Q*Q));}
    const auto actual=total(tr);
    require(frobenius(U-actual.first)+frobenius(B-actual.second)<1e-12,"first and last response half kicks and shared full interior kicks");
    // Exercise the general-spin exponential path as well.
    p.spin_float=1.;p.num_HilbertSpaceDimension=3;f::initialize_matrices(p);
    const auto spin_one=total(trajectory(p,true));
    require(frobenius(spin_one.second*spin_one.first-IDENTITY)<1e-12,"general-spin complex inverse closure");
    p.extra_interaction=Physics::ExtraInteraction{"quadrupolar",0.3};f::initialize_matrices(p);
    const auto with_rest=total(trajectory(p,true,true));
    const Operator exact=f::general_matrix_exponential(ComplexType{0.,-p.Tmax}*(hamiltonian(e)+H_REST));
    require(frobenius(with_rest.first-exact)<1e-12,"non-scalar extra interaction belongs entirely to precession");
    require(frobenius(with_rest.second*with_rest.first-IDENTITY)<1e-12,"general-spin inverse includes non-scalar extra interaction");
}

void check_ensemble_closure()
{
    auto p=parameters();p.num_RealTimeSteps=1;p.num_RealTimePoints=2;
    p.Tmax=0.3;p.delta_real_t=0.3;f::initialize_matrices(p);
    const c::StaggeredLayout grid{p.num_TimePoints,1};
    f::ComplexDynamicMatrix gamma(grid.dimension(),grid.dimension(),ComplexType{});
    const size_t e=grid.flat(Sector::Eta,0,1),v=grid.flat(Sector::Nu,0,0);
    gamma(e,e)=0.002;gamma(e,v)=gamma(v,e)=ComplexType{0.,0.001};
    for(size_t i=0;i<p.num_TimePoints;++i)
    {
        const size_t m=grid.flat(Sector::Matsubara,i,2);
        gamma(m,e)=gamma(e,m)=ComplexType{0.0003,0.0002};
        for(size_t j=0;j<p.num_TimePoints;++j)gamma(m,grid.flat(Sector::Matsubara,j,2))=0.002;
    }
    f::DenseComplexGaussianSampler sampler(gamma);
    require(sampler.latent_dimension()==3,"closure quadrature has three joint latents");
    if(sampler.latent_dimension()!=3)return;
    // Five-node standard-normal Gaussian quadrature, exact through degree 9.
    // Small fields make the remaining analytic-exponential tail negligible.
    const RealType r=std::sqrt(RealType{10.});
    const std::array<RealType,5> x{0.,std::sqrt(5-r),-std::sqrt(5-r),std::sqrt(5+r),-std::sqrt(5+r)};
    const std::array<RealType,5> w{8./15,(7+2*r)/60,(7+2*r)/60,(7-2*r)/60,(7-2*r)/60};
    ComplexType Z{},D{};RealType sample_defect{};
    for(size_t i=0;i<5;++i)for(size_t j=0;j<5;++j)for(size_t k=0;k<5;++k)
    {
        f::JointComplexGaussianSampler::LatentVector latent{x[i],x[j],x[k]};
        const auto tr=f::compute_contour_trajectory(p,sampler.field_from_latent(latent),
            f::MeanFieldTrajectory(p.num_RealTimePoints,FieldVector{}));
        Z+=w[i]*w[j]*w[k]*tr.partition_function;
        D+=w[i]*w[j]*w[k]*tr.final_closed_contour_trace;
        sample_defect=std::max(sample_defect,std::abs(tr.final_closed_contour_trace-tr.partition_function));
    }
    require(std::abs(D-Z)<1e-12,"joint Gaussian ensemble closes with noncommuting mixed thermal fields");
    require(sample_defect>1e-6,"closure is an ensemble identity, not samplewise normalization");
}
}
int main(){check_covariance();check_exponential_pair();check_skipped_final_trace();check_static_trajectory();check_prepared_mean();check_propagation();check_ensemble_closure();return failures?1:0;}
