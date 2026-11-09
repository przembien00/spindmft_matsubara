#include"../Functions/Staged_Solver.h"
#include<cmath>
#include<iostream>

namespace f=spinDMFT::Functions;
namespace c=spinDMFT::Contour;
namespace ps=spinDMFT::Parameter_Space;
int failures{};
void require(bool ok,const char* message)
{if(!ok){++failures;std::cerr<<"FAILED: "<<message<<'\n';}}

RealType difference(const Operator& a,const Operator& b)
{return f::matrix_difference(a,b);}

ps::ParameterSpace parameters(RealType spin)
{
    ps::ParameterSpace p;
    p.spin_float=spin;p.num_HilbertSpaceDimension=size_t(2*spin+1);
    p.spin_model=Physics::SpinModel{"ISO"};p.correlation_symmetry_type='D';
    p.time_translation_invariant=false;
    p.B=Physics::MagneticField{"z",1.3,0.,0.};
    p.noise=Physics::Noise{"none",0.};p.extra_interaction=Physics::ExtraInteraction{"none",0.};
    p.JQ=0.4;p.JL=-0.3;p.beta=0.8;
    p.num_TimeSteps=4;p.num_TimePoints=5;p.delta_t=0.2;
    p.num_RealTimeSteps=4;p.num_RealTimePoints=5;p.delta_real_t=0.3;p.Tmax=1.2;
    return p;
}

void check_rotation(RealType spin)
{
    auto p=parameters(spin);f::initialize_matrices(p);
    const RealType pi=std::acos(RealType{-1.});
    const auto R=f::pulse_rotation({0.,FieldVector{0.,0.,1.},pi/2});
    Operator expected=ZERO;
    for(size_t j=0;j<p.num_HilbertSpaceDimension;++j)
        expected(j,j)=std::exp(ComplexType{0.,pi/2*(spin-RealType(j))});
    require(difference(R,expected)<2e-13,"rotation has positive exponent and correct spin prefactor");
    require(difference(blaze::ctrans(R)*R,Operator(IDENTITY))<2e-13,"pulse rotation is unitary");
    require(difference(blaze::ctrans(R)*S_X*R,Operator(S_Y))<2e-13,"positive z pi/2 maps Heisenberg Sx to Sy");
    const auto Ry=f::pulse_rotation({0.,FieldVector{0.,1.,0.},pi/2});
    require(difference(blaze::ctrans(Ry)*S_X*Ry,Operator(-S_Z))<2e-13,"positive y pi/2 maps Heisenberg Sx to minus Sz");
    if(spin==0.5)
        require(difference(f::pulse_rotation({0.,FieldVector{1.,0.,0.},2*pi}),Operator(-IDENTITY))<2e-13,
            "spin-half 2pi rotation is minus identity");
}

void check_times(RealType spin,size_t q)
{
    auto p=parameters(spin);p.real_time_substeps=q;f::initialize_matrices(p);
    p.pulses={{0.,FieldVector{0.,1.,0.},0.7},
              {0.,FieldVector{0.,0.,1.},-0.2},
              {0.15,FieldVector{1.,0.,0.},0.3},
              {0.17,FieldVector{1.,0.,0.},0.9},
              {0.17,FieldVector{0.,1.,0.},-0.5}};
    const c::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps*q};
    f::JointComplexGaussianSampler::ContourFieldSample fields;
    fields.edge_field.resize(layout.dimension());fields.edge_field=ComplexType{};
    const f::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    f::ContourTrajectory tr;f::TrajectoryWorkspace workspace;
    const f::PreparedMeanField cached(p,mean);
    f::build_contour_trajectory(p,fields,cached,tr,workspace);
    const auto direct=f::compute_contour_trajectory(p,fields,mean);
    Operator actual=IDENTITY,backward=IDENTITY,exact=IDENTITY;
    size_t next{};RealType last{};
    const auto E=[&](RealType dt){return f::general_matrix_exponential(ComplexType{0.,-1.3*dt}*S_Z);};
    for(size_t t=0;t<p.num_RealTimePoints;++t)
    {
        const RealType time=RealType(t)*p.delta_real_t;
        while(next<p.pulses.size()&&p.pulses[next].time<=time+1e-14)
        {
            const auto& pulse=p.pulses[next++];
            const Operator generator=pulse.axis[0]*S_X+pulse.axis[1]*S_Y+pulse.axis[2]*S_Z;
            exact=f::general_matrix_exponential(ComplexType{0.,pulse.angle}*generator)*E(pulse.time-last)*exact;
            last=pulse.time;
        }
        exact=E(time-last)*exact;last=time;
        actual=tr.forward_steps[t]*actual;backward=backward*tr.backward_steps[t];
        require(difference(actual,exact)<3e-12,"physical-time pulse order matches exact piecewise propagation at every edge");
        require(difference(backward,Operator(blaze::ctrans(exact)))<3e-12,"backward pulse order is the inverse chronological order");
        require(difference(tr.forward_steps[t],direct.forward_steps[t])<2e-13,"cached and direct pulse propagation agree");
    }
    require(std::abs(tr.final_closed_contour_trace-tr.partition_function)<3e-12,
        "origin and off-grid pulses preserve the closed-contour trace");
}

void check_measurement_fields(RealType spin,size_t q)
{
    auto p=parameters(spin);p.real_time_substeps=q;f::initialize_matrices(p);
    p.pulses={{0.3,FieldVector{1.,0.,0.},0.7},
              {0.3,FieldVector{0.,1.,0.},-0.4},
              {1.2,FieldVector{0.,0.,1.},0.6}};
    const c::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps*q};
    f::JointComplexGaussianSampler::ContourFieldSample fields;
    fields.edge_field.resize(layout.dimension());fields.edge_field=ComplexType{};
    const f::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    const f::PreparedMeanField cached(p,mean);
    require(cached.pulses().empty(),"measurement-grid pulses do not create separate rotations");
    const RealType h=p.delta_real_t/RealType(q);
    require(std::abs(cached.pulse_fields()[q-1][0]+0.7/h)<2e-13
        &&std::abs(cached.pulse_fields()[q-1][1]-0.4/h)<2e-13,
        "coincident pulses add field vectors with minus angle over substep width");
    const auto tr=f::compute_contour_trajectory(p,fields,mean);
    Operator actual=IDENTITY,backward=IDENTITY,exact=IDENTITY;
    for(size_t t=0;t<p.num_RealTimePoints;++t)
    {
        if(t>0)for(size_t j=0;j<q;++j)
        {
            Operator generator=ComplexType{0.,-h*1.3}*S_Z;
            if(j+1==q&&t==1)generator+=ComplexType{0.,0.7}*S_X+ComplexType{0.,-0.4}*S_Y;
            if(j+1==q&&t==4)generator+=ComplexType{0.,0.6}*S_Z;
            exact=f::general_matrix_exponential(generator)*exact;
        }
        actual=tr.forward_steps[t]*actual;backward=backward*tr.backward_steps[t];
        require(difference(actual,exact)<3e-12,"measurement-edge pulses share the background Hamiltonian exponential");
        require(difference(backward,Operator(blaze::ctrans(exact)))<3e-12,
            "backward Hamiltonian field includes the same pulse with reversed evolution sign");
    }
    const Operator separated=f::pulse_rotation(p.pulses[1])*f::pulse_rotation(p.pulses[0])
        *f::general_matrix_exponential(ComplexType{0.,-1.3*p.delta_real_t}*S_Z);
    require(difference(tr.forward_steps[1],separated)>1e-3,
        "noncommuting measurement pulses use the combined generator rather than sequential rotations");
    require(std::abs(tr.final_closed_contour_trace-tr.partition_function)<3e-12,
        "Hamiltonian pulses preserve the zero-response contour trace");
}

void check_response()
{
    auto p=parameters(0.5);p.num_RealTimeSteps=1;p.num_RealTimePoints=2;p.Tmax=0.3;
    p.pulses={{0.1,FieldVector{1.,0.,0.},0.7},{0.3,FieldVector{0.,1.,0.},-0.4}};
    f::initialize_matrices(p);
    const c::StaggeredLayout layout{p.num_TimePoints,1};
    f::JointComplexGaussianSampler::ContourFieldSample fields;
    fields.edge_field.resize(layout.dimension());fields.edge_field=ComplexType{};
    const auto eta=layout.flat(c::StaggeredSector::Eta,0,0);
    fields.edge_field[eta]=ComplexType{0.2,0.1};
    fields.edge_field[layout.flat(c::StaggeredSector::Nu,0,1)]=ComplexType{-0.1,0.3};
    fields.edge_field[layout.flat(c::StaggeredSector::Nu,1,0)]=ComplexType{0.4,-0.2};
    const auto tr=f::compute_contour_trajectory(p,fields,f::MeanFieldTrajectory(2,FieldVector{}));
    const Operator H=ComplexType{0.2,0.1}*S_X+RealType{1.3}*S_Z+RealType{0.4/0.3}*S_Y;
    const Operator left=f::general_matrix_exponential(ComplexType{0.,-0.3/4}*ComplexType{-0.1,0.3}*S_Y);
    const Operator right=f::general_matrix_exponential(ComplexType{0.,-0.3/4}*ComplexType{0.4,-0.2}*S_X);
    const Operator R=f::general_matrix_exponential(ComplexType{0.,0.7}*S_X);
    const Operator expected=right*f::general_matrix_exponential(ComplexType{0.,-0.2}*H)
        *R*f::general_matrix_exponential(ComplexType{0.,-0.1}*H)*left;
    const Operator reverse=left*f::general_matrix_exponential(ComplexType{0.,0.1}*H)*blaze::ctrans(R)
        *f::general_matrix_exponential(ComplexType{0.,0.2}*H)*right;
    require(difference(tr.forward_steps[1],expected)<2e-13,"interior pulse splits precession containing the measurement-edge pulse field");
    require(difference(tr.backward_steps[1],reverse)<2e-13,"nonunitary backward branch keeps complex inverse precession and adjoint pulses");
}

void check_thermal_anchor()
{
    auto p=parameters(0.5);p.pulses={{0.,FieldVector{0.,1.,0.},std::acos(-1.)/2}};
    f::initialize_matrices(p);
    c::CorrelationSet prepared{'D',p.num_TimePoints,1};
    c::Magnetization pm{'D',1};pm.Re.front()[2]=-0.2;
    for(size_t pair=0;pair<prepared.Re.front().size();++pair)
        for(size_t tau=0;tau<p.num_TimePoints;++tau)
            prepared.Re.front()[pair][tau]=0.1+0.01*pair+0.002*tau;
    const auto before=f::self_consistent_equations(p,prepared,pm);
    c::CorrelationSet q;c::Magnetization m;
    f::seed_from_preparation(p,prepared,pm,q,m);
    require(std::abs(m.Re.front()[0]-0.2)<2e-13&&std::abs(m.Re.front()[2])<2e-13,
        "frozen real origin is post-pulse while preparation remains pre-pulse");
    const f::ThermalReference thermal{prepared,pm};
    const auto field=f::self_consistent_equations(p,q,m,&thermal);
    for(size_t i=0;i<3*p.num_TimePoints;++i)for(size_t j=0;j<3*p.num_TimePoints;++j)
        require(field.covariance(i,j)==before.covariance(i,j),"origin pulse leaves the entire frozen Matsubara covariance unchanged");
    for(size_t a=0;a<3;++a)require(field.imaginary_mean[a]==before.mean_time.front()[a],
        "imaginary propagation retains the unrotated thermal mean");
    require(std::abs(field.mean_time.front()[0])>0.01&&std::abs(field.imaginary_mean[0])<1e-14,
        "real and imaginary deterministic means differ after an origin pulse");
    const auto connected=f::connected_contour_primitive(q,m,&pm);
    const auto full=c::edge_value(q,0,0,2,0);
    require(std::abs(c::edge_value(connected,0,0,2,0)-(full-m.Re.front()[0]*pm.Re.front()[2]))<2e-13,
        "mixed disconnected term uses pre-pulse imaginary magnetization");
    const c::StaggeredLayout grid{p.num_TimePoints,p.num_RealTimeSteps};
    f::JointComplexGaussianSampler::ContourFieldSample fields;
    fields.edge_field.resize(grid.dimension());fields.edge_field=ComplexType{};
    f::ContourTrajectory tr,pre;f::TrajectoryWorkspace workspace,pre_workspace;
    const f::PreparedMeanField cache(p,field.mean_time,field.imaginary_mean);
    f::build_contour_trajectory(p,fields,cache,tr,workspace);
    f::JointComplexGaussianSampler::ContourFieldSample imaginary;
    imaginary.edge_field.resize(3*p.num_TimePoints+6);imaginary.edge_field=ComplexType{};
    f::build_matsubara_trajectory(p,imaginary,before.mean_time,pre,pre_workspace);
    require(difference(tr.imaginary_density_operator,pre.imaginary_density_operator)<2e-13,
        "real-time completion and pCN preparation use the same pre-pulse imaginary propagation");
}

int main()
{
    for(RealType spin:{RealType{0.5},RealType{1.}})
    {
        check_rotation(spin);
        for(size_t q:{size_t{1},size_t{2},size_t{3}})
        {check_times(spin,q);check_measurement_fields(spin,q);}
    }
    check_response();check_thermal_anchor();
    return failures?1:0;
}
