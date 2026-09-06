#include"../Functions/Functions.h"

#include<cmath>
#include<iostream>
#include<limits>
#include<Physics/CFET.h>

namespace func = spinDMFT::Functions;
namespace contour = spinDMFT::Contour;
namespace ps = spinDMFT::Parameter_Space;

namespace
{
int require( bool condition,const char* message )
{
    if( condition ) return 0;
    std::cerr<<"FAILED: "<<message<<'\n'; return 1;
}
RealType frobenius( const Operator& value )
{
    RealType result{};
    for( size_t i=0;i<value.rows();++i ) for( size_t j=0;j<value.columns();++j )
        result+=std::norm(value(i,j));
    return std::sqrt(result);
}

std::pair<Operator,Operator> cumulative_branch_propagators(
    const func::ContourTrajectory& trajectory, const size_t last )
{
    Operator forward=IDENTITY,backward=IDENTITY;
    for( size_t t=1;t<=last;++t )
    {
        forward=trajectory.forward_steps[t]*forward;
        backward=backward*trajectory.backward_steps[t];
    }
    return {forward,backward};
}

Operator evolved_density( const func::ContourTrajectory& trajectory,
                          const Operator& initial, const size_t last )
{
    Operator result=initial;
    for( size_t t=1;t<=last;++t )
        result=trajectory.forward_steps[t]*result*trajectory.backward_steps[t];
    return result;
}
}

int main()
{
    int failures{};
    ps::ParameterSpace p;
    p.spin_float=RealType{0.5}; p.num_HilbertSpaceDimension=2;
    p.spin_model=Physics::SpinModel{"ISO"};
    p.B=Physics::MagneticField{"z",RealType{0.7},RealType{0.},RealType{0.}};
    p.noise=Physics::Noise{"none",RealType{0.}};
    p.extra_interaction=Physics::ExtraInteraction{"none",RealType{0.}};
    p.beta=RealType{0.8}; p.num_TimeSteps=4; p.num_TimePoints=5;
    p.delta_t=p.beta/static_cast<RealType>(p.num_TimeSteps);
    p.num_RealTimeSteps=4; p.num_RealTimePoints=5; p.Tmax=RealType{0.6};
    p.delta_real_t=p.Tmax/static_cast<RealType>(p.num_RealTimeSteps);
    func::initialize_matrices(p);

    Operator diagonal(2,2,ComplexType{});
    diagonal(0,0)=ComplexType{1.,0.3}; diagonal(1,1)=ComplexType{-0.2,-0.4};
    const Operator exp_diagonal=func::general_matrix_exponential(diagonal);
    failures+=require(std::abs(exp_diagonal(0,0)-std::exp(diagonal(0,0)))<RealType{1e-12}
                     &&std::abs(exp_diagonal(1,1)-std::exp(diagonal(1,1)))<RealType{1e-12},
                      "general complex matrix exponential");

    const func::ComplexFieldVector complex_field{
        ComplexType{0.31,-0.17},ComplexType{-0.23,0.29},ComplexType{0.41,0.13}};
    const ComplexType scalar{0.19,-0.07};
    const ComplexType complex_step{-0.11,0.08};
    const Operator complex_hamiltonian=scalar*IDENTITY
        +complex_field[0]*S_X+complex_field[1]*S_Y+complex_field[2]*S_Z;
    const Operator direct_spin_half=func::spin_half_field_exponential(
        complex_field,scalar,complex_step);
    const Operator general_spin_half=func::general_matrix_exponential(
        complex_step*complex_hamiltonian);
    failures+=require(frobenius(direct_spin_half-general_spin_half)<RealType{2e-13},
                      "direct spin-half exponential matches general complex exponential");

    ps::ParameterSpace quadrupolar=p;
    const RealType quadrupolar_strength=RealType{0.37};
    quadrupolar.extra_interaction=Physics::ExtraInteraction{
        "quadrupolar",quadrupolar_strength};
    func::initialize_matrices(quadrupolar);
    const contour::ContourLayout quadrupolar_layout{
        quadrupolar.num_TimePoints,quadrupolar.num_RealTimePoints};
    func::DenseComplexGaussianSampler::FieldVector quadrupolar_zero_field(
        quadrupolar_layout.dimension(),ComplexType{});
    const func::MeanFieldTrajectory quadrupolar_mean(
        quadrupolar.num_RealTimePoints,FieldVector{});
    const auto quadrupolar_trajectory=func::compute_contour_trajectory(
        quadrupolar,quadrupolar_zero_field,quadrupolar_mean);
    const RealType expected_quadrupolar_Z=RealType{2.}
        *std::exp(-quadrupolar.beta*RealType{0.75}*quadrupolar_strength)
        *std::cosh(quadrupolar.beta*quadrupolar.B.m_h[2]/RealType{2.});
    failures+=require(std::abs(quadrupolar_trajectory.partition_function
                               -expected_quadrupolar_Z)<RealType{2e-12},
                      "spin-half direct exponential includes scalar quadrupolar term");
    func::initialize_matrices(p);

    const contour::ContourLayout layout{p.num_TimePoints,p.num_RealTimePoints};
    func::DenseComplexGaussianSampler::FieldVector zero_field(layout.dimension(),ComplexType{});
    const func::MeanFieldTrajectory mean_time(p.num_RealTimePoints,FieldVector{});
    const auto trajectory=func::compute_contour_trajectory(p,zero_field,mean_time);
    const RealType hz=p.B.m_h[2];
    const RealType expected_Z=RealType{2.}*std::cosh(p.beta*hz/RealType{2.});
    failures+=require(std::abs(trajectory.partition_function-expected_Z)<RealType{1e-12},
                      "uncoupled-spin partition function");
    const ComplexType mz_num=blaze::trace(trajectory.imaginary_density_operator*S_Z);
    const RealType expected_mz=-RealType{0.5}*std::tanh(p.beta*hz/RealType{2.});
    failures+=require(std::abs(mz_num/trajectory.partition_function-expected_mz)<RealType{1e-12},
                      "uncoupled-spin thermal magnetization sign");
    for( size_t t=0;t<p.num_RealTimePoints;++t )
        failures+=require(std::abs(blaze::trace(evolved_density(
                                      trajectory,trajectory.imaginary_density_operator,t))
                                  -trajectory.partition_function)<RealType{1e-11},
                          "equal forward/backward fields close the contour");

    auto varying_mean_time=mean_time;
    for( size_t t=0;t<varying_mean_time.size();++t )
    {
        varying_mean_time[t][0]=RealType{0.11}*static_cast<RealType>(t);
        varying_mean_time[t][2]=RealType{-0.07}*static_cast<RealType>(t*t);
    }
    const auto varying_mean_trajectory=func::compute_contour_trajectory(
        p,zero_field,varying_mean_time);
    for( size_t t=0;t<p.num_RealTimePoints;++t )
        failures+=require(std::abs(blaze::trace(evolved_density(
                                      varying_mean_trajectory,
                                      varying_mean_trajectory.imaginary_density_operator,t))
                                  -varying_mean_trajectory.partition_function)
                                  <RealType{1e-11},
                          "common time-dependent mean closes the real contour");

    const Operator H_old=p.B.m_h[2]*S_Z;
    const Operator H_new=varying_mean_time[1][0]*S_X
                        +(p.B.m_h[2]+varying_mean_time[1][2])*S_Z;
    std::array<Operator,3> expected_cfet4_exponentials{};
    std::transform(Physics::CFET::BETA_LIST.cbegin(),
                   Physics::CFET::BETA_LIST.cend(),
                   expected_cfet4_exponentials.begin(),
                   [&]( const Physics::CFET::BetaComponent& beta )
    {
        return func::general_matrix_exponential(
            ComplexType{RealType{0.},-p.delta_real_t}
            *(beta[0]*H_new+beta[1]*H_old));
    });
    const Operator expected_cfet4=expected_cfet4_exponentials[0]
                                 *expected_cfet4_exponentials[1]
                                 *expected_cfet4_exponentials[2];
    failures+=require(frobenius(
        varying_mean_trajectory.forward_steps[1]-expected_cfet4)
        <RealType{1e-12},"real-time step uses the spinDMFT CFET4 composition");
    const Operator cfet2_step=func::general_matrix_exponential(
        ComplexType{RealType{0.},-p.delta_real_t}
        *RealType{0.5}*(H_new+H_old));
    failures+=require(frobenius(expected_cfet4-cfet2_step)>RealType{1e-8},
                      "time-dependent noncommuting CFET4 step differs from CFET2");

    const size_t t_index=3;
    const auto [forward_t,backward_t]=cumulative_branch_propagators(
        trajectory,t_index);
    const Operator measured_x=backward_t*S_X*forward_t;
    const ComplexType gxx=blaze::trace(
        trajectory.imaginary_edge_insertions[0].back()*measured_x)
        /trajectory.partition_function;
    const RealType t=t_index*p.delta_real_t;
    const ComplexType expected_gxx{
        RealType{0.25}*std::cos(hz*t),expected_mz/RealType{2.}*std::sin(hz*t)};
    failures+=require(std::abs(gxx-expected_gxx)<RealType{1e-11},
                      "uncoupled transverse greater correlation");
    const Operator measured_z=backward_t*S_Z*forward_t;
    const ComplexType gzz=blaze::trace(
        trajectory.imaginary_edge_insertions[2].back()*measured_z)
        /trajectory.partition_function;
    failures+=require(std::abs(gzz-RealType{0.25})<RealType{1e-12},
                      "uncoupled longitudinal correlation");

    auto imaginary_edge_field=zero_field;
    RealType trapezoidal_edge_field_sum{};
    for( size_t k=0;k<p.num_TimePoints;++k )
    {
        const RealType value=RealType{0.1}*static_cast<RealType>(k+1);
        imaginary_edge_field[layout.flat(contour::Branch::Matsubara,k,2)]=value;
        trapezoidal_edge_field_sum+=(k==0||k==p.num_TimeSteps)
            ?RealType{0.5}*value:value;
    }
    const auto edge_trajectory=func::compute_contour_trajectory(
        p,imaginary_edge_field,mean_time);
    const RealType average_edge_field=trapezoidal_edge_field_sum/p.num_TimeSteps;
    const RealType expected_edge_Z=RealType{2.}*std::cosh(
        p.beta*(hz+average_edge_field)/RealType{2.});
    failures+=require(std::abs(edge_trajectory.partition_function-expected_edge_Z)
                              <RealType{1e-12},
                      "CFET4 uses distinct one-sided fields in the final imaginary interval");

    auto distinct_field=zero_field;
    for( size_t n=0;n<p.num_RealTimePoints;++n )
    {
        distinct_field[layout.flat(contour::Branch::Forward,n,0)]=ComplexType{0.3,0.1};
        distinct_field[layout.flat(contour::Branch::Backward,n,0)]=ComplexType{-0.2,0.05};
    }
    const auto distinct=func::compute_contour_trajectory(p,distinct_field,mean_time);
    const Operator closure=evolved_density(
        distinct,IDENTITY,p.num_RealTimePoints-1)-IDENTITY;
    failures+=require(frobenius(closure)>RealType{1e-4},
                      "distinct backward field is not replaced by the forward inverse");

    contour::CorrelationSet full{'D',4,3};
    func::ComplexMagnetizationTrajectory magnetization(3,func::ComplexFieldVector{});
    magnetization[0]={ComplexType{0.2,0.1},ComplexType{-0.1,0.05},
                      ComplexType{0.3,-0.04}};
    magnetization[1]={ComplexType{0.25,-0.03},ComplexType{-0.08,0.02},
                      ComplexType{0.24,0.07}};
    magnetization[2]={ComplexType{0.15,0.06},ComplexType{-0.04,-0.01},
                      ComplexType{0.18,0.09}};
    for( size_t ti=0;ti<full.Re.size();++ti )
        for( size_t pair=0;pair<full.Re[ti].size();++pair )
        {
            const auto ab=full.Re[ti].get_direction_pair(pair);
            const ComplexType disconnected=magnetization[ti][ab[0]]
                                          *magnetization[0][ab[1]];
            for( size_t tau=0;tau<full.Re[ti][pair].size();++tau )
            {
                const ComplexType target{
                    RealType{0.01}*(RealType{1.}+ti+pair+tau),
                    RealType{-0.005}*(RealType{1.}+RealType{2.}*ti+pair+tau)};
                const ComplexType value=target+disconnected;
                full.Re[ti][pair][tau]=std::real(value);
                full.Im[ti][pair][tau]=std::imag(value);
            }
        }
    const auto connected=func::connected_contour_primitive(full,magnetization);
    for( size_t ti=0;ti<connected.Re.size();++ti )
        for( size_t pair=0;pair<connected.Re[ti].size();++pair )
            for( size_t tau=0;tau<connected.Re[ti][pair].size();++tau )
            {
                const ComplexType expected{
                    RealType{0.01}*(RealType{1.}+ti+pair+tau),
                    RealType{-0.005}*(RealType{1.}+RealType{2.}*ti+pair+tau)};
                const ComplexType actual{connected.Re[ti][pair][tau],
                                         connected.Im[ti][pair][tau]};
                failures+=require(std::abs(actual-expected)<RealType{1e-13},
                    "connected primitive subtracts unconjugated m_a(t)m_b(0)");
            }

    const auto stationary=func::project_constant_magnetization(magnetization);
    failures+=require(stationary.size()==magnetization.size(),
                      "constant magnetization projection preserves the time grid");
    for( const auto& value:stationary )
        for( size_t c=0;c<3;++c )
            failures+=require(std::abs(value[c]-stationary.front()[c])<RealType{1e-13},
                "constant magnetization projection uses m(0) at every time");
    const auto stationary_connected=func::connected_contour_primitive(full,stationary);
    for( size_t pair=0;pair<stationary_connected.Re[0].size();++pair )
    {
        const auto ab=stationary_connected.Re[0].get_direction_pair(pair);
        const ComplexType original{full.Re[0][pair][0],full.Im[0][pair][0]};
        const ComplexType expected=original-stationary.front()[ab[0]]
                                           *stationary.front()[ab[1]];
        const ComplexType actual{stationary_connected.Re[0][pair][0],
                                 stationary_connected.Im[0][pair][0]};
        failures+=require(std::abs(actual-expected)<RealType{1e-13},
                          "stationary connected primitive reduces to m_a m_b subtraction");
    }

    ps::ParameterSpace mean_pspace=p;
    mean_pspace.num_TimeSteps=3;
    mean_pspace.num_TimePoints=4;
    mean_pspace.num_RealTimeSteps=2;
    mean_pspace.num_RealTimePoints=3;
    mean_pspace.JQ=RealType{0.};
    mean_pspace.JL=RealType{-2.};
    const auto self_consistent=func::self_consistent_equations(
        mean_pspace,full,magnetization);
    for( size_t ti=0;ti<magnetization.size();++ti )
        for( size_t c=0;c<3;++c )
            failures+=require(std::abs(self_consistent.mean_time[ti][c]
                -mean_pspace.JL*std::real(magnetization[ti][c]))<RealType{1e-13},
                "first-moment closure uses JL D Re m(t)");

    ps::ParameterSpace covariance_pspace=mean_pspace;
    covariance_pspace.JQ=RealType{1.};
    covariance_pspace.JL=RealType{0.};
    const auto endpoint_covariance=func::self_consistent_equations(
        covariance_pspace,full,magnetization);
    const contour::ContourLayout endpoint_layout{
        covariance_pspace.num_TimePoints,covariance_pspace.num_RealTimePoints};
    const size_t zero_row=endpoint_layout.flat(contour::Branch::Matsubara,0,0);
    const size_t beta_row=endpoint_layout.flat(
        contour::Branch::Matsubara,covariance_pspace.num_TimeSteps,0);
    const size_t real_column=endpoint_layout.flat(contour::Branch::Forward,1,1);
    failures+=require(std::abs(endpoint_covariance.covariance(
                                  zero_row,real_column)
                              -endpoint_covariance.covariance(
                                  beta_row,real_column))>RealType{1e-6},
        "0+ and beta- retain distinct lesser and greater mixed covariances");

    ps::ParameterSpace bath_pspace;
    bath_pspace.bath="harmonic";
    bath_pspace.bath_frequency=RealType{0.7};
    bath_pspace.bath_coupling=RealType{0.3};
    bath_pspace.bath_component='x';
    bath_pspace.correlation_symmetry_type='D';
    bath_pspace.beta=RealType{2.};
    bath_pspace.num_TimeSteps=4;
    bath_pspace.num_TimePoints=5;
    bath_pspace.delta_t=bath_pspace.beta/bath_pspace.num_TimeSteps;
    bath_pspace.num_RealTimeSteps=3;
    bath_pspace.num_RealTimePoints=4;
    bath_pspace.Tmax=RealType{0.6};
    bath_pspace.delta_real_t=bath_pspace.Tmax/bath_pspace.num_RealTimeSteps;
    const auto bath_primitive=func::harmonic_bath_primitive(bath_pspace);
    const RealType q=std::exp(-bath_pspace.beta*bath_pspace.bath_frequency);
    const RealType n= q/(RealType{1.}-q);
    for( size_t ti=0;ti<bath_pspace.num_RealTimePoints;++ti )
    {
        const RealType bath_t=ti*bath_pspace.delta_real_t;
        const ComplexType expected_greater=bath_pspace.bath_coupling
            *bath_pspace.bath_coupling
            *((n+RealType{1.})*std::exp(ComplexType{0.,-bath_pspace.bath_frequency*bath_t})
              +n*std::exp(ComplexType{0.,bath_pspace.bath_frequency*bath_t}));
        const ComplexType expected_lesser=bath_pspace.bath_coupling
            *bath_pspace.bath_coupling
            *(n*std::exp(ComplexType{0.,-bath_pspace.bath_frequency*bath_t})
              +(n+RealType{1.})*std::exp(ComplexType{0.,bath_pspace.bath_frequency*bath_t}));
        const ComplexType greater=contour::edge_value(
            bath_primitive,ti,0,0,bath_pspace.num_TimeSteps);
        const ComplexType lesser=contour::edge_value(bath_primitive,ti,0,0,0);
        failures+=require(std::abs(greater-expected_greater)<RealType{1e-13},
            "harmonic bath beta endpoint is the analytic greater function");
        failures+=require(std::abs(lesser-expected_lesser)<RealType{1e-13},
            "harmonic bath zero endpoint is the analytic lesser function");
        failures+=require(std::abs(std::conj(greater)-lesser)<RealType{1e-13},
            "harmonic bath endpoints obey KMS Hermiticity");
        failures+=require(std::abs(contour::edge_value(
            bath_primitive,ti,1,1,0))<RealType{1e-15},
            "single-component harmonic bath leaves orthogonal fields zero");
    }
    const auto bath_field=func::prescribed_harmonic_bath_field(bath_pspace);
    const contour::ContourLayout bath_layout{
        bath_pspace.num_TimePoints,bath_pspace.num_RealTimePoints};
    failures+=require(bath_field.covariance.rows()==bath_layout.dimension(),
        "harmonic bath covariance uses the complete M,+,- layout");
    failures+=require(bath_field.covariance_symmetry_error<RealType{1e-14}
                     &&bath_field.branch_identity_error<RealType{1e-14},
        "harmonic bath covariance is complex symmetric before and after mirroring");
    for( const auto& mean:bath_field.mean_time )
        failures+=require(blaze::length(mean)<RealType{1e-15},
            "prescribed harmonic bath has zero stochastic-field mean");
    const size_t plus_late=bath_layout.flat(contour::Branch::Forward,1,0);
    const size_t plus_zero=bath_layout.flat(contour::Branch::Forward,0,0);
    const size_t minus_late=bath_layout.flat(contour::Branch::Backward,1,0);
    const size_t minus_zero=bath_layout.flat(contour::Branch::Backward,0,0);
    const ComplexType bath_greater=contour::edge_value(
        bath_primitive,1,0,0,bath_pspace.num_TimeSteps);
    const ComplexType bath_lesser=contour::edge_value(bath_primitive,1,0,0,0);
    failures+=require(std::abs(bath_field.covariance(plus_late,minus_zero)
                              -bath_lesser)<RealType{1e-13},
        "harmonic bath +- block is the lesser function");
    failures+=require(std::abs(bath_field.covariance(minus_late,plus_zero)
                              -bath_greater)<RealType{1e-13},
        "harmonic bath -+ block is the greater function");
    failures+=require(std::abs(bath_greater-bath_lesser)>RealType{1e-6},
        "harmonic bath exercises a nonzero response covariance");

    contour::CorrelationSet old_iteration{'A',2,1};
    contour::CorrelationSet raw_iteration{'A',2,1};
    contour::CorrelationSet iteration_errors{'A',2,1};
    raw_iteration.Re[0][0][0]=RealType{0.2};
    iteration_errors.Re[0][0][0]=RealType{0.1};
    raw_iteration.Im[0][0][1]=RealType{-0.3};
    iteration_errors.Im[0][0][1]=RealType{0.1};
    func::ComplexMagnetizationTrajectory old_mag(1),raw_mag(1);
    std::vector<FieldVector> mag_Re_errors(1),mag_Im_errors(1);
    raw_mag[0][0]=ComplexType{RealType{0.3},RealType{0.4}};
    mag_Re_errors[0][0]=RealType{0.06};
    mag_Im_errors[0][0]=RealType{0.2};
    const auto residual=func::iteration_residual(
        old_iteration,raw_iteration,iteration_errors,old_mag,raw_mag,
        mag_Re_errors,mag_Im_errors);
    failures+=require(std::abs(residual.absolute-RealType{0.5})<RealType{1e-14},
        "iteration residual retains the largest absolute complex difference");
    failures+=require(std::abs(residual.standardized-RealType{5.})<RealType{1e-14},
        "iteration residual standardizes real and imaginary magnetization separately");
    iteration_errors.Re[0][0][0]=RealType{};
    const auto unresolved=func::iteration_residual(
        old_iteration,raw_iteration,iteration_errors,old_mag,raw_mag,
        mag_Re_errors,mag_Im_errors);
    failures+=require(std::isinf(unresolved.standardized),
        "a nonzero residual with zero statistical error cannot converge");
    contour::CorrelationSet exact_old{'A',2,1};
    contour::CorrelationSet exact_raw{'A',2,1};
    contour::CorrelationSet exact_errors{'A',2,1};
    exact_old.Re[0][0][0]=RealType{0.25};
    exact_raw.Re[0][0][0]=std::nextafter(
        RealType{0.25},RealType{1.});
    exact_old.Re[0][0][1]=RealType{0.25};
    exact_raw.Re[0][0][1]=std::nextafter(
        RealType{0.25},RealType{});
    const auto endpoint_roundoff=func::iteration_residual(
        exact_old,exact_raw,exact_errors,old_mag,old_mag,
        mag_Re_errors,mag_Im_errors);
    failures+=require(endpoint_roundoff.standardized==RealType{},
        "roundoff at zero-variance t=0 tau endpoints is resolved");
    exact_raw.Re[0][0][1]=RealType{0.24};
    const auto endpoint_mismatch=func::iteration_residual(
        exact_old,exact_raw,exact_errors,old_mag,old_mag,
        mag_Re_errors,mag_Im_errors);
    failures+=require(std::isinf(endpoint_mismatch.standardized),
        "a physical zero-variance t=0 tau endpoint mismatch remains unresolved");
    iteration_errors.Re[0][0][0]=std::numeric_limits<RealType>::quiet_NaN();
    const auto nonfinite=func::iteration_residual(
        old_iteration,raw_iteration,iteration_errors,old_mag,raw_mag,
        mag_Re_errors,mag_Im_errors);
    failures+=require(std::isinf(nonfinite.absolute)
                      &&std::isinf(nonfinite.standardized),
        "nonfinite iteration statistics cannot converge");
    return failures==0?0:1;
}
