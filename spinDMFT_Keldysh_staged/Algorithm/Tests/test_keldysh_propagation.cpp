#include"../Functions/Functions.h"

#include<cmath>
#include<iostream>
#include<limits>

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

func::ComplexFieldVector smooth_noncommuting_field( const RealType time )
{
    return {
        ComplexType{RealType{0.35}*std::cos(RealType{0.7}*time),0.},
        ComplexType{RealType{0.22}*std::sin(RealType{0.4}*time),0.},
        ComplexType{RealType{-0.18}*std::cos(RealType{1.1}*time),0.}};
}

}

// Integrate actual trajectories over three Gaussian latent coordinates using
// deterministic Gauss--Hermite quadrature. For commuting z fields the exact
// ensemble density follows from E exp(q.V)=exp(q.Gamma.q/2). This checks the
// nonlinear thermal + real-time propagation without Monte-Carlo tolerances.
int check_weighted_gaussian_propagation(ps::ParameterSpace p)
{
    int failures{};
    constexpr size_t order=10;
    using RealMatrix=blaze::DynamicMatrix<RealType,blaze::rowMajor>;
    blaze::SymmetricMatrix<RealMatrix> jacobi(order);
    for(size_t i=1;i<order;++i)jacobi(i-1,i)=std::sqrt(static_cast<RealType>(i));
    blaze::DynamicVector<RealType> nodes(order);
    RealMatrix vectors(order,order);
    blaze::eigen(jacobi,nodes,vectors);
    std::array<RealType,order> weights{};
    for(size_t i=0;i<order;++i)weights[i]=vectors(i,0)*vectors(i,0);

    // A synthetic rank-three joint complex covariance, including unequal
    // mixed blocks. No bath approximation enters the analytic reference.
    func::ComplexDynamicMatrix K{{{0.20,0.03},{0.03,0.02},{0.02,-0.01}},
                                  {{0.03,0.02},{0.14,0.05},{0.08,0.04}},
                                  {{0.02,-0.01},{0.08,0.04},{0.10,-0.02}}};
    const contour::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps};
    func::ComplexDynamicMatrix gamma(layout.dimension(),layout.dimension(),ComplexType{});
    for(size_t i=0;i<layout.dimension();++i)for(size_t j=0;j<layout.dimension();++j)
    {
        const auto a=layout.decode(i),b=layout.decode(j);
        if(a.component==2&&b.component==2)
            gamma(i,j)=K(static_cast<size_t>(a.sector),static_cast<size_t>(b.sector));
    }
    const func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    for(const std::array<RealType,3> noise_weights:
        {std::array<RealType,3>{1.,2.,2.},std::array<RealType,3>{1.,2.,8.},std::array<RealType,3>{1.,8.,2.}})
    {
        std::vector<RealType> diagonal_weights(layout.dimension());
        for(size_t i=0;i<diagonal_weights.size();++i)
        {
            const auto sector=layout.decode(i).sector;
            diagonal_weights[i]=sector==contour::StaggeredSector::Matsubara?noise_weights[0]
                :sector==contour::StaggeredSector::Eta?noise_weights[1]:noise_weights[2]/4;
        }
        func::DenseComplexGaussianSampler sampler(gamma,diagonal_weights);
        failures+=require(sampler.latent_dimension()==3,"commuting Gaussian reference has three latent modes");
        if(sampler.latent_dimension()!=3)continue;
        std::vector<Operator> average(p.num_RealTimePoints,Operator(2,2,ComplexType{}));
        for(size_t i=0;i<order;++i)for(size_t j=0;j<order;++j)for(size_t k=0;k<order;++k)
        {
            func::JointComplexGaussianSampler::LatentVector r{nodes[i],nodes[j],nodes[k]};
            const auto tr=func::compute_contour_trajectory(p,sampler.contour_field_from_latent(r),mean);
            Operator density=tr.imaginary_density_operator;
            const RealType weight=weights[i]*weights[j]*weights[k];
            for(size_t t=0;t<p.num_RealTimePoints;++t)
            {
                if(t>0)density=tr.forward_steps[t]*density*tr.backward_steps[t];
                average[t]+=weight*density;
            }
        }
        for(size_t t=0;t<p.num_RealTimePoints;++t)
        {
            const RealType time=t*p.delta_real_t;
            Operator exact(2,2,ComplexType{});
            for(size_t diagonal=0;diagonal<2;++diagonal)
            {
                const RealType spin=std::real(ComplexType(S_Z(diagonal,diagonal)));
                const std::array<ComplexType,3> q{
                    ComplexType{-p.beta*spin,0.},ComplexType{},ComplexType{0.,-time*spin}};
                ComplexType exponent=-p.beta*spin*p.B.m_h[2];
                for(size_t a=0;a<3;++a)for(size_t b=0;b<3;++b)exponent+=RealType{0.5}*q[a]*K(a,b)*q[b];
                exact(diagonal,diagonal)=std::exp(exponent);
            }
            failures+=require(frobenius(average[t]-exact)<RealType{2e-11},
                "weighted thermal and real-time propagation agrees with analytic Gaussian average");
        }
    }
    return failures;
}

int main()
{
    int failures{};
    ps::ParameterSpace p;
    p.real_time_substeps=1;
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
    failures+=check_weighted_gaussian_propagation(p);

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
    const contour::StaggeredLayout quadrupolar_layout{
        quadrupolar.num_TimePoints,quadrupolar.num_RealTimeSteps};
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

    const contour::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps};
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

    auto cubic_field=zero_field;
    for(size_t k=0;k<p.num_TimePoints;++k)
    {
        const RealType tau=k*p.delta_t;
        cubic_field[layout.flat(contour::StaggeredSector::Matsubara,k,2)]=
            0.1+0.2*tau-0.3*tau*tau+0.4*tau*tau*tau;
    }
    const auto cubic_trajectory=func::compute_contour_trajectory(p,cubic_field,mean_time);
    const RealType integrated=0.1*p.beta+0.1*p.beta*p.beta
        -0.1*std::pow(p.beta,3)+0.1*std::pow(p.beta,4);
    failures+=require(std::abs(cubic_trajectory.partition_function
        -2*std::cosh((p.beta*hz+integrated)/2))<2e-12,
        "Matsubara preparation integrates a cubic commuting field exactly");

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
        imaginary_edge_field[layout.flat(contour::StaggeredSector::Matsubara,k,2)]=value;
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
                      "Matsubara preparation uses distinct one-sided fields in the final imaginary interval");

    auto distinct_field=zero_field;
    for(size_t n=0;n<p.num_RealTimeSteps;++n)
        distinct_field[layout.flat(contour::StaggeredSector::Eta,n,0)]=ComplexType{0.05,0.075};
    for(size_t n=0;n<=p.num_RealTimeSteps;++n)
        distinct_field[layout.flat(contour::StaggeredSector::Nu,n,0)]=ComplexType{0.5,0.05};
    const auto distinct=func::compute_contour_trajectory(p,distinct_field,mean_time);
    const Operator closure=evolved_density(
        distinct,IDENTITY,p.num_RealTimePoints-1)-IDENTITY;
    failures+=require(frobenius(closure)>RealType{1e-4},
                      "distinct backward field is not replaced by the forward inverse");

    contour::CorrelationSet full{'D',4,3};
    func::MagTen magnetization_Re{'D',3},magnetization_Im{'D',3};
    const std::array<FieldVector,3> magnetization_Re_full{
        FieldVector{0.2,-0.1,0.3},FieldVector{0.25,-0.08,0.24},
        FieldVector{0.15,-0.04,0.18}};
    const std::array<FieldVector,3> magnetization_Im_full{
        FieldVector{0.1,0.05,-0.04},FieldVector{-0.03,0.02,0.07},
        FieldVector{0.06,-0.01,0.09}};
    for( size_t t=0;t<3;++t )
    {
        magnetization_Re[t]=func::MagVec{'D',magnetization_Re_full[t]};
        magnetization_Im[t]=func::MagVec{'D',magnetization_Im_full[t]};
    }
    for( size_t ti=0;ti<full.Re.size();++ti )
        for( size_t pair=0;pair<full.Re[ti].size();++pair )
        {
            const auto ab=full.Re[ti].get_direction_pair(pair);
            const ComplexType disconnected=
                ComplexType{magnetization_Re_full[ti][ab[0]],
                            magnetization_Im_full[ti][ab[0]]}
               *ComplexType{magnetization_Re_full[0][ab[1]],
                            magnetization_Im_full[0][ab[1]]};
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
    const auto connected=func::connected_contour_primitive(full, spinDMFT::Contour::Magnetization{magnetization_Re,magnetization_Im});
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

    const auto stationary_Re=func::project_constant_magnetization(magnetization_Re);
    const auto stationary_Im=func::project_constant_magnetization(magnetization_Im);
    failures+=require(stationary_Re.size()==magnetization_Re.size(),
                      "constant magnetization projection preserves the time grid");
    for( size_t t=0;t<stationary_Re.size();++t )
        for( size_t c=0;c<stationary_Re.num_components();++c )
            failures+=require(
                std::abs(stationary_Re[t][c]-stationary_Re.front()[c])<RealType{1e-13}
                &&std::abs(stationary_Im[t][c]-stationary_Im.front()[c])<RealType{1e-13},
                "constant magnetization projection uses m(0) at every time");
    const auto stationary_connected=func::connected_contour_primitive(full, spinDMFT::Contour::Magnetization{stationary_Re,stationary_Im});
    for( size_t pair=0;pair<stationary_connected.Re[0].size();++pair )
    {
        const auto ab=stationary_connected.Re[0].get_direction_pair(pair);
        const ComplexType original{full.Re[0][pair][0],full.Im[0][pair][0]};
        const ComplexType stationary_a{stationary_Re.front()[ab[0]],
                                       stationary_Im.front()[ab[0]]};
        const ComplexType stationary_b{stationary_Re.front()[ab[1]],
                                       stationary_Im.front()[ab[1]]};
        const ComplexType expected=original-stationary_a*stationary_b;
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
        mean_pspace,
        full,
        spinDMFT::Contour::Magnetization{magnetization_Re,magnetization_Im});
    for( size_t ti=0;ti<magnetization_Re.size();++ti )
        for( size_t c=0;c<3;++c )
            failures+=require(std::abs(self_consistent.mean_time[ti][c]
                -mean_pspace.JL*magnetization_Re[ti][c])<RealType{1e-13},
                "first-moment closure uses JL D Re m(t)");

    ps::ParameterSpace covariance_pspace=mean_pspace;
    covariance_pspace.JQ=RealType{1.};
    covariance_pspace.JL=RealType{0.};
    const auto endpoint_covariance=func::self_consistent_equations(
        covariance_pspace,
        full,
        spinDMFT::Contour::Magnetization{magnetization_Re,magnetization_Im});
    const contour::StaggeredLayout endpoint_layout{
        covariance_pspace.num_TimePoints,covariance_pspace.num_RealTimeSteps};
    const size_t zero_row=endpoint_layout.flat(contour::StaggeredSector::Matsubara,0,0);
    const size_t beta_row=endpoint_layout.flat(
        contour::StaggeredSector::Matsubara,covariance_pspace.num_TimeSteps,0);
    const size_t real_column=endpoint_layout.flat(contour::StaggeredSector::Eta,1,1);
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
    const auto prepared_bath=func::prescribed_harmonic_bath_field(bath_pspace,true);
    const size_t bath_nm=3*bath_pspace.num_TimePoints;
    failures+=require(prepared_bath.covariance.rows()==bath_nm
        &&prepared_bath.covariance.columns()==bath_nm,
        "harmonic preparation constructs only the Matsubara covariance");
    for(size_t i=0;i<bath_nm;++i)for(size_t j=0;j<bath_nm;++j)
        failures+=require(prepared_bath.covariance(i,j)==bath_field.covariance(i,j),
            "harmonic preparation matches the joint Matsubara block");
    const contour::StaggeredLayout bath_layout{
        bath_pspace.num_TimePoints,bath_pspace.num_RealTimeSteps};
    failures+=require(bath_field.covariance.rows()==bath_layout.dimension(),
        "harmonic bath covariance uses the complete staggered layout");
    failures+=require(bath_field.covariance==blaze::trans(bath_field.covariance),
        "harmonic bath covariance is complex symmetric");
    for( const auto& mean:bath_field.mean_time )
        failures+=require(blaze::length(mean)<RealType{1e-15},
            "prescribed harmonic bath has zero stochastic-field mean");
    const size_t eta=bath_layout.flat(contour::StaggeredSector::Eta,1,0);
    const size_t nu_past=bath_layout.flat(contour::StaggeredSector::Nu,0,0);
    const size_t nu_future=bath_layout.flat(contour::StaggeredSector::Nu,2,0);
    failures+=require(std::abs(bath_field.covariance(eta,nu_past))>1e-6,
        "harmonic bath exercises a nonzero retarded response covariance");
    failures+=require(bath_field.covariance(eta,nu_future)==ComplexType{}
        &&bath_field.covariance(nu_past,nu_future)==ComplexType{},
        "future response and response-response covariance are exactly zero");

    contour::CorrelationSet old_iteration{'A',2,1};
    contour::CorrelationSet raw_iteration{'A',2,1};
    contour::CorrelationSet iteration_errors{'A',2,1};
    raw_iteration.Re[0][0][0]=RealType{0.2};
    iteration_errors.Re[0][0][0]=RealType{0.1};
    raw_iteration.Im[0][0][1]=RealType{-0.3};
    iteration_errors.Im[0][0][1]=RealType{0.1};
    func::MagTen old_mag_Re{'D',1},old_mag_Im{'D',1};
    func::MagTen raw_mag_Re{'D',1},raw_mag_Im{'D',1};
    Observables::Magnetization::MagnetizationTensor<func::MagVec>
        mag_Re_errors{'D',1},mag_Im_errors{'D',1};
    raw_mag_Re[0][0]=RealType{0.3};
    raw_mag_Im[0][0]=RealType{0.4};
    mag_Re_errors[0][0]=RealType{0.06};
    mag_Im_errors[0][0]=RealType{0.2};
    auto scaled_mag_errors=RealType{2.}*mag_Re_errors;
    scaled_mag_errors+=mag_Re_errors;
    failures+=require(
        scaled_mag_errors.size()==1
        &&scaled_mag_errors.num_components()==3
        &&std::abs(scaled_mag_errors[0][0]-RealType{0.18})<RealType{1e-14},
        "magnetization tensor supports scalar multiplication and addition");
    Observables::Magnetization::MagnetizationTensor<func::MagVec> axial_mag{'C',2};
    axial_mag[1][0]=RealType{0.25};
    const auto full_axial_mag=axial_mag.expand();
    failures+=require(
        axial_mag[0].get_direction(0)==2
        &&full_axial_mag.size()==2
        &&full_axial_mag[1][0]==RealType{}
        &&full_axial_mag[1][1]==RealType{}
        &&full_axial_mag[1][2]==RealType{0.25},
        "magnetization tensor preserves symmetry-reduced MagVec storage");
    ps::ParameterSpace initial_mag_pspace;
    initial_mag_pspace.correlation_symmetry_type='C';
    initial_mag_pspace.num_RealTimePoints=2;
    initial_mag_pspace.initial_magnetization_linearized={1.,2.,3.};
    auto [zero_mag_Re,zero_mag_Im]
        =func::generate_initial_magnetization(initial_mag_pspace);
    failures+=require(
        zero_mag_Re.size()==2&&zero_mag_Re.num_components()==1
        &&zero_mag_Re[0][0]==RealType{}&&zero_mag_Im[1][0]==RealType{},
        "fresh magnetization initialization is zero");
    initial_mag_pspace.load_initial_spin_correlations=true;
    auto [imported_mag_Re,imported_mag_Im]
        =func::generate_initial_magnetization(initial_mag_pspace);
    failures+=require(
        imported_mag_Re[0][0]==RealType{3.}
        &&imported_mag_Re[1][0]==RealType{3.}
        &&imported_mag_Im[0][0]==RealType{},
        "loaded magnetization initializes every real-time point from imported m(0)");
    const auto residual=func::iteration_residual(
        old_iteration,
        raw_iteration,
        iteration_errors,
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{raw_mag_Re,raw_mag_Im},
        spinDMFT::Contour::Magnetization{mag_Re_errors,mag_Im_errors});
    failures+=require(std::abs(residual.absolute-RealType{0.5})<RealType{1e-14},
        "iteration residual retains the largest absolute complex difference");
    failures+=require(std::abs(residual.standardized-RealType{5.})<RealType{1e-14},
        "iteration residual standardizes real and imaginary magnetization separately");
    iteration_errors.Re[0][0][0]=RealType{};
    const auto unresolved=func::iteration_residual(
        old_iteration,
        raw_iteration,
        iteration_errors,
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{raw_mag_Re,raw_mag_Im},
        spinDMFT::Contour::Magnetization{mag_Re_errors,mag_Im_errors});
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
        exact_old,
        exact_raw,
        exact_errors,
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{mag_Re_errors,mag_Im_errors});
    failures+=require(endpoint_roundoff.standardized==RealType{},
        "roundoff at zero-variance t=0 tau endpoints is resolved");
    exact_raw.Re[0][0][1]=RealType{0.24};
    const auto endpoint_mismatch=func::iteration_residual(
        exact_old,
        exact_raw,
        exact_errors,
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{mag_Re_errors,mag_Im_errors});
    failures+=require(std::isinf(endpoint_mismatch.standardized),
        "a physical zero-variance t=0 tau endpoint mismatch remains unresolved");
    iteration_errors.Re[0][0][0]=std::numeric_limits<RealType>::quiet_NaN();
    const auto nonfinite=func::iteration_residual(
        old_iteration,
        raw_iteration,
        iteration_errors,
        spinDMFT::Contour::Magnetization{old_mag_Re,old_mag_Im},
        spinDMFT::Contour::Magnetization{raw_mag_Re,raw_mag_Im},
        spinDMFT::Contour::Magnetization{mag_Re_errors,mag_Im_errors});
    failures+=require(std::isinf(nonfinite.absolute)
                      &&std::isinf(nonfinite.standardized),
        "nonfinite iteration statistics cannot converge");
    contour::CorrelationSet constrained_old{'A',2,2},constrained_new{'A',2,2},constrained_errors{'A',2,2};
    func::MagTen constrained_old_mag{'D',2},constrained_new_mag{'D',2},constrained_mag_errors{'D',2};
    constrained_new.Re[0][0][0]=RealType{1.};
    constrained_new_mag[0][0]=RealType{1.};
    constrained_new.Re[1][0][0]=RealType{0.2};
    constrained_new_mag[1][0]=RealType{0.3};
    constrained_errors.Re[1][0][0]=RealType{0.1};
    constrained_mag_errors[1][0]=RealType{0.1};
    const auto frozen_residual=func::iteration_residual(
        constrained_old,
        constrained_new,
        constrained_errors,
        spinDMFT::Contour::Magnetization{constrained_old_mag,constrained_old_mag},
        spinDMFT::Contour::Magnetization{constrained_new_mag,constrained_old_mag},
        spinDMFT::Contour::Magnetization{constrained_mag_errors,constrained_mag_errors},
        true,
        false);
    failures+=require(std::abs(frozen_residual.absolute-RealType{0.3})<RealType{1e-14}
        &&std::abs(frozen_residual.standardized-RealType{3.})<RealType{1e-14},
        "frozen-origin residual excludes origin differences but retains real-time updates");
    const auto constant_residual=func::iteration_residual(
        constrained_old,
        constrained_new,
        constrained_errors,
        spinDMFT::Contour::Magnetization{constrained_old_mag,constrained_old_mag},
        spinDMFT::Contour::Magnetization{constrained_new_mag,constrained_old_mag},
        spinDMFT::Contour::Magnetization{constrained_mag_errors,constrained_mag_errors},
        true,
        true);
    failures+=require(std::abs(constant_residual.absolute-RealType{0.2})<RealType{1e-14}
        &&std::abs(constant_residual.standardized-RealType{2.})<RealType{1e-14},
        "constant-magnetization residual excludes constrained magnetizations");
    return failures==0?0:1;
}
