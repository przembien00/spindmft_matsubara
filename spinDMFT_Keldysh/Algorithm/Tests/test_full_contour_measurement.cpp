#include"../Functions/Functions.h"
#include"../Functions/PCN_Chain.h"

#include<algorithm>
#include<cmath>
#include<iostream>
#include<mpi.h>

namespace func = spinDMFT::Functions;
namespace contour = spinDMFT::Contour;
namespace ps = spinDMFT::Parameter_Space;
namespace rtd = spinDMFT::Run_Time_Data;

namespace
{
int require( bool condition, const char* message )
{
    if( condition ) return 0;
    std::cerr<<"FAILED: "<<message<<'\n';
    return 1;
}

// Independent reference: traverse every forward step, insert S at t, then
// traverse every backward step. No suffix/prefix factorization is used here.
ComplexType full_numerator( const func::ContourTrajectory& trajectory,
                            const Operator& initial, const Observable& spin,
                            size_t t )
{
    Operator value=initial;
    for( size_t k=0;k<trajectory.forward_steps.size();++k )
    {
        if( k>0 ) value=trajectory.forward_steps[k]*value;
        if( k==t ) value=spin*value;
    }
    for( size_t k=trajectory.backward_steps.size()-1;k>0;--k )
        value=trajectory.backward_steps[k]*value;
    return blaze::trace(value);
}

ComplexType prefix_numerator( const func::ContourTrajectory& trajectory,
                              const Operator& initial, const Observable& spin,
                              size_t t )
{
    Operator forward=IDENTITY,backward=IDENTITY;
    for( size_t k=1;k<=t;++k )
    {
        forward=trajectory.forward_steps[k]*forward;
        backward=backward*trajectory.backward_steps[k];
    }
    return blaze::trace(initial*backward*spin*forward);
}

int check_measurement( const ps::ParameterSpace& p,
                       const func::ContourTrajectory& trajectory, int rank,
                       const std::string& insertion_strategy="full-contour" )
{
    rtd::RunTimeData runtime(p,rank);
    for( size_t sample=0;sample<p.num_SamplesPerCore;++sample )
        func::compute_contour_correlations(
            runtime,trajectory,RealType{1.},insertion_strategy);
    contour::CorrelationSet correlations,errors;
    runtime.mpi_reduce_and_finalize(correlations,errors);
    const std::array<const Observable*,3> spins{&S_X,&S_Y,&S_Z};
    RealType correlation_error{},magnetization_error{},closure_error{};
    Operator density=trajectory.imaginary_density_operator;
    for( size_t t=0;t<p.num_RealTimePoints;++t )
    {
        if( t>0 ) density=trajectory.forward_steps[t]*density
                         *trajectory.backward_steps[t];
        const ComplexType closure{runtime.closed_contour_ratio_Re[t],
                                  runtime.closed_contour_ratio_Im[t]};
        closure_error=std::max(closure_error,std::abs(closure
            -blaze::trace(density)/trajectory.partition_function));
        for( size_t c=0;c<runtime.num_magnetization_components();++c )
        {
            const size_t a=runtime.magnetization_direction(c);
            const ComplexType expected=(insertion_strategy=="prefix"
                ?prefix_numerator(trajectory,
                    trajectory.imaginary_density_operator,*spins[a],t)
                :full_numerator(trajectory,
                    trajectory.imaginary_density_operator,*spins[a],t))
                /trajectory.partition_function;
            const ComplexType actual{runtime.magnetization_time_Re[t][a],
                                     runtime.magnetization_time_Im[t][a]};
            magnetization_error=std::max(magnetization_error,std::abs(actual-expected));
        }
        for( size_t pair=0;pair<runtime.num_correlation_components();++pair )
        {
            const auto ab=runtime.correlation_direction(pair);
            for( size_t tau=0;tau<p.num_TimePoints;++tau )
            {
                const ComplexType expected=(insertion_strategy=="prefix"
                    ?prefix_numerator(trajectory,
                        trajectory.imaginary_edge_insertions[ab[1]][tau],
                        *spins[ab[0]],t)
                    :full_numerator(trajectory,
                        trajectory.imaginary_edge_insertions[ab[1]][tau],
                        *spins[ab[0]],t))
                    /trajectory.partition_function;
                const ComplexType actual{correlations.Re[t][pair][tau],
                                         correlations.Im[t][pair][tau]};
                correlation_error=std::max(correlation_error,std::abs(actual-expected));
            }
        }
    }
    return require(correlation_error<RealType{2e-12},
                   "all correlation components and times include the full contour")
         + require(magnetization_error<RealType{2e-12},
                   "magnetization uses the same full-contour spin insertion")
         + require(closure_error<RealType{2e-12},
                   "prefix contour-closure diagnostic remains unchanged");
}
}

int main( int argc, char** argv )
{
    MPI_Init(&argc,&argv);
    int rank{},world_size{};
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    ps::ParameterSpace p(argc,argv,world_size,rank);
    p.sampling_strategy="independent";
    p.correlation_symmetry_type='D';
    p.num_SamplesPerCore=2; p.num_blocks=2;
    p.beta=RealType{0.8}; p.num_TimeSteps=4; p.num_TimePoints=5;
    p.delta_t=p.beta/p.num_TimeSteps;
    p.num_RealTimeSteps=4; p.num_RealTimePoints=5;
    p.Tmax=RealType{0.6}; p.delta_real_t=p.Tmax/p.num_RealTimeSteps;
    int failures{};
    {
        p.antithetic_pairs=true;
        p.num_SamplesPerCore=12; p.num_blocks=4;
        rtd::RunTimeData pair_blocks(p,rank);
        failures+=require(pair_blocks.num_blocks()==3
                          &&pair_blocks.samples_per_block()==4,
            "antithetic jackknife blocks are adjusted to contain complete pairs");
        p.antithetic_pairs=false;
        p.num_SamplesPerCore=2; p.num_blocks=2;
    }
    func::initialize_matrices(p);
    const contour::ContourLayout layout{p.num_TimePoints,p.num_RealTimePoints};
    func::DenseComplexGaussianSampler::FieldVector fields(layout.dimension(),ComplexType{});
    const func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    for( size_t k=0;k<p.num_TimePoints;++k )
        for( size_t c=0;c<3;++c )
            fields[layout.flat(contour::Branch::Matsubara,k,c)]={
                RealType{0.09}*(k+1)*(c+1),RealType{0.03}*(k+c+1)};
    for( size_t k=0;k<p.num_RealTimePoints;++k )
        for( size_t c=0;c<3;++c )
        {
            fields[layout.flat(contour::Branch::Forward,k,c)]={
                RealType{0.13}*(k+1)+RealType{0.07}*(c+1)*(k*k+1),
                RealType{-0.06}*(k+c+1)};
            fields[layout.flat(contour::Branch::Backward,k,c)]={
                RealType{-0.17}*(k+1)+RealType{0.04}*(c+1)*(k*k+1),
                RealType{0.08}*(2*k+c+1)};
        }
    const auto distinct=func::compute_contour_trajectory(p,fields,mean);
    failures+=check_measurement(p,distinct,rank);
    failures+=check_measurement(p,distinct,rank,"prefix");

    // Earlier insertions must depend on both future forward and backward steps.
    for( const bool change_forward : {false,true} )
    {
        auto changed=distinct;
        auto& steps=change_forward?changed.forward_steps:changed.backward_steps;
        steps.back()=func::general_matrix_exponential(
            ComplexType{0.11,-0.19}*S_X+ComplexType{-0.07,0.23}*S_Y)*steps.back();
        failures+=check_measurement(p,changed,rank);
        for( const size_t t : {size_t{0},size_t{2}} )
            failures+=require(std::abs(full_numerator(changed,
                changed.imaginary_edge_insertions[0][1],S_X,t)
                -full_numerator(distinct,distinct.imaginary_edge_insertions[0][1],S_X,t))
                >RealType{1e-5},"future branch steps affect earlier correlations");
    }

    for( size_t k=0;k<p.num_RealTimePoints;++k )
        for( size_t c=0;c<3;++c )
            fields[layout.flat(contour::Branch::Backward,k,c)]=
                fields[layout.flat(contour::Branch::Forward,k,c)];
    const auto equal=func::compute_contour_trajectory(p,fields,mean);
    failures+=check_measurement(p,equal,rank);
    Operator forward=IDENTITY,backward=IDENTITY;
    for( size_t t=0;t<p.num_RealTimePoints;++t )
    {
        if( t>0 )
        {
            forward=equal.forward_steps[t]*forward;
            backward=backward*equal.backward_steps[t];
        }
        failures+=require(std::abs(full_numerator(equal,
            equal.imaginary_edge_insertions[0][1],S_X,t)
            -blaze::trace(equal.imaginary_edge_insertions[0][1]*backward*S_X*forward))
            <RealType{2e-12},"equal complex noncommuting branch fields recover the prefix estimator");
    }

    // Degenerate trajectory: no real-time steps, only the t=0 insertion.
    auto zero=distinct;
    zero.forward_steps.resize(1); zero.backward_steps.resize(1);
    p.num_RealTimeSteps=0; p.num_RealTimePoints=1;
    failures+=check_measurement(p,zero,rank);

    p.self_consistency=true;
    p.Iteration_Limit=20;
    p.iteration_error_sigma_threshold=RealType{5.};
    auto prime_diagnostics=[]( rtd::RunTimeData& runtime )
    {
        runtime.covariance_symmetry_errors.push_back(RealType{});
        runtime.branch_identity_errors.push_back(RealType{});
        runtime.gaussian_factor_reconstruction_errors.push_back(RealType{});
        runtime.average_phase_magnitudes.push_back(RealType{1.});
        runtime.denominator_constancy_residuals.push_back(RealType{});
        runtime.magnetization_time_Im.assign(1,FieldVector{});
        runtime.magnetization_time_Im_stds.assign(1,FieldVector{});
        runtime.num_Iterations=1;
    };
    rtd::RunTimeData accepted(p,rank);
    prime_diagnostics(accepted);
    accepted.record_iteration_error(RealType{0.2},RealType{4.99});
    failures+=require(accepted.terminate(),
        "standardized residual below q passes the fixed-point stopping rule");
    rtd::RunTimeData rejected(p,rank);
    prime_diagnostics(rejected);
    rejected.record_iteration_error(RealType{1e-8},RealType{5.});
    failures+=require(!rejected.terminate(),
        "standardized residual equal to q does not pass the strict stopping rule");

    // A zero-rank Gaussian has constant field and partition function.  Every
    // pCN proposal must therefore be accepted, including the zero-dimensional
    // latent-state edge case.
    p.sampling_strategy="pcn";
    p.antithetic_pairs=false;
    p.num_RealTimeSteps=1; p.num_RealTimePoints=2;
    const contour::ContourLayout pcn_layout{p.num_TimePoints,p.num_RealTimePoints};
    func::ComplexDynamicMatrix zero_covariance(
        pcn_layout.dimension(),pcn_layout.dimension(),ComplexType{});
    func::DenseComplexGaussianSampler zero_sampler(zero_covariance);
    const func::MeanFieldTrajectory zero_mean(p.num_RealTimePoints,FieldVector{});
    std::mt19937 pcn_engine{12345};
    func::PCNChain constant_chain(p,zero_sampler,zero_mean,RealType{0.4},pcn_engine);
    for( size_t step=0;step<6;++step ) constant_chain.step();
    failures+=require(constant_chain.proposed()==6&&constant_chain.accepted()==6,
        "constant-likelihood pCN accepts every proposal");
    failures+=require(constant_chain.real_partition()>RealType{},
        "pCN current state has positive real partition function");

    // Antithetic pCN keeps both signs in one state and uses their summed real
    // partition function as the Metropolis likelihood.  For the zero-rank
    // Gaussian both signs coincide, so the paired estimator must reduce
    // exactly to the ordinary pCN estimator while still counting one state.
    p.antithetic_pairs=true;
    std::mt19937 paired_engine{54321};
    func::PCNChain paired_chain(
        p,zero_sampler,zero_mean,RealType{0.4},paired_engine);
    for( size_t step=0;step<6;++step ) paired_chain.step();
    const RealType paired_partition=
        std::real(paired_chain.trajectory().partition_function)
       +std::real(paired_chain.antithetic_trajectory().partition_function);
    failures+=require(paired_chain.uses_antithetic_pairs()
                      &&std::abs(paired_chain.real_partition()-paired_partition)
                         <RealType{1e-13},
        "antithetic pCN likelihood is the sum over both latent signs");
    failures+=require(paired_chain.proposed()==6&&paired_chain.accepted()==6,
        "constant sign-symmetrized likelihood accepts every pCN proposal");

    rtd::RunTimeData paired_runtime(p,rank);
    for( size_t sample=0;sample<p.num_SamplesPerCore;++sample )
        func::compute_contour_pair_correlations(
            paired_runtime,paired_chain.trajectory(),
            paired_chain.antithetic_trajectory(),
            RealType{1.}/paired_chain.real_partition());
    contour::CorrelationSet paired_mean,paired_error;
    paired_runtime.mpi_reduce_and_finalize(paired_mean,paired_error);

    p.antithetic_pairs=false;
    rtd::RunTimeData ordinary_runtime(p,rank);
    for( size_t sample=0;sample<p.num_SamplesPerCore;++sample )
        func::compute_contour_correlations(
            ordinary_runtime,constant_chain.trajectory(),
            RealType{1.}/constant_chain.real_partition());
    contour::CorrelationSet ordinary_mean,ordinary_error;
    ordinary_runtime.mpi_reduce_and_finalize(ordinary_mean,ordinary_error);
    failures+=require(
        func::max_contour_difference(paired_mean,ordinary_mean)
            <RealType{1e-13},
        "sign-symmetrized pair measurement is accumulated as one pCN state");

    // Deliberately correlated block means: merging adjacent equal blocks must
    // increase the reported uncertainty, demonstrating that the pCN path does
    // not treat chain states as independent jackknife samples.
    p.num_RealTimeSteps=0; p.num_RealTimePoints=1;
    p.num_SamplesPerCore=16; p.num_blocks=8;
    rtd::RunTimeData correlated(p,rank);
    const std::array<RealType,16> correlated_values{
        0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1};
    for( const RealType value:correlated_values )
    {
        correlated.begin_sample(ComplexType{1.,0.},RealType{1.});
        correlated.accumulate_closed_contour_trace(0,ComplexType{1.,0.});
        correlated.accumulate_edge_correlation(0,0,0,ComplexType{value,0.});
        correlated.end_sample();
    }
    contour::CorrelationSet correlated_mean,correlated_error;
    correlated.mpi_reduce_and_finalize(correlated_mean,correlated_error);
    failures+=require(correlated.blocking_curve_mean_errors.size()>=3,
        "pCN statistics expose a multi-scale blocking curve");
    failures+=require(correlated.blocking_curve_max_errors[1]
                     >correlated.blocking_curve_max_errors[0],
        "pCN batch merging detects autocorrelation variance inflation");
    failures+=require(correlated.contour_tau_int.Re[0][0][0]>RealType{0.5},
        "pCN statistics report an autocorrelation time above the iid value");
    MPI_Finalize();
    return failures==0?0:1;
}
