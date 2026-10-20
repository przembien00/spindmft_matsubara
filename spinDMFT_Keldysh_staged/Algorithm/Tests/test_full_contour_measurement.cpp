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
                       const std::string& insertion_strategy="closed-contour" )
{
    rtd::RunTimeData runtime(p,rank);
    for( size_t sample=0;sample<p.num_SamplesPerCore;++sample )
        func::compute_contour_correlations(
            runtime,trajectory,RealType{1.},insertion_strategy);
    contour::CorrelationSet correlations;
    rtd::MagTen magnetization_Re,magnetization_Im;
    runtime.mpi_reduce_and_finalize(
        correlations,magnetization_Re,magnetization_Im);
    const std::array<const Observable*,3> spins{&S_X,&S_Y,&S_Z};
    RealType correlation_error{},magnetization_error{},closure_error{};
    Operator final_density=trajectory.imaginary_density_operator;
    for( size_t t=1;t<p.num_RealTimePoints;++t )
        final_density=trajectory.forward_steps[t]*final_density
                     *trajectory.backward_steps[t];
    const ComplexType observable_denominator=
        p.correlation_normalization=="closed-contour"
        ?blaze::trace(final_density):trajectory.partition_function;
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
                /observable_denominator;
            const ComplexType actual{magnetization_Re[t][c],
                                     magnetization_Im[t][c]};
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
                    /observable_denominator;
                const ComplexType actual{correlations.Re[t][pair][tau],
                                         correlations.Im[t][pair][tau]};
                correlation_error=std::max(correlation_error,std::abs(actual-expected));
            }
        }
    }
    return require(correlation_error<RealType{2e-12},
                   "all correlation components and times include the closed contour")
         + require(magnetization_error<RealType{2e-12},
                   "magnetization uses the same closed-contour spin insertion")
         + require(closure_error<RealType{2e-12},
                   "prefix contour-closure diagnostic remains unchanged");
}
}

template<typename Matrix>
int check_fused_accumulation(ps::ParameterSpace p,int rank)
{
    int failures{};
    p.sampling_strategy="independent";p.num_SamplesPerCore=45;p.num_blocks=5;
    p.num_Samples=45*p.get_num_Cores();
    const size_t size=p.uses_staggered_strang()
        ?contour::StaggeredLayout{p.num_TimePoints,p.num_RealTimeSteps*p.real_time_steps_per_interval()}.dimension()
        :contour::ContourLayout{p.num_TimePoints,p.num_RealTimePoints}.dimension();
    func::JointComplexGaussianSampler::ContourFieldSample field;
    field.edge_field.resize(size);
    func::BasicContourTrajectory<Matrix> trajectory;
    func::BasicTrajectoryWorkspace<Matrix> propagation;
    func::MeasurementWorkspace buffered_workspace,fused_workspace;
    rtd::MeasuredSample sample;
    const func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{0.03,-0.02,0.05});
    for(const std::string normalization:{"partition-function","closed-contour"})
    {
        p.correlation_normalization=normalization;
        rtd::RunTimeData buffered(p,rank),fused(p,rank);
        // Vary complex weights and numerators across statistical blocks.
        for(size_t s=0;s<p.num_SamplesPerCore;++s)
        {
            for(size_t i=0;i<size;++i)
                field.edge_field[i]={0.12*std::sin(0.3*(i+2*s+rank)),0.04*std::cos(0.2*(2*i+s))};
            func::build_contour_trajectory(p,field,mean,trajectory,propagation,false,func::FinalTraceMode::Skip);
            func::measure_contour_observables(buffered,trajectory,sample,buffered_workspace,"prefix");
            buffered.accumulate_sample(sample);
            func::measure_and_accumulate_independent(fused,trajectory,fused_workspace,"prefix");
        }
        contour::CorrelationSet a,b;
        rtd::MagTen ar,ai,br,bi;
        buffered.mpi_reduce_and_finalize(a,ar,ai);
        fused.mpi_reduce_and_finalize(b,br,bi);
        failures+=require(func::max_contour_difference(a,b)==RealType{},
            "fused measurement preserves complex-ratio means exactly");
        failures+=require(func::max_contour_difference(buffered.contour_sample_stds,fused.contour_sample_stds)==RealType{},
            "fused measurement preserves statistical-block errors exactly");
        for(size_t t=0;t<ar.size();++t)for(size_t c=0;c<ar[t].size();++c)
            failures+=require(ar[t][c]==br[t][c]&&ai[t][c]==bi[t][c]
                &&buffered.magnetization_time_Re_stds[t][c]==fused.magnetization_time_Re_stds[t][c]
                &&buffered.magnetization_time_Im_stds[t][c]==fused.magnetization_time_Im_stds[t][c],
                "fused magnetization means and errors agree exactly");
        failures+=require(buffered.closed_contour_ratio_Re==fused.closed_contour_ratio_Re
            &&buffered.closed_contour_ratio_Im==fused.closed_contour_ratio_Im
            &&buffered.closed_contour_residual_Re_sample_stds==fused.closed_contour_residual_Re_sample_stds
            &&buffered.closed_contour_residual_Im_sample_stds==fused.closed_contour_residual_Im_sample_stds
            &&buffered.closed_contour_residual_abs_sample_stds==fused.closed_contour_residual_abs_sample_stds
            &&buffered.partition_sum_Re==fused.partition_sum_Re
            &&buffered.partition_sum_Im==fused.partition_sum_Im
            &&buffered.partition_abs_sums==fused.partition_abs_sums,
            "fused closure and partition diagnostics agree exactly");
    }
    // The shared estimator also retains the full-contour insertion strategy.
    p.num_SamplesPerCore=1;p.num_blocks=1;
    rtd::RunTimeData full(p,rank);
    func::measure_and_accumulate_independent(full,trajectory,fused_workspace,"closed-contour");
    bool rejected{};
    trajectory.imaginary_edge_insertions[0][0](0,0)=std::numeric_limits<RealType>::infinity();
    rtd::RunTimeData invalid(p,rank);
    try{func::measure_and_accumulate_independent(invalid,trajectory,fused_workspace,"prefix");}
    catch(const std::runtime_error&){rejected=true;}
    failures+=require(rejected,"fused measurement rejects nonfinite correlation numerators");
    p.sampling_strategy="pcn";
    rtd::RunTimeData pcn(p,rank);
    rejected=false;
    try{func::measure_and_accumulate_independent(pcn,trajectory,fused_workspace,"prefix");}
    catch(const std::invalid_argument&){rejected=true;}
    failures+=require(rejected,"independent fused path rejects pCN accumulators");
    return failures;
}

int check_pcn_cache(ps::ParameterSpace p,int rank)
{
    int failures{};
    p.sampling_strategy="pcn";p.num_SamplesPerCore=96;p.num_blocks=16;
    p.beta=RealType{2.};p.delta_t=p.beta/p.num_TimeSteps;
    const contour::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps*p.real_time_steps_per_interval()};
    func::ComplexDynamicMatrix covariance(layout.dimension(),layout.dimension(),ComplexType{});
    for(size_t i=0;i<layout.dimension();++i)
        covariance(i,i)=layout.decode(i).sector==contour::StaggeredSector::Matsubara?RealType{1.4}:RealType{0.03};
    const func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    for(const std::string normalization:{"partition-function","closed-contour"})
    {
        p.correlation_normalization=normalization;
        func::DenseComplexGaussianSampler sampler(covariance),reference_sampler(covariance);
        std::mt19937 engine{924},reference_engine{924};
        const RealType step=RealType{0.95},retention=std::sqrt(RealType{1.}-step*step);
        func::PCNChain chain(p,sampler,mean,step,engine);
        const auto weight=[&](const func::ContourTrajectory& tr)
        {return normalization=="closed-contour"?tr.final_closed_contour_trace:tr.partition_function;};
        const auto positive=[](ComplexType z){return std::isfinite(std::real(z))&&std::isfinite(std::imag(z))&&std::real(z)>0.;};
        func::JointComplexGaussianSampler::LatentVector latent;
        func::ContourTrajectory reference;
        for(size_t attempt=0;attempt<128;++attempt)
        {
            latent=reference_sampler.draw_latent(reference_engine);
            reference=func::compute_contour_trajectory(p,reference_sampler.contour_field_from_latent(latent,false),mean);
            if(positive(weight(reference)))break;
        }
        std::uniform_real_distribution<RealType> uniform(0.,1.);
        rtd::RunTimeData cached(p,rank),eager(p,rank);
        rtd::MeasuredSample measured;
        func::MeasurementWorkspace workspace;
        size_t rejections{};bool valid=false;
        for(size_t i=0;i<p.num_SamplesPerCore;++i)
        {
            const auto innovation=reference_sampler.draw_latent(reference_engine);
            auto proposal=latent;
            for(size_t j=0;j<proposal.size();++j)proposal[j]=retention*latent[j]+step*innovation[j];
            const auto proposed=func::compute_contour_trajectory(p,reference_sampler.contour_field_from_latent(proposal,false),mean);
            bool accepted=false;
            if(positive(weight(proposed)))
            {
                const RealType log_alpha=std::log(std::real(weight(proposed)))-std::log(std::real(weight(reference)));
                accepted=std::log(uniform(reference_engine))<std::min(RealType{},log_alpha);
            }
            if(accepted){latent=proposal;reference=proposed;}
            else ++rejections;
            failures+=require(chain.step()==accepted,"Matsubara-first pCN preserves eager acceptance decisions");
            failures+=require(std::abs(chain.real_sampling_weight()-std::real(weight(reference)))<RealType{1e-12},
                              "pCN likelihood matches eager full-contour reference");
            if(accepted||!valid)
            {
                func::measure_contour_observables(cached,chain.trajectory(),measured,workspace,p.spin_insertion_strategy);
                valid=true;
            }
            cached.accumulate_sample(measured,RealType{1.}/chain.real_sampling_weight());
            func::compute_contour_correlations(eager,reference,RealType{1.}/std::real(weight(reference)),p.spin_insertion_strategy);
        }
        failures+=require(rejections>0,"cache regression exercises rejected proposals");
        failures+=require(engine==reference_engine,"delayed propagation preserves pCN RNG consumption");
        contour::CorrelationSet cached_mean,eager_mean;
        rtd::MagTen cached_mag_Re,cached_mag_Im,eager_mag_Re,eager_mag_Im;
        cached.mpi_reduce_and_finalize(
            cached_mean,cached_mag_Re,cached_mag_Im);
        eager.mpi_reduce_and_finalize(
            eager_mean,eager_mag_Re,eager_mag_Im);
        failures+=require(func::max_contour_difference(cached_mean,eager_mean)<RealType{1e-12},"cached rejected states preserve complex-ratio means");
        failures+=require(func::max_contour_difference(cached.contour_sample_stds,eager.contour_sample_stds)<RealType{1e-12},"cached rejected states preserve pCN block errors");
        // Ratio-influence variances include subtractive cancellation. Compare
        // tau relatively: these nearly constant components can have tau >400,
        // and subtractive cancellation amplifies accumulation roundoff to about
        // 1e-9 relative. Means and block errors remain checked to 1e-12.
        RealType tau_error{};
        for(const auto pair:{std::make_pair(&cached.contour_tau_int.Re,&eager.contour_tau_int.Re),
                             std::make_pair(&cached.contour_tau_int.Im,&eager.contour_tau_int.Im)})
            for(size_t t=0;t<pair.first->size();++t)
                for(size_t c=0;c<(*pair.first)[t].size();++c)
                    for(size_t k=0;k<(*pair.first)[t][c].size();++k)
                    {
                        const auto a=(*pair.first)[t][c][k],b=(*pair.second)[t][c][k];
                        const RealType relative=std::abs(a-b)/(RealType{1.}+std::max(std::abs(a),std::abs(b)));
                        tau_error=std::max(tau_error,relative);
                    }
        failures+=require(tau_error<RealType{1e-8},"cached rejected states preserve autocorrelation statistics");
    }
    return failures;
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
    func::initialize_matrices(p);
    {
        auto stopping=p;stopping.Iteration_Limit=10;
        rtd::RunTimeData data(stopping,rank);
        data.record_complex_field_diagnostics(0.,0.,0.,0,0);
        data.average_phase_magnitudes.push_back(1.);
        rtd::MagTen zero_im{'D',p.num_RealTimePoints};
        data.magnetization_time_Im_stds=rtd::MagTen{'D',p.num_RealTimePoints};
        // Recorded Strang failures had enormous absolute updates but less
        // than five sigma because their sampling variance also exploded.
        data.record_iteration_error(7.225e11,4.897);
        data.num_Iterations=4;
        failures+=require(!data.terminate(zero_im),
            "inflated error bars cannot declare an exploding real-time iteration converged");
        data.record_iteration_error(0.05,4.);
        failures+=require(data.terminate(zero_im),
            "resolved sub-spin-scale residuals retain the existing standardized stopping rule");
    }
    {
        auto strang=p;
        strang.real_time_evolution="staggered-strang";strang.real_time_substeps=2;
        strang.gaussian_factorization="dense";
        const contour::StaggeredLayout grid{p.num_TimePoints,2*p.num_RealTimeSteps};
        func::DenseComplexGaussianSampler::FieldVector z(grid.dimension());
        for(size_t i=0;i<z.size();++i)z[i]=ComplexType{0.003*RealType(i+1),-0.002*RealType(i+1)};
        const auto tr=func::compute_contour_trajectory(strang,z,
            func::MeanFieldTrajectory(p.num_RealTimePoints,FieldVector{}));
        // Compare inserted observables with independent explicit traversals
        // of nonunitary Strang branch maps, including mixed thermal insertions.
        failures+=check_measurement(strang,tr,rank,"prefix");
        failures+=check_measurement(strang,tr,rank,"closed-contour");
        failures+=check_fused_accumulation<func::SpinHalfOperator>(strang,rank);
    }
    failures+=check_pcn_cache(p,rank);
    failures+=check_fused_accumulation<Operator>(p,rank);
    failures+=check_fused_accumulation<func::SpinHalfOperator>(p,rank);
    const contour::StaggeredLayout layout{p.num_TimePoints,p.num_RealTimeSteps*p.real_time_steps_per_interval()};
    func::DenseComplexGaussianSampler::FieldVector fields(layout.dimension(),ComplexType{});
    const func::MeanFieldTrajectory mean(p.num_RealTimePoints,FieldVector{});
    for(size_t i=0;i<fields.size();++i)
        fields[i]=ComplexType{0.03*std::sin(RealType(i)),0.02*std::cos(RealType(i))};
    const auto distinct=func::compute_contour_trajectory(p,fields,mean);
    Operator distinct_final_density=distinct.imaginary_density_operator;
    for( size_t t=1;t<p.num_RealTimePoints;++t )
        distinct_final_density=distinct.forward_steps[t]*distinct_final_density
                             *distinct.backward_steps[t];
    failures+=require(std::abs(distinct.final_closed_contour_trace
                              -blaze::trace(distinct_final_density))
                      <RealType{2e-12},
        "trajectory retains the final closed-contour sampling weight");
    failures+=check_measurement(p,distinct,rank);
    failures+=check_measurement(p,distinct,rank,"prefix");
    p.correlation_normalization="closed-contour";
    failures+=check_measurement(p,distinct,rank);
    failures+=check_measurement(p,distinct,rank,"prefix");
    p.correlation_normalization="partition-function";

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
                >RealType{1e-7},"future branch steps affect earlier correlations");
    }

    // Zero response gives identical complex branch Hamiltonians.
    for(size_t k=0;k<=p.num_RealTimeSteps*p.real_time_steps_per_interval();++k)
        for(size_t c=0;c<3;++c)
            fields[layout.flat(contour::StaggeredSector::Nu,k,c)]=ComplexType{};
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

    // Closed-contour normalization is a ratio of accumulated means, not an
    // average of trajectorywise N/D values, for either sampling path.
    p.correlation_normalization="closed-contour";
    p.num_SamplesPerCore=2; p.num_blocks=2;
    for( const std::string strategy : {"independent","pcn"} )
    {
        p.sampling_strategy=strategy;
        rtd::RunTimeData normalized(p,rank);
        for( const auto sample : std::array<std::array<RealType,2>,2>{
                 std::array<RealType,2>{RealType{1.},RealType{2.}},
                 std::array<RealType,2>{RealType{3.},RealType{12.}}} )
        {
            normalized.begin_sample(ComplexType{1.,0.});
            normalized.accumulate_closed_contour_trace(
                0,ComplexType{sample[0],0.});
            normalized.accumulate_edge_correlation(
                0,0,0,ComplexType{sample[1],0.});
            normalized.accumulate_magnetization(
                0,0,ComplexType{sample[1],0.});
            normalized.end_sample();
        }
        contour::CorrelationSet normalized_mean;
        rtd::MagTen normalized_mag_Re,normalized_mag_Im;
        normalized.mpi_reduce_and_finalize(
            normalized_mean,normalized_mag_Re,normalized_mag_Im);
        failures+=require(
            std::abs(normalized_mean.Re[0][0][0]-RealType{3.5})
                <RealType{1e-13},
            "closed-contour correlation normalization uses sum N over sum D(T)");
        failures+=require(
            std::abs(normalized_mag_Re[0][0]-RealType{3.5})
                <RealType{1e-13},
            "closed-contour magnetization normalization uses sum M over sum D(T)");
    }
    p.correlation_normalization="partition-function";
    p.sampling_strategy="independent";

    // A pCN importance estimator must retain the correspondingly reweighted
    // complex denominator. Otherwise the finite-chain arithmetic mean of
    // N/Re Z violates the exact spin-half identity N=Z/4 whenever Z has a
    // fluctuating phase.
    p.sampling_strategy="pcn";
    p.num_SamplesPerCore=2; p.num_blocks=2;
    rtd::RunTimeData pcn_ratio(p,rank);
    for( const ComplexType Z : std::array<ComplexType,2>{
             ComplexType{1.,0.4},ComplexType{2.,-0.7}} )
    {
        const RealType inverse_likelihood=RealType{1.}/std::real(Z);
        pcn_ratio.begin_sample(Z,inverse_likelihood);
        pcn_ratio.accumulate_edge_correlation(
            0,0,0,RealType{0.25}*Z);
        pcn_ratio.end_sample();
    }
    contour::CorrelationSet pcn_ratio_mean;
    rtd::MagTen pcn_ratio_mag_Re,pcn_ratio_mag_Im;
    pcn_ratio.mpi_reduce_and_finalize(
        pcn_ratio_mean,pcn_ratio_mag_Re,pcn_ratio_mag_Im);
    failures+=require(
        std::abs(ComplexType{pcn_ratio_mean.Re[0][0][0],
                             pcn_ratio_mean.Im[0][0][0]}-RealType{0.25})
            <RealType{1e-13},
        "pCN retains its reweighted denominator for exact spin identities");
    failures+=require(
        pcn_ratio.contour_sample_stds.Re[0][0][0]<RealType{1e-13}
            &&pcn_ratio.contour_sample_stds.Im[0][0][0]<RealType{1e-13},
        "pCN ratio blocking reports zero error for an exact spin identity");
    p.sampling_strategy="independent";

    p.self_consistency=true;
    p.Iteration_Limit=20;
    p.iteration_error_sigma_threshold=RealType{5.};
    auto prime_diagnostics=[]( rtd::RunTimeData& runtime )
    {
        runtime.covariance_symmetry_errors.push_back(RealType{});
        runtime.gaussian_factor_reconstruction_errors.push_back(RealType{});
        runtime.average_phase_magnitudes.push_back(RealType{1.});
        runtime.denominator_constancy_residuals.push_back(RealType{});
        runtime.magnetization_time_Im_stds=rtd::MagTen{'D',1};
        runtime.num_Iterations=1;
    };
    const rtd::MagTen zero_magnetization{'D',1};
    rtd::RunTimeData accepted(p,rank);
    prime_diagnostics(accepted);
    accepted.record_iteration_error(RealType{0.2},RealType{4.99});
    failures+=require(accepted.terminate(zero_magnetization),
        "standardized residual below q passes the fixed-point stopping rule");
    rtd::RunTimeData rejected(p,rank);
    prime_diagnostics(rejected);
    rejected.record_iteration_error(RealType{1e-8},RealType{5.});
    failures+=require(!rejected.terminate(zero_magnetization),
        "standardized residual equal to q does not pass the strict stopping rule");

    // A non-default q must also control the imaginary-magnetization gate.
    auto shared_threshold_p=p;
    shared_threshold_p.iteration_error_sigma_threshold=RealType{2.};
    rtd::MagTen imaginary_magnetization{'D',1};
    imaginary_magnetization.front()[0]=RealType{0.03};
    rtd::RunTimeData imaginary_rejected(shared_threshold_p,rank);
    prime_diagnostics(imaginary_rejected);
    imaginary_rejected.magnetization_time_Im_stds.front()[0]=RealType{0.01};
    imaginary_rejected.record_iteration_error(RealType{0.01},RealType{0.1});
    failures+=require(!imaginary_rejected.terminate(imaginary_magnetization),
        "shared q=2 rejects imaginary magnetization at three standard errors");
    shared_threshold_p.iteration_error_sigma_threshold=RealType{4.};
    rtd::RunTimeData imaginary_accepted(shared_threshold_p,rank);
    prime_diagnostics(imaginary_accepted);
    imaginary_accepted.magnetization_time_Im_stds.front()[0]=RealType{0.01};
    imaginary_accepted.record_iteration_error(RealType{0.01},RealType{0.1});
    failures+=require(imaginary_accepted.terminate(imaginary_magnetization),
        "shared q=4 accepts imaginary magnetization at three standard errors");

    // Preparation follows the absolute Matsubara rule even when measured
    // noise makes the standardized residual look acceptable (or infinite).
    auto prep_p=p;
    prep_p.num_SamplesPerCore=10000;
    prep_p.num_Samples=10000*prep_p.get_num_Cores();
    prep_p.spin_float=RealType{0.5};
    rtd::RunTimeData prep_rejected(prep_p,rank,true);
    prime_diagnostics(prep_rejected);
    const RealType prep_threshold=RealType{5.}*RealType{0.25}
        /std::sqrt(RealType{3.}*static_cast<RealType>(prep_p.num_Samples));
    failures+=require(std::abs(prep_rejected.matsubara_iteration_error_threshold()
        -prep_threshold)<RealType{1e-14},
        "preparation threshold uses spin and total preparation samples");
    prep_rejected.record_iteration_error(RealType{0.02},RealType{0.1});
    failures+=require(!prep_rejected.terminate(zero_magnetization),
        "large sampling errors cannot relax Matsubara convergence");
    rtd::RunTimeData prep_boundary(prep_p,rank,true);
    prime_diagnostics(prep_boundary);
    prep_boundary.record_iteration_error(prep_threshold,RealType{});
    failures+=require(!prep_boundary.terminate(zero_magnetization),
        "absolute Matsubara threshold is strict");
    rtd::RunTimeData prep_accepted(prep_p,rank,true);
    prime_diagnostics(prep_accepted);
    prep_accepted.record_iteration_error(prep_threshold/RealType{2.},
        std::numeric_limits<RealType>::infinity());
    failures+=require(prep_accepted.terminate(zero_magnetization),
        "standardized residual is only diagnostic during preparation");

    contour::CorrelationSet prep_old{'D',4,1},prep_raw{'D',4,1};
    rtd::MagTen prep_old_mag{'D',1},prep_raw_mag{'D',1};
    prep_raw.Re[0][0][0]=RealType{0.12}; // mean absolute change = 0.03
    prep_raw.Im[0][1][1]=RealType{-0.2}; // mean absolute change = 0.05
    prep_raw_mag[0][0]=RealType{0.04};
    rtd::RunTimeData prep_metric(prep_p,rank,true);
    prep_metric.record_matsubara_iteration_error(
        prep_old,prep_raw,prep_old_mag,prep_raw_mag,RealType{0.1});
    failures+=require(std::abs(prep_metric.absolute_iteration_errors.back()
        -RealType{0.05})<RealType{1e-14},
        "Matsubara metric maximizes component mean absolute Re/Im changes");
    prep_raw_mag[0][0]=RealType{-0.08};
    prep_metric.record_matsubara_iteration_error(
        prep_old,prep_raw,prep_old_mag,prep_raw_mag,RealType{0.1});
    failures+=require(std::abs(prep_metric.absolute_iteration_errors.back()
        -RealType{0.08})<RealType{1e-14},
        "Matsubara metric includes the real magnetization change");
    prep_raw.Im[0][0][0]=std::numeric_limits<RealType>::quiet_NaN();
    prep_metric.record_matsubara_iteration_error(
        prep_old,prep_raw,prep_old_mag,prep_raw_mag,RealType{});
    failures+=require(std::isinf(prep_metric.absolute_iteration_errors.back()),
        "nonfinite Matsubara residual cannot converge");

    // A zero-rank Gaussian has constant field and partition function.  Every
    // pCN proposal must therefore be accepted, including the zero-dimensional
    // latent-state edge case.
    p.sampling_strategy="pcn";
    p.num_RealTimeSteps=1; p.num_RealTimePoints=2;
    const contour::StaggeredLayout pcn_layout{p.num_TimePoints,p.num_RealTimeSteps*p.real_time_steps_per_interval()};
    func::ComplexDynamicMatrix zero_covariance(
        pcn_layout.dimension(),pcn_layout.dimension(),ComplexType{});
    func::DenseComplexGaussianSampler zero_sampler(zero_covariance);
    const func::MeanFieldTrajectory zero_mean(p.num_RealTimePoints,FieldVector{});
    std::mt19937 pcn_engine{12345};
    func::PCNChain constant_chain(p,zero_sampler,zero_mean,RealType{0.4},pcn_engine);
    for( size_t step=0;step<6;++step ) constant_chain.step();
    failures+=require(constant_chain.proposed()==6&&constant_chain.accepted()==6,
        "constant-likelihood pCN accepts every proposal");
    failures+=require(constant_chain.real_sampling_weight()>RealType{},
        "pCN current state has a positive real sampling weight");

    // With fixed closed-contour normalization, the pCN likelihood is Re D(T),
    // not Re Z_M. A single fluctuating response mode makes the two
    // quantities distinct while retaining a positive likelihood.
    p.correlation_normalization="closed-contour";
    func::ComplexDynamicMatrix forward_covariance(
        pcn_layout.dimension(),pcn_layout.dimension(),ComplexType{});
    const size_t forward_mode=pcn_layout.flat(
        contour::StaggeredSector::Nu,1,0);
    forward_covariance(forward_mode,forward_mode)=RealType{0.02};
    func::DenseComplexGaussianSampler forward_sampler(forward_covariance);
    std::mt19937 closed_contour_engine{24680};
    func::PCNChain closed_contour_chain(
        p,forward_sampler,zero_mean,RealType{0.4},closed_contour_engine);
    failures+=require(
        std::abs(closed_contour_chain.real_sampling_weight()
                 -std::real(closed_contour_chain.trajectory()
                                .final_closed_contour_trace))
            <RealType{1e-13},
        "closed-contour pCN uses Re D(T) as its likelihood");
    failures+=require(
        std::abs(closed_contour_chain.trajectory().final_closed_contour_trace
                 -closed_contour_chain.trajectory().partition_function)
            >RealType{1e-8},
        "closed-contour pCN test distinguishes D(T) from Z_M");
    p.correlation_normalization="partition-function";

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
    contour::CorrelationSet correlated_mean;
    rtd::MagTen correlated_mag_Re,correlated_mag_Im;
    correlated.mpi_reduce_and_finalize(
        correlated_mean,correlated_mag_Re,correlated_mag_Im);
    failures+=require(correlated.blocking_curve_mean_errors.size()>=3,
        "pCN statistics expose a multi-scale blocking curve");
    failures+=require(correlated.blocking_curve_max_errors[1]
                     >correlated.blocking_curve_max_errors[0],
        "pCN batch merging detects autocorrelation variance inflation");
    failures+=require(correlated.contour_tau_int.Re[0][0][0]>RealType{0.5},
        "pCN statistics report an autocorrelation time above the iid value");
    // Exercise the retained general-matrix measurement path independently
    // of the specialized spin-1/2 representation.
    p.sampling_strategy="independent";p.correlation_normalization="partition-function";
    p.spin_float=RealType{1.};p.num_HilbertSpaceDimension=3;
    p.num_TimeSteps=4;p.num_TimePoints=5;p.delta_t=p.beta/4;
    p.num_RealTimeSteps=4;p.num_RealTimePoints=5;p.delta_real_t=RealType{0.1};
    p.num_SamplesPerCore=2;p.num_blocks=2;
    func::initialize_matrices(p);
    const contour::StaggeredLayout spin_one_layout{5,4*p.real_time_steps_per_interval()};
    func::JointComplexGaussianSampler::FieldVector spin_one_fields(spin_one_layout.dimension());
    for(size_t i=0;i<spin_one_fields.size();++i)
        spin_one_fields[i]={0.1*std::sin(RealType(i)),0.02*std::cos(RealType(i))};
    const auto spin_one=func::compute_contour_trajectory(p,spin_one_fields,
        func::MeanFieldTrajectory(5,FieldVector{}));
    failures+=check_measurement(p,spin_one,rank,"closed-contour");
    failures+=check_measurement(p,spin_one,rank,"prefix");
    failures+=check_fused_accumulation<Operator>(p,rank);
    MPI_Finalize();
    return failures==0?0:1;
}
