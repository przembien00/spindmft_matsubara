#include"main_header.h"
#include"Functions/Staged_Solver.h"
#include"Functions/Joint_Matsubara_PCN.h"
#include<algorithm>
#include<cmath>
#include<limits>
#include<memory>
#include<stdexcept>

namespace
{
struct StageState
{
    CorrelationSet input,raw;
    MagTen re,im,raw_re,raw_im;
    explicit StageState(const ps::ParameterSpace& p)
        :raw(p.correlation_symmetry_type,p.num_TimePoints,p.num_RealTimePoints),
         raw_re(p.correlation_symmetry_type,p.num_RealTimePoints),
         raw_im(p.correlation_symmetry_type,p.num_RealTimePoints)
    {
        std::tie(re,im)=func::generate_initial_magnetization(p);
        input=func::generate_initial_correlations(p,re);
    }
};

RealType matrix_difference(const func::ComplexDynamicMatrix& a,const func::ComplexDynamicMatrix& b)
{
    RealType difference{},scale{};
    for(size_t i=0;i<a.rows();++i)for(size_t j=0;j<a.columns();++j)
    {difference+=std::norm(a(i,j)-b(i,j));scale+=std::norm(b(i,j));}
    return std::sqrt(difference)/std::max(std::sqrt(scale),RealType{1.});
}

// The two stages share propagation, estimators, MPI reduction, and stopping
// diagnostics. Preparation uses the spinDMFT_cplx_field absolute stopping rule;
// real-time feedback uses the standardized contour residual.
void run_stage(const ps::ParameterSpace& p,int rank,StageState& state,rtd::RunTimeData& data,
    tmm::DerivedTimeMeasure& clock,const StageState* preparation,
    const func::ComplexDynamicMatrix* frozen,
    std::vector<RealType>& frozen_errors)
{
    const bool real=preparation!=nullptr;
    const size_t nm=3*p.num_TimePoints;
    clock.enter_loop();
    do
    {
        clock.new_iteration();
        print::print_R0(rank,std::string{"\n--- "}+(real?"Mixed/real-time":"Matsubara")
            +" iteration "+std::to_string(data.num_Iterations+1)+" ---\n");
        auto field=p.uses_harmonic_bath()?func::prescribed_harmonic_bath_field(p,true)
            :func::self_consistent_equations(p,state.input,state.re,state.im,true);
        std::unique_ptr<func::JointComplexGaussianSampler> sampler;
        if(real)
        {
            const auto A=func::matsubara_block(field,p.num_TimePoints);
            const RealType mismatch=matrix_difference(A,*frozen);
            frozen_errors.push_back(mismatch);
            if(!std::isfinite(mismatch)||mismatch>p.covariance_tolerance)
                throw std::runtime_error("MM feedback changed despite the frozen origin");
            // Exact copy also prevents different roundoff paths from changing
            // fixed pseudo-covariance at the factorization boundary.
            for(size_t i=0;i<nm;++i)for(size_t j=0;j<nm;++j)
                field.covariance(i,j)=(*frozen)(i,j);
            // Preserve the branch-field sampling norm even when the sampled
            // coordinates are staggered eta/nu. Every mixed block is retained.
            sampler=func::make_staged_real_sampler(p,field.covariance);
        }
        else
            sampler=std::make_unique<func::DenseComplexGaussianSampler>(
                func::matsubara_block(field,p.num_TimePoints));
        if(!std::isfinite(sampler->reconstruction_error())
            ||sampler->reconstruction_error()>p.takagi_tolerance)
            throw std::runtime_error("joint Gaussian factor failed pseudo-covariance reconstruction");
        data.record_complex_field_diagnostics(field.covariance_symmetry_error,
            real?field.branch_identity_error:RealType{},sampler->reconstruction_error(),
            RealType{},sampler->latent_dimension(),sampler->largest_factorization_dimension());
        // Separate stage streams, also for the same seed and iteration index.
        std::mt19937 engine{static_cast<unsigned>(rd::throw_seed(data.generated_seed,
            2*data.num_Iterations+(real?1:0),rank))};
        clock.measure(real?"joint field distribution":"Matsubara field distribution",true);
        data.begin_iteration_accumulation();
        func::ContourTrajectory trajectory;
        func::TrajectoryWorkspace trajectory_workspace;
        func::MeasurementWorkspace measurement;
        rtd::MeasuredSample sample;
        const bool pcn=p.sampling_strategy=="pcn";
        std::unique_ptr<func::JointMatsubaraSampler> chain_sampler;
        std::unique_ptr<func::MatsubaraPCNChain> chain;
        size_t burn_accepted{},burn_proposed{},burn_zero{};
        if(pcn)
        {
            if(real)
            {
                const auto* joint=dynamic_cast<const func::DenseComplexGaussianSampler*>(sampler.get());
                if(!joint)throw std::logic_error("joint Matsubara pCN requires dense joint Takagi");
                chain_sampler=std::make_unique<func::JointMatsubaraSampler>(*joint,nm);
            }
            chain=std::make_unique<func::MatsubaraPCNChain>(p,
                real?*chain_sampler:*sampler,field.mean_time,engine);
            for(size_t i=0;i<p.mh_burn_in;++i)chain->step();
            burn_accepted=chain->accepted();burn_proposed=chain->proposed();
            burn_zero=chain->rejected_zero_weight();
        }
        if(!pcn)
        {
            // One immutable deterministic trajectory per iteration, shared by
            // all samples. Only its real-time interpolation nodes are cached.
            const func::PreparedMeanField prepared_mean(p,field.mean_time);
            const auto sample_independent=[&](auto& independent_trajectory,auto& independent_workspace)
            {
                if(real)
                {
                    // Preserve sample-major RNG order and statistical blocks.
                    constexpr size_t batch_size=32;
                    for(size_t first=0;first<p.num_SamplesPerCore;)
                    {
                        const size_t count=std::min(batch_size,p.num_SamplesPerCore-first);
                        const auto fields=sampler->draw_contour_batch(engine,count,p.uses_cf4());
                        for(const auto& field_sample:fields)
                        {
                            func::build_contour_trajectory(p,field_sample,prepared_mean,
                                independent_trajectory,independent_workspace,false,func::FinalTraceMode::Skip);
                            func::measure_contour_observables(data,independent_trajectory,sample,measurement,"prefix");
                            data.accumulate_sample(sample);
                        }
                        first+=count;
                    }
                }
                else
                    for(size_t index=0;index<p.num_SamplesPerCore;++index)
                    {
                        func::measure_matsubara_sample(p,*sampler,engine,field.mean_time,data,
                            independent_trajectory,independent_workspace,sample,measurement);
                        data.accumulate_sample(sample);
                    }
            };
            if(p.num_HilbertSpaceDimension==2)
            {
                func::SpinHalfContourTrajectory spin_half_trajectory;
                func::SpinHalfTrajectoryWorkspace spin_half_workspace;
                sample_independent(spin_half_trajectory,spin_half_workspace);
            }
            else sample_independent(trajectory,trajectory_workspace);
        }
        else for(size_t index=0;index<p.num_SamplesPerCore;++index)
        {
            chain->step();
            func::JointComplexGaussianSampler::ContourFieldSample fields;
            if(real)
            {
                // Refresh the joint factor null-space after EVERY step,
                // including rejection, retaining the M-active coordinates.
                const auto latent=chain_sampler->draw_joint_latent(chain->latent(),engine);
                fields=sampler->contour_field_from_latent(latent,p.uses_cf4());
            }
            else
            {
                const auto M=sampler->field_from_latent(chain->latent());
                fields.edge_field.resize(M.size()+6);fields.edge_field=ComplexType{};
                for(size_t i=0;i<M.size();++i)fields.edge_field[i]=M[i];
            }
            func::build_contour_trajectory(p,fields,field.mean_time,trajectory,trajectory_workspace);
            func::measure_contour_observables(data,trajectory,sample,measurement,"prefix");
            data.accumulate_sample(sample,RealType{1.}/chain->weight());
        }
        if(pcn)data.record_pcn_diagnostics(chain->accepted()-burn_accepted,
            chain->proposed()-burn_proposed,chain->rejected_zero_weight()-burn_zero,RealType{});
        clock.measure(real?"mixed real-time sampling":"Matsubara sampling",true);
        data.mpi_reduce_and_finalize(state.raw,state.raw_re,state.raw_im);
        clock.measure("MPI communication",true);
        auto active=state.raw;
        auto active_re=state.raw_re,active_im=state.raw_im;
        if(real)
        {
            // Mask frozen coordinates for residuals and updates, never for
            // output or uncertainty diagnostics.
            func::freeze_origin(active,active_re,active_im,
                preparation->input,preparation->re,preparation->im);
            if(p.constant_magnetization_time)
                for(size_t t=1;t<p.num_RealTimePoints;++t)
                {active_re[t]=preparation->re.front();active_im[t]=preparation->im.front();}
        }
        const auto residual=func::iteration_residual(state.input,active,data.contour_sample_stds,
            state.re,state.im,active_re,active_im,
            data.magnetization_time_Re_stds,data.magnetization_time_Im_stds);
        if(real)
            data.record_iteration_error(residual.absolute,residual.standardized);
        else
            data.record_matsubara_iteration_error(state.input,active,state.re,active_re,
                residual.standardized);
        if(p.uses_harmonic_bath())
        {state.input=active;state.re=active_re;state.im=active_im;}
        else
        {
            state.input=func::mix_correlations(state.input,active,p.mixing_alpha);
            state.re=func::mix_magnetization_tensor(state.re,active_re,p.mixing_alpha);
            state.im=func::mix_magnetization_tensor(state.im,active_im,p.mixing_alpha);
        }
        if(real)func::freeze_origin(state.input,state.re,state.im,
            preparation->input,preparation->re,preparation->im);
        data.finalize_iteration_step();
        clock.measure("iteration-step finalization",true);
        // Assess imaginary-magnetization diagnostics on the measured result.
        // Substituting the frozen origin here could hide a changed ensemble.
    }while(!data.terminate(state.raw_im));
    clock.leave_loop();
    clock.measure(real?"real-time self consistency":"Matsubara self consistency",true);
}
}

int main(int argc,char* const argv[])
{
    MPI_Init(nullptr,nullptr);
    int world{},rank{};MPI_Comm_size(MPI_COMM_WORLD,&world);MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    try
    {
        tmm::DerivedTimeMeasure clock(rank,world);
        const ps::ParameterSpace parameters(argc,argv,world,rank);
        if(parameters.gaussian_factorization!="dense"
            ||parameters.spin_insertion_strategy!="prefix"
            ||parameters.correlation_normalization!="partition-function")
            throw std::invalid_argument("staged solver requires dense factorization, prefix insertions, and partition-function normalization");
        if(!parameters.self_consistency&&!parameters.uses_harmonic_bath())
            throw std::invalid_argument("noselfcons cannot prepare a converged self-consistent Matsubara bath; use bath=harmonic for prescribed-bath validation");
        if(parameters.fft_cross_frequency_cutoff>=RealType{})
            throw std::invalid_argument("the staged dense sampler does not support an FFT cutoff");
        print::print_R0(rank,parameters.create_essentials_string());
        func::initialize_matrices(parameters);
        auto pM=parameters;
        pM.num_RealTimeSteps=0;pM.num_RealTimePoints=1;pM.Tmax=RealType{};
        pM.num_SamplesPerCore=parameters.preparation_samples?parameters.preparation_samples:parameters.num_SamplesPerCore;
        if(pM.num_SamplesPerCore>std::numeric_limits<size_t>::max()/parameters.get_num_Cores())
            throw std::invalid_argument("preparation sample count overflows size_t");
        pM.num_Samples=pM.num_SamplesPerCore*pM.get_num_Cores();
        pM.Iteration_Limit=parameters.preparation_iteration_limit?parameters.preparation_iteration_limit:parameters.Iteration_Limit;
        StageState preparation(pM);
        rtd::RunTimeData prep_data(pM,rank,true);
        prep_data.generated_seed=rd::generate_seed(prep_data.get_seed_str(),rank);
        clock.measure("initialization",true);
        std::vector<RealType> frozen_errors;
        run_stage(pM,rank,preparation,prep_data,clock,nullptr,nullptr,frozen_errors);
        if(prep_data.termination!="by convergence"&&!parameters.uses_harmonic_bath())
        {
            auto failed=pM;failed.filename_extension+="__preparation_failed";
            stoc::HDF5_Storage storage(rank,failed,"by iteration limit");
            storage.store_main(failed,prep_data,preparation.raw,preparation.raw_re,preparation.raw_im);
            storage.store_preparation(pM,prep_data,preparation.input,preparation.re,preparation.im,
                preparation.raw,preparation.raw_re,preparation.raw_im,{},{},frozen_errors,false);
            clock.stop();storage.store_time(clock);storage.finalize();
            print::print_R0(rank,"Matsubara preparation did not converge; real-time stage was not started.\n");
            MPI_Finalize();return 2;
        }
        const auto prepared_field=parameters.uses_harmonic_bath()
            ?func::prescribed_harmonic_bath_field(pM,true)
            :func::self_consistent_equations(pM,preparation.input,preparation.re,preparation.im,true);
        const auto frozen=func::matsubara_block(prepared_field,pM.num_TimePoints);
        // Imaginary-only reference; the real loop always uses the current joint law.
        const auto reference=func::autonne_takagi(frozen);
        if(!std::isfinite(reference.reconstruction_error)
            ||reference.reconstruction_error>parameters.takagi_tolerance)
            throw std::runtime_error("prepared Gaussian factor failed pseudo-covariance reconstruction");
        const auto& preparation_factor=reference.L;
        StageState real(parameters);
        func::seed_from_preparation(parameters,preparation.input,preparation.re,preparation.im,
            real.input,real.re,real.im);
        rtd::RunTimeData real_data(parameters,rank);real_data.generated_seed=prep_data.generated_seed;
        print::print_R0(rank,std::string{"\nMatsubara preparation frozen. Starting "}
            +(parameters.sampling_strategy=="pcn"?"joint-Takagi Matsubara pCN":"joint-Takagi")
            +" mixed/real-time iteration.\n");
        run_stage(parameters,rank,real,real_data,clock,&preparation,&frozen,frozen_errors);
        stoc::HDF5_Storage storage(rank,parameters,real_data.termination);
        storage.store_main(parameters,real_data,real.raw,real.raw_re,real.raw_im);
        storage.store_preparation(pM,prep_data,preparation.input,preparation.re,preparation.im,
            preparation.raw,preparation.raw_re,preparation.raw_im,frozen,preparation_factor,frozen_errors,true);
        clock.measure("storing",true);clock.stop();storage.store_time(clock);storage.finalize();
        MPI_Finalize();return real_data.termination=="by iteration limit"?3:0;
    }
    catch(const std::exception& exception)
    {
        std::cerr<<"Rank "<<rank<<": staged solver failed: "<<exception.what()<<'\n';
        MPI_Abort(MPI_COMM_WORLD,1);return 1;
    }
}
