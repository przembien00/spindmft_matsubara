#include"main_header.h"
#include"Functions/Staged_Solver.h"
#include"Functions/Joint_Matsubara_PCN.h"
#include<algorithm>
#include<cmath>
#include<memory>
#include<stdexcept>

int main(int argc,char* const argv[])
{
    MPI_Init(nullptr,nullptr);
    int world{},rank{};MPI_Comm_size(MPI_COMM_WORLD,&world);MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    tmm::DerivedTimeMeasure clock(rank,world);
    const ps::ParameterSpace parameters(argc,argv,world,rank);
    print::print_R0(rank,parameters.create_essentials_string());
    func::initialize_matrices(parameters);
    auto preparation_magnetization=func::generate_initial_magnetization(parameters,true);
    auto preparation_correlations=func::generate_initial_correlations(parameters,preparation_magnetization.Re);
    rtd::RunTimeData runtime(parameters,rank);
    runtime.generated_seed=rd::generate_seed(runtime.get_seed_str(),rank);
    clock.measure("initialization",true);

    // Self-consistency loop for imaginary-time correlations.
    clock.enter_loop();
    do
    {
        clock.new_iteration();
        const auto previous_correlations=preparation_correlations;
        const auto previous_magnetization=preparation_magnetization;
        print::print_R0(rank,"\n--- Matsubara iteration "
            +std::to_string(runtime.preparation.num_Iterations+1)+" ---\n");
        const auto field=parameters.uses_harmonic_bath()?func::prescribed_harmonic_bath_field(parameters,true)
            :func::self_consistent_equations(parameters,preparation_correlations,
                preparation_magnetization);
        func::DenseComplexGaussianSampler sampler(
            field.covariance,parameters.takagi_tolerance);
        runtime.begin_iteration_accumulation(preparation_correlations,preparation_magnetization,true);
        runtime.record_complex_field_diagnostics(field.covariance_symmetry_error,
            sampler.reconstruction_error(),
            RealType{},sampler.latent_dimension(),sampler.largest_factorization_dimension());
        // Even stream indices are reserved for Matsubara preparation.
        std::mt19937 engine{static_cast<unsigned>(rd::throw_seed(runtime.generated_seed,
            2*runtime.preparation.num_Iterations,rank))};
        clock.measure("Matsubara field distribution",true);
        func::MeasurementWorkspace measurement;
        func::MeasuredSample sample;
        const bool pcn=parameters.sampling_strategy=="pcn";
        std::unique_ptr<func::MatsubaraPCNChain> chain;
        size_t burn_accepted{},burn_proposed{},burn_zero{};
        if(pcn)
        {
            chain=std::make_unique<func::MatsubaraPCNChain>(parameters,sampler,field.mean_time,engine);
            for(size_t i=0;i<parameters.mh_burn_in;++i)chain->step();
            burn_accepted=chain->accepted();burn_proposed=chain->proposed();
            burn_zero=chain->rejected_zero_weight();
        }
        if(!pcn)
        {
            const auto sample_independent=[&](auto& independent_trajectory,auto& independent_workspace)
            {
                for(size_t index=0;index<parameters.num_PreparationSamplesPerCore;++index)
                {
                    func::measure_matsubara_sample(parameters,sampler,engine,field.mean_time,runtime,
                        independent_trajectory,independent_workspace,sample,measurement);
                    runtime.accumulate_sample(sample);
                }
            };
            if(parameters.num_HilbertSpaceDimension==2)
            {
                func::SpinHalfContourTrajectory spin_half_trajectory;
                func::SpinHalfTrajectoryWorkspace spin_half_workspace;
                sample_independent(spin_half_trajectory,spin_half_workspace);
            }
            else
            {
                func::ContourTrajectory trajectory;
                func::TrajectoryWorkspace trajectory_workspace;
                sample_independent(trajectory,trajectory_workspace);
            }
        }
        else
        {
            func::ContourTrajectory trajectory;
            func::TrajectoryWorkspace trajectory_workspace;
            for(size_t index=0;index<parameters.num_PreparationSamplesPerCore;++index)
            {
                chain->step();
                const auto M=sampler.field_from_latent(chain->latent());
                func::JointComplexGaussianSampler::ContourFieldSample fields;
                fields.edge_field.resize(M.size()+6);fields.edge_field=ComplexType{};
                for(size_t i=0;i<M.size();++i)fields.edge_field[i]=M[i];
                func::build_matsubara_trajectory(parameters,fields,field.mean_time,trajectory,trajectory_workspace);
                func::measure_contour_observables(runtime,trajectory,sample,measurement);
                runtime.accumulate_sample(sample,RealType{1.}/chain->weight());
            }
        }
        if(pcn)runtime.record_pcn_diagnostics(chain->accepted()-burn_accepted,
            chain->proposed()-burn_proposed,chain->rejected_zero_weight()-burn_zero);
        clock.measure("Matsubara sampling",true);
        runtime.mpi_reduce_and_normalize();
        clock.measure("MPI communication",true);
        const auto residual=func::iteration_residual(previous_correlations,preparation_correlations,
            runtime.preparation.contour_sample_stds,previous_magnetization,
            preparation_magnetization,
            runtime.preparation.magnetization_time_stds);
        runtime.record_matsubara_iteration_error(previous_correlations,preparation_correlations,
            previous_magnetization.Re,preparation_magnetization.Re,residual.standardized);
        runtime.finalize_iteration_step();
        clock.measure("iteration-step finalization",true);
    }while(!runtime.terminate(preparation_magnetization.Im));
    clock.leave_loop();
    clock.measure("Matsubara self consistency",true);

    // Start real-time evolution only after successful preparation.
    if(runtime.preparation.termination!="by convergence"&&!parameters.uses_harmonic_bath())
    {
        stoc::HDF5_Storage storage(rank,parameters,"by iteration limit","__preparation_failed");
        storage.store_main(parameters,runtime,preparation_correlations,
            preparation_magnetization,true);
        storage.store_preparation(parameters,runtime,preparation_correlations,
            preparation_magnetization,false);
        clock.stop();storage.store_time(clock);storage.finalize();
        print::print_R0(rank,"Matsubara preparation did not converge; real-time stage was not started.\n");
        MPI_Finalize();return 2;
    }
    const auto prepared_field=parameters.uses_harmonic_bath()
        ?func::prescribed_harmonic_bath_field(parameters,true)
        :func::self_consistent_equations(parameters,preparation_correlations,
            preparation_magnetization);
    const auto& frozen=prepared_field.covariance;
    CorrelationSet correlations;
    func::Magnetization magnetization;
    func::seed_from_preparation(parameters,preparation_correlations,
        preparation_magnetization,
        correlations,magnetization);
    print::print_R0(rank,std::string{"\nMatsubara preparation frozen. Starting "}
        +(parameters.sampling_strategy=="pcn"?"joint-Takagi Matsubara pCN":"joint-Takagi")
        +" mixed/real-time iteration.\n");

    // Self-consistency loop for mixed and real-time correlations.
    const size_t nm=3*parameters.num_TimePoints;
    clock.enter_loop();
    do
    {
        clock.new_iteration();
        print::print_R0(rank,"\n--- Mixed/real-time iteration "
            +std::to_string(runtime.real_time.num_Iterations+1)+" ---\n");
        // Keep the measured quantities for output; constrain only the feedback.
        auto previous_correlations=correlations;
        auto previous_magnetization=magnetization;
        func::freeze_origin(previous_correlations,previous_magnetization,
            preparation_correlations,preparation_magnetization);
        if(parameters.constant_magnetization_time)
            for(size_t t=1;t<parameters.num_RealTimePoints;++t)
            {
                previous_magnetization.Re[t]=preparation_magnetization.Re.front();
                previous_magnetization.Im[t]=preparation_magnetization.Im.front();
            }
        auto field=parameters.uses_harmonic_bath()?func::prescribed_harmonic_bath_field(parameters)
            :func::self_consistent_equations(parameters,previous_correlations,
                previous_magnetization);
        const auto A=func::matsubara_block(field,parameters.num_TimePoints);
        const RealType mismatch=func::matrix_difference(A,frozen);
        runtime.frozen_errors.push_back(mismatch);
        if(!std::isfinite(mismatch)||mismatch>parameters.covariance_tolerance)
            throw std::runtime_error("MM feedback changed despite the frozen origin");
        // Preserve the branch-field sampling norm even when the sampled
        // coordinates are staggered eta/nu. Every mixed block is retained.
        const auto sampler=func::make_staged_real_sampler(parameters,field.covariance);
        runtime.begin_iteration_accumulation(correlations,magnetization);
        runtime.record_complex_field_diagnostics(field.covariance_symmetry_error,
            sampler->reconstruction_error(),
            RealType{},sampler->latent_dimension(),sampler->largest_factorization_dimension());
        // Odd stream indices keep real-time sampling separate from preparation.
        std::mt19937 engine{static_cast<unsigned>(rd::throw_seed(runtime.generated_seed,
            2*runtime.real_time.num_Iterations+1,rank))};
        clock.measure("joint field distribution",true);
        func::MeasurementWorkspace measurement;
        func::MeasuredSample sample;
        const bool pcn=parameters.sampling_strategy=="pcn";
        std::unique_ptr<func::JointMatsubaraSampler> chain_sampler;
        std::unique_ptr<func::MatsubaraPCNChain> chain;
        size_t burn_accepted{},burn_proposed{},burn_zero{};
        if(pcn)
        {
            const auto* joint=dynamic_cast<const func::DenseComplexGaussianSampler*>(sampler.get());
            if(!joint)throw std::logic_error("joint Matsubara pCN requires dense joint Takagi");
            chain_sampler=std::make_unique<func::JointMatsubaraSampler>(*joint,nm);
            chain=std::make_unique<func::MatsubaraPCNChain>(parameters,*chain_sampler,field.mean_time,engine);
            for(size_t i=0;i<parameters.mh_burn_in;++i)chain->step();
            burn_accepted=chain->accepted();burn_proposed=chain->proposed();
            burn_zero=chain->rejected_zero_weight();
        }
        if(!pcn)
        {
            // One immutable deterministic trajectory per iteration, shared by
            // all samples. Only its real-time interpolation nodes are cached.
            const func::PreparedMeanField prepared_mean(parameters,field.mean_time);
            const auto sample_independent=[&](auto& independent_trajectory,auto& independent_workspace)
            {
                // Preserve sample-major RNG order and statistical blocks.
                constexpr size_t batch_size=32;
                for(size_t first=0;first<parameters.num_SamplesPerCore;)
                {
                    const size_t count=std::min(batch_size,parameters.num_SamplesPerCore-first);
                    const auto fields=sampler->draw_contour_batch(engine,count);
                    for(const auto& field_sample:fields)
                    {
                        func::build_contour_trajectory(parameters,field_sample,prepared_mean,
                            independent_trajectory,independent_workspace,false,func::FinalTraceMode::Skip);
                        func::measure_and_accumulate_independent(runtime,independent_trajectory,measurement);
                    }
                    first+=count;
                }
            };
            if(parameters.num_HilbertSpaceDimension==2)
            {
                func::SpinHalfContourTrajectory spin_half_trajectory;
                func::SpinHalfTrajectoryWorkspace spin_half_workspace;
                sample_independent(spin_half_trajectory,spin_half_workspace);
            }
            else
            {
                func::ContourTrajectory trajectory;
                func::TrajectoryWorkspace trajectory_workspace;
                sample_independent(trajectory,trajectory_workspace);
            }
        }
        else
        {
            const func::PreparedMeanField prepared_mean(parameters,field.mean_time);
            func::ContourTrajectory trajectory;
            func::TrajectoryWorkspace trajectory_workspace;
            for(size_t index=0;index<parameters.num_SamplesPerCore;++index)
            {
                chain->step();
                // Refresh the joint factor null-space after EVERY step,
                // including rejection, retaining the M-active coordinates.
                const auto latent=chain_sampler->draw_joint_latent(chain->latent(),engine);
                const auto fields=sampler->contour_field_from_latent(latent);
                func::build_contour_trajectory(parameters,fields,prepared_mean,trajectory,trajectory_workspace);
                func::measure_contour_observables(runtime,trajectory,sample,measurement);
                runtime.accumulate_sample(sample,RealType{1.}/chain->weight());
            }
        }
        if(pcn)runtime.record_pcn_diagnostics(chain->accepted()-burn_accepted,
            chain->proposed()-burn_proposed,chain->rejected_zero_weight()-burn_zero);
        clock.measure("mixed real-time sampling",true);
        runtime.mpi_reduce_and_normalize();
        clock.measure("MPI communication",true);
        // Frozen coordinates do not contribute to the iteration residual.
        const auto residual=func::iteration_residual(previous_correlations,correlations,
            runtime.real_time.contour_sample_stds,previous_magnetization,
            magnetization,
            runtime.real_time.magnetization_time_stds,
            true,parameters.constant_magnetization_time);
        runtime.record_iteration_error(residual.absolute,residual.standardized);
        runtime.finalize_iteration_step();
        clock.measure("iteration-step finalization",true);
        // Assess imaginary-magnetization diagnostics on the measured result.
        // Substituting the frozen origin here could hide a changed ensemble.
    }while(!runtime.terminate(magnetization.Im));
    clock.leave_loop();
    clock.measure("real-time self consistency",true);

    stoc::HDF5_Storage storage(rank,parameters,runtime.real_time.termination);
    storage.store_main(parameters,runtime,correlations,magnetization);
    storage.store_preparation(parameters,runtime,preparation_correlations,
        preparation_magnetization,true);
    clock.measure("storing",true);clock.stop();storage.store_time(clock);storage.finalize();
    MPI_Finalize();return runtime.real_time.termination=="by iteration limit"?3:0;
}
