#include "main_header.h"
#include <mpi.h>


#include <sys/resource.h>
#include <iomanip>
extern std::array<double,4> audit_profile;
static void audit_probe(const Sampler& original,const Model& model,const ps::ParameterSpace& p,const std::vector<C>& mean) {
    std::cout << std::setprecision(12);
    Sampler warm=original; std::mt19937 warm_rng(419);
    warm.draw_contour_batch(warm_rng,64); // Construct shared batch-factor storage before timing.
    std::vector<Field> reference;
    for(size_t batch_size:{size_t(1),size_t(8),size_t(32),size_t(64)}) {
        Sampler sampler=original;
        std::mt19937 rng(613);
        double elapsed=0.,difference=0.;
        for(size_t offset=0;offset<128;offset+=batch_size) {
            const double start=MPI_Wtime();
            auto batch=sampler.draw_contour_batch(rng,batch_size);
            elapsed+=MPI_Wtime()-start;
            for(size_t j=0;j<batch_size;++j) {
                if(batch_size==1) reference.push_back(std::move(batch[j].edge_field));
                else for(size_t k=0;k<reference[offset+j].size();++k)
                    difference=std::max(difference,std::abs(reference[offset+j][k]-batch[j].edge_field[k]));
            }
        }
        std::cout << "AUDIT batch="<<batch_size<<" draw_seconds_per_sample="<<elapsed/128<<" max_field_difference="<<difference<<"\n";
    }
    audit_profile.fill(0.);
    double checksum=0.;
    const double start=MPI_Wtime();
    for(size_t i=0;i<16;++i) {
        const auto sample=func::measure<ContourCorrelations>(model,p,p.real_time_layout,reference[i],mean);
        checksum+=std::abs(sample.partition)+std::abs(sample.correlations.values.back());
    }
    std::cout<<"AUDIT measure_seconds_per_sample="<<(MPI_Wtime()-start)/16<<" checksum="<<checksum<<"\n";
    std::cout<<"AUDIT imaginary_ms="<<1000*audit_profile[0]/16<<" insertions_ms="<<1000*audit_profile[1]/16<<" real_propagation_ms="<<1000*audit_profile[2]/16<<" observables_ms="<<1000*audit_profile[3]/16<<"\n";
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::cout<<"AUDIT dimension="<<original.size()<<" latent="<<original.latent_dimension()<<" largest_block="<<original.largest_factorization_dimension()<<" peak_rss_platform_units="<<usage.ru_maxrss<<"\n";
}

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0,ranks=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    const ps::ParameterSpace parameters(argc,argv,ranks,rank);
    tmm::DerivedTimeMeasure clock(rank);
    const Model model=mfm::read_model(parameters);
    print::print_R0(rank,parameters.create_essentials_string(model));
    rtd::RunTimeData runtime(parameters,model);
    MatsubaraCorrelations preparation_correlations=func::generate_initial_correlations(model,parameters);
    ClusterMagnetization preparation_magnetization=func::generate_initial_magnetization(model,parameters);
    auto next_preparation_correlations=preparation_correlations;
    auto next_preparation_magnetization=preparation_magnetization;
    const bool pcn=parameters.sampling_strategy=="pcn";
    clock.measure("initialization");
    while(runtime.preparation.termination.empty()) {
        clock.new_iteration();
        print::print_R0(rank,"\n--- Matsubara iteration "
            +std::to_string(runtime.preparation.num_Iterations+1)+" ---\n");
        const size_t samples=parameters.num_PreparationSamplesPerCore;
        // Construct the current field law before drawing any trajectories.
        const auto covariance=mfm::self_consistent_equations(model,parameters,
            preparation_correlations,preparation_magnetization,parameters.preparation_layout);
        Sampler sampler(covariance,parameters.takagi_tolerance);
        const auto mean=parameters.uses_harmonic_bath()?std::vector<C>(model.operators.size()):mfm::mean_field(model,preparation_magnetization);
        std::seed_seq seed{static_cast<unsigned>(parameters.generated_seed),static_cast<unsigned>(parameters.generated_seed>>32),
            static_cast<unsigned>(rank),static_cast<unsigned>(runtime.preparation.num_Iterations),0u};
        std::mt19937 engine(seed);
        clock.measure("Matsubara field distribution");
        std::unique_ptr<func::MatsubaraPCNChain> chain;
        if(pcn) {
            chain=std::make_unique<func::MatsubaraPCNChain>(model,parameters,sampler,mean,engine);
            for(size_t burn=0;burn<parameters.mh_burn_in;++burn) chain->step();
            chain->reset_statistics();
        }
        func::initialize_sampling(next_preparation_correlations,next_preparation_magnetization,runtime.preparation,
            samples,parameters.num_blocks,pcn);
        timing::DurationEstimator progress(rank,pcn?"pCN simulation":"Monte-Carlo simulation",
            {samples},{"sampling"},true);
        std::vector<Sampler::ContourFieldSample> batch;
        for(size_t index=0;index<samples;++index) {
            const double sample_start=MPI_Wtime();
            Field field;
            if(pcn) {
                chain->step();
                field=chain->draw_field();
            } else {
                if(index%parameters.num_SamplesPerSet==0)
                    batch=sampler.draw_contour_batch(engine,std::min(parameters.num_SamplesPerSet,samples-index));
                field=std::move(batch[index%parameters.num_SamplesPerSet].edge_field);
            }
            const auto sample=func::measure<MatsubaraCorrelations>(model,parameters,parameters.preparation_layout,field,mean);
            func::accumulate_sample(next_preparation_correlations,next_preparation_magnetization,sample,pcn?1./chain->weight():1.);
            func::record_sample_statistics(runtime.preparation,sample,pcn?1./chain->weight():1.);
            progress.obtain(timing::DurationQuantity{"sampling",MPI_Wtime()-sample_start});
        }
        clock.measure("Matsubara sampling");
        // Reduce paired numerator/partition sums and estimate their errors.
        func::reduce_and_normalize(next_preparation_correlations,next_preparation_magnetization,runtime.preparation,ranks);
        runtime.preparation.factorization_error=sampler.reconstruction_error();
        runtime.preparation.latent_dimension=sampler.latent_dimension();
        runtime.preparation.largest_factorization_dimension=sampler.largest_factorization_dimension();
        double rates[2]{pcn?chain->acceptance_rate():0.,pcn?chain->zero_modulus_rejection_rate():0.};
        MPI_Allreduce(MPI_IN_PLACE,rates,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        runtime.preparation.acceptance_rate=rates[0]/double(ranks);
        runtime.preparation.zero_modulus_rejection_rate=rates[1]/double(ranks);
        clock.measure("MPI communication");
        runtime.record(parameters,model,preparation_correlations,preparation_magnetization,
            next_preparation_correlations,next_preparation_magnetization,rank);
        std::swap(preparation_correlations,next_preparation_correlations);
        std::swap(preparation_magnetization,next_preparation_magnetization);
        clock.measure("iteration-step finalization");
        clock.finish_iteration(true);
    }
    clock.finish_stage(true);
    if(runtime.preparation.termination=="by iteration limit" || runtime.preparation.termination=="by hand") {
        stoc::HDF5_Storage storage(rank,parameters,runtime.preparation.termination,"__preparation_failed");
        storage.store_main(parameters,model,runtime,preparation_correlations,preparation_magnetization);
        storage.store_preparation(parameters,model,runtime,preparation_correlations,preparation_magnetization);
        clock.measure("storing"); clock.stop(); storage.store_time(clock); storage.finalize();
        print::print_R0(rank,"Matsubara preparation did not converge; real-time stage was not started.\n");
        MPI_Finalize(); return 2;
    }
    ContourCorrelations correlations(parameters.num_RealTimeSteps+1,model.sites,parameters.num_TimeSteps+1,parameters.correlation_symmetry_type,model.correlation_categories);
    TimeMagnetization magnetization(parameters.num_RealTimeSteps+1,model.sites,parameters.correlation_symmetry_type);
    func::seed_from_preparation(preparation_correlations,preparation_magnetization,correlations,magnetization);
    auto next_correlations=correlations;
    auto next_magnetization=magnetization;
    print::print_R0(rank,std::string{"\nMatsubara preparation frozen. Starting "}
        +(parameters.sampling_strategy=="pcn"?"joint-Takagi Matsubara pCN":"joint-Takagi")
        +" mixed/real-time iteration.\n");
    while(runtime.real_time.termination.empty()) {
        clock.new_iteration();
        print::print_R0(rank,"\n--- Mixed/real-time iteration "
            +std::to_string(runtime.real_time.num_Iterations+1)+" ---\n");
        func::constrain_feedback(correlations,magnetization,preparation_correlations,preparation_magnetization,parameters);
        const size_t samples=parameters.num_SamplesPerCore;
        // Construct the current field law before drawing any trajectories.
        const auto covariance=mfm::self_consistent_equations(model,parameters,correlations,magnetization,parameters.real_time_layout);
        Sampler sampler(covariance,parameters.real_time_layout.ni*parameters.real_time_layout.components,
            parameters.real_time_layout.intervals*parameters.real_time_layout.components,parameters.takagi_tolerance);
        const auto mean=parameters.uses_harmonic_bath()?std::vector<C>((parameters.num_RealTimeSteps+1)*model.operators.size()):mfm::mean_field(model,magnetization);
        std::seed_seq seed{static_cast<unsigned>(parameters.generated_seed),static_cast<unsigned>(parameters.generated_seed>>32),
            static_cast<unsigned>(rank),static_cast<unsigned>(runtime.real_time.num_Iterations),1u};
        std::mt19937 engine(seed);
        clock.measure("joint field distribution");
        audit_probe(sampler,model,parameters,mean);
        std::unique_ptr<func::MatsubaraPCNChain> chain;
        if(pcn) {
            chain=std::make_unique<func::MatsubaraPCNChain>(model,parameters,sampler,mean,engine);
            for(size_t burn=0;burn<parameters.mh_burn_in;++burn) chain->step();
            chain->reset_statistics();
        }
        func::initialize_sampling(next_correlations,next_magnetization,runtime.real_time,samples,parameters.num_blocks,pcn);
        timing::DurationEstimator progress(rank,pcn?"pCN simulation":"Monte-Carlo simulation",
            {samples},{"sampling"},true);
        std::vector<Sampler::ContourFieldSample> batch;
        for(size_t index=0;index<samples;++index) {
            const double sample_start=MPI_Wtime();
            Field field;
            if(pcn) {
                chain->step();
                field=chain->draw_field();
            } else {
                if(index%parameters.num_SamplesPerSet==0)
                    batch=sampler.draw_contour_batch(engine,std::min(parameters.num_SamplesPerSet,samples-index));
                field=std::move(batch[index%parameters.num_SamplesPerSet].edge_field);
            }
            const auto sample=func::measure<ContourCorrelations>(model,parameters,parameters.real_time_layout,field,mean);
            func::accumulate_sample(next_correlations,next_magnetization,sample,pcn?1./chain->weight():1.);
            func::record_sample_statistics(runtime.real_time,sample,pcn?1./chain->weight():1.);
            progress.obtain(timing::DurationQuantity{"sampling",MPI_Wtime()-sample_start});
        }
        clock.measure("mixed real-time sampling");
        // Reduce paired numerator/partition sums and estimate their errors.
        func::reduce_and_normalize(next_correlations,next_magnetization,runtime.real_time,ranks);
        runtime.real_time.factorization_error=sampler.reconstruction_error();
        runtime.real_time.latent_dimension=sampler.latent_dimension();
        runtime.real_time.largest_factorization_dimension=sampler.largest_factorization_dimension();
        double rates[2]{pcn?chain->acceptance_rate():0.,pcn?chain->zero_modulus_rejection_rate():0.};
        MPI_Allreduce(MPI_IN_PLACE,rates,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        runtime.real_time.acceptance_rate=rates[0]/double(ranks);
        runtime.real_time.zero_modulus_rejection_rate=rates[1]/double(ranks);
        clock.measure("MPI communication");
        runtime.record(parameters,model,correlations,magnetization,next_correlations,next_magnetization,rank);
        std::swap(correlations,next_correlations);
        std::swap(magnetization,next_magnetization);
        clock.measure("iteration-step finalization");
        clock.finish_iteration(false);
    }
    clock.finish_stage(false);
    stoc::HDF5_Storage storage(rank,parameters,runtime.real_time.termination);
    storage.store_main(parameters,model,runtime,correlations,magnetization);
    storage.store_preparation(parameters,model,runtime,preparation_correlations,preparation_magnetization);
    clock.measure("storing"); clock.stop(); storage.store_time(clock); storage.finalize();
    const int result=runtime.real_time.termination=="by iteration limit"?3:0;
    MPI_Finalize(); return result;
}
