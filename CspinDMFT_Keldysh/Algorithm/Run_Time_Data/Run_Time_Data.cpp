#include "../main_header.h"
#include <filesystem>
#include <iostream>
#include <mpi.h>
namespace CspinDMFT::Run_Time_Data {
template<class Correlations>
IterationResidual iteration_residual(const Correlations& old,const MagnetizationFor<Correlations>& old_moments,
    const Correlations& next,const MagnetizationFor<Correlations>& next_moments,
    const BasicStageData<Correlations>& statistics,bool constant) {
    constexpr bool frozen=std::is_same_v<Correlations,ContourCorrelations>;
    IterationResidual out;
    const auto update=[&](C previous,C current,C re_error,C im_error) {
        const C diff=current-previous;
        out.absolute=std::max(out.absolute,std::abs(diff));
        const double floor=128.*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(previous),std::abs(current)});
        for(auto pair:{std::pair<double,double>{std::abs(diff.real()),re_error.real()},
                       std::pair<double,double>{std::abs(diff.imag()),im_error.real()}}) {
            const double ratio=pair.first<=floor?0.:pair.second>floor?pair.first/pair.second:std::numeric_limits<double>::infinity();
            out.standardized=std::max(out.standardized,ratio);
        }
    };
    for(size_t t=frozen?1:0;t<time_points(old);++t) {
        for(size_t i=0;i<old.sites;++i) if(!constant || !frozen) for(size_t a:old_moments.directions)
            update(magnetization(old_moments,t,3*i+a),magnetization(next_moments,t,3*i+a),
                magnetization(statistics.magnetization_error_re,t,3*i+a),magnetization(statistics.magnetization_error_im,t,3*i+a));
        for(auto ij:old.categories) for(auto ab:old.directions) for(size_t k=0;k<old.imaginary_points;++k) {
            const auto a=3*ij[0]+ab[0],b=3*ij[1]+ab[1];
            update(correlation(old,t,a,b,k),correlation(next,t,a,b,k),
                correlation(statistics.correlation_error_re,t,a,b,k),correlation(statistics.correlation_error_im,t,a,b,k));
        }
    }
    return out;
}
double RunTimeData::preparation_threshold(const ParameterSpace& p,const Model& model) const {
    if(p.absolute_iteration_error_option!="conservative") return p.absolute_iteration_error_tolerance;
    double scale=0.;
    for(double spin:model.spins) scale=std::max(scale,spin*(spin+1.)/3.);
    return p.iteration_error_sigma_threshold*scale/std::sqrt(3.*double(p.num_PreparationSamplesPerCore*p.get_num_Cores()));
}
template<class Correlations>
bool RunTimeData::record(const ParameterSpace& p,const Model& model,const Correlations& old,
    const MagnetizationFor<Correlations>& old_moments,const Correlations& current,
    const MagnetizationFor<Correlations>& current_moments,int rank) {
    constexpr bool prep=std::is_same_v<Correlations,MatsubaraCorrelations>;
    auto& statistics=stage<Correlations>();
    auto error=iteration_residual(old,old_moments,current,current_moments,statistics,p.constant_magnetization_time);
    statistics.absolute_iteration_errors.push_back(error.absolute); statistics.standardized_iteration_errors.push_back(error.standardized);
    statistics.average_phase_magnitudes.push_back(statistics.average_phase);
    statistics.gaussian_factor_reconstruction_errors.push_back(statistics.factorization_error);
    statistics.covariance_symmetry_errors.push_back(statistics.covariance_error);
    statistics.partition_sum_Re.push_back(statistics.partition.real()*double(statistics.samples_per_core*p.get_num_Cores()));
    statistics.partition_sum_Im.push_back(statistics.partition.imag()*double(statistics.samples_per_core*p.get_num_Cores()));
    statistics.mh_acceptance_rates.push_back(statistics.acceptance_rate);
    statistics.effective_sample_sizes.push_back(statistics.effective_sample_size);
    statistics.sample_counts.push_back(double(statistics.samples_per_core));
    ++statistics.num_Iterations;
    double ceiling=0.;
    for(double s:model.spins) ceiling=std::max(ceiling,std::max(s,s*s));
    bool diagnostics=statistics.average_phase>=p.minimum_phase_magnitude
        && (!p.check_gaussian_reconstruction || statistics.factorization_error<=p.takagi_tolerance)
        && statistics.covariance_error<=p.covariance_tolerance;
    for(size_t i=0;i<model.sites;++i) for(auto a:magnetization_directions(p.correlation_symmetry_type)) {
        double floor=128.*std::numeric_limits<double>::epsilon();
        if(std::abs(magnetization(current_moments,0,3*i+a).imag())>p.iteration_error_sigma_threshold*magnetization(statistics.magnetization_error_im,0,3*i+a).real()+floor) diagnostics=false;
    }
    int stop_requested=rank==0 && std::filesystem::exists("terminate_"+std::to_string(p.Run_ID));
    MPI_Bcast(&stop_requested,1,MPI_INT,0,MPI_COMM_WORLD);
    bool fixed=prep || p.iteration_error_mode=="absolute"
        ?error.absolute<=preparation_threshold(p,model)
        :error.standardized<=p.iteration_error_sigma_threshold && error.absolute<ceiling;
    if(!p.self_consistency) statistics.termination=p.uses_harmonic_bath()?"prescribed harmonic bath":"no self consistency";
    else if(fixed && diagnostics) statistics.termination="by convergence";
    else if(stop_requested) statistics.termination="by hand";
    else if(statistics.num_Iterations>=(prep?p.Preparation_Iteration_Limit:p.Iteration_Limit)) statistics.termination="by iteration limit";
    double closure_residual=0.;
    for(C value:statistics.closure) closure_residual=std::max(closure_residual,std::abs(value-C(1.)));
    const auto format=[&](double value) { return print::round_value_to_string(value,p.num_PrintDigits); };
    print::print_R0(rank,std::string{"Iteration step finished:\n"}
        +print::quantity_to_output_line(36,p.uses_harmonic_bath()?"initial-guess difference (unused)":prep
            ?"raw Matsubara absolute residual":"raw complex fixed-point residual",format(error.absolute))
        +print::quantity_to_output_line(36,"largest standardized residual",format(error.standardized))
        +print::quantity_to_output_line(36,prep?"absolute Matsubara threshold":"standardized residual threshold",
            format(prep?preparation_threshold(p,model):p.iteration_error_sigma_threshold))
        +print::quantity_to_output_line(36,p.sampling_strategy=="pcn"?"sampled-Z phase magnitude":"average phase magnitude",format(statistics.average_phase))
        +print::quantity_to_output_line(36,"complex-weight N_eff",format(statistics.effective_sample_size))
        +print::quantity_to_output_line(36,"largest Gaussian factor block",std::to_string(statistics.largest_factorization_dimension))
        +print::quantity_to_output_line(36,"closed-contour trace residual",format(closure_residual)));
    if(p.sampling_strategy=="pcn")
        print::print_R0(rank,print::quantity_to_output_line(36,"MH acceptance rate",format(statistics.acceptance_rate))
            +print::quantity_to_output_line(36,"zero-modulus rejection rate",format(statistics.zero_modulus_rejection_rate)));
    if(statistics.termination=="by convergence")
        print::print_R0(rank,"\033[1;32mTerminating: fixed point and contour diagnostics converged.\033[0m\n");
    else if(statistics.termination=="by iteration limit")
        print::print_R0(rank,"\033[1;31mTerminating because the iteration limit was reached.\033[0m\n");
    return !statistics.termination.empty();
}
}

namespace CspinDMFT::Run_Time_Data {
template IterationResidual iteration_residual(const MatsubaraCorrelations&,const ClusterMagnetization&,const MatsubaraCorrelations&,const ClusterMagnetization&,const PreparationData&,bool);
template IterationResidual iteration_residual(const ContourCorrelations&,const TimeMagnetization&,const ContourCorrelations&,const TimeMagnetization&,const RealTimeData&,bool);
template bool RunTimeData::record(const ParameterSpace&,const Model&,const MatsubaraCorrelations&,const ClusterMagnetization&,const MatsubaraCorrelations&,const ClusterMagnetization&,int);
template bool RunTimeData::record(const ParameterSpace&,const Model&,const ContourCorrelations&,const TimeMagnetization&,const ContourCorrelations&,const TimeMagnetization&,int);
}
