#pragma once
#include "../Functions/Functions.h"
namespace CspinDMFT::Run_Time_Data {
struct IterationResidual { double absolute=0.,standardized=0.; };
struct StageData {
    std::string termination;
    std::vector<double> absolute_iteration_errors,standardized_iteration_errors;
    std::vector<double> average_phase_magnitudes,gaussian_factor_reconstruction_errors,covariance_symmetry_errors;
    std::vector<double> partition_sum_Re,partition_sum_Im,effective_sample_sizes,mh_acceptance_rates,sample_counts;
    size_t num_Iterations=0;
};
// Runtime owns errors and diagnostics only. Observable values stay in main's
// two correlation and magnetization buffers.
template<class Correlations>
struct BasicStageData:StageData {
    Correlations correlation_error_re,correlation_error_im;
    MagnetizationFor<Correlations> magnetization_error_re,magnetization_error_im;
    std::vector<C> closure;
    C partition{};
    double average_phase=0.,factorization_error=0.,covariance_error=0.,acceptance_rate=0.,effective_sample_size=0.;
    double zero_modulus_rejection_rate=0.;
    size_t latent_dimension=0,largest_factorization_dimension=0,actual_num_blocks=0,samples_per_block=0,samples_per_core=0;
    std::vector<double> blocking_lengths,blocking_errors;
    // Per-iteration paired blocks for ratio errors and sampling diagnostics.
    std::vector<C> sample_blocks;
    size_t sample_count=0;
    double absolute_partition_sum=0.,squared_partition_sum=0.;
    int nonfinite_sample=0;
    bool correlated_samples=false;
    BasicStageData(size_t nt,size_t sites,size_t ni,char symmetry,const std::vector<IndexPair>& categories={})
        :correlation_error_re(make_correlations<Correlations>(nt,sites,ni,symmetry,categories)),
         correlation_error_im(make_correlations<Correlations>(nt,sites,ni,symmetry,categories)),
         magnetization_error_re(make_magnetization<MagnetizationFor<Correlations>>(nt,sites,symmetry)),
         magnetization_error_im(make_magnetization<MagnetizationFor<Correlations>>(nt,sites,symmetry)),closure(nt) {}
};
using PreparationData=BasicStageData<MatsubaraCorrelations>;
using RealTimeData=BasicStageData<ContourCorrelations>;
class RunTimeData {
public:
    PreparationData preparation;
    RealTimeData real_time;
    unsigned long long generated_seed=0;
    RunTimeData(const ParameterSpace& p,const Model& model)
        :preparation(1,model.sites,p.num_TimeSteps+1,p.correlation_symmetry_type,model.correlation_categories),
         real_time(p.num_RealTimeSteps+1,model.sites,p.num_TimeSteps+1,p.correlation_symmetry_type,model.correlation_categories),generated_seed(p.generated_seed) {}
    template<class Correlations> auto& stage() {
        if constexpr(std::is_same_v<Correlations,MatsubaraCorrelations>) return preparation;
        else return real_time;
    }
    template<class Correlations> const auto& stage() const {
        if constexpr(std::is_same_v<Correlations,MatsubaraCorrelations>) return preparation;
        else return real_time;
    }
    template<class Correlations>
    bool record(const ParameterSpace&,const Model&,const Correlations& previous,
        const MagnetizationFor<Correlations>& previous_magnetization,const Correlations& current,
        const MagnetizationFor<Correlations>& current_magnetization,int rank);
    double preparation_threshold(const ParameterSpace&,const Model&) const;
};
template<class Correlations>
IterationResidual iteration_residual(const Correlations&,const MagnetizationFor<Correlations>&,
    const Correlations&,const MagnetizationFor<Correlations>&,const BasicStageData<Correlations>&,bool constant_magnetization);
}
namespace CspinDMFT { using Run_Time_Data::BasicStageData; }
