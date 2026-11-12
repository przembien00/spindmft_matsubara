#pragma once
#include "../Mean_Field_Models/Mean_Field_Models.h"
namespace CspinDMFT::Run_Time_Data { template<class Correlations> struct BasicStageData; }
namespace CspinDMFT::Functions {
Matrix exponential(const Matrix&);
MatsubaraCorrelations generate_initial_correlations(const Model&,const ParameterSpace&);
ClusterMagnetization generate_initial_magnetization(const Model&,const ParameterSpace&);
ContourCorrelations isolated_correlations(const Model&,const ParameterSpace&,size_t times);
TimeMagnetization isolated_magnetization(const Model&,const ParameterSpace&,size_t times);
MatsubaraCorrelations import_initial_correlations(const Model&,const ParameterSpace&);
ClusterMagnetization import_initial_magnetization(const Model&,const ParameterSpace&);
template<class Correlations>
struct BasicSample {
    Correlations correlations;
    MagnetizationFor<Correlations> magnetization;
    std::vector<C> closure;
    C partition;
};
using PreparationSample=BasicSample<MatsubaraCorrelations>;
using RealTimeSample=BasicSample<ContourCorrelations>;
// partition_only returns just the partition; all observable containers are empty.
template<class Correlations>
BasicSample<Correlations> measure(const Model&,const ParameterSpace&,const Layout&,const Field&,const std::vector<C>& mean,bool partition_only=false);
template<class Correlations>
void initialize_sampling(Correlations&,MagnetizationFor<Correlations>&,
    Run_Time_Data::BasicStageData<Correlations>&,size_t samples,size_t requested_blocks,bool pcn);
template<class Correlations>
void accumulate_sample(Correlations&,MagnetizationFor<Correlations>&,
    const BasicSample<Correlations>&,double normalization=1.);
template<class Correlations>
void record_sample_statistics(Run_Time_Data::BasicStageData<Correlations>&,
    const BasicSample<Correlations>&,double normalization=1.);
template<class Correlations>
void reduce_and_normalize(Correlations&,MagnetizationFor<Correlations>&,
    Run_Time_Data::BasicStageData<Correlations>&,int ranks);
void seed_from_preparation(const MatsubaraCorrelations&,const ClusterMagnetization&,ContourCorrelations&,TimeMagnetization&);
void freeze_origin(ContourCorrelations&,TimeMagnetization&,const MatsubaraCorrelations&,const ClusterMagnetization&);
void constrain_feedback(ContourCorrelations&,TimeMagnetization&,const MatsubaraCorrelations&,const ClusterMagnetization&,const ParameterSpace&);
}
namespace CspinDMFT {
using Functions::exponential;
using Functions::measure;
using Functions::BasicSample;
}
