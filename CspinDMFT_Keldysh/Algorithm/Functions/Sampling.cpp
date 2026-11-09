#include "../Run_Time_Data/Run_Time_Data.h"
#include <mpi.h>
#include <climits>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace CspinDMFT::Functions {
using Run_Time_Data::BasicStageData;
namespace {
bool finite(C value) { return std::isfinite(value.real()) && std::isfinite(value.imag()); }
void sum(std::vector<C>& values) {
    for(size_t offset=0;offset<values.size();) {
        int count=static_cast<int>(std::min(values.size()-offset,size_t(INT_MAX)));
        MPI_Allreduce(MPI_IN_PLACE,values.data()+offset,count,MPI_C_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD); offset+=size_t(count);
    }
}
}
template<class Correlations>
void initialize_sampling(Correlations& correlations,MagnetizationFor<Correlations>& moments,
    BasicStageData<Correlations>& statistics,size_t samples,size_t requested_blocks,bool pcn) {
    size_t blocks=std::min(requested_blocks,samples/2);
    const size_t width=1+correlations.values.size()+moments.values.size()+statistics.closure.size();
    if(pcn) {
        size_t power=1;
        while(power<=blocks/2) power*=2;
        blocks=power;
        while(blocks>2 && samples%blocks) blocks/=2;
    } else while(blocks>2 && samples%blocks) --blocks;
    if(blocks && samples%blocks) blocks=1;
    if(blocks<2) throw std::invalid_argument("sample count must admit at least two equal blocks (use an even count)");
    if(blocks>size_t(INT_MAX)/width) throw std::runtime_error("MPI block buffer exceeds supported collective count");
    std::fill(correlations.values.begin(),correlations.values.end(),C{});
    std::fill(moments.values.begin(),moments.values.end(),C{});
    std::fill(statistics.closure.begin(),statistics.closure.end(),C{});
    statistics.sample_blocks.assign(blocks*width,C{});
    statistics.sample_count=0;
    statistics.absolute_partition_sum=statistics.squared_partition_sum=0.;
    statistics.nonfinite_sample=0;
    statistics.partition=C{};
    statistics.actual_num_blocks=blocks; statistics.samples_per_block=samples/blocks; statistics.samples_per_core=samples;
    statistics.correlated_samples=pcn;
}
template<class Correlations>
void accumulate_sample(Correlations& correlations,MagnetizationFor<Correlations>& moments,
    const BasicSample<Correlations>& sample,double normalization) {
    if(correlations.values.size()!=sample.correlations.values.size() || moments.values.size()!=sample.magnetization.values.size())
        throw std::invalid_argument("sample/observable shape mismatch");
    auto correlation=correlations.values.begin();
    for(C value:sample.correlations.values) *correlation+++=normalization*value;
    auto magnetization=moments.values.begin();
    for(C value:sample.magnetization.values) *magnetization+++=normalization*value;
}
template<class Correlations>
void record_sample_statistics(BasicStageData<Correlations>& statistics,const BasicSample<Correlations>& sample,double normalization) {
    const size_t nx=statistics.correlation_error_re.values.size(),nm=statistics.magnetization_error_re.values.size(),nt=statistics.closure.size();
    const size_t width=1+nx+nm+nt;
    if(statistics.sample_count>=statistics.samples_per_core || sample.correlations.values.size()!=nx
        || sample.magnetization.values.size()!=nm || sample.closure.size()!=nt)
        throw std::invalid_argument("sample/statistics shape mismatch");
    auto block=statistics.sample_blocks.begin()+(statistics.sample_count/statistics.samples_per_block)*width;
    const C z=sample.partition*normalization;
    *block+++=z; statistics.partition+=z;
    statistics.absolute_partition_sum+=std::abs(z); statistics.squared_partition_sum+=std::norm(z);
    if(!finite(z)) statistics.nonfinite_sample=1;
    for(const auto* array:{&sample.correlations.values,&sample.magnetization.values,&sample.closure})
        for(C value:*array) { value*=normalization; *block+++=value; if(!finite(value)) statistics.nonfinite_sample=1; }
    auto closure=statistics.closure.begin();
    for(C value:sample.closure) *closure+++=normalization*value;
    ++statistics.sample_count;
}
template<class Correlations>
void reduce_and_normalize(Correlations& next_correlations,MagnetizationFor<Correlations>& next_moments,
    BasicStageData<Correlations>& statistics,int ranks) {
    if(statistics.sample_count!=statistics.samples_per_core) throw std::logic_error("incomplete sampling before MPI reduction");
    const size_t samples=statistics.samples_per_core,blocks=statistics.actual_num_blocks;
    const size_t nx=next_correlations.values.size(),nm=next_moments.values.size(),nt=statistics.closure.size(),width=1+nx+nm+nt;
    const bool pcn=statistics.correlated_samples;
    auto& local=statistics.sample_blocks; auto& bad=statistics.nonfinite_sample;
    auto& abs_z=statistics.absolute_partition_sum; auto& abs_zsq=statistics.squared_partition_sum;
    MPI_Allreduce(MPI_IN_PLACE,&bad,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(bad) throw std::runtime_error("nonfinite sample: reduce time/beta or refine the grid");
    // Point estimates come from the observable containers; paired blocks are
    // retained only to estimate errors with the same complex denominator.
    std::vector<C> total;
    total.reserve(width); total.push_back(statistics.partition);
    total.insert(total.end(),next_correlations.values.begin(),next_correlations.values.end());
    total.insert(total.end(),next_moments.values.begin(),next_moments.values.end());
    total.insert(total.end(),statistics.closure.begin(),statistics.closure.end());
    sum(total);
    std::vector<C> pooled(local.size()*size_t(ranks));
    MPI_Allgather(local.data(),int(local.size()),MPI_C_DOUBLE_COMPLEX,pooled.data(),int(local.size()),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&abs_z,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&abs_zsq,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    const C z=total[0];
    if(!finite(z) || !std::isfinite(abs_z) || std::abs(z)<=1e-12*std::max(abs_z,1e-300)) throw std::runtime_error("unresolved complex partition denominator");
    statistics.blocking_lengths.clear();
    statistics.blocking_errors.clear();
    statistics.partition=z/double(samples*size_t(ranks)); statistics.average_phase=std::abs(z)/abs_z;
    statistics.actual_num_blocks=blocks; statistics.samples_per_block=samples/blocks; statistics.samples_per_core=samples;
    statistics.effective_sample_size=std::norm(z)/abs_zsq;
    std::vector<double> error_re(width),error_im(width);
    size_t per_rank=blocks,block_length=samples/blocks;
    // Independent draws use paired delete-block ratio jackknife. For pCN,
    // retain the largest error under successive coarsening within each chain.
    for(;;) {
        const size_t count=per_rank*size_t(ranks);
        double largest=0.;
        for(size_t j=1;j<width;++j) {
            if(!finite(total[j]/z)) throw std::runtime_error("nonfinite observable ratio");
            C average{}; std::vector<C> leave(count);
            for(size_t b=0;b<count;++b) {
                C denominator=z-pooled[b*width];
                if(std::abs(denominator)<=1e-12*std::max(abs_z,1e-300)) throw std::runtime_error("unresolved jackknife denominator");
                leave[b]=(total[j]-pooled[b*width+j])/denominator;
                if(!finite(leave[b])) throw std::runtime_error("nonfinite jackknife ratio");
                average+=leave[b]/double(count);
            }
            double vr=0.,vi=0.;
            for(C value:leave) { vr+=std::pow((value-average).real(),2); vi+=std::pow((value-average).imag(),2); }
            double re=std::sqrt(double(count-1)/double(count)*vr),im=std::sqrt(double(count-1)/double(count)*vi);
            error_re[j]=std::max(error_re[j],re); error_im[j]=std::max(error_im[j],im);
            largest=std::max({largest,re,im});
        }
        statistics.blocking_lengths.push_back(double(block_length)); statistics.blocking_errors.push_back(largest);
        if(!pcn || per_rank<4 || per_rank%2) break;
        const size_t next=per_rank/2;
        std::vector<C> coarse(size_t(ranks)*next*width,C{});
        for(size_t rank_i=0;rank_i<size_t(ranks);++rank_i) for(size_t b=0;b<next;++b) for(size_t j=0;j<width;++j)
            coarse[(rank_i*next+b)*width+j]=pooled[(rank_i*per_rank+2*b)*width+j]+pooled[(rank_i*per_rank+2*b+1)*width+j];
        pooled=std::move(coarse); per_rank=next; block_length*=2;
    }
    for(size_t j=0;j<nx;++j) { next_correlations.values[j]=total[1+j]/z; statistics.correlation_error_re.values[j]=error_re[1+j]; statistics.correlation_error_im.values[j]=error_im[1+j]; }
    for(size_t j=0;j<nm;++j) { next_moments.values[j]=total[1+nx+j]/z; statistics.magnetization_error_re.values[j]=error_re[1+nx+j]; statistics.magnetization_error_im.values[j]=error_im[1+nx+j]; }
    for(size_t j=0;j<nt;++j) statistics.closure[j]=total[1+nx+nm+j]/z;
}
template void initialize_sampling(MatsubaraCorrelations&,ClusterMagnetization&,BasicStageData<MatsubaraCorrelations>&,size_t,size_t,bool);
template void initialize_sampling(ContourCorrelations&,TimeMagnetization&,BasicStageData<ContourCorrelations>&,size_t,size_t,bool);
template void accumulate_sample(MatsubaraCorrelations&,ClusterMagnetization&,const BasicSample<MatsubaraCorrelations>&,double);
template void accumulate_sample(ContourCorrelations&,TimeMagnetization&,const BasicSample<ContourCorrelations>&,double);
template void record_sample_statistics(BasicStageData<MatsubaraCorrelations>&,const BasicSample<MatsubaraCorrelations>&,double);
template void record_sample_statistics(BasicStageData<ContourCorrelations>&,const BasicSample<ContourCorrelations>&,double);
template void reduce_and_normalize(MatsubaraCorrelations&,ClusterMagnetization&,BasicStageData<MatsubaraCorrelations>&,int);
template void reduce_and_normalize(ContourCorrelations&,TimeMagnetization&,BasicStageData<ContourCorrelations>&,int);
}
