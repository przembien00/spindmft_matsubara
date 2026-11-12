#include "main_header.h"
#include <mpi.h>
#include <iostream>
#include <stdexcept>

using namespace CspinDMFT;
namespace {
void near(C actual,C expected,const char* message) {
    if(std::abs(actual-expected)>2e-12*std::max(1.,std::abs(expected)))
        throw std::runtime_error(message);
}
template<class Correlations>
void check(int rank,int ranks,bool pcn,char symmetry) {
    constexpr size_t samples=16,blocks=8;
    const size_t times=std::is_same_v<Correlations,MatsubaraCorrelations>?1:3;
    auto x=make_correlations<Correlations>(times,2,5,symmetry,{{0,1},{0,0}});
    auto m=make_magnetization<MagnetizationFor<Correlations>>(times,2,symmetry);
    Run_Time_Data::BasicStageData<Correlations> statistics(times,2,5,symmetry,{{0,1},{0,0}});
    func::initialize_sampling(x,m,statistics,samples,blocks,pcn);
    BasicSample<Correlations> sample{x,m,std::vector<C>(times),C{}};
    for(size_t s=0;s<samples;++s) {
        sample.partition=C(1.+.03*s+.2*rank,.07*(rank+1));
        for(size_t j=0;j<x.values.size();++j)
            sample.correlations.values[j]=C(.1*s+.02*j-.03*rank,.04*s-.1*rank+.001*j);
        for(size_t j=0;j<m.values.size();++j)
            sample.magnetization.values[j]=C(.03*s-.02*j+.1*rank,.01*s+.03*rank);
        for(size_t t=0;t<times;++t) sample.closure[t]=sample.partition+C(.01*t*s,-.02*t*rank);
        const double weight=pcn?1./(1.+.1*s+.05*rank):1.;
        func::accumulate_sample(x,m,sample,weight);
        func::record_sample_statistics(statistics,sample,weight);
    }
    // Small test-only oracle reproduces the original replicated estimator.
    // Production must never gather these block arrays.
    const size_t width=1+x.values.size()+m.values.size()+times;
    std::vector<C> pooled(statistics.sample_blocks.size()*size_t(ranks));
    MPI_Allgather(statistics.sample_blocks.data(),int(statistics.sample_blocks.size()),MPI_C_DOUBLE_COMPLEX,
        pooled.data(),int(statistics.sample_blocks.size()),MPI_C_DOUBLE_COMPLEX,MPI_COMM_WORLD);
    std::vector<C> total(width);
    for(size_t b=0;b<blocks*size_t(ranks);++b)for(size_t j=0;j<width;++j) total[j]+=pooled[b*width+j];
    std::vector<double> re(width),im(width),lengths,errors;
    size_t per_rank=blocks,length=samples/blocks;
    for(;;) {
        const size_t count=per_rank*size_t(ranks);
        double largest=0.;
        for(size_t j=1;j<width;++j) {
            std::vector<C> leave(count);
            C average{};
            for(size_t b=0;b<count;++b) {
                leave[b]=(total[j]-pooled[b*width+j])/(total[0]-pooled[b*width]);
                average+=leave[b]/double(count);
            }
            double vr=0.,vi=0.;
            for(C value:leave) { vr+=std::pow((value-average).real(),2); vi+=std::pow((value-average).imag(),2); }
            const double er=std::sqrt(double(count-1)/double(count)*vr),ei=std::sqrt(double(count-1)/double(count)*vi);
            re[j]=std::max(re[j],er); im[j]=std::max(im[j],ei);
            largest=std::max({largest,er,ei});
        }
        lengths.push_back(double(length)); errors.push_back(largest);
        if(!pcn || per_rank<4) break;
        const size_t next=per_rank/2;
        std::vector<C> coarse(size_t(ranks)*next*width);
        for(size_t r=0;r<size_t(ranks);++r)for(size_t b=0;b<next;++b)for(size_t j=0;j<width;++j)
            coarse[(r*next+b)*width+j]=pooled[(r*per_rank+2*b)*width+j]+pooled[(r*per_rank+2*b+1)*width+j];
        pooled=std::move(coarse); per_rank=next; length*=2;
    }
    func::reduce_and_normalize(x,m,statistics,ranks);
    for(size_t j=0;j<x.values.size();++j) {
        near(x.values[j],total[1+j]/total[0],"correlation ratio");
        near(statistics.correlation_error_re.values[j],re[1+j],"real correlation jackknife");
        near(statistics.correlation_error_im.values[j],im[1+j],"imaginary correlation jackknife");
    }
    for(size_t j=0;j<m.values.size();++j) {
        const size_t k=1+x.values.size()+j;
        near(m.values[j],total[k]/total[0],"magnetization ratio");
        near(statistics.magnetization_error_re.values[j],re[k],"real magnetization jackknife");
        near(statistics.magnetization_error_im.values[j],im[k],"imaginary magnetization jackknife");
    }
    for(size_t t=0;t<times;++t) near(statistics.closure[t],total[width-times+t]/total[0],"closure ratio");
    if(statistics.blocking_lengths!=lengths) throw std::runtime_error("blocking lengths");
    for(size_t l=0;l<errors.size();++l) near(statistics.blocking_errors[l],errors[l],"blocking envelope");
    near(statistics.partition,total[0]/double(samples*size_t(ranks)),"mean partition");
}
void check_invalid_denominator(int rank,int ranks) {
    MatsubaraCorrelations x(1,5,'A'); ClusterMagnetization m(1,'A');
    Run_Time_Data::PreparationData statistics(1,1,5,'A');
    func::initialize_sampling(x,m,statistics,8,4,false);
    BasicSample<MatsubaraCorrelations> sample{x,m,std::vector<C>(1),C{}};
    for(size_t i=0;i<8;++i) {
        sample.partition=rank==0 && i==0?1.:0.;
        func::accumulate_sample(x,m,sample);
        func::record_sample_statistics(statistics,sample);
    }
    bool caught=false;
    try { func::reduce_and_normalize(x,m,statistics,ranks); }
    catch(const std::runtime_error& e) { caught=std::string(e.what())=="unresolved jackknife denominator"; }
    if(!caught) throw std::runtime_error("invalid rank-local denominator must fail on every rank");
}
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank,ranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
        for(char symmetry:std::string("ABCD"))for(bool pcn:{false,true}) {
            check<MatsubaraCorrelations>(rank,ranks,pcn,symmetry);
            check<ContourCorrelations>(rank,ranks,pcn,symmetry);
        }
        check_invalid_denominator(rank,ranks);
        if(rank==0) std::cout<<"distributed jackknife agrees with pooled reference on "<<ranks<<" ranks\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
