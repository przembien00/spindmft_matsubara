from pathlib import Path
root=Path(__file__).resolve().parents[2]
out=Path(__file__).resolve().parent
source=(root/'Algorithm/main.cpp').read_text()
helper=r'''
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
'''
source=source.replace('int main(int argc,char** argv)',helper+'\nint main(int argc,char** argv)')
source=source.replace('clock.measure("joint field distribution");','clock.measure("joint field distribution");\n        audit_probe(sampler,model,parameters,mean);')
(out/'probe_main.cpp').write_text(source)

functions=(root/'Algorithm/Functions/Functions.cpp').read_text().replace('#include "../main_header.h"','#include "main_header.h"\n#include <mpi.h>\nstd::array<double,4> audit_profile{};')
functions=functions.replace('    const size_t nc=l.components','    double audit_start=MPI_Wtime();\n    const size_t nc=l.components')
functions=functions.replace('    std::vector<Matrix> insertions(nc*ni);','    audit_profile[0]+=MPI_Wtime()-audit_start; audit_start=MPI_Wtime();\n    std::vector<Matrix> insertions(nc*ni);')
functions=functions.replace('    Matrix forward=model.identity,backward=model.identity;','    audit_profile[1]+=MPI_Wtime()-audit_start;\n    Matrix forward=model.identity,backward=model.identity;')
functions=functions.replace('    for(size_t t=0;t<nt;++t) {\n        if(t)', '    for(size_t t=0;t<nt;++t) {\n        audit_start=MPI_Wtime();\n        if(t)')
functions=functions.replace('        Matrix closed=forward*prefix.back()*backward;', '        audit_profile[2]+=MPI_Wtime()-audit_start; audit_start=MPI_Wtime();\n        Matrix closed=forward*prefix.back()*backward;')
functions=functions.replace('    }\n    return sample;','        audit_profile[3]+=MPI_Wtime()-audit_start;\n    }\n    return sample;')
(out/'probe_functions.cpp').write_text(functions)
