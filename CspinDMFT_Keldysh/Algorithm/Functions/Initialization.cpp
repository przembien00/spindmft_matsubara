#include "../main_header.h"
#include "../Storage_Concept/HDF5_Utilities.h"
namespace CspinDMFT::Functions {
namespace {
using namespace HDF5_Utilities;
struct InitialFile {
    Handle file,parameters;
    char symmetry=0;
    InitialFile(const Model& model,const ParameterSpace& p)
        :file(H5Fopen(p.initial_correlations_src_file.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose),
         parameters(H5Gopen2(file,"parameters",H5P_DEFAULT),H5Gclose) {
        if(read(parameters,"num_Spins",{},true)[0]!=double(model.sites)) throw std::invalid_argument("initial file cluster size differs");
        Handle attr(H5Aopen(parameters,"correlation_symmetry_type",H5P_DEFAULT),H5Aclose);
        Handle type(H5Aget_type(attr),H5Tclose);
        if(H5Tget_class(type)==H5T_STRING && H5Tis_variable_str(type)<=0) {
            std::vector<char> text(H5Tget_size(type)+1,0); check(H5Aread(attr,type,text.data())); symmetry=text[0];
        } else if(H5Tget_class(type)==H5T_INTEGER) check(H5Aread(attr,H5T_NATIVE_CHAR,&symmetry));
        else throw std::invalid_argument("unsupported initial symmetry metadata");
        direction_pairs(symmetry);
    }
};
}
MatsubaraCorrelations import_initial_correlations(const Model& model,const ParameterSpace& p) {
    const InitialFile source(model,p);
    const double old_beta=read(source.parameters,"beta",{},true)[0];
    const size_t ni=static_cast<size_t>(read(source.parameters,"num_ImagTimePoints",{},true)[0]);
    MatsubaraCorrelations correlations(model.sites,p.num_TimeSteps+1,p.correlation_symmetry_type,model.correlation_categories);
    if(ni<2) throw std::invalid_argument("invalid initial imaginary grid");
    if(!p.extrapolate_initial_spin_correlations && (ni!=correlations.imaginary_points || std::abs(old_beta-p.beta)>1e-12))
        throw std::invalid_argument("initial imaginary grid differs; use extrapolate");
    if(p.beta>old_beta+1e-12) throw std::invalid_argument("initial beta interval is too short for interpolation");
    MatsubaraCorrelations loaded(model.sites,ni,source.symmetry,model.correlation_categories);
    for(auto ij:loaded.categories) {
        const auto i=ij[0],j=ij[1];
        std::string pair=std::to_string(i+1)+"-"+std::to_string(j+1);
        auto re=read(source.file,("preparation/Re_correlation/"+pair).c_str(),{loaded.directions.size(),ni});
        auto im=read(source.file,("preparation/Im_correlation/"+pair).c_str(),{loaded.directions.size(),ni});
        for(size_t d=0;d<loaded.directions.size();++d) for(size_t k=0;k<ni;++k)
            loaded.set(3*i+loaded.directions[d][0],3*j+loaded.directions[d][1],k,C(re[d*ni+k],im[d*ni+k]));
    }
    for(auto ij:correlations.categories) for(auto ab:correlations.directions)
        for(size_t k=0;k<correlations.imaginary_points;++k) {
            const auto i=ij[0],j=ij[1];
            double coordinate=p.beta*double(k)/double(correlations.imaginary_points-1)/old_beta*double(ni-1);
            size_t left=std::min(static_cast<size_t>(std::floor(coordinate)),ni-2);
            double fraction=coordinate-double(left);
            correlations.set(3*i+ab[0],3*j+ab[1],k,(1.-fraction)*loaded.x(3*i+ab[0],3*j+ab[1],left)
                +fraction*loaded.x(3*i+ab[0],3*j+ab[1],left+1));
        }
    return correlations;
}
ClusterMagnetization import_initial_magnetization(const Model& model,const ParameterSpace& p) {
    const InitialFile source(model,p);
    ClusterMagnetization loaded(model.sites,source.symmetry),moments(model.sites,p.correlation_symmetry_type);
    const auto re=read(source.file,"preparation/Re_magnetization",{1,model.sites,3});
    const auto im=read(source.file,"preparation/Im_magnetization",{1,model.sites,3});
    for(size_t i=0;i<model.sites;++i) for(size_t a:loaded.directions)
        loaded.set(3*i+a,C(re[3*i+a],im[3*i+a]));
    for(size_t i=0;i<model.sites;++i) for(size_t a:moments.directions) moments.set(3*i+a,loaded.m(3*i+a));
    return moments;
}
}
