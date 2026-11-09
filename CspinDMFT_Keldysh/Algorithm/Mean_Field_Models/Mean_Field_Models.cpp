#include "../main_header.h"
#include <Physics/Spin.h>

namespace CspinDMFT::Mean_Field_Models {
double Model::weight(size_t i,size_t j,size_t k,size_t l) const {
    return correlation_weights.at(((i*sites+j)*sites+k)*sites+l);
}
void Model::initialize(const ParameterSpace& p) {
    coupling=p.spinmf_cmodel.coupling_matrix;
    size_t dim=1;
    for(double s:spins) {
        if(!(s>0) || !std::isfinite(s) || std::abs(2*s-std::round(2*s))>1e-12)
            throw std::invalid_argument("spin lengths must be positive integer or half integer");
        const size_t local=static_cast<size_t>(std::llround(2*s+1));
        if(local>4096 || dim>4096/local) throw std::invalid_argument("dense cluster Hilbert space exceeds 4096 states");
        dim*=local;
    }
    identity=blaze::IdentityMatrix<C>(dim);
    operators.clear();
    for(size_t i=0;i<sites;++i) {
        Physics::Spin::Observable sx,sy,sz;
        Physics::Spin::write_spin_matrices(spins[i],sx,sy,sz);
        for(const auto* s:{&sx,&sy,&sz})
            operators.emplace_back(Physics::Spin::create_linear_spin_term(spins,i,*s));
    }
    hamiltonian=Matrix(dim,dim,C{});
    for(size_t i=0;i<sites;++i) {
        if(p.chemical_shift!="none") {
            auto fields=Print_Routines::split_string_at_delimiter(p.chemical_shift,',');
            if(fields.size()!=sites) throw std::invalid_argument("chemshift requires one value per site");
            size_t end=0; double shift=std::stod(fields[i],&end);
            if(end!=fields[i].size() || !std::isfinite(shift)) throw std::invalid_argument("invalid chemical shift");
            hamiltonian+=shift*operators[3*i+2];
        }
        if(p.extra_interaction=="quadrupolar") hamiltonian+=3.*p.interaction_strength*(operators[3*i+2]*operators[3*i+2]);
        for(size_t j=i+1;j<sites;++j) for(size_t a=0;a<3;++a) for(size_t b=0;b<3;++b)
            hamiltonian+=p.rescale_spinspin*J[i*sites+j]*p.spinspin_cmodel.coupling_matrix(a,b)*(operators[3*i+a]*operators[3*j+b]);
    }
    for(auto& w:expectation_weights) w*=p.rescale_meanfield;
    for(auto& w:correlation_weights) w*=p.rescale_meanfield*p.rescale_meanfield;
}
}

namespace CspinDMFT::Mean_Field_Models {
template<class Correlations>
Correlations connected(const Correlations& correlations,const MagnetizationFor<Correlations>& moments) {
    Correlations out=make_correlations_like(correlations,'D');
    for(size_t t=0;t<time_points(out);++t) for(auto ij:out.categories) for(auto ab:out.directions)
        for(size_t k=0;k<out.imaginary_points;++k) {
            const auto a=3*ij[0]+ab[0],b=3*ij[1]+ab[1];
            set_correlation(out,t,a,b,k,correlation(correlations,t,a,b,k)-magnetization(moments,t,a)*magnetization(moments,0,b));
        }
    return out;
}
template<class Correlations>
Correlations bath_primitive(const Model& m,const Correlations& correlations,const MagnetizationFor<Correlations>& moments) {
    Correlations source=connected(correlations,moments);
    auto out=make_correlations<Correlations>(time_points(correlations),m.sites,correlations.imaginary_points,'D');
    std::vector<bool> included(m.sites*m.sites,false);
    for(auto kl:source.categories) included[kl[0]*m.sites+kl[1]]=true;
    for(size_t i=0;i<m.sites;++i) for(size_t j=0;j<m.sites;++j)
        for(auto kl:source.categories) {
            const auto k=kl[0],l=kl[1];
            // A representative stands for both orientations unless each is measured explicitly.
            double w=m.weight(i,j,k,l);
            if(k!=l && !included[l*m.sites+k]) w*=2.;
            if(w==0.) continue;
            for(size_t a=0;a<3;++a) for(size_t b=0;b<3;++b) for(size_t c=0;c<3;++c) for(size_t d=0;d<3;++d) {
                const double factor=w*m.coupling(a,c)*m.coupling(b,d); if(factor==0.) continue;
                for(size_t t=0;t<time_points(correlations);++t) for(size_t tau=0;tau<correlations.imaginary_points;++tau)
                    set_correlation(out,t,3*i+a,3*j+b,tau,correlation(out,t,3*i+a,3*j+b,tau)+factor*correlation(source,t,3*k+c,3*l+d,tau));
            }
        }
    return out;
}
template<class Magnetization>
std::vector<C> mean_field(const Model& m,const Magnetization& moments) {
    std::vector<C> mean(time_points(moments)*3*m.sites);
    for(size_t t=0;t<time_points(moments);++t) for(size_t i=0;i<m.sites;++i)
        for(size_t k=0;k<m.sites;++k) for(size_t a=0;a<3;++a) for(size_t b=0;b<3;++b)
            mean[t*3*m.sites+3*i+a]+=m.expectation_weights[i*m.sites+k]*m.coupling(a,b)*magnetization(moments,t,3*k+b);
    return mean;
}
template<class Correlations>
Correlations prescribed_bath(const Model& model,const ParameterSpace& p,const Correlations& grid) {
    auto out=make_correlations<Correlations>(time_points(grid),model.sites,grid.imaginary_points,'D');
    const size_t a=std::string("xyz").find(p.bath_component);
    const double occupation=1./std::expm1(p.beta*p.bath_frequency);
    for(size_t t=0;t<time_points(out);++t) for(size_t tau=0;tau<out.imaginary_points;++tau) {
        const double imaginary=p.beta*double(tau)/double(p.num_TimeSteps);
        const C phase(0.,p.bath_frequency*p.delta_real_t*double(t));
        C value=p.bath_coupling*p.bath_coupling*(occupation+1.)*(std::exp(-p.bath_frequency*imaginary+phase)
            +std::exp(-p.bath_frequency*(p.beta-imaginary)-phase));
        for(size_t i=0;i<model.sites;++i) set_correlation(out,t,3*i+a,3*i+a,tau,value);
    }
    return out;
}
template<class Correlations>
Matrix self_consistent_equations(const Model& model,const ParameterSpace& p,
    const Correlations& correlations,const MagnetizationFor<Correlations>& moments,const Contour::Layout& layout) {
    auto field_correlations=p.uses_harmonic_bath()?prescribed_bath(model,p,correlations)
        :bath_primitive(model,correlations,moments);
    if(p.noise.m_name!="none") for(size_t i=0;i<model.sites;++i)
        for(size_t t=0;t<time_points(field_correlations);++t)
            for(size_t k=0;k<field_correlations.imaginary_points;++k)
                set_correlation(field_correlations,t,3*i+2,3*i+2,k,
                    correlation(field_correlations,t,3*i+2,3*i+2,k)+p.noise.m_variance_in(2,2));
    return Contour::covariance(field_correlations,layout,p.real_time_substeps);
}
template Matrix self_consistent_equations(const Model&,const ParameterSpace&,const MatsubaraCorrelations&,const ClusterMagnetization&,const Contour::Layout&);
template Matrix self_consistent_equations(const Model&,const ParameterSpace&,const ContourCorrelations&,const TimeMagnetization&,const Contour::Layout&);
template MatsubaraCorrelations connected(const MatsubaraCorrelations&,const ClusterMagnetization&);
template ContourCorrelations connected(const ContourCorrelations&,const TimeMagnetization&);
template MatsubaraCorrelations bath_primitive(const Model&,const MatsubaraCorrelations&,const ClusterMagnetization&);
template ContourCorrelations bath_primitive(const Model&,const ContourCorrelations&,const TimeMagnetization&);
template std::vector<C> mean_field(const Model&,const ClusterMagnetization&);
template std::vector<C> mean_field(const Model&,const TimeMagnetization&);
template MatsubaraCorrelations prescribed_bath(const Model&,const ParameterSpace&,const MatsubaraCorrelations&);
template ContourCorrelations prescribed_bath(const Model&,const ParameterSpace&,const ContourCorrelations&);
}

#include "../Storage_Concept/HDF5_Utilities.h"
namespace CspinDMFT::Mean_Field_Models {
using namespace HDF5_Utilities;
Model read_model(const ParameterSpace& p) {
    Handle f(H5Fopen(p.configuration_path().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
    Model m;
    m.sites=p.num_Spins;
    m.spins=H5Aexists(f,"spin lengths")>0?read(f,"spin lengths",{m.sites},true):std::vector<double>(m.sites,.5);
    m.J=read(f,"spin-spin couplings",{m.sites,m.sites});
    m.expectation_weights=read(f,"expectation weights",{m.sites,m.sites});
    m.correlation_weights=read(f,"correlation weights",{m.sites,m.sites,m.sites,m.sites});
    for(size_t i=0;i<m.sites;++i) for(size_t j=0;j<m.sites;++j) {
        if(std::abs(m.J[i*m.sites+j]-m.J[j*m.sites+i])>1e-12) throw std::invalid_argument("nonsymmetric spin-spin couplings");
        if(i==j && std::abs(m.J[i*m.sites+i])>1e-12) throw std::invalid_argument("spin-spin coupling diagonal must be zero");
        for(size_t k=0;k<m.sites;++k) for(size_t l=0;l<m.sites;++l)
            if(std::abs(m.weight(i,j,k,l)-m.weight(j,i,k,l))>1e-12)
                throw std::invalid_argument("correlation weights must be symmetric in field sites");
    }
    // Read the representative pairs used for measurement and feedback.
    if(H5Lexists(f,"correlation_categories",H5P_DEFAULT)>0) {
        Handle dataset(H5Dopen2(f,"correlation_categories",H5P_DEFAULT),H5Dclose);
        Handle space(H5Dget_space(dataset),H5Sclose);
        hsize_t dims[2]{};
        if(H5Sget_simple_extent_ndims(space)!=2) throw std::invalid_argument("invalid correlation_categories rank");
        H5Sget_simple_extent_dims(space,dims,nullptr);
        if(dims[1]!=2 || !dims[0] || dims[0]>m.sites*m.sites) throw std::invalid_argument("invalid category shape");
        std::vector<unsigned long long> data(2*dims[0]);
        check(H5Dread(dataset,H5T_NATIVE_ULLONG,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data()));
        for(size_t i=0;i<dims[0];++i) {
            if(data[2*i]>=m.sites || data[2*i+1]>=m.sites) throw std::invalid_argument("invalid category site");
            m.correlation_categories.push_back({size_t(data[2*i]),size_t(data[2*i+1])});
        }
    } else for(size_t i=0;i<m.sites;++i) for(size_t j=0;j<m.sites;++j) m.correlation_categories.push_back({i,j});
    if(m.correlation_categories.empty()) throw std::invalid_argument("configuration must have at least one correlation category");
    std::vector<bool> included(m.sites*m.sites,false);
    for(auto ij:m.correlation_categories) {
        const auto index=ij[0]*m.sites+ij[1];
        if(included[index]) throw std::invalid_argument("duplicate correlation category");
        included[index]=true;
    }
    // The configuration folds equivalent source pairs into these representatives.
    for(size_t i=0;i<m.sites;++i) for(size_t j=0;j<m.sites;++j)
        for(size_t k=0;k<m.sites;++k) for(size_t l=0;l<m.sites;++l)
            if(!included[k*m.sites+l] && m.weight(i,j,k,l)!=0.
                && !(included[l*m.sites+k] && std::abs(m.weight(i,j,k,l)-m.weight(i,j,l,k))<=1e-12))
                throw std::invalid_argument("correlation weights reference an unlisted correlation category");
    m.initialize(p);
    return m;
}
}
