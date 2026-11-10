#include "main_header.h"
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <chrono>
#include "Storage_Concept/HDF5_Utilities.h"
#include <Keldysh_Contour/Contour_Kernel.h>
using namespace CspinDMFT;
namespace {
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
void near(C a,C b,double tolerance,const char* message) { require(std::abs(a-b)<tolerance,message); }
template<class T,class=void> struct has_times:std::false_type {};
template<class T> struct has_times<T,std::void_t<decltype(std::declval<T>().times)>>:std::true_type {};
static_assert(!has_times<MatsubaraCorrelations>::value);
static_assert(!has_times<ClusterMagnetization>::value);
static_assert(has_times<ContourCorrelations>::value);
static_assert(has_times<TimeMagnetization>::value);
double norm(const Matrix& a) {
    double result=0;
    for(size_t i=0;i<a.rows();++i) for(size_t j=0;j<a.columns();++j) result+=std::norm(a(i,j));
    return std::sqrt(result);
}
}
namespace {
template<class Correlations>
void check_sample_accumulation(size_t times,char symmetry,bool pcn,const std::vector<IndexPair>& categories={}) {
    auto correlations=make_correlations<Correlations>(times,2,5,symmetry,categories);
    auto moments=make_magnetization<MagnetizationFor<Correlations>>(times,2,symmetry);
    Run_Time_Data::BasicStageData<Correlations> statistics(times,2,5,symmetry,categories);
    BasicSample<Correlations> sample{correlations,moments,std::vector<C>(times),C{}};
    std::fill(correlations.values.begin(),correlations.values.end(),C(99.));
    std::fill(moments.values.begin(),moments.values.end(),C(99.));
    func::initialize_sampling(correlations,moments,statistics,8,4,pcn);
    C expected_c{},expected_m{},expected_z{},expected_closure{};
    for(size_t i=0;i<8;++i) {
        const C c(.2+double(i),.1*double(i)),m(.4*double(i),-.3),z(1.+.1*double(i),.2),closure(.7,.02*double(i));
        const double normalization=pcn?1./(1.+double(i)):1.;
        std::fill(sample.correlations.values.begin(),sample.correlations.values.end(),c);
        std::fill(sample.magnetization.values.begin(),sample.magnetization.values.end(),m);
        std::fill(sample.closure.begin(),sample.closure.end(),closure); sample.partition=z;
        func::accumulate_sample(correlations,moments,sample,normalization);
        func::record_sample_statistics(statistics,sample,normalization);
        expected_c+=normalization*c; expected_m+=normalization*m;
        expected_z+=normalization*z; expected_closure+=normalization*closure;
    }
    for(C value:correlations.values)near(value,expected_c,1e-13,"direct correlation accumulation");
    for(C value:moments.values)near(value,expected_m,1e-13,"direct magnetization accumulation");
    near(statistics.partition,expected_z,1e-13,"paired partition accumulation");
    for(C value:statistics.closure)near(value,expected_closure,1e-13,"closure statistics accumulation");
    const size_t width=1+correlations.values.size()+moments.values.size()+times;
    C block_c{},block_z{};
    for(size_t i=0;i<statistics.actual_num_blocks;++i) {
        block_z+=statistics.sample_blocks[i*width]; block_c+=statistics.sample_blocks[i*width+1];
    }
    near(block_c,expected_c,1e-13,"paired correlation blocks"); near(block_z,expected_z,1e-13,"paired partition blocks");
    require(statistics.sample_count==8 && statistics.nonfinite_sample==0,"sample diagnostics");
    func::initialize_sampling(correlations,moments,statistics,8,4,pcn);
    for(C value:correlations.values)near(value,0.,1e-13,"iteration clears correlations");
    for(C value:moments.values)near(value,0.,1e-13,"iteration clears magnetization");
    require(statistics.sample_count==0 && statistics.partition==C{},"iteration clears diagnostics");
}
void check_initial_import(const Model& model,ParameterSpace p) {
    using namespace CspinDMFT::HDF5_Utilities;
    const auto path=std::filesystem::temp_directory_path()/("cluster-initial-import-"
        +std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".hdf5");
    const auto attribute=[](hid_t group,const char* name,double value) {
        Handle space(H5Screate(H5S_SCALAR),H5Sclose);
        Handle attr(H5Acreate2(group,name,H5T_NATIVE_DOUBLE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
        check(H5Awrite(attr,H5T_NATIVE_DOUBLE,&value));
    };
    const auto dataset=[](hid_t group,const std::string& name,const std::vector<hsize_t>& shape,const std::vector<double>& values) {
        Handle space(H5Screate_simple(int(shape.size()),shape.data(),nullptr),H5Sclose);
        Handle data(H5Dcreate2(group,name.c_str(),H5T_NATIVE_DOUBLE,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
        check(H5Dwrite(data,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()));
    };
    const auto corr=[](size_t i,size_t j,size_t d,size_t k) { return C(.31+.01*i+.02*j+.001*d+.0001*k,-.07-.003*d); };
    const auto moment=[](size_t a) { return C(.11+.01*a,-.021-.002*a); };
    {
        Handle file(H5Fcreate(path.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
        Handle parameters(H5Gcreate2(file,"parameters",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
        attribute(parameters,"num_Spins",double(model.sites)); attribute(parameters,"beta",p.beta);
        attribute(parameters,"num_ImagTimePoints",double(p.num_TimeSteps+1));
        Handle space(H5Screate(H5S_SCALAR),H5Sclose);
        Handle symmetry(H5Acreate2(parameters,"correlation_symmetry_type",H5T_NATIVE_CHAR,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
        const char type='D'; check(H5Awrite(symmetry,H5T_NATIVE_CHAR,&type));
        Handle prep(H5Gcreate2(file,"preparation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
        Handle re(H5Gcreate2(prep,"Re_correlation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
        Handle im(H5Gcreate2(prep,"Im_correlation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose);
        for(size_t i=0;i<model.sites;++i) for(size_t j=0;j<model.sites;++j) {
            std::vector<double> r,v;
            for(size_t d=0;d<9;++d) for(size_t k=0;k<=p.num_TimeSteps;++k) { C value=corr(i,j,d,k); r.push_back(value.real()); v.push_back(value.imag()); }
            const auto pair=std::to_string(i+1)+"-"+std::to_string(j+1);
            dataset(re,pair,{9,p.num_TimeSteps+1},r); dataset(im,pair,{9,p.num_TimeSteps+1},v);
        }
        std::vector<double> r,v;
        for(size_t a=0;a<3*model.sites;++a) { r.push_back(moment(a).real()); v.push_back(moment(a).imag()); }
        dataset(prep,"Re_magnetization",{1,model.sites,3},r); dataset(prep,"Im_magnetization",{1,model.sites,3},v);
    }
    p.init_corr_mode="import"; p.initial_correlations_src_file=path.string();
    for(char symmetry:std::string("DB")) {
        p.correlation_symmetry_type=symmetry;
        const auto correlations=func::generate_initial_correlations(model,p);
        const auto magnetization=func::generate_initial_magnetization(model,p);
        for(size_t i=0;i<model.sites;++i) {
            for(size_t a:magnetization.directions) near(magnetization.m(3*i+a),moment(3*i+a),1e-14,"generator imports complex magnetization");
            for(auto ij:correlations.categories) if(ij[0]==i) for(auto ab:correlations.directions) for(size_t k=0;k<=p.num_TimeSteps;++k)
                near(correlations.x(3*i+ab[0],3*ij[1]+ab[1],k),corr(i,ij[1],3*ab[0]+ab[1],k),1e-14,"generator imports correlations");
        }
    }
    std::filesystem::remove(path);
}
}
int main() {
    try {
        ParameterSpace p; p.init_nondiag_corr="isolated"; p.init_diag_corr="isolated"; p.num_TimeSteps=4; p.num_RealTimeSteps=3; p.real_time_substeps=2; p.chemical_shift=".4,.4"; p.beta=1.7;
        Model model; model.sites=2; model.spins={.5,1.}; model.J={0,.7,.7,0};
        model.expectation_weights.resize(4); model.correlation_weights.resize(16); model.initialize(p);
        require(model.identity.rows()==6,"mixed spin Hilbert dimension");
        Matrix comm=model.operators[0]*model.operators[1]-model.operators[1]*model.operators[0];
        require(norm(Matrix(comm-C(0,1)*model.operators[2]))<1e-12,"spin algebra");
        Matrix a=model.hamiltonian+C(0,.3)*model.operators[0];
        Matrix inverse_check=exponential(a)*exponential(-a);
        require(norm(Matrix(inverse_check-model.identity))<1e-11,"non-Hermitian exponential inverse");
        Layout layout{5,6,6};
        for(size_t i=0;i<layout.dimension();++i) {
            auto n=layout.decode(i); require(layout.flat(n.sector,n.point,n.component)==i,"layout roundtrip");
        }
        auto exact=func::isolated_correlations(model,p,4);
        auto exactmag=func::isolated_magnetization(model,p,4);
        Field zero(layout.dimension(),C{});
        auto measured=measure<ContourCorrelations>(model,p,layout,zero,mean_field(model,exactmag));
        for(size_t i=0;i<exact.values.size();++i)
            near(measured.correlations.values[i]/measured.partition,exact.values[i],1e-11,"isolated cluster correlation");
        for(size_t i=0;i<exactmag.values.size();++i)
            near(measured.magnetization.values[i]/measured.partition,exactmag.values[i],1e-11,"isolated cluster magnetization");
        for(C z:measured.closure) near(z,measured.partition,1e-11,"isolated closure");
        // A static complex field on every site/direction must reproduce the exact
        // non-Hermitian cluster Hamiltonian, including its mixed insertions.
        Field static_fields(layout.dimension(),C{});
        Model shifted=model;
        for(size_t a=0;a<6;++a) {
            C v(.03*double(a+1),-.02*double(a+2));
            shifted.hamiltonian+=v*model.operators[a];
            for(size_t k=0;k<5;++k) static_fields[layout.flat(Sector::M,k,a)]=v;
            for(size_t k=0;k<6;++k) static_fields[layout.flat(Sector::Eta,k,a)]=v;
        }
        auto complex_exact=func::isolated_correlations(shifted,p,4);
        auto complex_sample=measure<ContourCorrelations>(model,p,layout,static_fields,std::vector<C>(24));
        for(size_t i=0;i<complex_exact.values.size();++i)
            near(complex_sample.correlations.values[i]/complex_sample.partition,
                 complex_exact.values[i],2e-11,"static complex cluster fields");
        Model sparse_model=model; sparse_model.correlation_categories={{0,1},{0,0}};
        const auto sparse_sample=measure<ContourCorrelations>(sparse_model,p,layout,static_fields,std::vector<C>(24));
        require(sparse_sample.correlations.categories==sparse_model.correlation_categories,"configured category order");
        require(sparse_sample.correlations.values.size()==4*2*9*5,"only configured real-time pairs allocated");
        for(auto ij:sparse_sample.correlations.categories) for(auto ab:sparse_sample.correlations.directions)
            for(size_t t=0;t<4;++t) for(size_t k=0;k<5;++k)
                near(sparse_sample.correlations.x(t,3*ij[0]+ab[0],3*ij[1]+ab[1],k),
                    complex_sample.correlations.x(t,3*ij[0]+ab[0],3*ij[1]+ab[1],k),1e-13,"representative measurement equals full reference");
        const auto sparse_prepared=func::generate_initial_correlations(sparse_model,p);
        require(sparse_prepared.values.size()==2*9*5,"only configured preparation pairs allocated");
        auto sparse_real=sparse_sample.correlations;
        auto sparse_moments=sparse_sample.magnetization;
        func::seed_from_preparation(sparse_prepared,func::generate_initial_magnetization(sparse_model,p),sparse_real,sparse_moments);
        for(auto ij:sparse_prepared.categories) for(auto ab:sparse_prepared.directions)
            near(sparse_real.x(3,3*ij[0]+ab[0],3*ij[1]+ab[1],2),sparse_prepared.x(3*ij[0]+ab[0],3*ij[1]+ab[1],2),1e-13,"representative origin seeding");
        bool missing_pair_rejected=false;
        try { sparse_prepared.x(3,3,0); } catch(const std::out_of_range&) { missing_pair_rejected=true; }
        require(missing_pair_rejected,"unlisted pairs have no invented trajectory symmetry");
        // Uniform z response commutes with this ISO Hamiltonian. The full
        // branch maps have a closed form independent of the Strang construction.
        Field response(layout.dimension(),C{});
        C nu(.2,.1);
        for(size_t k=0;k<7;++k) for(size_t a:{size_t{2},size_t{5}})
            response[layout.flat(Sector::Nu,k,a)]=nu;
        auto response_sample=measure<ContourCorrelations>(model,p,layout,response,std::vector<C>(24));
        Matrix rho=exponential(-p.beta*model.hamiltonian);
        Matrix nufield=nu*(model.operators[2]+model.operators[5]);
        for(size_t t=0;t<4;++t) {
            Matrix f=exponential(C(0,-p.delta_real_t*double(t))*(model.hamiltonian+.5*nufield));
            Matrix b=exponential(C(0,p.delta_real_t*double(t))*(model.hamiltonian-.5*nufield));
            Matrix closed=f*rho*b;
            near(response_sample.closure[t],blaze::trace(closed),2e-11,"response endpoint factors");
            for(size_t a=0;a<6;++a) {
                Matrix product=rho*b*model.operators[a]*f;
                near(response_sample.magnetization.m(t,a),blaze::trace(product),2e-11,"response prefix insertion");
            }
        }
        ContourCorrelations x(4,2,5,'D');
        for(size_t t=0;t<4;++t) for(size_t i=0;i<6;++i) for(size_t j=0;j<6;++j) for(size_t k=0;k<5;++k)
            x.set(t,i,j,k,C(.1*t+.03*i+.02*j+.01*k,.04*t-.01*i+.03*j-.02*k));
        near(kernel(x,{Sector::Eta,1,0},{Sector::Nu,2,5},2),0.,1e-14,"retarded support");
        near(kernel(x,{Sector::Eta,1,0},{Sector::Nu,1,5},2),x.x(0,0,5,4)-x.x(0,0,5,0),1e-14,"response sign");
        near(kernel(x,{Sector::M,4,1},{Sector::M,1,5},2),x.x(0,5,1,3),1e-14,"Matsubara site exchange");
        near(kernel(x,{Sector::Eta,1,1},{Sector::Eta,4,5},1),.5*(x.x(3,5,1,0)+x.x(3,5,1,4)),1e-14,"negative lag site exchange");
        Matrix cov=covariance(x,layout,2);
        require(norm(Matrix(cov-blaze::trans(cov)))<1e-13,"unconjugated covariance symmetry");
        // Fixed sector multipliers preserve the physical law and pCN factor rows.
        Sampler sector_sampler(cov,layout.ni*layout.components,layout.intervals*layout.components,1e-8);
        std::vector<double> reference_weights(layout.dimension(),1.);
        for(size_t i=0;i<reference_weights.size();++i) {
            const auto sector=layout.decode(i).sector;
            if(sector==Sector::Eta)reference_weights[i]=2.;
            else if(sector==Sector::Nu)reference_weights[i]=.5;
        }
        Sampler reference_sampler(cov,reference_weights,1e-8);
        const auto physical_factor=sector_sampler.leading_factor_rows(layout.dimension());
        const auto reference_factor=reference_sampler.leading_factor_rows(layout.dimension());
        require(norm(Matrix(physical_factor*blaze::trans(physical_factor)-cov))<1e-8,"sector scaling reconstructs physical covariance");
        require(norm(Matrix(physical_factor*blaze::ctrans(physical_factor)
            -reference_factor*blaze::ctrans(reference_factor)))<1e-8,"sector scaling preserves sampling law");
        std::mt19937 sampler_engine(42);
        const auto latent=sector_sampler.draw_latent(sampler_engine);
        const Field expected_field=physical_factor*latent;
        const auto field=sector_sampler.field_from_latent(latent);
        for(size_t i=0;i<field.size();++i)near(field[i],expected_field[i],1e-12,"sector field inverse scaling");
        for(const auto& sample:sector_sampler.draw_contour_batch(sampler_engine,3))
            require(sample.edge_field.size()==layout.dimension(),"sector batch field size");
        // Match the shared single-site implementation on every one-site staggered block.
        namespace ct=spinDMFT::Contour;
        ct::CorrelationSet single('D',5,4);
        ContourCorrelations scalar(4,1,5,'D');
        for(size_t t=0;t<4;++t) for(size_t i=0;i<3;++i) for(size_t j=0;j<3;++j) for(size_t k=0;k<5;++k) {
            C v=x.x(t,i,j,k); scalar.set(t,i,j,k,v);
            for(size_t d=0;d<9;++d) {
                auto pair=single.Re[t].get_direction_pair(d);
                if(pair[0]==i && pair[1]==j) { single.Re[t][d][k]=v.real(); single.Im[t][d][k]=v.imag(); }
            }
        }
        Layout sl{5,6,3};
        for(size_t i=0;i<sl.dimension();++i) for(size_t j=0;j<sl.dimension();++j) {
            auto f=sl.decode(i),s=sl.decode(j);
            auto convert=[](Node n) { return ct::StaggeredIndex{static_cast<ct::StaggeredSector>(n.sector),n.point,n.component}; };
            near(kernel(scalar,f,s,2),ct::staggered_correlation(single,convert(f),convert(s),2),1e-13,"staged kernel equivalence");
        }
        Model weights; weights.sites=2; weights.coupling=Physics::FieldMatrix{{1.,0.,0.},{0.,-.5,0.},{0.,0.,1.}};
        weights.expectation_weights={1,2,3,4}; weights.correlation_weights.resize(16);
        weights.correlation_weights[1]=2.; // W[0,0,0,1]
        TimeMagnetization xmag(4,2,'D'); xmag.set(0,0,.2); xmag.set(0,4,.3);
        auto bath=bath_primitive(weights,x,xmag);
        near(bath.x(0,0,1,2),-1.*(x.x(0,0,4,2)-xmag.m(0,0)*xmag.m(0,4)),1e-14,"connected four-index feedback");
        ContourCorrelations sparse_feedback(4,2,5,'D',{{0,1}});
        for(auto ab:sparse_feedback.directions) for(size_t t=0;t<4;++t) for(size_t k=0;k<5;++k)
            sparse_feedback.set(t,ab[0],3+ab[1],k,x.x(t,ab[0],3+ab[1],k));
        Model representative_weights=weights;
        representative_weights.correlation_weights[1]=1.; // Per orientation, previously 2.
        const auto representative_bath=bath_primitive(representative_weights,sparse_feedback,xmag);
        require(representative_bath.values.size()==bath.values.size(),"field kernel retains all cluster site pairs");
        for(size_t i=0;i<bath.values.size();++i)near(representative_bath.values[i],bath.values[i],1e-13,"representative feedback includes both source orientations");
        const auto representative_covariance=covariance(representative_bath,layout,2);
        const auto full_covariance=covariance(bath,layout,2);
        require(norm(Matrix(representative_covariance-full_covariance))<1e-13,"representative feedback matches full contour covariance");
        Model ordered_weights=representative_weights;
        ordered_weights.correlation_weights[2]=1.; // Explicit reverse orientation.
        const auto ordered_bath=bath_primitive(ordered_weights,x,xmag);
        near(ordered_bath.x(0,0,1,2),-.5*(x.x(0,0,4,2)-xmag.m(0,0)*xmag.m(0,4)
            +x.x(0,3,1,2)-xmag.m(0,3)*xmag.m(0,1)),1e-14,"ordered source pairs each counted once");
        ordered_weights.correlation_weights[0]=3.;
        const auto diagonal_bath=bath_primitive(ordered_weights,x,xmag);
        near(diagonal_bath.x(0,0,1,2)-ordered_bath.x(0,0,1,2),
            -1.5*(x.x(0,0,1,2)-xmag.m(0,0)*xmag.m(0,1)),1e-14,"diagonal source pair counted once");
        const auto prepared=func::generate_initial_correlations(model,p);
        const auto prepared_moments=func::generate_initial_magnetization(model,p);
        func::freeze_origin(exact,exactmag,prepared,prepared_moments);
        near(exact.x(0,2,4,3),prepared.x(2,4,3),1e-14,"frozen origin");
        const Layout ml{5,0,6};
        const auto measured_preparation=measure<MatsubaraCorrelations>(model,p,ml,Field(ml.dimension(),C{}),mean_field(model,prepared_moments));
        for(size_t j=0;j<prepared.values.size();++j)
            near(measured_preparation.correlations.values[j]/measured_preparation.partition,prepared.values[j],1e-11,"imaginary-only preparation measurement");
        for(size_t i=0;i<ml.dimension();++i) for(size_t j=0;j<ml.dimension();++j)
            near(kernel(prepared,ml.decode(i),ml.decode(j),1),kernel(exact,ml.decode(i),ml.decode(j),1),1e-13,"Matsubara and real-time origin kernels agree");
        // Allocation and access enforce the native A/B/C/D schemes, rather than
        // storing discarded full tensors or a one-point real-time axis.
        for(char symmetry:std::string("ABCD")) {
            MatsubaraCorrelations c(2,5,symmetry);
            ClusterMagnetization m(2,symmetry);
            require(c.values.size()==4*direction_pairs(symmetry).size()*5,"symmetry-reduced Matsubara allocation");
            require(m.values.size()==2*magnetization_directions(symmetry).size(),"symmetry-reduced first moments");
            c.set(0,0,2,.25);
            near(c.x(1,1,2),symmetry=='D'?0.:.25,1e-14,"diagonal alias under spin symmetry");
            if(symmetry!='D') near(c.x(0,2,2),0.,1e-14,"forbidden correlation has exact zero");
            if(symmetry=='B' || symmetry=='C') { m.set(2,.1); near(m.m(2),.1,1e-14,"native axial magnetization"); }
            ContourCorrelations real(4,2,5,symmetry); TimeMagnetization real_m(4,2,symmetry);
            require(real.values.size()==4*c.values.size(),"real-time correlation allocation");
            func::seed_from_preparation(c,m,real,real_m);
            for(size_t t=0;t<4;++t) near(real.x(t,0,0,2),.25,1e-14,"seeding distinct container types");
        }
        for(char symmetry:std::string("ABCD")) for(bool pcn:{false,true}) {
            check_sample_accumulation<MatsubaraCorrelations>(1,symmetry,pcn);
            check_sample_accumulation<ContourCorrelations>(4,symmetry,pcn);
            check_sample_accumulation<MatsubaraCorrelations>(1,symmetry,pcn,{{0,1},{0,0}});
            check_sample_accumulation<ContourCorrelations>(4,symmetry,pcn,{{0,1},{0,0}});
        }
        check_initial_import(model,p);
        check_initial_import(sparse_model,p);
        std::cout<<"cluster core checks passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
