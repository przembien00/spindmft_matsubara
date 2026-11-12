#include "main_header.h"
#include <mpi.h>
std::array<double,4> audit_profile{};

namespace CspinDMFT::Functions {
Matrix exponential( const Matrix& A )
{
    if( A.rows()!=A.columns() )
        throw std::invalid_argument("matrix exponential requires a square matrix");
    const size_t n=A.rows();
    if( n==0 ) return Matrix{};

    RealType norm_one{};
    for( size_t j=0;j<n;++j )
    {
        RealType column_sum{};
        for( size_t i=0;i<n;++i ) column_sum+=std::abs(A(i,j));
        norm_one=std::max(norm_one,column_sum);
    }
    constexpr RealType theta13=RealType{5.371920351148152};
    const int squarings=norm_one>theta13
        ?std::max(0,static_cast<int>(std::ceil(std::log2(norm_one/theta13)))):0;
    const RealType scale=std::ldexp(RealType{1.},squarings);
    const Matrix B=A/scale;
    const Matrix B2=B*B,B4=B2*B2,B6=B4*B2;
    const Matrix I=blaze::IdentityMatrix<ComplexType,blaze::rowMajor>(n);
    const std::array<RealType,14> c{
        RealType{64764752532480000.},RealType{32382376266240000.},
        RealType{7771770303897600.}, RealType{1187353796428800.},
        RealType{129060195264000.},  RealType{10559470521600.},
        RealType{670442572800.},     RealType{33522128640.},
        RealType{1323241920.},      RealType{40840800.},
        RealType{960960.},          RealType{16380.},
        RealType{182.},             RealType{1.}
    };
    const Matrix U=B*(B6*(c[13]*B6+c[11]*B4+c[9]*B2)
                       +c[7]*B6+c[5]*B4+c[3]*B2+c[1]*I);
    const Matrix V=B6*(c[12]*B6+c[10]*B4+c[8]*B2)
                    +c[6]*B6+c[4]*B4+c[2]*B2+c[0]*I;
    Matrix result=blaze::inv(V-U)*(V+U);
    for( int s=0;s<squarings;++s ) result=result*result;
    return result;
}

namespace {
C contract(const Matrix& a,const Matrix& b) {
    C sum{};
    for(size_t i=0;i<a.rows();++i) for(size_t j=0;j<a.columns();++j) sum+=a(i,j)*b(j,i);
    return sum;
}
}
template<class Correlations>
Correlations isolated_correlations_impl(const Model& model,const ParameterSpace& p,size_t nt) {
    Correlations result=make_correlations<Correlations>(nt,model.sites,p.num_TimeSteps+1,p.correlation_symmetry_type,model.correlation_categories);
    Matrix rho=exponential(-p.beta*model.hamiltonian);
    C z=blaze::trace(rho);
    for(size_t t=0;t<nt;++t) {
        Matrix u=exponential(C(0,-p.delta_real_t*double(t))*model.hamiltonian);
        Matrix v=exponential(C(0,p.delta_real_t*double(t))*model.hamiltonian);
        for(auto ij:result.categories) for(auto ab:result.directions) {
            const auto i=ij[0],j=ij[1];
            Matrix spin=v*model.operators[3*i+ab[0]]*u;
            for(size_t k=0;k<result.imaginary_points;++k) {
                double tau=p.beta*double(k)/double(p.num_TimeSteps);
                Matrix left=exponential(-(p.beta-tau)*model.hamiltonian),right=exponential(-tau*model.hamiltonian);
                Matrix insertion=left*model.operators[3*j+ab[1]]*right;
                set_correlation(result,t,3*i+ab[0],3*j+ab[1],k,contract(insertion,spin)/z);
            }
        }
    }
    return result;
}
MatsubaraCorrelations generate_initial_correlations(const Model& model,const ParameterSpace& p) {
    if(p.init_corr_mode=="import") return import_initial_correlations(model,p);
    auto result=isolated_correlations_impl<MatsubaraCorrelations>(model,p,1);
    std::function<double(double)> diagonal=[](double) { return 0.; };
    std::function<double(double)> offdiagonal=[](double) { return 0.; };
    if(p.init_diag_corr!="zero" && p.init_diag_corr!="isolated")
        diagonal=Physics::DiagonalSpinCorrelation(p.init_diag_corr,p.beta,p.beta/p.correlation_periods).f;
    if(p.init_nondiag_corr=="sin") offdiagonal=Physics::NonDiagonalSpinCorrelation("sin",p.beta,p.beta/p.correlation_periods).f;
    for(auto ij:result.categories) for(auto ab:result.directions)
        for(size_t k=0;k<result.imaginary_points;++k) {
            const auto i=ij[0],j=ij[1];
            const double tau=p.beta*double(k)/double(p.num_TimeSteps),spin=model.spins[i];
            if(ab[0]==ab[1] && p.init_diag_corr!="isolated") result.set(3*i+ab[0],3*j+ab[1],k,i==j?spin*(spin+1.)/3.*diagonal(tau):0.);
            if(ab[0]!=ab[1] && p.init_nondiag_corr!="isolated") result.set(3*i+ab[0],3*j+ab[1],k,i==j?spin*(spin+1.)/3.*offdiagonal(tau):0.);
        }
    return result;
}
ClusterMagnetization isolated_magnetization_impl(const Model& model,const ParameterSpace& p) {
    ClusterMagnetization result(model.sites,p.correlation_symmetry_type);
    Matrix rho=exponential(-p.beta*model.hamiltonian);
    const C z=blaze::trace(rho);
    for(size_t i=0;i<model.sites;++i) for(size_t a:result.directions) result.set(3*i+a,contract(rho,model.operators[3*i+a])/z);
    return result;
}
ClusterMagnetization generate_initial_magnetization(const Model& model,const ParameterSpace& p) {
    if(p.init_corr_mode=="import") return import_initial_magnetization(model,p);
    return isolated_magnetization_impl(model,p);
}
ContourCorrelations isolated_correlations(const Model& model,const ParameterSpace& p,size_t nt) {
    return isolated_correlations_impl<ContourCorrelations>(model,p,nt);
}
TimeMagnetization isolated_magnetization(const Model& model,const ParameterSpace& p,size_t nt) {
    const auto thermal=isolated_magnetization_impl(model,p);
    TimeMagnetization result(nt,model.sites,p.correlation_symmetry_type);
    for(size_t t=0;t<nt;++t) for(size_t i=0;i<model.sites;++i) for(size_t a:result.directions) result.set(t,3*i+a,thermal.m(3*i+a));
    return result;
}
template<class Correlations>
BasicSample<Correlations> measure(const Model& model,const ParameterSpace& p,const Layout& l,const Field& field,const std::vector<C>& mean,bool partition_only) {
    double audit_start=MPI_Wtime();
    const size_t nc=l.components,ni=l.ni,nt=l.intervals?p.num_RealTimeSteps+1:1;
    if(field.size()!=l.dimension() || mean.size()!=nt*nc) throw std::invalid_argument("measurement grid mismatch");
    auto node_h=[&](Sector sector,size_t point,bool rest) {
        Matrix h=rest?model.hamiltonian:Matrix(model.identity.rows(),model.identity.columns(),C{});
        for(size_t a=0;a<nc;++a) h+=field[l.flat(sector,point,a)]*model.operators[a];
        return h;
    };
    auto imag_h=[&](double coordinate) {
        // Local four-knot interpolation of sampled Matsubara edge fields, as in staged CF4.
        size_t interval=static_cast<size_t>(std::floor(coordinate));
        size_t start=std::min(interval?interval-1:0,ni-4);
        Matrix h=model.hamiltonian;
        for(size_t a=0;a<nc;++a) {
            C v=mean[a];
            for(size_t j=0;j<4;++j) {
                double w=1;
                for(size_t k=0;k<4;++k) if(k!=j) w*=(coordinate-double(start+k))/(double(j)-double(k));
                v+=w*field[l.flat(Sector::M,start+j,a)];
            }
            h+=v*model.operators[a];
        }
        return h;
    };
    std::vector<Matrix> steps(ni-1),prefix(ni),suffix(ni);
    prefix[0]=model.identity; suffix.back()=model.identity;
    double root=std::sqrt(3.),early=(3.-2.*root)/12.,late=(3.+2.*root)/12.,h=p.beta/double(ni-1);
    for(size_t k=0;k+1<ni;++k) {
        Matrix a=imag_h(double(k)+.5-root/6.),b=imag_h(double(k)+.5+root/6.);
        steps[k]=exponential(-h*(early*a+late*b))*exponential(-h*(late*a+early*b));
        prefix[k+1]=steps[k]*prefix[k];
    }
    for(size_t k=ni-1;k-->0;) suffix[k]=suffix[k+1]*steps[k];
    auto correlations=make_correlations<Correlations>(nt,model.sites,ni,p.correlation_symmetry_type,model.correlation_categories);
    auto magnetization=make_magnetization<MagnetizationFor<Correlations>>(nt,model.sites,p.correlation_symmetry_type);
    if(partition_only) return BasicSample<Correlations>{std::move(correlations),std::move(magnetization),{},blaze::trace(prefix.back())};
    audit_profile[0]+=MPI_Wtime()-audit_start; audit_start=MPI_Wtime();
    std::vector<Matrix> insertions(nc*ni);
    for(size_t b=0;b<nc;++b) for(size_t k=0;k<ni;++k) insertions[b*ni+k]=suffix[k]*model.operators[b]*prefix[k];
    BasicSample<Correlations> sample{std::move(correlations),std::move(magnetization),std::vector<C>(nt),blaze::trace(prefix.back())};
    audit_profile[1]+=MPI_Wtime()-audit_start;
    Matrix forward=model.identity,backward=model.identity;
    const double dt=p.delta_real_t/double(p.real_time_substeps);
    for(size_t t=0;t<nt;++t) {
        audit_start=MPI_Wtime();
        if(t) for(size_t j=0;j<p.real_time_substeps;++j) {
            size_t micro=(t-1)*p.real_time_substeps+j;
            Matrix eta=node_h(Sector::Eta,micro,true);
            const double fraction=(double(j)+.5)/double(p.real_time_substeps);
            for(size_t a=0;a<nc;++a)
                eta+=((1.-fraction)*mean[(t-1)*nc+a]+fraction*mean[t*nc+a])*model.operators[a];
            Matrix left=exponential(C(0,-dt/4.)*node_h(Sector::Nu,micro,false));
            Matrix right=exponential(C(0,-dt/4.)*node_h(Sector::Nu,micro+1,false));
            Matrix u=right*exponential(C(0,-dt)*eta)*left;
            Matrix b=left*exponential(C(0,dt)*eta)*right;
            forward=u*forward; backward=backward*b;
        }
        audit_profile[2]+=MPI_Wtime()-audit_start; audit_start=MPI_Wtime();
        Matrix closed=forward*prefix.back()*backward;
        sample.closure[t]=blaze::trace(closed);
        for(size_t i=0;i<model.sites;++i) for(size_t a:sample.magnetization.directions) {
            Matrix insertion=backward*model.operators[3*i+a]*forward;
            set_magnetization(sample.magnetization,t,3*i+a,contract(prefix.back(),insertion));
        }
        for(auto ij:sample.correlations.categories) for(auto ab:sample.correlations.directions) {
            const auto i=ij[0],j=ij[1];
            Matrix insertion=backward*model.operators[3*i+ab[0]]*forward;
            for(size_t k=0;k<ni;++k)
                set_correlation(sample.correlations,t,3*i+ab[0],3*j+ab[1],k,contract(insertions[(3*j+ab[1])*ni+k],insertion));
        }
        audit_profile[3]+=MPI_Wtime()-audit_start;
    }
    return sample;
}
}

namespace CspinDMFT::Functions {
template PreparationSample measure<MatsubaraCorrelations>(const Model&,const ParameterSpace&,const Layout&,const Field&,const std::vector<C>&,bool);
template RealTimeSample measure<ContourCorrelations>(const Model&,const ParameterSpace&,const Layout&,const Field&,const std::vector<C>&,bool);
}
