#include "../main_header.h"

namespace CspinDMFT::Contour {
namespace {
C interpolate(const ContourCorrelations& x,double t,size_t a,size_t b,size_t k) {
    if(t<0 || t>double(x.times-1)) throw std::out_of_range("kernel lag");
    size_t interval=static_cast<size_t>(std::floor(t));
    if(t==double(interval)) return x.x(interval,a,b,k);
    size_t count=std::min(size_t{4},x.times),start=std::min(interval?interval-1:0,x.times-count);
    C value{};
    for(size_t j=0;j<count;++j) {
        double w=1;
        for(size_t l=0;l<count;++l) if(l!=j) w*=(t-double(start+l))/(double(j)-double(l));
        value+=w*x.x(start+j,a,b,k);
    }
    return value;
}
}
C kernel(const ContourCorrelations& x,Node f,Node s,size_t q) {
    using S=Sector;
    if(q==0) throw std::invalid_argument("zero substeps");
    if(f.sector==S::Nu) return s.sector==S::Eta?kernel(x,s,f,q):C{};
    if(s.sector==S::Nu) {
        if(f.sector!=S::Eta || s.point>f.point) return {};
        double lag=(double(f.point-s.point)+.5)/double(q);
        return interpolate(x,lag,f.component,s.component,x.imaginary_points-1)-interpolate(x,lag,f.component,s.component,0);
    }
    if(f.sector==S::M && s.sector==S::M) {
        if(f.point>s.point) return x.x(0,s.component,f.component,f.point-s.point);
        if(f.point<s.point) return x.x(0,f.component,s.component,s.point-f.point);
        if(f.component>s.component) return kernel(x,s,f,q);
        return .5*(x.x(0,f.component,s.component,0)+x.x(0,f.component,s.component,x.imaginary_points-1));
    }
    if(f.sector==S::M) return kernel(x,s,f,q);
    if(s.sector==S::M) return interpolate(x,(double(f.point)+.5)/double(q),f.component,s.component,s.point);
    if(f.point<s.point || (f.point==s.point && f.component>s.component)) return kernel(x,s,f,q);
    double lag=double(f.point-s.point)/double(q);
    return .5*(interpolate(x,lag,f.component,s.component,0)+interpolate(x,lag,f.component,s.component,x.imaginary_points-1));
}
template<class Correlations>
Matrix covariance(const Correlations& bath,const Layout& l,size_t q) {
    Matrix cov(l.dimension(),l.dimension(),C{});
    for(size_t i=0;i<l.dimension();++i) for(size_t j=i;j<l.dimension();++j)
        cov(i,j)=cov(j,i)=kernel(bath,l.decode(i),l.decode(j),q);
    return cov;
}
}

namespace CspinDMFT::Contour {
C kernel(const MatsubaraCorrelations& x,Node f,Node s,size_t q) {
    if(!q || f.sector!=Sector::M || s.sector!=Sector::M) throw std::invalid_argument("Matsubara covariance only has imaginary-time nodes");
    if(f.point>s.point) return x.x(s.component,f.component,f.point-s.point);
    if(f.point<s.point) return x.x(f.component,s.component,s.point-f.point);
    if(f.component>s.component) return kernel(x,s,f,q);
    return .5*(x.x(f.component,s.component,0)+x.x(f.component,s.component,x.imaginary_points-1));
}
template Matrix covariance(const MatsubaraCorrelations&,const Layout&,size_t);
template Matrix covariance(const ContourCorrelations&,const Layout&,size_t);
}
