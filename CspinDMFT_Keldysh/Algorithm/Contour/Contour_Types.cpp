#include "Contour_Types.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <Observables/Correlations.h>
#include <Observables/Tensors.h>
#include <Observables/Magnetization.h>
namespace CspinDMFT::Contour {
namespace {
size_t product(std::initializer_list<size_t> factors) {
    size_t value=1;
    for(size_t factor:factors) {
        if(factor && value>std::numeric_limits<size_t>::max()/factor) throw std::invalid_argument("observable allocation overflows size_t");
        value*=factor;
    }
    return value;
}
size_t direction(char symmetry,const std::vector<IndexPair>& pairs,size_t a,size_t b) {
    a%=3; b%=3;
    if(a==b && ((symmetry=='A') || (a==1 && symmetry!='D'))) a=b=0;
    auto found=std::find(pairs.begin(),pairs.end(),IndexPair{a,b});
    return found==pairs.end()?pairs.size():size_t(found-pairs.begin());
}
std::vector<IndexPair> correlation_categories(size_t sites,const std::vector<IndexPair>& requested) {
    if(!requested.empty()) return requested;
    std::vector<IndexPair> categories;
    for(size_t i=0;i<sites;++i) for(size_t j=0;j<sites;++j) categories.push_back({i,j});
    return categories;
}
std::vector<size_t> category_lookup(size_t sites,const std::vector<IndexPair>& categories) {
    std::vector<size_t> indices(product({sites,sites}),categories.size());
    for(size_t n=0;n<categories.size();++n) {
        const auto ij=categories[n];
        if(ij[0]>=sites || ij[1]>=sites) throw std::invalid_argument("correlation category outside cluster");
        auto& index=indices[ij[0]*sites+ij[1]];
        if(index!=categories.size()) throw std::invalid_argument("duplicate correlation category");
        index=n;
    }
    return indices;
}
size_t corr_index(size_t t,size_t sites,size_t ni,size_t nd,size_t a,size_t b,size_t k,size_t d,
    const std::vector<size_t>& indices,size_t categories) {
    if(a>=3*sites || b>=3*sites || k>=ni) throw std::out_of_range("cluster correlation index");
    const size_t category=indices[a/3*sites+b/3];
    if(category==categories) throw std::out_of_range("site pair is not a configured correlation category");
    return ((t*categories+category)*nd+d)*ni+k;
}
size_t mag_direction(const std::vector<size_t>& dirs,size_t a) {
    auto found=std::find(dirs.begin(),dirs.end(),a%3);
    return found==dirs.end()?dirs.size():size_t(found-dirs.begin());
}
void require_origin(size_t t) { if(t!=0) throw std::out_of_range("preparation has no real-time axis"); }
}
std::vector<IndexPair> direction_pairs(char symmetry) {
    if(std::string("ABCD").find(symmetry)==std::string::npos) throw std::invalid_argument("invalid correlation symmetry");
    using Tensor=Observables::Tensors::CorrelationTensor<Observables::Correlations::CorrelationVector>;
    return Tensor(symmetry,size_t{0}).get_direction_pairs();
}
std::vector<size_t> magnetization_directions(char symmetry) {
    if(std::string("ABCD").find(symmetry)==std::string::npos) throw std::invalid_argument("invalid magnetization symmetry");
    return Observables::Magnetization::determine_magnetization_directions(symmetry);
}
MatsubaraCorrelations::MatsubaraCorrelations(size_t n,size_t ni,char s,const std::vector<IndexPair>& pairs)
    :sites(n),imaginary_points(ni),symmetry(s),directions(direction_pairs(s)),
     categories(correlation_categories(n,pairs)),values(product({categories.size(),directions.size(),ni})),
     category_indices(category_lookup(n,categories)) {}
ContourCorrelations::ContourCorrelations(size_t nt,size_t n,size_t ni,char s,const std::vector<IndexPair>& pairs)
    :times(nt),sites(n),imaginary_points(ni),symmetry(s),directions(direction_pairs(s)),
     categories(correlation_categories(n,pairs)),values(product({nt,categories.size(),directions.size(),ni})),
     category_indices(category_lookup(n,categories)) {}
C MatsubaraCorrelations::x(size_t a,size_t b,size_t k) const {
    const size_t d=direction(symmetry,directions,a,b);
    size_t index=corr_index(0,sites,imaginary_points,directions.size(),a,b,k,d,category_indices,categories.size());
    return d==directions.size()?C{}:values.at(index);
}
void MatsubaraCorrelations::set(size_t a,size_t b,size_t k,C v) {
    const size_t d=direction(symmetry,directions,a,b);
    size_t index=corr_index(0,sites,imaginary_points,directions.size(),a,b,k,d,category_indices,categories.size());
    if(d==directions.size()) { if(v!=C{}) throw std::invalid_argument("correlation direction forbidden by symmetry"); return; }
    values.at(index)=v;
}
C ContourCorrelations::x(size_t t,size_t a,size_t b,size_t k) const {
    if(t>=times) throw std::out_of_range("correlation time index");
    const size_t d=direction(symmetry,directions,a,b);
    size_t index=corr_index(t,sites,imaginary_points,directions.size(),a,b,k,d,category_indices,categories.size());
    return d==directions.size()?C{}:values.at(index);
}
void ContourCorrelations::set(size_t t,size_t a,size_t b,size_t k,C v) {
    if(t>=times) throw std::out_of_range("correlation time index");
    const size_t d=direction(symmetry,directions,a,b);
    size_t index=corr_index(t,sites,imaginary_points,directions.size(),a,b,k,d,category_indices,categories.size());
    if(d==directions.size()) { if(v!=C{}) throw std::invalid_argument("correlation direction forbidden by symmetry"); return; }
    values.at(index)=v;
}
ClusterMagnetization::ClusterMagnetization(size_t n,char s):sites(n),symmetry(s),directions(magnetization_directions(s)),values(product({n,directions.size()})) {}
TimeMagnetization::TimeMagnetization(size_t nt,size_t n,char s):times(nt),sites(n),symmetry(s),directions(magnetization_directions(s)),values(product({nt,n,directions.size()})) {}
C ClusterMagnetization::m(size_t a) const {
    if(a>=3*sites) throw std::out_of_range("magnetization index");
    size_t d=mag_direction(directions,a); return d==directions.size()?C{}:values.at(a/3*directions.size()+d);
}
void ClusterMagnetization::set(size_t a,C v) {
    if(a>=3*sites) throw std::out_of_range("magnetization index");
    size_t d=mag_direction(directions,a);
    if(d==directions.size()) { if(v!=C{}) throw std::invalid_argument("magnetization direction forbidden by symmetry"); return; }
    values.at(a/3*directions.size()+d)=v;
}
C TimeMagnetization::m(size_t t,size_t a) const {
    if(t>=times || a>=3*sites) throw std::out_of_range("magnetization index");
    size_t d=mag_direction(directions,a); return d==directions.size()?C{}:values.at((t*sites+a/3)*directions.size()+d);
}
void TimeMagnetization::set(size_t t,size_t a,C v) {
    if(t>=times || a>=3*sites) throw std::out_of_range("magnetization index");
    size_t d=mag_direction(directions,a);
    if(d==directions.size()) { if(v!=C{}) throw std::invalid_argument("magnetization direction forbidden by symmetry"); return; }
    values.at((t*sites+a/3)*directions.size()+d)=v;
}
C correlation(const MatsubaraCorrelations& c,size_t t,size_t a,size_t b,size_t k) { require_origin(t); return c.x(a,b,k); }
C correlation(const ContourCorrelations& c,size_t t,size_t a,size_t b,size_t k) { return c.x(t,a,b,k); }
void set_correlation(MatsubaraCorrelations& c,size_t t,size_t a,size_t b,size_t k,C v) { require_origin(t); c.set(a,b,k,v); }
void set_correlation(ContourCorrelations& c,size_t t,size_t a,size_t b,size_t k,C v) { c.set(t,a,b,k,v); }
C magnetization(const ClusterMagnetization& m,size_t t,size_t a) { require_origin(t); return m.m(a); }
C magnetization(const TimeMagnetization& m,size_t t,size_t a) { return m.m(t,a); }
void set_magnetization(ClusterMagnetization& m,size_t t,size_t a,C v) { require_origin(t); m.set(a,v); }
void set_magnetization(TimeMagnetization& m,size_t t,size_t a,C v) { m.set(t,a,v); }
MatsubaraCorrelations make_correlations_like(const MatsubaraCorrelations& c,char s) { return {c.sites,c.imaginary_points,s,c.categories}; }
ContourCorrelations make_correlations_like(const ContourCorrelations& c,char s) { return {c.times,c.sites,c.imaginary_points,s,c.categories}; }
ClusterMagnetization make_magnetization_like(const ClusterMagnetization& m) { return {m.sites,m.symmetry}; }
TimeMagnetization make_magnetization_like(const TimeMagnetization& m) { return {m.times,m.sites,m.symmetry}; }
size_t Layout::dimension() const {
    const size_t limit=std::numeric_limits<size_t>::max();
    if(ni==0 || components==0 || ni==limit || (intervals && intervals>(limit-ni-1)/2)) throw std::invalid_argument("invalid or overflowing contour grid");
    return product({components,ni+(intervals?2*intervals+1:0)});
}
size_t Layout::flat(Sector s,size_t p,size_t a) const {
    const size_t count=s==Sector::M?ni:s==Sector::Eta?intervals:intervals+1;
    if(a>=components || p>=count || (s!=Sector::M && !intervals)) throw std::out_of_range("field node");
    return ((s==Sector::M?0:s==Sector::Eta?ni:ni+intervals)+p)*components+a;
}
Node Layout::decode(size_t i) const {
    if(i>=dimension()) throw std::out_of_range("field index");
    size_t p=i/components,a=i%components;
    if(p<ni) return {Sector::M,p,a};
    if(p<ni+intervals) return {Sector::Eta,p-ni,a};
    return {Sector::Nu,p-ni-intervals,a};
}
}
