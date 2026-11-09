#pragma once
#include <vector>
#include <array>
#include <type_traits>
#include <stdexcept>
#include "../matrices.h"
namespace CspinDMFT::Contour {
using IndexPair=std::array<size_t,2>;
std::vector<IndexPair> direction_pairs(char symmetry);
std::vector<size_t> magnetization_directions(char symmetry);

// Preparation has imaginary time only. Both classes store only the selected
// spin-direction basis and the configured representative site pairs.
class MatsubaraCorrelations {
public:
    size_t sites,imaginary_points;
    char symmetry;
    std::vector<IndexPair> directions,categories;
    std::vector<C> values;
    MatsubaraCorrelations(size_t sites,size_t imaginary_points,char symmetry,const std::vector<IndexPair>& categories={});
    C x(size_t a,size_t b,size_t tau) const;
    void set(size_t a,size_t b,size_t tau,C value);
private:
    std::vector<size_t> category_indices;
};
class ContourCorrelations {
public:
    size_t times,sites,imaginary_points;
    char symmetry;
    std::vector<IndexPair> directions,categories;
    std::vector<C> values;
    ContourCorrelations(size_t times,size_t sites,size_t imaginary_points,char symmetry,const std::vector<IndexPair>& categories={});
    C x(size_t t,size_t a,size_t b,size_t tau) const;
    void set(size_t t,size_t a,size_t b,size_t tau,C value);
private:
    std::vector<size_t> category_indices;
};
class ClusterMagnetization {
public:
    size_t sites;
    char symmetry;
    std::vector<size_t> directions;
    std::vector<C> values;
    ClusterMagnetization(size_t sites,char symmetry);
    C m(size_t a) const;
    void set(size_t a,C value);
};
class TimeMagnetization {
public:
    size_t times,sites;
    char symmetry;
    std::vector<size_t> directions;
    std::vector<C> values;
    TimeMagnetization(size_t times,size_t sites,char symmetry);
    C m(size_t t,size_t a) const;
    void set(size_t t,size_t a,C value);
};
// Shared kernels use these overloads; the preparation containers themselves
// expose no real-time index, size, or placeholder real-time slice.
inline size_t time_points(const MatsubaraCorrelations&) { return 1; }
inline size_t time_points(const ClusterMagnetization&) { return 1; }
inline size_t time_points(const ContourCorrelations& x) { return x.times; }
inline size_t time_points(const TimeMagnetization& x) { return x.times; }
C correlation(const MatsubaraCorrelations&,size_t,size_t,size_t,size_t);
C correlation(const ContourCorrelations&,size_t,size_t,size_t,size_t);
void set_correlation(MatsubaraCorrelations&,size_t,size_t,size_t,size_t,C);
void set_correlation(ContourCorrelations&,size_t,size_t,size_t,size_t,C);
C magnetization(const ClusterMagnetization&,size_t,size_t);
C magnetization(const TimeMagnetization&,size_t,size_t);
void set_magnetization(ClusterMagnetization&,size_t,size_t,C);
void set_magnetization(TimeMagnetization&,size_t,size_t,C);
MatsubaraCorrelations make_correlations_like(const MatsubaraCorrelations&,char symmetry);
ContourCorrelations make_correlations_like(const ContourCorrelations&,char symmetry);
ClusterMagnetization make_magnetization_like(const ClusterMagnetization&);
TimeMagnetization make_magnetization_like(const TimeMagnetization&);

enum class Sector { M,Eta,Nu };
struct Node { Sector sector; size_t point,component; };
struct Layout {
    size_t ni,intervals,components;
    size_t dimension() const;
    size_t flat(Sector,size_t,size_t) const;
    Node decode(size_t) const;
};
}
namespace CspinDMFT {
using Contour::MatsubaraCorrelations;
using Contour::ContourCorrelations;
using Contour::ClusterMagnetization;
using Contour::TimeMagnetization;
using Contour::Sector;
using Contour::Node;
using Contour::Layout;
using Contour::IndexPair;
using Contour::direction_pairs;
using Contour::magnetization_directions;
using Contour::time_points;
using Contour::correlation;
using Contour::set_correlation;
using Contour::magnetization;
using Contour::set_magnetization;
using Contour::make_correlations_like;
using Contour::make_magnetization_like;
}

namespace CspinDMFT::Contour {
template<class Correlations>
using MagnetizationFor=std::conditional_t<std::is_same_v<Correlations,MatsubaraCorrelations>,ClusterMagnetization,TimeMagnetization>;
template<class Correlations>
Correlations make_correlations(size_t nt,size_t sites,size_t ni,char symmetry,const std::vector<IndexPair>& categories={}) {
    if constexpr(std::is_same_v<Correlations,MatsubaraCorrelations>) {
        if(nt!=1) throw std::invalid_argument("Matsubara correlations have no real-time axis");
        return {sites,ni,symmetry,categories};
    } else return {nt,sites,ni,symmetry,categories};
}
template<class Magnetization>
Magnetization make_magnetization(size_t nt,size_t sites,char symmetry) {
    if constexpr(std::is_same_v<Magnetization,ClusterMagnetization>) {
        if(nt!=1) throw std::invalid_argument("thermal magnetization has no real-time axis");
        return {sites,symmetry};
    } else return {nt,sites,symmetry};
}
}
namespace CspinDMFT {
using Contour::MagnetizationFor;
using Contour::make_correlations;
using Contour::make_magnetization;
}
