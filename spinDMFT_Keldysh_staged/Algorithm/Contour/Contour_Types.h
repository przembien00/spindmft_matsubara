#pragma once

#include<array>
#include<limits>
#include<stdexcept>
#include<vector>
#include<utility>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Magnetization.h>

namespace spinDMFT::Contour
{

namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;

using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;

enum class Branch : size_t
{
    Matsubara = 0,
    Forward = 1,
    Backward = 2
};

struct ContourIndex
{
    Branch branch{};
    size_t point{};
    size_t component{};
};

// Checked layout of the physical field values
// [ V_M(tau_k), V_+(t_n), V_-(t_n) ], with distinct one-sided Matsubara
// endpoints k=0,...,N_tau.
class ContourLayout
{
 public:
    ContourLayout() = default;
    ContourLayout( size_t num_imaginary_points, size_t num_real_points )
        : m_num_imaginary_points(num_imaginary_points),
          m_num_real_points(num_real_points)
    {
        if( m_num_imaginary_points == 0 )
            throw std::invalid_argument("the Keldysh contour needs at least one imaginary point");
        if( m_num_real_points == 0 )
            throw std::invalid_argument("the Keldysh contour needs at least one real-time point");
    }

    size_t num_imaginary_points() const { return m_num_imaginary_points; }
    size_t num_real_points() const { return m_num_real_points; }
    size_t num_contour_points() const
    { return m_num_imaginary_points + 2*m_num_real_points; }
    size_t dimension() const { return 3*num_contour_points(); }

    size_t branch_size( Branch branch ) const
    {
        return branch == Branch::Matsubara
            ? m_num_imaginary_points : m_num_real_points;
    }
    size_t point_offset( Branch branch ) const
    {
        switch( branch )
        {
            case Branch::Matsubara: return 0;
            case Branch::Forward: return m_num_imaginary_points;
            case Branch::Backward: return m_num_imaginary_points+m_num_real_points;
        }
        throw std::logic_error("unknown contour branch");
    }
    size_t flat( Branch branch, size_t point, size_t component ) const
    {
        if( point >= branch_size(branch) || component >= 3 )
            throw std::out_of_range("Keldysh contour index out of range");
        return 3*(point_offset(branch)+point)+component;
    }
    ContourIndex decode( size_t flat_index ) const
    {
        if( flat_index >= dimension() )
            throw std::out_of_range("flat Keldysh contour index out of range");
        const size_t component = flat_index%3;
        const size_t point = flat_index/3;
        if( point < m_num_imaginary_points )
            return {Branch::Matsubara,point,component};
        if( point < m_num_imaginary_points+m_num_real_points )
            return {Branch::Forward,point-m_num_imaginary_points,component};
        return {Branch::Backward,
                point-m_num_imaginary_points-m_num_real_points,component};
    }

 private:
    size_t m_num_imaginary_points{};
    size_t m_num_real_points{};
};

// [real_time_index](real_spin, imaginary_spin)[imaginary_time_edge].
using ContourCorrelation = std::vector<CorrTen>;

enum class StaggeredSector { Matsubara, Eta, Nu };
struct StaggeredIndex
{
    StaggeredSector sector;
    size_t point,component;
};

// One joint draw [V_M(tau_k), eta(t_{n+1/2}), nu(t_m)]. Interior
// response boundaries occur only once and are reused by adjacent steps.
class StaggeredLayout
{
 public:
    StaggeredLayout(size_t imaginary_points,size_t intervals)
        : m_imaginary_points(imaginary_points),m_intervals(intervals)
    {
        const size_t limit=std::numeric_limits<size_t>::max()/3;
        if(imaginary_points==0 || intervals==0 || imaginary_points>=limit
            || intervals>(limit-imaginary_points-1)/2)
            throw std::invalid_argument("invalid staggered field grid");
    }
    size_t dimension() const { return 3*(m_imaginary_points+2*m_intervals+1); }
    size_t sector_size(StaggeredSector sector) const
    {
        return sector==StaggeredSector::Matsubara?m_imaginary_points
            :sector==StaggeredSector::Eta?m_intervals:m_intervals+1;
    }
    size_t flat(StaggeredSector sector,size_t point,size_t component) const
    {
        if(point>=sector_size(sector)||component>=3)
            throw std::out_of_range("staggered field index out of range");
        const size_t offset=sector==StaggeredSector::Matsubara?0
            :sector==StaggeredSector::Eta?m_imaginary_points:m_imaginary_points+m_intervals;
        return 3*(offset+point)+component;
    }
    StaggeredIndex decode(size_t index) const
    {
        if(index>=dimension())throw std::out_of_range("staggered field index out of range");
        const size_t point=index/3,c=index%3;
        if(point<m_imaginary_points)return {StaggeredSector::Matsubara,point,c};
        if(point<m_imaginary_points+m_intervals)
            return {StaggeredSector::Eta,point-m_imaginary_points,c};
        return {StaggeredSector::Nu,point-m_imaginary_points-m_intervals,c};
    }
 private:
    size_t m_imaginary_points,m_intervals;
};

// Complex, symmetry-reduced first moments on the real-time grid.
struct Magnetization
{
    using Tensor=Observables::Magnetization::MagnetizationTensor<
        Observables::Magnetization::MagnetizationVector>;
    Magnetization() = default;
    Magnetization(char symmetry,size_t num_real_points)
        : Re(symmetry,num_real_points),Im(symmetry,num_real_points) {}
    Magnetization(Tensor re,Tensor im) : Re(std::move(re)),Im(std::move(im)) {}
    Tensor Re{},Im{};
};

struct CorrelationSet
{
    CorrelationSet() = default;
    CorrelationSet( const char symmetry, const size_t num_imaginary_points,
                    const size_t num_real_points )
        : Re( num_real_points, CorrTen{ symmetry, num_imaginary_points } ),
          Im( num_real_points, CorrTen{ symmetry, num_imaginary_points } )
    {}

    ContourCorrelation Re{};
    ContourCorrelation Im{};
};

}
