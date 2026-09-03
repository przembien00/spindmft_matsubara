#pragma once

#include<array>
#include<stdexcept>
#include<vector>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>

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
