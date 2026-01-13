#pragma once

#include<vector>
#include<functional>
#include<algorithm>
#include<iostream>
#include"../Globals/Types.h"
#include"O_Error_Handling.h"

namespace Observables::Correlations
{
// ===================== NAMESPACES =======================
namespace error = Observables::Error_Handling;


// ===================== USING STATEMENTS =================
using Vector = std::vector<RealType>;


// ========================================================
// ============= CORRELATION VECTOR HEADER ================
// ========================================================
/* contains the correlations between two time-dependent observables 
due to time translation invariance only a vector is saved (wrapper around std::vector) */
class CorrelationVector
{
 public:
    // CONSTRUCTORS
    CorrelationVector() = default;
    explicit CorrelationVector( Vector&& correlation ) : m_correlation(std::move(correlation)), m_num_TimePoints(m_correlation.size()) {}
    explicit CorrelationVector( const size_t num_TimePoints, const RealType x = RealType() ) : m_correlation(num_TimePoints,x), m_num_TimePoints(num_TimePoints) {}
    CorrelationVector( const size_t num_TimePoints, std::function<RealType(RealType)> func, RealType delta_t );
    
    // PUBLIC METHODS
    RealType& at( const size_t time ){ return m_correlation[time]; }
    const RealType& at( const size_t time ) const { return m_correlation[time]; }
    const size_t get_num_TimePoints() const { return m_num_TimePoints; }
    const size_t size() const { return m_num_TimePoints; }
    RealType* data() { return m_correlation.data(); }
    void print( const size_t my_rank = 0 ) const;
    const Vector& get_correlation_vector() const { return m_correlation; }
    CorrelationVector zero_according_to() const;
    
    // ITERATORS
    auto begin() { return m_correlation.begin(); }
    auto end() { return m_correlation.end(); }
    auto cbegin() const { return m_correlation.cbegin(); }
    auto cend() const { return m_correlation.cend(); }
    auto rbegin() { return m_correlation.rbegin(); }
    auto rend() { return m_correlation.rend(); }
    auto crbegin() const { return m_correlation.crbegin(); }
    auto crend() const { return m_correlation.crend(); }

    // OPERATORS 
    RealType& operator[]( const size_t index ) { return m_correlation[index]; };
    const RealType& operator[]( const size_t index ) const { return m_correlation[index]; };

    CorrelationVector& operator=(const Vector& other);
    CorrelationVector& operator=(const CorrelationVector& other);
    CorrelationVector& operator*=(const RealType& factor);
    CorrelationVector& operator+=(const CorrelationVector& other);

    // FRIENDS
    friend std::ostream& operator<<( std::ostream& os, const CorrelationVector& CV );
    friend CorrelationVector operator*( const RealType& factor, const CorrelationVector& v );
    friend CorrelationVector operator+( const CorrelationVector& lhs, const CorrelationVector& rhs );
    friend CorrelationVector operator-( const CorrelationVector& lhs, const CorrelationVector& rhs );
    friend CorrelationVector operator+( const CorrelationVector& C, const RealType& x );

 private:
    // PRIVATE MEMBERS
    Vector m_correlation{};
    size_t m_num_TimePoints{};
};


// ========================================================
// ========== CORRELATION VECTOR IMPLEMENTATION ===========
// ========================================================
// constructor that resizes and fills the elements from a function that depends solely on the time difference
inline CorrelationVector::CorrelationVector( const size_t num_TimePoints, std::function<RealType(RealType)> func, RealType delta_t ):
    m_num_TimePoints( num_TimePoints )
{
    m_correlation.resize( m_num_TimePoints );
    for( size_t t = 0; t < m_num_TimePoints; ++t )
    {
        m_correlation[t] = func( (RealType) t * delta_t  );
    }
}
 
// print
inline void CorrelationVector::print( const size_t my_rank ) const
{
    if( my_rank == 0 )
    {
        for( auto v : m_correlation )
        {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
}

// return a CorrelationVector of zero's with the same size
inline CorrelationVector CorrelationVector::zero_according_to() const
{
    return CorrelationVector(m_num_TimePoints);
}

// Operators:
// copy construct from correlation vector
inline CorrelationVector& CorrelationVector::operator=(const CorrelationVector& other)
{
    m_correlation = other.get_correlation_vector();
    m_num_TimePoints = other.size();
    return *this;
}

// copy construct from normal vector
inline CorrelationVector& CorrelationVector::operator=(const Vector& other)
{
    m_correlation = other;
    m_num_TimePoints = other.size();
    return *this;
}

// multiply assign with factor
inline CorrelationVector& CorrelationVector::operator*=(const RealType& factor)
{
    std::transform( m_correlation.begin(), m_correlation.end(), m_correlation.begin(), 
        [&factor](RealType& x){ return std::move(x)*factor; } );
    return *this;
}

// add another correlation vector
inline CorrelationVector& CorrelationVector::operator+=(const CorrelationVector& other)
{
    std::transform( m_correlation.begin(), m_correlation.end(), other.cbegin(), m_correlation.begin(), 
        [](RealType& x, const RealType& y){ return x + y; } );
    return *this;
}


// ==================== FRIENDS ====================
// outstream operator
inline std::ostream& operator<<( std::ostream& os, const CorrelationVector& CV )
{
    for( auto v : CV.m_correlation )
    {
        os << v << "\n";
    }
    return os;
}

// multiply with factor
inline CorrelationVector operator*(const RealType& factor, const CorrelationVector& v )
{
    CorrelationVector w(v.size());
    std::transform( v.cbegin(), v.cend(), w.begin(), [&factor](RealType x){ return std::move(x)*factor; } );
    return w;
}

// add two correlation vectors
inline CorrelationVector operator+( const CorrelationVector& lhs, const CorrelationVector& rhs )
{
    if( lhs.size() != rhs.size() )
    {
        error::SIZE_MISMATCH( __PRETTY_FUNCTION__ );
    }
    CorrelationVector add(lhs.size());
    std::transform( lhs.cbegin(), lhs.cend(), rhs.cbegin(), add.begin(), []( const RealType& l, const RealType& r )
    {
        return l + r;
    } );
    return add;
}

// subtract two correlation vectors
inline CorrelationVector operator-( const CorrelationVector& lhs, const CorrelationVector& rhs )
{
    if( lhs.size() != rhs.size() )
    {
        error::SIZE_MISMATCH( __PRETTY_FUNCTION__ );
    }
    CorrelationVector diff(lhs.size());
    std::transform( lhs.cbegin(), lhs.cend(), rhs.cbegin(), diff.begin(), []( const RealType& l, const RealType& r )
    {
        return l - r;
    } );
    return diff;
}

// add a scalar to a correlation vector
inline CorrelationVector operator+(const CorrelationVector& C, const RealType& x)
{
    CorrelationVector add(C.size());
    std::transform( C.cbegin(), C.cend(), add.begin(), [&x](const RealType& c){ return c + x; } );
    return add;
}



};