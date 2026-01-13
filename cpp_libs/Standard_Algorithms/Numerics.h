#pragma once

#include<vector>
#include<cmath>
#include<complex>
#include"Error_Handling.h"

namespace Standard_Algorithms::Numerics
{

namespace error = Standard_Algorithms::Error_Handling;

using uint = unsigned int;

// =========================================================
// ================ NUMERIC EXTRAPOLATION ==================
// =========================================================
/* extrapolates an std::vector linearly to a new equidistant discretization 
assumes that both discretizations start at 0 
sets everything behind the last data point to zero */
template<typename RType>
std::vector<RType> extrapolate( const std::vector<RType>& old_v, const uint new_size, const RType old_dx, const RType new_dx )
{
    const uint old_size = old_v.size();
    std::vector<RType> new_v{};
    for( uint new_n = 0; new_n < new_size; ++new_n )
    {
        RType new_x = static_cast<RType>(new_n) * new_dx;
        uint old_lower_n = std::floor( new_x / old_dx );
        uint old_upper_n = old_lower_n + 1;

        if( old_upper_n >= old_size ) // behind the last data point
        {
            new_v.emplace_back( RType{0.} );
        }
        else // within the data
        {
            RType old_lower_y = old_v[old_lower_n];
            RType old_upper_y = old_v[old_upper_n];
            RType m = (old_upper_y - old_lower_y) / old_dx;
            RType b = old_upper_y - old_dx * static_cast<RType>(old_upper_n) * m;
            new_v.emplace_back( m * new_x + b );
        }
    }
    return new_v;
}


// =========================================================
// ================= NUMERIC AVERAGING =====================
// =========================================================
template<typename T>
double mean( const T& v )
{
    return std::accumulate(v.cbegin(),v.cend(),double{})/static_cast<double>(v.size());
}

template<typename T>
double stddev( const T& v )
{
    double mean = std::accumulate(v.cbegin(),v.cend(),double{})/static_cast<double>(v.size());
    return std::sqrt( std::accumulate(v.cbegin(),v.cend(),double{},[&mean](auto sum, auto& c){return sum+std::pow(c-mean,2);})/static_cast<double>(v.size()) ); // careful: this stddev considers division by N not N-1 
}


// =========================================================
// ======== FLOATING POINT NUMERICS AND CASTING ============
// =========================================================
template<typename T>
bool is_equal(const T& value1, const T& value2)
{
    if constexpr (std::is_floating_point_v<T>)
    {
        return std::abs(value1 - value2) < std::numeric_limits<T>::epsilon() * std::max(std::abs(value1), std::abs(value2));
    }
    else
    {
        return value1 == value2;
    }
}

template<typename T>
bool is_equal_soft(const T& value1, const T& value2)
{
    if constexpr (std::is_floating_point_v<T>)
    {
        return std::abs(value1 - value2) < T{1000.} * std::numeric_limits<T>::epsilon() * std::max(std::abs(value1), std::abs(value2));
    }
    else
    {
        return value1 == value2;
    }
}

template<typename T1, typename T2>
bool is_equal(const T1& value1, const T1& value2, const T2& epsilon)
{
    if constexpr (std::is_floating_point_v<T1>)
    {
        return std::abs(value1 - value2) < std::abs(epsilon);
    }
    else
    {
        return value1 == value2;
    }
}

/* checks, whether a value is numerically zero
yields true if smaller than epsilon from std::numerical_limits, e.g., epsilon = 2.22045e-16 for double */
template<typename T>
bool is_zero(const T& value)
{
    if constexpr (std::is_floating_point_v<T>)
    {
        return std::abs(value) < std::numeric_limits<T>::epsilon();
    }
    else
    {
        return value == T{};
    }
}

/* checks, whether a value is numerically zero
yields true if smaller than 10^3*epsilon from std::numerical_limits, e.g., epsilon = 2.22045e-16 for double */
template<typename T>
bool is_zero_soft(const T& value)
{
    if constexpr (std::is_floating_point_v<T>)
    {
        return std::abs(value) < T{1000.} * std::numeric_limits<T>::epsilon();
    }
    else
    {
        return value == T{};
    }
}

/* checks, whether a value is numerically zero
yields true if smaller than a given tolerance */
template<typename T1, typename T2>
bool is_zero(const T1& value, const T2& epsilon)
{
    return std::abs(value) < epsilon;
}

/* cast a complex value to a real value if the imaginary part is smaller than epsilon from std::numerical_limits, e.g., epsilon = 2.22045e-16 for double */
template<typename T>
T cast_if_real( std::complex<T>&& value )
{
    if( !is_zero(std::imag(value)) )
    {
        error::IMAGINARY_PART( __PRETTY_FUNCTION__, std::imag(value) );
    }
    return std::real( value );
}

/* cast a complex value to a real value if the imaginary part is smaller than 10^3*epsilon from std::numerical_limits, e.g., epsilon = 2.22045e-16 for double */
template<typename T>
T cast_if_real_soft( std::complex<T>&& value )
{
    if( !is_zero_soft(std::imag(value)) )
    {
        error::IMAGINARY_PART( __PRETTY_FUNCTION__, std::imag(value) );
    }
    return std::real( value );
}

/* cast a complex value to a real value if the imaginary part is below a given tolerance
USE_TYPE needs to be defined for this function to work */
template<typename T>
T cast_if_real( std::complex<T>&& value, T&& epsilon )
{
    if( !is_zero(std::imag(value),epsilon) )
    {
        error::IMAGINARY_PART( __PRETTY_FUNCTION__, std::imag(value) );
    }
    return std::real( value );
}


};
