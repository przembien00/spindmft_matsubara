#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Multivariate_Gaussian::Error_Handling
{

// rounds a value to a string with custom precision
template<typename T>
std::string round_value_to_string( const T& val, const size_t num_Digits )
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision( num_Digits ) << val;
  return ss.str();
}

inline void UNDEFINED_TYPE( const std::string function_name )
{
    throw std::runtime_error( "caught undefined type in " + function_name );
}

template<typename T>
void NEGATIVE_VARIANCE( const T& variance, const std::string function_name )
{
    throw std::runtime_error( "caught negative variance sigma^2 = \'" + round_value_to_string( variance, 2 ) + "\' in " + function_name );
}

inline void INVALID_SIZE( const size_t size, const std::string function_name )
{
    throw std::runtime_error( "caught invalid container size \'" + std::to_string( size ) + "\' in " + function_name );
}

inline void SYMMETRY_TYPE( const char symmetry_type, const std::string function_name )
{
    throw std::runtime_error( "symmetry type \'" + std::string{symmetry_type} + "\' not implemented in " + function_name );
}

}