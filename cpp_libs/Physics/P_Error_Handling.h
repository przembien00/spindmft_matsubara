#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Physics::Error_Handling
{

using uint = unsigned int;

// rounds a value to a string with custom precision
template<typename T>
std::string round_value_to_string( const T& val, const uint num_Digits )
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision( num_Digits ) << val;
  return ss.str();
}

inline void SPIN_MODEL( const std::string& name, const std::string function_name )
{
    throw std::runtime_error( "invalid spin model name \'" + name + "\' in " + function_name );
}

inline void MAGNETIC_FIELD( const std::string& name, const std::string function_name )
{
    throw std::runtime_error( "invalid magnetic field name \'" + name + "\' in " + function_name );
}

inline void CORRELATION( const std::string& name, const std::string function_name )
{
    throw std::runtime_error( "invalid correlation name \'" + name + "\' in " + function_name );
}

inline void NOISE( const std::string name, const std::string function_name )
{
    throw std::runtime_error( "noise with name " + name + " not implemented in " + function_name );
}

inline void EXTRAINTERACTION( const std::string name, const std::string function_name )
{
    throw std::runtime_error( "extra interaction with name " + name + " not implemented in " + function_name );
}

inline void MATRIX_COMPUTATION( const std::string function_name, const std::string matrix_exponential_computation )
{
    throw std::runtime_error( "matrix exponential computation \"" + matrix_exponential_computation + "\" unknown in " + function_name );
}

//void INVALID_SIZE( const uint size, const std::string function_name )
//{
//    throw std::runtime_error( "caught invalid container size \'" + std::to_string( size ) + "\' in " + function_name );
//}


};