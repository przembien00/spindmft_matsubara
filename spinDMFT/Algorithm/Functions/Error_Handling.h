#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace spinDMFT::Functions::Error_Handling
{

// rounds a value to a string with custom precision
template<typename T>
std::string round_value_to_string( const T& val, const uint num_Digits )
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision( num_Digits ) << val;
  return ss.str();
}

inline void SYMMETRY_TYPE( const char symmetry_type, const std::string function_name )
{
    throw std::runtime_error( "symmetry type \'" + std::string{symmetry_type} + "\' invalid in " + function_name );
}

inline void IMAGINARY_TRACE( const RealType imaginary_part, const std::string function_name )
{
    throw std::runtime_error( "caught large imaginary part \'" + round_value_to_string( imaginary_part, 10 ) + "\' in a trace in " + function_name );
}

inline void SELF_CONSISTENCY_MODEL( const std::string sc_model, const std::string function_name )
{
    throw std::runtime_error( "self-consistency model \'" + std::string{sc_model} + "\' not implemented in " + function_name );
}


}