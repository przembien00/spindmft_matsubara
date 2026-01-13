#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace SpinED::Parameter_Space::Error_Handling
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

void SPIN_MODEL( const std::string& name, const std::string function_name )
{
    throw std::runtime_error( "invalid spin model name \'" + name + "\' in " + function_name );
}

}