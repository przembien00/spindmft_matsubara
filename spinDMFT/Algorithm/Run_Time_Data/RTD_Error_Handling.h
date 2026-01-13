#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace spinDMFT::Run_Time_Data::Error_Handling
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

void EIGENVALUE_TRUNCATION( const std::string scheme, const std::string function_name )
{
    throw std::runtime_error( "eigenvalue truncation scheme \'" + scheme + "\' invalid in " + function_name );
}

void CRITICAL_EIGENVALUE_RATIO( const RealType ratio, const RealType threshold_ratio, const std::string function_name )
{
    throw std::runtime_error( "terminated the iteration because the eigenvalue ratio is above the threshold : \'" + round_value_to_string(ratio, 3) + "\' > " + round_value_to_string(threshold_ratio, 3) + " in " + function_name );
}


}