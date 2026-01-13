#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace File_Management::Error_Handling
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

//void INVALID_SIZE( const uint size, const std::string function_name )
//{
//    throw std::runtime_error( "caught invalid container size \'" + std::to_string( size ) + "\' in " + function_name );
//}


}