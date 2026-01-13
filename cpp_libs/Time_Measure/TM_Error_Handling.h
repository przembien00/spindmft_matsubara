#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Time_Measure::Error_Handling
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

inline void TIME_MEASURE( const std::string what, const std::string function_name )
{
    throw std::runtime_error( "time measurement \'" + what + "\' undefined in " + function_name );
}

inline void TIME_STOPPED( const std::string function_name )
{
    throw std::runtime_error( "clock was already stopped in " + function_name );
}

inline void TIME_STOPPED_IN_LOOP( const size_t loop_level, const std::string function_name )
{
    throw std::runtime_error( "clock was stopped in loop level " + std::to_string(loop_level) + " in " + function_name );
}

};