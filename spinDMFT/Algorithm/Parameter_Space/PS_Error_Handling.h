#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace spinDMFT::Parameter_Space::Error_Handling
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

void IMPORT_FILE_NOT_FOUND( const std::string function_name, const std::string filename )
{
    throw std::runtime_error( "import file \"" + filename + "\" not found" );
}

void INIT_CORRELATIONS_MISMATCH( const std::string function_name )
{
    throw std::runtime_error( "data of read initial data do not agree with the new data in " + function_name );
}

void INIT_CORRELATIONS_CONTAIN_NANS( const std::string function_name )
{
    throw std::runtime_error( "initial correlations contain NaN's in " + function_name );
}

void SRC_DIRECTORY( const std::string src_dir, const std::string function_name )
{
    throw std::runtime_error( "source directory " + src_dir + " invalid in " + function_name );
}

};