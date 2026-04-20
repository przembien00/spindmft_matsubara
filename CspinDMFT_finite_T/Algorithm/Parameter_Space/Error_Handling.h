#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace DMFT_parameter_space::Error_Handling
{

void CONFIG_NOT_SET( const std::string function_name )
{
    throw std::runtime_error( "Configuration file is not set in " + function_name );
}

void IMPORT_FILE_NOT_FOUND( const std::string function_name, const std::string filename )
{
    throw std::runtime_error( "import file \"" + filename + "\" not found" );
}

void INIT_CORRELATIONS_UNKNOWN( const std::string function_name, const std::string correlation_name )
{
    throw std::runtime_error( "correlation name \'" + correlation_name + "\' unknown in " + function_name );
}

void INIT_CORRELATIONS_SRC_FOLDER( const std::string function_name, const std::string foldername )
{
    throw std::runtime_error( "initial correlations src folder \'" + foldername + "\' invalid in " + function_name );
}

void INIT_CORRELATIONS_MISMATCH( const std::string function_name )
{
    throw std::runtime_error( "data of read initial data do not agree with the new data in " + function_name );
}

void INIT_CORRELATIONS_CONTAIN_NANS( const std::string function_name )
{
    throw std::runtime_error( "initial correlations contain NaN's in " + function_name );
}

}