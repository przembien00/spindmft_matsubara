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

void CONFLICTING_OPTIONS( const std::string function_name, const std::string opt1, const std::string opt2 )
{
    throw std::runtime_error( "options '" + opt1 + "' and '" + opt2 + "' cannot be used together in " + function_name );
}

void NOT_BIPARTITE( const std::string function_name )
{
    throw std::runtime_error( "stagfield requires a bipartite lattice, but J graph is not 2-colorable in " + function_name );
}

void DISCONNECTED_GRAPH( const std::string function_name )
{
    throw std::runtime_error( "stagfield requires a connected lattice, but J graph has disconnected components in " + function_name );
}

}