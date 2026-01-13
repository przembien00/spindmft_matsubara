#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Mean_Field_Models::Error_Handling
{

void CONFIG_FILE_NOT_FOUND( const std::string function_name, const std::string filename )
{
    throw std::runtime_error( "config file \"" + filename + "\" not found" );
}

}