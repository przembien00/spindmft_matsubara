#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace spinDMFT::Storage_Concept::Error_Handling
{

using uint = unsigned int;

void CREATE_FILE( const std::string& filename, const std::string function_name )
{
    throw std::runtime_error( "creating file with name \'" + filename + "\' failed in " + function_name );
}


}