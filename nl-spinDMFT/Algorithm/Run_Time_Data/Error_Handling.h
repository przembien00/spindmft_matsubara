#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Run_Time_Data::Error_Handling
{

void EIGENVALUE_TRUNCATION( const std::string scheme, const std::string function_name )
{
    throw std::runtime_error( "eigenvalue truncation scheme \'" + scheme + "\' invalid in " + function_name );
}

}