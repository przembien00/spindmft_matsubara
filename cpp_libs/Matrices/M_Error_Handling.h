#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Matrices::Error_Handling
{

inline void UNDEFINED_TYPE( const std::string function_name )
{
    throw std::runtime_error( "caught undefined type in " + function_name );
}

};