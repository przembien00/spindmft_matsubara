#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Standard_Algorithms::Error_Handling
{

using uint = unsigned int;

template<typename RType>
void IMAGINARY_PART( const std::string function_name, const RType imag )
{
    throw std::runtime_error( "imaginary part is " + std::to_string(imag) + " and therefore to large in " + function_name );
}

}