#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Functions::Error_Handling
{

void MATRIX_COMPUTATION( const std::string function_name, const std::string matrix_exponential_computation )
{
    throw std::runtime_error( "matrix exponential computation \"" + matrix_exponential_computation + "\" unknown in " + function_name );
}

void INIT_CORRELATIONS_UNKNOWN( const std::string function_name, const std::string correlation_name )
{
    throw std::runtime_error( "initial correlation option \"" + correlation_name + "\" unknown in " + function_name );
}

void IMAGINARY_PART( const std::string function_name, const RealType imag )
{
    throw std::runtime_error( "imaginary part is " + std::to_string(imag) + " and therefore to large in " + function_name );
}

}
