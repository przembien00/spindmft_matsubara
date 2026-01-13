#pragma once

#include<stdexcept>
#include<string>
#include<sstream>

namespace Observables::Error_Handling
{

using uint = unsigned int;

inline void SIZE_MISMATCH( const std::string function_name )
{
    throw std::runtime_error( "two correlation vectors have different sizes in " + function_name );
}

inline void SYMMETRY_TYPE( const char symmetry_type, const std::string function_name )
{
    throw std::runtime_error( "symmetry type \'" + std::string(1,symmetry_type) + "\' invalid in " + function_name );
}

inline void SYMMETRY_TYPE( const size_t number_of_correlations, const std::string function_name )
{
    throw std::runtime_error( "symmetry type with \'" + std::to_string(number_of_correlations) + "\' correlations does not exist in " + function_name );
}

inline void INVALID_INDEX( const size_t index, const std::string function_name )
{
    throw std::runtime_error( "correlation tensor index \'" + std::to_string(index) + "\' invalid in " + function_name );
}

inline void IS_ZERO( const std::string object, const std::string function_name )
{
    throw std::runtime_error( "requested object \'" + object + "\' is zero in " + function_name );
}

inline void IS_NOT_STORED( const std::string object, const std::string function_name )
{
    throw std::runtime_error( "requested object \'" + object + "\' is not stored in " + function_name );
}

inline void NOT_SAVED( const std::string object, const std::string function_name )
{
    throw std::runtime_error( "requested object \'" + object + "\' is not explicitly saved in " + function_name );
}

inline void INVALID_SIZE( const size_t size, const std::string function_name )
{
    throw std::runtime_error( "container size \'" + std::to_string(size) + "\' invalid in " + function_name );
}

}