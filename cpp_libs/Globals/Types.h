#pragma once

typedef unsigned int uint;

#ifdef USE_FLOAT
typedef float RealType;
#else // default to double
typedef double RealType; // default
#endif

#include<complex>
typedef std::complex<RealType> ComplexType;

