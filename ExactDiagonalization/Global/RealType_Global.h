// Define RealType and uint
// If RealType is changed MPI_REALTYPE has to be changed as well
#pragma once

typedef unsigned int uint;

#ifdef USE_DOUBLE
typedef double RealType;
#endif 

#ifdef USE_FLOAT 
typedef float RealType;
#endif
