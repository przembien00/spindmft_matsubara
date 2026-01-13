#pragma once
#include<mpi.h>

#ifdef USE_FLOAT
inline MPI_Datatype MPI_REALTYPE = MPI_FLOAT;
#else // default to double
inline MPI_Datatype MPI_REALTYPE = MPI_DOUBLE;
#endif