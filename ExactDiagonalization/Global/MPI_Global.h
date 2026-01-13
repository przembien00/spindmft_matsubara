// Define MPI Macros
#pragma once
#include<mpi.h>
#include"RealType_Global.h"

#ifdef USE_DOUBLE
inline MPI_Datatype MPI_REALTYPE = MPI_DOUBLE;
#endif

#ifdef USE_FLOAT
inline MPI_Datatype MPI_REALTYPE = MPI_FLOAT;
#endif
