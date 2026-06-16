#pragma once

#include<Globals/Matrix_Types.h>

// Matrices on the Hilbert Space
using Observable = HermitianMatrix;
using Operator = GeneralMatrix;

// Specific matrices built external at run time
inline Observable ZERO{};
inline Observable IDENTITY{};
inline Observable S_X{};
inline Observable S_Y{};
inline Observable S_Z{};
inline Observable H_REST{};