#pragma once

#include"Types.h"
#include<blaze/Math.h>

// Matrices on Real Space
using Rotation = blaze::StaticMatrix<RealType,3UL,3UL>;
using FieldMatrix = Rotation;
using DiagonalRotation = blaze::DiagonalMatrix<blaze::StaticMatrix<RealType,3UL,3UL>>;
using GeneralRealMatrix = blaze::DynamicMatrix<RealType,blaze::rowMajor>;

// Vectors on Real Space
using RealSpaceVector = blaze::StaticVector<RealType,3UL,blaze::columnVector>;
using FieldVector = RealSpaceVector;

// Hilbert Space Matrices : Dense
using GeneralMatrix = blaze::DynamicMatrix<ComplexType,blaze::rowMajor>;
using HermitianMatrix = blaze::HermitianMatrix<GeneralMatrix>;
using DiagonalMatrix = blaze::DiagonalMatrix<GeneralMatrix>;

// Hilbert Space Matrices : Sparse
using SparseGeneralMatrix = blaze::CompressedMatrix<ComplexType,blaze::rowMajor>;
using SparseHermitianMatrix = blaze::HermitianMatrix<SparseGeneralMatrix>;