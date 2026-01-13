#pragma once

#include<Globals/Matrix_Types.h>

// Dense Matrix Types
using Observable = blaze::HermitianMatrix<blaze::DynamicMatrix<ComplexType,blaze::rowMajor>>;
using Operator = blaze::DynamicMatrix<ComplexType, blaze::rowMajor>;
using DiagonalOperator = blaze::DiagonalMatrix<blaze::DynamicMatrix<ComplexType, blaze::rowMajor>>;
using Matrix = blaze::DynamicMatrix<RealType, blaze::rowMajor>;
using SymmMatrix = blaze::SymmetricMatrix<Matrix>;
using RVector = blaze::DynamicVector<RealType, blaze::columnVector>;
using SymmMatrixOfVector = blaze::SymmetricMatrix<blaze::DynamicMatrix<RVector, blaze::rowMajor>>;

// Sparse Matrix Types
using SparseObservable = blaze::HermitianMatrix<blaze::CompressedMatrix<ComplexType, blaze::rowMajor>>;
using SparseOperator = blaze::CompressedMatrix<ComplexType, blaze::rowMajor>;

// Specific matrices built external at run time
inline SparseObservable ZERO{};
inline SparseObservable IDENTITY{};
inline std::vector<SparseObservable> S_X_LIST{};
inline std::vector<SparseObservable> S_Y_LIST{};
inline std::vector<SparseObservable> S_Z_LIST{};
inline SparseObservable H_CLUSTER{};