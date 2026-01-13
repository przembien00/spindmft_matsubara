// Define Complex Type and Global Matrix Types
#pragma once

#include"RealType_Global.h"
#include<complex>
#include<blaze/Math.h>

typedef std::complex<RealType> ComplexType;

// Quantum Operators
typedef blaze::HermitianMatrix<blaze::DynamicMatrix<ComplexType,blaze::rowMajor>> Observable;
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<RealType,blaze::rowMajor>> RealObservable;
typedef blaze::DynamicMatrix<ComplexType,blaze::rowMajor> Operator;
typedef blaze::DynamicMatrix<RealType,blaze::rowMajor> RealOperator;

// Sparse Quantum Operators
typedef blaze::DiagonalMatrix<Operator> DiagonalOperator;
typedef blaze::DiagonalMatrix<RealOperator> RealDiagonalOperator;
typedef blaze::HermitianMatrix<blaze::CompressedMatrix<ComplexType,blaze::rowMajor>> SparseObservable;
typedef blaze::SymmetricMatrix<blaze::CompressedMatrix<RealType,blaze::rowMajor>> SparseRealObservable;
typedef blaze::CompressedMatrix<ComplexType,blaze::rowMajor> SparseOperator;
typedef blaze::CompressedMatrix<RealType,blaze::rowMajor> SparseRealOperator;

// Others
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<RealType,blaze::rowMajor>> SpinSpinCouplings;
typedef blaze::DynamicVector<RealType,blaze::columnVector> EigenValues;
