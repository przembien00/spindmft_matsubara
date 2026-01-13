/* if the Macro "EIGEN" is defined, the diagonalization will be done with eigen3 if the corresponding "USE_<TYPE>" is defined */
#pragma once

#include<iostream>
#include<vector>
#include<blaze/Math.h>
#include"M_Error_Handling.h"

#ifdef EIGEN
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Eigenvalues>
#ifdef USE_FLOAT
    using EigenMatrixType = Eigen::MatrixXf;
    using ComplexEigenMatrixType = Eigen::MatrixXcf;
#else // default to double
    using EigenMatrixType = Eigen::MatrixXd;
    using ComplexEigenMatrixType = Eigen::MatrixXcd;
#endif
#endif

namespace Matrices::Diagonalization
{

template<typename Matrix, typename EigVal, typename OrthoTrafo>
void diagonalize_real( const Matrix& matrix, EigVal& eigvals, OrthoTrafo& O )
{
#ifdef EIGEN // diagonalization using eigen
    // cast the real matrix from blaze to eigen
    EigenMatrixType EIGEN_matrix( matrix.rows(), matrix.rows() );
    for( size_t row = 0; row < matrix.rows(); ++row ) // copy element wise
    {
        for( size_t col = 0; col < matrix.rows(); ++col )
        {
            EIGEN_matrix( row, col ) = matrix( row, col );
        }
    }

    // diagonalize the eigen covariance matrix
    Eigen::SelfAdjointEigenSolver<EigenMatrixType> solver;
    solver.compute( EIGEN_matrix );

    // initialize eigvals and O properly
    eigvals.resize( matrix.rows() );
    O.resize( matrix.rows(), matrix.rows() );

    // copy results to eigvals and O
    std::copy( solver.eigenvalues().cbegin(), solver.eigenvalues().cend(), eigvals.begin() ); // copy eigenvalues
    for( size_t row = 0; row < matrix.rows(); ++row ) // copy eigenvectors
    {
        std::copy( solver.eigenvectors().col(row).cbegin(), solver.eigenvectors().col(row).cend(), O.begin(row) );
    }
#else // diagonalization using blaze and llapack
    blaze::eigen( matrix, eigvals, O ); // D = O M O^T
#endif
}


template<typename Matrix, typename EigVal, typename UnitaryTrafo>
void diagonalize_cplx( const Matrix& matrix, EigVal& eigvals, UnitaryTrafo& U )
{
#ifdef EIGEN // diagonalization using eigen
    // cast the real matrix from blaze to eigen
    ComplexEigenMatrixType EIGEN_matrix( matrix.rows(), matrix.rows() );
    for( size_t row = 0; row < matrix.rows(); ++row ) // copy element wise
    {
        for( size_t col = 0; col < matrix.rows(); ++col )
        {
            EIGEN_matrix( row, col ) = matrix( row, col );
        }
    }

    // diagonalize the eigen covariance matrix
    Eigen::SelfAdjointEigenSolver<ComplexEigenMatrixType> solver;
    solver.compute( EIGEN_matrix );

    // initialize eigvals and U properly
    eigvals.resize( matrix.rows() );
    U.resize( matrix.rows(), matrix.rows() );

    // copy results to eigvals and U
    std::copy( solver.eigenvalues().cbegin(), solver.eigenvalues().cend(), eigvals.begin() ); // copy eigenvalues
    for( size_t row = 0; row < matrix.rows(); ++row ) // copy eigenvectors
    {
        std::copy( solver.eigenvectors().col(row).cbegin(), solver.eigenvectors().col(row).cend(), U.begin(row) );
    }
#else // diagonalization using blaze and llapack
    blaze::eigen( matrix, eigvals, U ); // D = U M U^+
#endif
}



};