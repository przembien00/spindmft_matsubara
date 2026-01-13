#pragma once

#include"../Global/RealType_Global.h"
#include"../Global/Matrices_Global.h"
#include"Diagonalization/Diagonalization.h"
#include<iostream>

namespace SpinED::Functions
{

namespace diag = Diagonalization;

// diagonalize the real Hamiltonian
inline void diagonalize_real( const RealObservable& H, EigenValues& lambda, RealOperator& U )
{
    diag::diagonalize_real( H, lambda, U );
    //blaze::eigen( H, lambda, U );
    //U = blaze::trans( U ); // after this the diagonalized matrix is given by D = Udagg M U
}


// exponentiate a diagonal matrix
template<typename DiagOp>
DiagOp exponentiate( const DiagOp& D )
{
    DiagOp E{};
    E.resize( D.rows() );
    for( uint row = 0; row < D.rows(); ++row )
    {
        E(row,row) = std::exp( D(row,row) );
    }
    return E;
}

// exponentiate a diagonal matrix times some prefactor
template<typename DiagOp, typename T>
DiagOp exponentiate( const T& prefactor, const DiagOp& D )
{
    DiagOp E{};
    E.resize( D.rows() );
    for( uint row = 0; row < D.rows(); ++row )
    {
        E(row,row) = std::exp( prefactor * D(row,row) );
    }
    return E;
}

// exponentiate a diagonal matrix of general type times some prefactor
template<typename GenOp, typename T>
GenOp exponentiate_gen( const T& prefactor, const GenOp& D )
{
    GenOp E(D.rows(),D.rows());
    for( uint row = 0; row < D.rows(); ++row )
    {
        E(row,row) = std::exp( prefactor * D(row,row) );
    }
    return E;
}

// write real eigenvalues to the diagonal of a matrix
template<typename DiagOp>
void write_to_diagonal( DiagOp& result, const EigenValues& lambda )
{
    result.resize( lambda.size() );
    for( uint row = 0; row < lambda.size(); ++row )
    {
        result(row,row) = lambda[row];
    }
}

// truncate compressed matrix
template<typename SparseOpe, typename T>
void truncate_sparse( SparseOpe& O, const T truncation )
{
    std::vector<std::tuple<uint,uint>> truncate_elements{};
    for( uint row = 0; row < O.rows(); ++row )
    {
        std::for_each( O.cbegin(row), O.cend(row),
            [&truncate_elements,&truncation,&row]( const auto& element )
        {
            if( std::abs(element.value()) < truncation )
            {
                truncate_elements.emplace_back( std::tuple<uint,uint>{row, element.index()} );
            }
        } );
    }
    std::for_each( truncate_elements.cbegin(), truncate_elements.cend(),
        [&O]( const auto& index )
    {
        O.erase( std::get<0>(index), std::get<1>(index) );
    } );
}

};