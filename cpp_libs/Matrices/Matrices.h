#pragma once

#include<blaze/Math.h>
#include"../Standard_Algorithms/Numerics.h"

namespace Matrices 
{
    
// ===================== NAMESPACES ==============
namespace num = Standard_Algorithms::Numerics;


// ===============================================
// ========== CASTING AND TRUNCATION =============
// ===============================================
// count numerical zeros of a matrix by std::numerical_limits epsilon
template<typename M>
size_t count_zeros( const M& matrix )
{
    size_t counter{};
    for( size_t i=0; i<matrix.rows(); ++i )
    {
        std::for_each( matrix.cbegin(i), matrix.cend(i), [&]( const auto& x )
        {
            if( num::is_zero(x) )
            {
                counter++;
            }
        } );
    }
    return counter;
}

// count numerical zeros of a matrix by given tolerance
template<typename M, typename T>
size_t count_zeros( const M& matrix, const T& epsilon)
{
    size_t counter{};
    for( size_t i=0; i<matrix.rows(); ++i )
    {
        std::for_each( matrix.cbegin(i), matrix.cend(i), [&]( const auto& x )
        {
            if( num::is_zero(x,epsilon) )
            {
                counter++;
            }
        } );
    }
    return counter;
}

/* cast blaze dense matrix to sparse
THIS FUNCTION HAS NOT YET BEEN PROPERLY TESTED
IT MAY NOT WORK FOR COLUMN MAJOR MATRICES */
template<typename SM, typename M>
SM cast_to_sparse( const M& matrix )
{
    SM sparse_matrix{};
    sparse_matrix.resize(matrix.rows(),matrix.columns());
    for( size_t i=0; i<matrix.rows(); ++i )
    {
        std::for_each( matrix.cbegin(i), matrix.cend(i), [&,j=0]( const auto& x ) mutable
        {
            if( !num::is_zero(x) )
            {
                sparse_matrix(i,j++) = x;
            }
        } );
    }

    sparse_matrix.shrinkToFit();
    return sparse_matrix;
}

// truncate a sparse matrix via std::numerical_limits epsilon
template<typename SM>
inline void truncate_sparse( SM& matrix )
{
    matrix.erase( []( auto x ){ return num::is_zero(x); } );
    matrix.shrinkToFit();
}

// truncate a sparse matrix via given tolerance
template<typename SM, typename T>
inline void truncate_sparse( SM& matrix, const T& epsilon )
{
    matrix.erase( [&epsilon]( auto x ){ return num::is_zero(x,epsilon); } );
    matrix.shrinkToFit();
}


// ===============================================
// ============= TENSOR PRODUCTS =================
// ===============================================
/* compute a tensor product of two square matrices
the method rows() must be available */
template<typename M>
M tensor_product( const M& A, const M& B )
{
    uint dimA = A.rows();
    uint dimB = B.rows();
    M C{};
    C.resize( dimA*dimB );
    for( uint rowA=0; rowA<dimA; rowA++ )
    {
        for( uint colA=0; colA<dimA; colA++ )
        {
            for( uint rowB=0; rowB<dimB; rowB++ )
            {
                for( uint colB=0; colB<dimB; colB++ )
                {
                    C( dimB*rowA + rowB, dimB*colA + colB ) = A( rowA, colA ) * B( rowB, colB );
                }
            }
        }
    }
    return C;
}

// compute a tensor product of N matrices
template<typename M>
M Nth_order_tensor_product( const std::vector<M>& local_matrices ) // copy intentional
{
    M start{local_matrices[0]};
    return std::accumulate(local_matrices.cbegin()+1, local_matrices.cend(), start, [](const M& prod, const M& matrix)
    {
        return tensor_product(prod,matrix);
    });
}

};