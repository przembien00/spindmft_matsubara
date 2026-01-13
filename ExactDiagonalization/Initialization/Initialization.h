#pragma once

#include"../Global/MPI_Global.h"
#include"../Global/RealType_Global.h"
#include"../Global/Matrices_Global.h"
#include"../Global/InlineMatrices_Global.h"

namespace SpinED::Initialization
{

// initialization function : initialize ZERO and IDENTITY
// deprecated: allocation errors in case of large Hilbert spaces
/*
inline void build_trivial_spin_matrices( const uint num_Spins, const uint num_HilbertSpaceDimension )
{
  ZERO = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( num_HilbertSpaceDimension, num_HilbertSpaceDimension );
  ZERO_REAL = blaze::ZeroMatrix<RealType,blaze::rowMajor>( num_HilbertSpaceDimension, num_HilbertSpaceDimension );
  IDENTITY = blaze::IdentityMatrix<ComplexType,blaze::rowMajor>( num_HilbertSpaceDimension );
  IDENTITY_REAL = blaze::IdentityMatrix<RealType,blaze::rowMajor>( num_HilbertSpaceDimension );
}
*/

// ========================= FORWARD DECLARATIONS =========================
template<typename Obs, typename HPauli>
void build_linear_spin_operator_hermit( Obs& result, const uint num_Spins, const uint index, const HPauli& SIGMA_ALPHA );
template<typename Obs1, typename Obs2>
void compute_multi_tensor_product_hermit( Obs1& result, std::vector<Obs2>& local_operators );
template<typename Obs>
Obs compute_tensor_product_hermit( const Obs& A, const Obs& B );

template<typename Ope, typename NPauli>
void build_linear_spin_operator_nonhermit( Ope& result, const uint num_Spins, const uint index, const NPauli& SIGMA_ALPHA );
template<typename Ope1, typename Ope2>
void compute_multi_tensor_product_nonhermit( Ope1& result, std::vector<Ope2>& local_operators );
template<typename Ope>
Ope compute_tensor_product_nonhermit( const Ope& A, const Ope& B );

template<typename DiagOp, typename DPauli>
void build_linear_spin_operator_diag( DiagOp& result, const uint num_Spins, const uint index, const DPauli& SIGMA_ALPHA );
template<typename DiagOp1, typename DiagOp2>
void compute_multi_tensor_product_diag( DiagOp1& result, std::vector<DiagOp2>& local_operators );
template<typename DiagOp>
DiagOp compute_tensor_product_diag( const DiagOp& A, const DiagOp& B );


// ========================= TEMPLATE FUNCTION IMPLEMENTATIONS =========================
// ---------- HERMITIAN ----------
// initialization function : build a linear hermitian spin operator S^alpha_i
template<typename Obs, typename HPauli>
void build_linear_spin_operator_hermit( Obs& result, const uint num_Spins, const uint index, const HPauli& SIGMA_ALPHA )
{
    std::vector<Obs> v{ num_Spins-1, SIGMA_0_REAL };
    v.insert( v.begin()+index, SIGMA_ALPHA );
    compute_multi_tensor_product_hermit( result, v );
}

// function : compute a tensor product of N hermitian operators
template<typename Obs1, typename Obs2>
void compute_multi_tensor_product_hermit( Obs1& result, std::vector<Obs2>& local_operators )
{
  if( local_operators.size() == 1 ) // only one operator given -> return it
  {
    result = local_operators[0];
  }
  else // vector of operators given -> compute tensor product and return it
  {
    local_operators[0] = compute_tensor_product_hermit( local_operators[0], local_operators[1] ); // replace first operator by the tensor product
    local_operators.erase( local_operators.begin()+1 ); // erase second operator
    compute_multi_tensor_product_hermit( result, local_operators );
  }
}

// function : compute a tensor product of 2 hermitian operators
template<typename Obs>
Obs compute_tensor_product_hermit( const Obs& A, const Obs& B )
{
    uint dimA = A.rows();
    uint dimB = B.rows();
    Obs C;
    C.resize( dimA*dimB );
    for( uint rowA=0; rowA<dimA; rowA++ ){
      for( uint colA=0; colA<dimA; colA++ ){
        for( uint rowB=0; rowB<dimB; rowB++ ){
          for( uint colB=0; colB<dimB; colB++ ){
            C( dimB*rowA + rowB, dimB*colA + colB ) = A( rowA, colA ) * B( rowB, colB ); 
          }
        }
      }
    }
    return C;
}


// ---------- NON-HERMITIAN ----------
// initialization function : build a linear non-hermitian spin operator S^alpha_i (such as 1/i sigma_y)
template<typename Ope, typename NPauli>
void build_linear_spin_operator_nonhermit( Ope& result, const uint num_Spins, const uint index, const NPauli& SIGMA_ALPHA )
{
    std::vector<Ope> v{ num_Spins-1, SIGMA_0_REAL };
    v.insert( v.begin()+index, SIGMA_ALPHA );
    compute_multi_tensor_product_nonhermit( result, v );
}

// function : compute a tensor product of N non-hermitian operators
template<typename Ope1, typename Ope2>
void compute_multi_tensor_product_nonhermit( Ope1& result, std::vector<Ope2>& local_operators )
{
  if( local_operators.size() == 1 ) // only one operator given -> return it
  {
    result = local_operators[0];
  }
  else // vector of operators given -> compute tensor product and return it
  {
    local_operators[0] = compute_tensor_product_nonhermit( local_operators[0], local_operators[1] ); // replace first operator by the tensor product
    local_operators.erase( local_operators.begin()+1 ); // erase second operator
    compute_multi_tensor_product_nonhermit( result, local_operators );
  }
}

// function : compute a tensor product of 2 square non-hermitian operators
template<typename Ope>
Ope compute_tensor_product_nonhermit( const Ope& A, const Ope& B )
{
    uint dimA = A.rows();
    uint dimB = B.rows();
    Ope C;
    C.resize( dimA*dimB, dimA*dimB );
    for( uint rowA=0; rowA<dimA; rowA++ ){
      for( uint colA=0; colA<dimA; colA++ ){
        for( uint rowB=0; rowB<dimB; rowB++ ){
          for( uint colB=0; colB<dimB; colB++ ){
            C( dimB*rowA + rowB, dimB*colA + colB ) = A( rowA, colA ) * B( rowB, colB ); 
          }
        }
      }
    }
    return C;
}

// ---------- DIAGONAL ----------
// initialization function : build a linear diagonal spin operator S^alpha_i (such as sigma_z)
template<typename DiagOp, typename DPauli>
void build_linear_spin_operator_diag( DiagOp& result, const uint num_Spins, const uint index, const DPauli& SIGMA_ALPHA )
{
    std::vector<DiagOp> v{ num_Spins-1, SIGMA_0_REAL };
    v.insert( v.begin()+index, SIGMA_ALPHA );
    compute_multi_tensor_product_diag( result, v );
}

// function : compute a tensor product of N diagonal operators
template<typename DiagOp1, typename DiagOp2>
void compute_multi_tensor_product_diag( DiagOp1& result, std::vector<DiagOp2>& local_operators )
{
  if( local_operators.size() == 1 ) // only one operator given -> return it
  {
    result = local_operators[0];
  }
  else // vector of operators given -> compute tensor product and return it
  {
    local_operators[0] = compute_tensor_product_diag( local_operators[0], local_operators[1] ); // replace first operator by the tensor product
    local_operators.erase( local_operators.begin()+1 ); // erase second operator
    compute_multi_tensor_product_diag( result, local_operators );
  }
}

// function : compute a tensor product of 2 diagonal operators
template<typename DiagOp>
DiagOp compute_tensor_product_diag( const DiagOp& A, const DiagOp& B )
{
    uint dimA = A.rows();
    uint dimB = B.rows();
    DiagOp C;
    C.resize( dimA*dimB );
    for( uint rowA=0; rowA<dimA; rowA++ ){
        for( uint rowB=0; rowB<dimB; rowB++ ){
            C( dimB*rowA + rowB, dimB*rowA + rowB ) = A( rowA, rowA ) * B( rowB, rowB ); 
        }
    }
    return C;
}

};