#pragma once

#include<string>
#include"../Globals/Types.h"
#include"../Globals/Matrix_Types.h"

#include"../Matrices/Matrices.h"
namespace mat = Matrices;

namespace Physics::Spin
{

// ================ USING STATEMENTS ================
using Observable = HermitianMatrix;
using Operator = GeneralMatrix;
using SparseObservable = SparseHermitianMatrix;
using SparseOperator = SparseGeneralMatrix;


// ========= FUNCTIONS FOR A SINGLE SPIN ============
inline RealType convert_spin_string_to_float( const std::string& spin_string )
{
    auto pos = spin_string.find("/");
    if( pos != std::string::npos ) // interpret 1/2, 3/2, ...
    {
        size_t pos = spin_string.find('/');
        size_t numerator = std::stoi(spin_string.substr(0, pos));
        size_t denominator = std::stoi(spin_string.substr(pos+1, std::string::npos));
        return static_cast<RealType>(numerator) / static_cast<RealType>(denominator);
    }
    else // interpret 0.5, 1.0, 1.5, ...
    {
        return static_cast<RealType>(std::stod(spin_string));
    }
}

template<typename T>
std::string convert_spin_float_to_string( const T& spin_float )
{
    size_t two_spin = std::lround(2*spin_float);
    if( two_spin % 2 == 0 ) // integer spin
    {
        return std::to_string(std::lround(spin_float)); 
    }
    else // half integer spin
    {
        return std::to_string(two_spin) + "/2"; 
    }
}

// write the spin matrices in dependence of the spin length
inline void write_spin_matrices( const RealType& S, Observable& S_X, Observable& S_Y, Observable& S_Z )
{
    size_t dim = std::lround(2*S+1);
    S_X = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( dim, dim );
    S_Y = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( dim, dim );
    S_Z = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( dim, dim );
    for( size_t row = 1; row < dim; ++row )
    {
        S_X(row-1,row) = ComplexType( 0.5 * std::sqrt((S+1)*2*row - row*(row+1)),  0.  );
        S_Y(row-1,row) = ComplexType( 0.,  -0.5 * std::sqrt((S+1)*2*row - row*(row+1)) );
    }
    for( size_t row = 1; row <= dim; ++row )
    {
        S_Z(row-1,row-1) = ComplexType( S+1.-row,  0. );
    }
}

// create zero matrix for a single spin
inline Observable create_local_zero( const RealType& S )
{
    size_t dim = std::lround(2*S+1);
    return blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( dim, dim );
}

// create identity matrix for a single spin
inline Observable create_local_identity( const RealType& S )
{
    size_t dim = std::lround(2*S+1);
    return blaze::IdentityMatrix<ComplexType,blaze::rowMajor>( dim );
}


// ========= FUNCTIONS FOR SEVERAL SPINS ============
// return Hilbert space dimension from list of spin values
inline size_t return_Hilbert_space_dimension( const std::vector<RealType>& spin_float_list )
{
    return std::accumulate( spin_float_list.cbegin(), spin_float_list.cend(), size_t{1}, []( size_t product, const double& spin )
    {
        return product * (std::lround(2*spin)+1);
    } );
}

// return ZERO and IDENTITY matrix in dependence of Hilbert space dimension
inline std::tuple<SparseObservable,SparseObservable> create_zero_and_identity( const size_t num_HilbertSpaceDimension )
{
    SparseObservable ZERO = blaze::ZeroMatrix<ComplexType,blaze::rowMajor>( num_HilbertSpaceDimension, num_HilbertSpaceDimension );
    SparseObservable IDENTITY = blaze::IdentityMatrix<ComplexType,blaze::rowMajor>( num_HilbertSpaceDimension );  
    return std::make_tuple(ZERO,IDENTITY);
}

// create a linear spin term S^alpha_i
inline Observable create_linear_spin_term( const std::vector<RealType>& spin_float_list, size_t& index, const Observable& S_ALPHA )
{
    std::vector<Observable> local_spin_list{};
    std::for_each( spin_float_list.cbegin(), spin_float_list.cend(), [&local_spin_list,index,comp_index=0](const RealType& S) mutable
    {
        if( comp_index++ != index )
        {
            local_spin_list.emplace_back( create_local_identity(S) );
        }
    } );
    local_spin_list.insert( local_spin_list.begin()+index, S_ALPHA );
    return mat::Nth_order_tensor_product( local_spin_list );
}

};