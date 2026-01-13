#pragma once
#include<string>
#include<array>
#include"../Globals/Matrix_Types.h"
#include"../Matrices/Diagonalization.h"
#include"P_Error_Handling.h"
#include"../Standard_Algorithms/Numerics.h"

namespace Physics::CFET
{

// ===================== NAMESPACES ========================
namespace error = Physics::Error_Handling;
namespace num = Standard_Algorithms::Numerics;
namespace diag = Matrices::Diagonalization;

// ===================== USING STATEMENTS ==================
using BetaComponent = blaze::StaticVector<RealType,2UL,blaze::columnVector>;
using BetaComponentList = std::array<BetaComponent,3>;
using Operator = GeneralMatrix;
using Observable = HermitianMatrix;
using DiagonalOperator = DiagonalMatrix;

// ===================== INLINE VARIABLES ==================
// CFET4-opt beta-vectors
inline const BetaComponent BETA_1{ 
    RealType{0.5} * ( RealType{11.}/RealType{40.} + RealType{3.}*RealType{20.}/RealType{87.} + RealType{5.}*RealType{7.}/RealType{50.}), 
    RealType{0.5} * ( RealType{11.}/RealType{40.} - RealType{3.}*RealType{20.}/RealType{87.} + RealType{5.}*RealType{7.}/RealType{50.}) };

inline const BetaComponent BETA_2{ 
    RealType{0.5} * ( RealType{9. }/RealType{20.} - RealType{5.}*RealType{7. }/RealType{25.}),
    RealType{0.5} * ( RealType{9. }/RealType{20.} - RealType{5.}*RealType{7. }/RealType{25.}) };

inline const BetaComponent BETA_3{ 
    RealType{0.5} * ( RealType{11.}/RealType{40.} - RealType{3.}*RealType{20.}/RealType{87.} + RealType{5.}*RealType{7.}/RealType{50.}),
    RealType{0.5} * ( RealType{11.}/RealType{40.} + RealType{3.}*RealType{20.}/RealType{87.} + RealType{5.}*RealType{7.}/RealType{50.}) };

inline const BetaComponentList BETA_LIST{ BETA_1, BETA_2, BETA_3 };


// ===================== FUNCTIONS =========================
/* computes the sinus cardinal for any x */
inline RealType sinc( const RealType& x )
{
    return num::is_zero(x) ? RealType{1.} : std::sinh(x) / x;
}


/* computes the short step propagator for a single spin 1/2 according to the Hamiltonian H(t) = Q(t) * S, considering CFET 4 opt 
the propagator is computed by U_step = exp(-i F_1 * sigma) exp(-i F_2 * sigma) exp(-i F_3 * sigma)
the F-vectors are computed by F_i = dt/2 * T * betai 
the matrix exponentials are evaluated using exp(- F*sigma) = cosh(|F|) sigma^0 -  sinc(|F|) F * sigma */
inline Operator CFET4opt_for_single_spin_one_half( const FieldVector& Q_new, const FieldVector& Q_old, const RealType& dt, const Observable& S_X, const Observable& S_Y, const Observable& S_Z, const Observable& IDENTITY )
{
    // A.) Construct the T Matrix
    blaze::StaticMatrix<RealType,3UL,2UL> T{}; // T = (Q_new Q_old)
    T( 0, 0 ) = Q_new[0], T( 0, 1 ) = Q_old[0]; // T_x 
    T( 1, 0 ) = Q_new[1], T( 1, 1 ) = Q_old[1]; // T_y
    T( 2, 0 ) = Q_new[2], T( 2, 1 ) = Q_old[2]; // T_z

    // B.) Compute the Matrix Exponentials
    std::array<Operator,3> Ms{}; // matrix exponent vectors
    std::transform( BETA_LIST.cbegin(), BETA_LIST.cend(), Ms.begin(), [&]( const BetaComponent& beta )
    {
        RealSpaceVector F = RealType{0.5} * dt * T * beta; // matrix exponent vector, 1/2 is due to conversion S -> sigma
        RealType Fnorm = blaze::norm( F );
        return std::cosh(Fnorm) * IDENTITY  -  RealType{2.0} * sinc( Fnorm ) * ( F[0]*S_X + F[1]*S_Y + F[2]*S_Z ); // factor 2 due to conversion sigma -> S
    } );
    
    // C.) Compute and return the Short-Step Propagator 
    return Ms[0] * Ms[1] * Ms[2];
}


/* computes the short step propagator considering CFET 2
the propagator is computed by U_step = exp( -i * dt * 0.5*( H_new + H_old ) ) */
inline Operator CFET2( const Observable& H_new, const Observable& H_old, const RealType& dt, const Observable& ZERO, const Observable& IDENTITY, const std::string matrix_exponential_computation = "diagonalization" )
{
    Observable Exponent = RealType{0.5} * ( H_new + H_old ); // CFET2

    if( matrix_exponential_computation == "taylor" ) // computation via truncated taylor expansion (not recommended, large numerical errors)
    {
        return IDENTITY  +  ComplexType( 0.0, -dt ) * Exponent  +  ComplexType( -dt*dt /2.0 , 0.0 ) * Exponent * Exponent; 
    }
    else if( matrix_exponential_computation == "diagonalization" ) // computation via diagonalization
    {
        // diagonalize:
        Operator UnitaryTransformation{};
        blaze::DynamicVector<RealType,blaze::columnVector> EigenValues{};
        Observable Exponent_DENSE{Exponent};
        diag::diagonalize_cplx( Exponent_DENSE, EigenValues, UnitaryTransformation );

        // compute exp(-i dt D) with D being the diagonalized exponent:
        DiagonalOperator DiagonalExponentMatrix{ ZERO };
        std::for_each( EigenValues.begin(), EigenValues.end(), [&,index=0]( const RealType& value ) mutable
        { 
            DiagonalExponentMatrix(index,index) = exp( ComplexType(RealType{0.0},RealType{-1.0}*dt) * value ); 
            index++;
        } );

        // transpose:
        return ctrans( UnitaryTransformation ) * DiagonalExponentMatrix * UnitaryTransformation;
    }
    else
    {
        error::MATRIX_COMPUTATION( __PRETTY_FUNCTION__, matrix_exponential_computation );
        return Operator{}; // suppress warnings
    }
}

};