#pragma once

#include<string>
#include<vector>
#include<iostream>
#include"../Global/Matrices_Global.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"spinDMFT_Tensors/Tensors.h"

namespace SpinED::Models
{


namespace ps = Parameter_Space;
namespace ten = spinDMFT::Tensors;

using Rotation = blaze::StaticMatrix<RealType,3UL,3UL>;
using MeanTen = ten::MeanTensor<std::vector<RealType>>;
using CorrTen = ten::CorrelationTensor<std::vector<RealType>>;
using state = std::vector<uint>;

// ======================== CLASS DEFINITIONS ========================
// REALTYPE PARAMETER CLASS
struct Parameter
{
    Parameter( const std::string& name, const RealType& value ) : m_name(name), m_value(value) {};
    const std::string m_name{};
    RealType m_value{};
};


// MODEL BASE CLASS
/* explain the difference between real and complex (e.g. for Sy)
*/
class Model 
{
 public:
    virtual ~Model() = default;

    // PUBLIC MEMBERS
    MeanTen m_means;
    CorrTen m_correlations;
    const std::string m_name;

    // PURELY VIRTUAL PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const = 0;
    virtual std::string init_state_info( const uint num_Digits ) const = 0;
    virtual std::vector<Parameter> return_params() const = 0;
    virtual void compute_Hamiltonian() = 0;
    virtual void compute_density_operator() = 0;
    virtual void diagonalize() = 0;
    virtual void transform_to_diagonal_basis() = 0;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) = 0;

 protected:
     // CONSTRUCTORS
    Model( const std::string& name, const char mean_symmetry_type, const char corr_symmetry_type,
        const ps::ParameterSpace& pspace, const bool H_is_real ):
        m_means( MeanTen{mean_symmetry_type} ),
        m_correlations( CorrTen{corr_symmetry_type} ),
        m_name( name ),
        num_Spins( pspace.num_Spins ),
        num_HilbertSpaceDimension( pspace.num_HilbertSpaceDimension ),
        J( pspace.J ),
        m_H_is_real( H_is_real )
    {}

    // PROTECTED MEMBERS
    const uint num_Spins;
    const uint num_HilbertSpaceDimension;
    const SpinSpinCouplings J;
    const bool m_H_is_real; // Hamiltonian is real-valued
};


// ISO DISORDERED MODEL DERIVED CLASS
class ISO_Disordered_Model : public Model
{
 public:
    // CONSTRUCTOR
    ISO_Disordered_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    RealObservable H{};         // (real-symmetric) Hamiltonian
    RealType rho{};             // rho -> 1/d 
    EigenValues omega{};       // (real) Eigenvalues of H
    RealOperator K{};           // orthogonal transformation that diagonalizes H
    RealDiagonalOperator H_D{}; // (real) diagonalized Hamiltonian
    RealObservable Sz_D{};      // (real) S1z in the diagonal basis
};

// ISO DISORDERED BLOCKWISE MODEL DERIVED CLASS
class ISO_Disordered_Blockwise_Model : public Model
{
 public:
    // CONSTRUCTOR
    ISO_Disordered_Blockwise_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    std::vector<RealObservable> H{};            // blockwise (real-symmetric) Hamiltonian
    const uint num_blocks{};
    const uint half_blocks{}; // due to symmetry only the blocks with zero or negativ total magnetization need to be considered
    std::vector<std::vector<state>> states{};   // blockwise basis states 
    RealType rho{};                             // rho -> 1/d 
    std::vector<EigenValues> omega{};           // blockwise (real) Eigenvalues of H
    std::vector<RealOperator> K{};              // blockwise orthogonal transformation that diagonalizes H
    std::vector<RealDiagonalOperator> H_D{};    // blockwise (real) diagonalized Hamiltonian
    std::vector<SparseRealObservable> Sz_D{};   // blockwise (real) S1z in the diagonal basis
};



// ISO POLARIZED MODEL DERIVED CLASS
class ISO_Polarized_Model : public Model
{
 public:
    // CONSTRUCTOR
    ISO_Polarized_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    RealObservable H{};         // (real-symmetric) Hamiltonian
    RealDiagonalOperator rho{}; // rho -> contains only Sz_tot which is diagonal
    EigenValues omega{};       // (real) Eigenvalues of H
    RealOperator K{};           // orthogonal transformation that diagonalizes H
    RealDiagonalOperator H_D{}; // (real) diagonalized Hamiltonian
    RealObservable rho_D{};     // (real) rho in the diagonal basis
    RealObservable Sx_D{};      // (real) S1x in the diagonal basis 
    RealOperator   Sy_D{};      // (real antisymmetric) S1y in the diagonal basis (i is left out)
    RealObservable Sz_D{};      // (real) S1z in the diagonal basis 
    RealOperator Sx_D_rho_D{};  // (real non-symmetric ) product of Sx and rho
    RealOperator Sz_D_rho_D{};  // (real non-symmetric ) product of Sz and rho
    const RealType h_z;         // polarization field
};

// XXZ DISORDERED MODEL DERIVED CLASS
class XXZ_Disordered_Model : public Model
{
 public:
    // CONSTRUCTOR
    XXZ_Disordered_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    RealObservable H{};         // (real-symmetric) Hamiltonian
    RealType rho{};             // rho -> 1/d 
    EigenValues omega{};       // (real) Eigenvalues of H
    RealOperator K{};           // orthogonal transformation that diagonalizes H
    RealDiagonalOperator H_D{}; // (real) diagonalized Hamiltonian
    RealObservable Sx_D{};      // (real) S1x in the diagonal basis 
    RealOperator   Sy_D{};      // (real antisymmetric) S1y in the diagonal basis (i is left out)
    RealObservable Sz_D{};      // (real) S1z in the diagonal basis
    const RealType lambda;      // anisotropy prefactor
};

// DRF DISORDERED MODEL DERIVED CLASS
class DRF_Disordered_Model : public Model
{
 public:
    // CONSTRUCTOR
    DRF_Disordered_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    RealObservable H{};         // (real-symmetric) Hamiltonian
    RealType rho{};             // rho -> 1/d 
    EigenValues omega{};       // (real) Eigenvalues of H
    RealOperator K{};           // orthogonal transformation that diagonalizes H
    RealDiagonalOperator H_D{}; // (real) diagonalized Hamiltonian
    RealObservable Sx1_D{};      // (real) Sx1 in the diagonal basis 
    RealObservable Sx2_D{};      // (real) Sx2 in the diagonal basis 
    RealObservable Sx3_D{};      // (real) Sx3 in the diagonal basis 
};

// DRF DISORDERED BLOCKWISE MODEL DERIVED CLASS
class DRF_Disordered_Blockwise_Model : public Model
{
 public:
    // CONSTRUCTOR
    DRF_Disordered_Blockwise_Model( const ps::ParameterSpace& pspace );

    // PUBLIC METHODS
    virtual std::string compact_info( const uint num_Digits ) const override;
    virtual std::string init_state_info( const uint num_Digits ) const override;
    virtual std::vector<Parameter> return_params() const override;
    virtual void compute_Hamiltonian() override;
    virtual void compute_density_operator() override;
    virtual void diagonalize() override;
    virtual void transform_to_diagonal_basis() override;
    virtual void compute_means_and_autocorrelations_at( const RealType& time ) override;
    
 private:
    // PRIVATE MEMBERS
    std::vector<RealObservable> H{};            // blockwise (real-symmetric) Hamiltonian
    const uint num_blocks{};
    const uint half_blocks{}; // due to symmetry only the blocks with 0 or (-) total magnetization required for some calculations
    std::vector<std::vector<state>> states{};   // blockwise basis states 
    RealType rho{};                             // rho -> 1/d 
    std::vector<EigenValues> omega{};           // blockwise (real) Eigenvalues of H
    std::vector<RealOperator> K{};              // blockwise orthogonal transformation that diagonalizes H
    std::vector<RealDiagonalOperator> H_D{};    // blockwise (real) diagonalized Hamiltonian
    std::vector<SparseRealObservable> Sz_D{};   // blockwise (real) S1z in the diagonal basis
    SparseRealOperator Sp_D{};                  // (real) S1+ in the diagonal basis
};



}