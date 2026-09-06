#pragma once

#include<array>
#include<random>
#include<tuple>
#include<vector>
#include<Globals/Types.h>
#include<Standard_Algorithms/Standard_Algorithms.h>
#include<Observables/Magnetization.h>
#include"Complex_Gaussian.h"
#include"../Contour/Contour_Types.h"
#include"../Contour/Contour_Kernel.h"
#include"../matrices.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace spinDMFT::Functions
{

namespace contour = spinDMFT::Contour;
namespace mag = Observables::Magnetization;
namespace ps = spinDMFT::Parameter_Space;
namespace rtd = spinDMFT::Run_Time_Data;

using Corr = contour::Corr;
using CorrTen = contour::CorrTen;
using CorrelationSet = contour::CorrelationSet;
using MagVec = mag::MagnetizationVector;
using OperatorTrajectories = std::array<std::vector<Operator>,3>;
using ComplexMagnetizationTrajectory = std::vector<ComplexFieldVector>;
using MeanFieldTrajectory = std::vector<FieldVector>;

struct SelfConsistentField
{
    MeanFieldTrajectory mean_time{};
    ComplexDynamicMatrix covariance{};
    RealType covariance_symmetry_error{};
    RealType branch_identity_error{};
};

void initialize_matrices( const ps::ParameterSpace& pspace );

CorrelationSet generate_initial_correlations( const ps::ParameterSpace& pspace,
                                              const FieldVector& spin_expectation );

// The stored mixed primitive is X^{ab}(t,tau)=<S_b(-i tau) S_a(t)>.
// Its connected counterpart therefore subtracts the unconjugated product
// m_a(t)m_b(0). The imaginary-branch one-point function is represented by
// its equilibrium value at the contour origin.
CorrelationSet connected_contour_primitive(
    const CorrelationSet& correlations,
    const ComplexMagnetizationTrajectory& magnetization_time );

// Build E[V V^T] for [V_M(tau_k),V_+(t_n),V_-(t_n)], k=0,...,N_tau,
// retaining distinct 0+ and beta- Matsubara variables.
// The covariance is
// complex symmetric (unconjugated), and its branch blocks are reconstructed
// from the edge-grid mixed primitive before the upper triangle is mirrored.
SelfConsistentField self_consistent_equations(
    const ps::ParameterSpace& pspace,
    const CorrelationSet& correlations,
    const ComplexMagnetizationTrajectory& magnetization_time );

// Analytic mixed primitive of one thermal oscillator
// X=g(a+a^dagger), coupled to bathComponent.  Its endpoints are K lesser at
// tau=0 and K greater at tau=beta, and obey K greater(t)^*=K lesser(t).
CorrelationSet harmonic_bath_primitive( const ps::ParameterSpace& pspace );

// Fixed zero-mean field distribution obtained from the analytic harmonic-bath
// contour correlation.  It bypasses spinDMFT first- and second-moment feedback.
SelfConsistentField prescribed_harmonic_bath_field(
    const ps::ParameterSpace& pspace );

struct ContourTrajectory
{
    ComplexType partition_function{};
    ComplexType final_closed_contour_trace{};
    Operator imaginary_density_operator{};
    OperatorTrajectories imaginary_edge_insertions{};
    // One-step real-branch maps. Measurements build B_-(T,0) U_+(t,T)
    // suffixes and sweep U_+(t,0), inserting S between them at each t.
    std::vector<Operator> forward_steps{};
    std::vector<Operator> backward_steps{};
};

Operator general_matrix_exponential( const Operator& matrix );

// Exact exponential for a spin-1/2 Hamiltonian
// H=scalar*I+field dot S, S=sigma/2.  The complex dot product entering the
// Pauli identity is unconjugated, so this also covers non-Hermitian fields.
Operator spin_half_field_exponential(
    const ComplexFieldVector& field, ComplexType scalar,
    ComplexType contour_step );

// Construct the unnormalized trajectory. Sampling strategy is handled by the
// caller. pCN uses the real part of the selected observable denominator as its
// likelihood: Z_M for partition-function normalization and D(T) for fixed
// closed-contour normalization.
ContourTrajectory compute_contour_trajectory(
    const ps::ParameterSpace& pspace,
    const DenseComplexGaussianSampler::FieldVector& joint_field,
    const MeanFieldTrajectory& mean_field_time );

// Correlations and magnetization use either the closed-contour insertion
// B_-(T,0) U_+(t,T) S U_+(t,0), or the prefix insertion
// U_+(0,t) S B_-(t,0), selected by insertion_strategy.
void compute_contour_correlations( rtd::RunTimeData& rtdata,
                                   const ContourTrajectory& trajectory,
                                   RealType observable_normalization=RealType{1.},
                                   const std::string& insertion_strategy="closed-contour" );

// Treat two sign-related trajectories as one Markov observation. Their raw
// numerators are added before the shared normalization and before sample
// squares/block statistics are finalized.
void compute_contour_pair_correlations(
    rtd::RunTimeData& rtdata,
    const ContourTrajectory& positive_trajectory,
    const ContourTrajectory& negative_trajectory,
    RealType observable_normalization,
    const std::string& insertion_strategy="closed-contour" );

// Pure correlations are boundary views of the contour tensor, not independent traces:
// G_imag^{ab}(tau)=C^{ba}(0,tau), G_real^{ab}(t)=C^{ab}(t,beta).
CorrTen imaginary_time_slice( const contour::ContourCorrelation& correlations );
CorrTen real_time_slice( const contour::ContourCorrelation& correlations );

RealType max_contour_difference( const CorrelationSet& old_correlations,
                                 const CorrelationSet& new_correlations );

CorrelationSet mix_correlations( const CorrelationSet& old_correlations,
                                 const CorrelationSet& raw_correlations,
                                 RealType alpha );

ComplexMagnetizationTrajectory mix_magnetization_trajectory(
    const ComplexMagnetizationTrajectory& old_values,
    const ComplexMagnetizationTrajectory& raw_values,
    RealType alpha );

ComplexMagnetizationTrajectory project_constant_magnetization(
    const ComplexMagnetizationTrajectory& values );

RealType max_magnetization_difference(
    const ComplexMagnetizationTrajectory& old_values,
    const ComplexMagnetizationTrajectory& raw_values );

struct IterationResidual
{
    RealType absolute{};
    RealType standardized{};
};

IterationResidual iteration_residual(
    const CorrelationSet& old_correlations,
    const CorrelationSet& raw_correlations,
    const CorrelationSet& standard_errors,
    const ComplexMagnetizationTrajectory& old_magnetization,
    const ComplexMagnetizationTrajectory& raw_magnetization,
    const std::vector<FieldVector>& magnetization_Re_errors,
    const std::vector<FieldVector>& magnetization_Im_errors );

}
