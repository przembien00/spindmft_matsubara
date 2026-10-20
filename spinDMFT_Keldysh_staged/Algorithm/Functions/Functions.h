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
using MagTen = mag::MagnetizationTensor<MagVec>;
using OperatorTrajectories = std::array<std::vector<Operator>,3>;
using MeanFieldTrajectory = std::vector<FieldVector>;

struct SelfConsistentField
{
    MeanFieldTrajectory mean_time{};
    ComplexDynamicMatrix covariance{};
    CovarianceSource covariance_source{};
    RealType covariance_symmetry_error{};
};

void initialize_matrices( const ps::ParameterSpace& pspace );

std::pair<MagTen,MagTen> generate_initial_magnetization(
    const ps::ParameterSpace& pspace );

CorrelationSet generate_initial_correlations( const ps::ParameterSpace& pspace,
                                              const MagTen& magnetization_Re );

// The stored mixed primitive is X^{ab}(t,tau)=<S_b(-i tau) S_a(t)>.
// Its connected counterpart therefore subtracts the unconjugated product
// m_a(t)m_b(0). The imaginary-branch one-point function is represented by
// its equilibrium value at the contour origin.
CorrelationSet connected_contour_primitive(
    const CorrelationSet& correlations,
    const MagTen& magnetization_Re,
    const MagTen& magnetization_Im );

// Build E[V V^T] for [V_M(tau_k),V_+(t_n),V_-(t_n)], k=0,...,N_tau,
// retaining distinct 0+ and beta- Matsubara variables.
// The covariance is
// complex symmetric (unconjugated), and its branch blocks are reconstructed
// from the edge-grid mixed primitive before the upper triangle is mirrored.
SelfConsistentField self_consistent_equations(
    const ps::ParameterSpace& pspace,
    const CorrelationSet& correlations,
    const MagTen& magnetization_Re,
    const MagTen& magnetization_Im,
    bool materialize_covariance=true );

// Analytic mixed primitive of one thermal oscillator
// X=g(a+a^dagger), coupled to bathComponent.  Its endpoints are K lesser at
// tau=0 and K greater at tau=beta, and obey K greater(t)^*=K lesser(t).
CorrelationSet harmonic_bath_primitive( const ps::ParameterSpace& pspace );

// Fixed zero-mean field distribution obtained from the analytic harmonic-bath
// contour correlation.  It bypasses spinDMFT first- and second-moment feedback.
SelfConsistentField prescribed_harmonic_bath_field(
    const ps::ParameterSpace& pspace, bool materialize_covariance=true );

using SpinHalfOperator = blaze::StaticMatrix<ComplexType,2UL,2UL,blaze::rowMajor>;

// The propagation/measurement templates are instantiated in Functions.cpp for
// Operator and SpinHalfOperator. Static trajectories stay static through measurement.
template<typename Matrix>
struct BasicContourTrajectory
{
    ComplexType partition_function{};
    // NaN when final-trace evaluation is explicitly skipped.
    ComplexType final_closed_contour_trace{};
    Matrix imaginary_density_operator{};
    std::array<std::vector<Matrix>,3> imaginary_edge_insertions{};
    // One-step real-branch maps. Measurements build B_-(T,0) U_+(t,T)
    // suffixes and sweep U_+(t,0), inserting S between them at each t.
    std::vector<Matrix> forward_steps{};
    std::vector<Matrix> backward_steps{};
};

template<typename Matrix>
struct BasicTrajectoryWorkspace
{
    std::vector<ComplexFieldVector> imaginary_fields;
    std::vector<Matrix> imaginary_steps,prefix,suffix;
};

using ContourTrajectory = BasicContourTrajectory<Operator>;
using TrajectoryWorkspace = BasicTrajectoryWorkspace<Operator>;
using SpinHalfContourTrajectory = BasicContourTrajectory<SpinHalfOperator>;
using SpinHalfTrajectoryWorkspace = BasicTrajectoryWorkspace<SpinHalfOperator>;

// Immutable snapshot of one iteration's deterministic mean, including its
// propagation-node values. Own the edges as well so nodes cannot silently
// become stale when the caller updates the next iteration's mean trajectory.
class PreparedMeanField
{
 public:
    PreparedMeanField(const ps::ParameterSpace& pspace,const MeanFieldTrajectory& mean);
    void validate_grid(const ps::ParameterSpace& pspace) const;
    const MeanFieldTrajectory& edges() const { return m_edges; }
    const MeanFieldTrajectory& midpoints() const { return m_midpoints; }
 private:
    MeanFieldTrajectory m_edges,m_midpoints;
    size_t m_steps,m_substeps;
};

// An imaginary-only preparation is sufficient for partition-function pCN
// acceptance. Completion reuses those Matsubara maps after acceptance.
// Independent measurements already compute the closure trajectory, so callers
// that do not consume final_closed_contour_trace can skip its separate sweep.
enum class FinalTraceMode { Compute, Skip };
template<typename Matrix>
void build_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field,
    const MeanFieldTrajectory& mean, BasicContourTrajectory<Matrix>& result,
    BasicTrajectoryWorkspace<Matrix>& workspace, bool imaginary_only=false,
    FinalTraceMode final_trace=FinalTraceMode::Compute );
template<typename Matrix>
void complete_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field,
    const MeanFieldTrajectory& mean, BasicContourTrajectory<Matrix>& result,
    BasicTrajectoryWorkspace<Matrix>& workspace,
    FinalTraceMode final_trace=FinalTraceMode::Compute );

template<typename Matrix>
void build_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field,
    const PreparedMeanField& mean, BasicContourTrajectory<Matrix>& result,
    BasicTrajectoryWorkspace<Matrix>& workspace, bool imaginary_only=false,
    FinalTraceMode final_trace=FinalTraceMode::Compute );
template<typename Matrix>
void complete_contour_trajectory( const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field,
    const PreparedMeanField& mean, BasicContourTrajectory<Matrix>& result,
    BasicTrajectoryWorkspace<Matrix>& workspace,
    FinalTraceMode final_trace=FinalTraceMode::Compute );

Operator general_matrix_exponential( const Operator& matrix );

// Exact exponential for a spin-1/2 Hamiltonian
// H=scalar*I+field dot S, S=sigma/2.  The complex dot product entering the
// Pauli identity is unconjugated, so this also covers non-Hermitian fields.
Operator spin_half_field_exponential(
    const ComplexFieldVector& field, ComplexType scalar,
    ComplexType contour_step );

// Return exp(step*H) and exp(-step*H), sharing the Pauli coefficients.
// The second map is the actual inverse, including for complex fields/scalars.
std::pair<Operator,Operator> spin_half_field_exponential_pair(
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

ContourTrajectory compute_contour_trajectory(
    const ps::ParameterSpace& pspace,
    const JointComplexGaussianSampler::ContourFieldSample& field_sample,
    const MeanFieldTrajectory& mean_field_time );

struct MeasurementWorkspace
{
    std::vector<SpinHalfOperator> spin_half_left;
    std::vector<Operator> general_left;
    std::array<std::vector<std::array<ComplexType,4>>,3> spin_half_insertions;
};

template<typename Matrix>
void measure_contour_observables( const rtd::RunTimeData& layout,
    const BasicContourTrajectory<Matrix>& trajectory, rtd::MeasuredSample& sample,
    MeasurementWorkspace& workspace,
    const std::string& insertion_strategy="closed-contour" );

// Same estimators, streamed directly into independent total/block sums.
// Does not materialize a MeasuredSample; pCN callers must use the buffered API.
template<typename Matrix>
void measure_and_accumulate_independent(rtd::RunTimeData& data,
    const BasicContourTrajectory<Matrix>& trajectory,MeasurementWorkspace& workspace,
    const std::string& insertion_strategy="closed-contour");

// Correlations and magnetization use either the closed-contour insertion
// B_-(T,0) U_+(t,T) S U_+(t,0), or the prefix insertion
// U_+(0,t) S B_-(t,0), selected by insertion_strategy.
void compute_contour_correlations( rtd::RunTimeData& rtdata,
                                   const ContourTrajectory& trajectory,
                                   RealType observable_normalization=RealType{1.},
                                   const std::string& insertion_strategy="closed-contour" );

RealType max_contour_difference( const CorrelationSet& old_correlations,
                                 const CorrelationSet& new_correlations );

CorrelationSet mix_correlations( const CorrelationSet& old_correlations,
                                 const CorrelationSet& raw_correlations,
                                 RealType alpha );

MagTen mix_magnetization_tensor( const MagTen& old_values,
                                 const MagTen& raw_values,
                                 RealType alpha );

MagTen project_constant_magnetization( const MagTen& values );

struct IterationResidual
{
    RealType absolute{};
    RealType standardized{};
};

IterationResidual iteration_residual(
    const CorrelationSet& old_correlations,
    const CorrelationSet& raw_correlations,
    const CorrelationSet& standard_errors,
    const MagTen& old_magnetization_Re,
    const MagTen& old_magnetization_Im,
    const MagTen& raw_magnetization_Re,
    const MagTen& raw_magnetization_Im,
    const MagTen& magnetization_Re_errors,
    const MagTen& magnetization_Im_errors );

}
