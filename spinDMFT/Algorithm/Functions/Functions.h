#pragma once

#include<Globals/Types.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Magnetization.h>
#include<Frequency_Multivariate_Gaussian/Frequency_Covariance_Matrix.h>
#include<Frequency_Multivariate_Gaussian/Frequency_Noise_Vectors.h>
#include<FFT/RealFFT.h>
#include"../matrices.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace spinDMFT::Functions
{

namespace stda = Standard_Algorithms;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace mag = Observables::Magnetization;
namespace fmvg = Frequency_Multivariate_Gaussian;
namespace ps = spinDMFT::Parameter_Space;
namespace rtd = spinDMFT::Run_Time_Data;

using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using MagVec = mag::MagnetizationVector;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;

void initialize_matrices( const ps::ParameterSpace& );

CorrTen generate_initial_spin_correlations( const ps::ParameterSpace& );

std::tuple<FieldVector, fmvg::FrequencyCovarianceMatrix> self_consistent_equations( const ps::ParameterSpace& pspace, const CorrTen& spin_correlations, const FieldVector& spin_expval );

std::tuple<RealType, std::vector<Operator>, std::vector<Operator>> compute_propagators( const ps::ParameterSpace& pspace, std::vector<FieldVector>& noise_vectors, FieldVector& mf_mean );

void compute_S_of_t( const ps::ParameterSpace& pspace, const std::vector<Operator>& propagators, const std::vector<Operator>& propagators_inv, std::vector<Operator>& S_x_of_t, std::vector<Operator>& S_y_of_t, std::vector<Operator>& S_z_of_t );

// Accumulates the per-sample correlation O_raw / Z(V) into spin_expval, spin_correlations,
// and the running sum of squares used to estimate the standard error of the mean. Only the
// first-moment components the symmetry type permits are stored in the MagVec's (and measured).
void compute_spin_correlations( rtd::RunTimeData& rtdata,
    CorrTen& spin_correlations_R, CorrTen& spin_correlations_I,
    MagVec& spin_expval, MagVec& spin_expval_sqsum,
    const RealType Z,
    const std::vector<Operator>& S_x_of_t,
    const std::vector<Operator>& S_y_of_t,
    const std::vector<Operator>& S_z_of_t );

void MPI_share_results( rtd::RunTimeData& rtdata, CorrTen& spin_correlations_R, CorrTen& spin_correlations_I, MagVec& spin_expval, MagVec& spin_expval_sqsum );

// divide the summed-over-samples accumulator by N to obtain the Monte-Carlo mean
void normalize( CorrTen& spin_correlations, RealType N );

// reconstruct the upper-half (tau > beta/2) correlations and their std errors from the
// computed lower half using g^{ab}(beta - tau) = g^{ba}(tau)^*
void symmetrize_upper_half( CorrTen& Re, CorrTen& Im, CorrTen& Re_std, CorrTen& Im_std, CorrTen& tau_Re, CorrTen& tau_Im );

// Build a time-domain mean-field trajectory from V given in the diagonal frequency basis:
//   V_full[block] = ortho[block] * V_diag[block], then inverse-DFT to imaginary time.
// The convention matches fmvg::FrequencyNoiseVectors::fourier_back_transform; the inverse
// DFT is performed component-wise in O(N log N) via the supplied FFTW plan (fft).
std::vector<FieldVector> diag_freq_to_time( const std::vector<fmvg::Vector>& V_diag,
                                            const fmvg::OrthogonalTransformationList& ortho,
                                            FFT::InverseRealDFT& fft );

// A single preconditioned-Crank-Nicolson chain for sampling V ~ pi propto p(V) Z(V).
// State (V in the diagonal Matsubara basis, Z, and the cached propagators and S^a(t))
// lives inside the object; observables are read off via the accessors and folded into
// the global accumulators by compute_spin_correlations.
class PCNChain
{
 public:
    // Cold-start constructor: draws an initial V from the prior and retries up to 32
    // times if Z happens to be non-positive (numerically degenerate).
    PCNChain( const ps::ParameterSpace& pspace,
              const fmvg::EigenValuesList& eig,
              const fmvg::OrthogonalTransformationList& ortho,
              const FieldVector& mf_mean,
              RealType pcn_beta,
              std::mt19937& engine );

    // Warm-start constructor: seeds the chain with the supplied V_diag (already in the
    // current iteration's diagonal basis). Modes frozen by the prior (eig <= 0) are
    // forced to zero so the seed stays consistent with the proposal distribution.
    PCNChain( const ps::ParameterSpace& pspace,
              const fmvg::EigenValuesList& eig,
              const fmvg::OrthogonalTransformationList& ortho,
              const FieldVector& mf_mean,
              RealType pcn_beta,
              std::mt19937& engine,
              std::vector<fmvg::Vector> V_diag_initial );

    // One pCN proposal + Metropolis accept/reject. Returns true on acceptance.
    bool step();

    // Accessors used by compute_spin_correlations.
    RealType Z() const { return m_Z; }
    const std::vector<Operator>& S_x() const { return m_S_x; }
    const std::vector<Operator>& S_y() const { return m_S_y; }
    const std::vector<Operator>& S_z() const { return m_S_z; }

    // Current state in the basis-independent Matsubara representation:
    //   V_full[b] = ortho[b] * V_diag[b].
    // Use this to carry state across self-consistent iterations (where the per-block
    // diagonalization changes but the Matsubara frequencies do not).
    std::vector<fmvg::Vector> V_full() const;

 private:
    void draw_prior_into( std::vector<fmvg::Vector>& V_diag );
    void rebuild_propagators_and_S(); // refreshes Z, propagators, S_*(t) from current V_diag

    // borrowed references (lifetime owned by main)
    const ps::ParameterSpace& m_pspace;
    const fmvg::OrthogonalTransformationList& m_ortho;
    FieldVector m_mf_mean;
    std::mt19937& m_engine;

    // sampling configuration
    RealType m_pcn_beta;
    RealType m_pcn_root; // sqrt(1 - pcn_beta^2)
    std::vector<std::array<std::normal_distribution<RealType>, 3>> m_prior_dists;
    std::uniform_real_distribution<RealType> m_uniform01;

    // reusable FFTW plan for the diagonal-frequency -> imaginary-time transform
    FFT::InverseRealDFT m_fft;

    // chain state
    std::vector<fmvg::Vector> m_V_diag;
    RealType m_Z;
    std::vector<Operator> m_propagators;
    std::vector<Operator> m_propagators_inv;
    std::vector<Operator> m_S_x, m_S_y, m_S_z;
};

};