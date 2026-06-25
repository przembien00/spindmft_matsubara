#pragma once

#include<string>
#include<fstream>
#include<memory>
#include<random>
#include<array>
#include<Globals/Types.h>
#include"../matrices.h"

#include<Multivariate_Gaussian/Multivariate_Gaussian_Blocks.h>
#include<Multivariate_Gaussian/Symmetry_Schemes.h>
#include<FFT/RealFFT.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Clusters.h>
#include<Physics/Physics.h>

#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace Functions
{

namespace mvgb = Multivariate_Gaussian::Blocks;
namespace mss = mvgb::Symmetry_Schemes;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace clu = Observables::Clusters;
namespace ph = Physics;
namespace ps = DMFT_parameter_space;
namespace rtd = Run_Time_Data;

using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;
using SiteFields = std::vector<FieldVector>;
using TimeTrajectory = std::vector<SiteFields>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using CluCorrTen = clu::CorrelationCluster<CorrTen>;

class FrequencyCovarianceCluster
{
public:
    FrequencyCovarianceCluster() = default;
    FrequencyCovarianceCluster( const CluCorrTen& correlations, const char symmetry_type, const size_t num_spins );

    void fill( const CluCorrTen& correlations, const char symmetry_type, const size_t num_spins );
    void diagonalize( mvgb::EigenValuesBlocks& eig, mvgb::OrthogonalTransformationBlocks& ortho ) const;

    // Build a single imaginary-time mean-field trajectory from one frequency-space sample
    // given in the full (non-diagonal) block basis: full[flat_index(frequency,block)][row]
    // holds the same coefficients that sample_time_trajectories reads from the rotated
    // Gaussian noise. Used by the pCN sampler, which keeps its state in frequency space.
    TimeTrajectory trajectory_from_full_frequency( const mvgb::EigenValuesBlocks& full ) const;

private:
    char m_symmetry_type{};
    size_t m_num_spins{};
    size_t m_num_frequencies{};
    std::vector<size_t> m_block_sizes{};
    std::vector<Multivariate_Gaussian::SymmetricMatrix> m_covariances{};
    // reusable FFTW plan for the per-site, per-component freq -> imaginary-time transform;
    // built in fill() once m_num_frequencies is known (shared_ptr so the class stays copyable
    // and the plan is reused across all pCN steps of one self-consistent iteration).
    std::shared_ptr<FFT::InverseRealDFT> m_fft{};

    size_t flat_index( const size_t frequency, const size_t block ) const;
};

namespace Initialization
{
    void write_general_spin_matrices( const std::vector<RealType>& spin_float_list, const size_t num_HilbertSpaceDimension );
    void write_uncoupled_spin_matrices();

    void write_cluster_Hamiltonian( const size_t num_Spins, const size_t num_HilbertSpaceDimension, const SymmMatrix& J, const ph::SpinModel& spinspin_cmodel, const ph::ChemicalShift& chemical_shift, const ph::LocalExtraInteraction& local_extra_interaction );

	    CluCorrTen generate_initial_environment_spin_correlations( const ps::ParameterSpace& pspace, const std::string& component = "Re" );
}

Operator compute_shortstep_propagator( const SparseObservable& old_Hamiltonian, const SparseObservable& new_Hamiltonian, const ps::ParameterSpace& pspace );

std::tuple<RealType, std::vector<Operator>, std::vector<Operator>> compute_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace );

std::tuple<std::vector<RealType>, std::vector<std::vector<Operator>>, std::vector<std::vector<Operator>>> compute_uncoupled_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace );

// ======================== preconditioned Crank-Nicolson sampler ========================
// Instead of importance sampling V ~ p(V) and reweighting by Z(V), the pCN chain targets
// pi(V) propto p(V) Z(V) directly, where Z(V) is the cluster partition function for the
// sampled mean-field trajectory. The chain state lives in the diagonal frequency basis,
// where p(V) factorizes into independent N(0, eig) Gaussians per block component, so the
// pCN proposal V' = sqrt(1-beta^2) V + beta xi (xi ~ p(V)) preserves p(V) exactly and the
// Metropolis ratio collapses to min(1, Z(V')/Z(V)).
//
// In the uncoupled-spins mode Z(V) factorizes over sites, Z = prod_i Z_i, and per-sample
// observables are divided site-by-site (Z_i for on-site, Z_i Z_j for inter-site).

// pCN versions of the observable accumulators: they fold the per-sample value O_raw / Z
// into the running correlation sums and the sum of squares (used for the standard error).
// No covariance/partition-function accumulation is needed because V ~ pi already.
void compute_spin_observables_mh( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<Operator>& forward_propagators, const std::vector<Operator>& backward_propagators, const RealType Z, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace );
void compute_uncoupled_spin_observables_mh( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<std::vector<Operator>>& forward_propagators, const std::vector<std::vector<Operator>>& backward_propagators, const std::vector<RealType>& Z_i_list, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace );

// sum the per-core pCN accumulators (means and sums of squares) across all MPI cores
void MPI_share_results_mh( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, SiteFields& spin_expvals, rtd::RunTimeData& rtdata );

// divide the summed-over-samples accumulators by the global sample count N
void normalize_mh( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, const RealType N );

// Reconstruct the upper-half (tau > beta/2) correlations from the computed lower half via the
// per-component reflection symmetry g^{ab}_{ij}(beta-tau) = g^{ab}_{ij}(tau)^* (Re even, Im
// odd about beta/2). This follows from trace cyclicity g^{ab}_{ij}(beta-tau) = g^{ba}_{ji}(tau)
// together with Hermiticity g^{ab}_{ij}(tau)^* = g^{ba}_{ji}(tau), so it is a SAME-component
// relation (no site/direction transpose). The observable accumulators only fill [0, beta/2];
// this fills the rest (and mirrors the per-point standard errors, with no sign change).
void symmetrize_upper_half( CluCorrTen& Re, CluCorrTen& Im, CluCorrTen& Re_std, CluCorrTen& Im_std, CluCorrTen& Re_tau, CluCorrTen& Im_tau );

// A single preconditioned-Crank-Nicolson chain for sampling V ~ pi propto p(V) Z(V).
// The chain owns its frequency-space state and the cached imaginary-time propagators, and
// exposes them through accessors so the observable accumulators can read off the current
// state every production step. Both the coupled and uncoupled cluster paths are supported.
class PCNChainCluster
{
 public:
    // Cold-start constructor: draws V from the prior and retries on a non-positive Z.
    PCNChainCluster( const ps::ParameterSpace& pspace,
                     const FrequencyCovarianceCluster& cov,
                     const mvgb::EigenValuesBlocks& eig,
                     const mvgb::OrthogonalTransformationBlocks& ortho,
                     const SiteFields& mean_fields,
                     RealType pcn_beta,
                     std::mt19937& engine );

    // Warm-start constructor: seeds the chain from a frequency-space state carried over
    // from the previous self-consistent iteration (full_freq[flat] in the full block basis).
    PCNChainCluster( const ps::ParameterSpace& pspace,
                     const FrequencyCovarianceCluster& cov,
                     const mvgb::EigenValuesBlocks& eig,
                     const mvgb::OrthogonalTransformationBlocks& ortho,
                     const SiteFields& mean_fields,
                     RealType pcn_beta,
                     std::mt19937& engine,
                     const mvgb::EigenValuesBlocks& full_freq_initial );

    // One pCN proposal + Metropolis accept/reject. Returns true on acceptance.
    bool step();

    // Accessors for the coupled path.
    RealType Z() const { return m_Z; }
    const std::vector<Operator>& forward() const { return m_forward; }
    const std::vector<Operator>& backward() const { return m_backward; }

    // Accessors for the uncoupled path.
    const std::vector<RealType>& Z_i_list() const { return m_Z_i_list; }
    const std::vector<std::vector<Operator>>& forward_uncoupled() const { return m_forward_unc; }
    const std::vector<std::vector<Operator>>& backward_uncoupled() const { return m_backward_unc; }

    // Current state in the basis-independent full frequency representation,
    //   full_freq[flat] = ortho[flat] * V_diag[flat],
    // used to carry state across self-consistent iterations.
    mvgb::EigenValuesBlocks full_frequency() const;

 private:
    void build_prior();
    void draw_prior_into( mvgb::EigenValuesBlocks& V_diag );
    mvgb::EigenValuesBlocks rotate_to_full( const mvgb::EigenValuesBlocks& V_diag ) const;
    // evaluate Z (and log Z), the forward propagators and the short-step propagators for a
    // given diagonal-basis state. Backward propagators are NOT built here: the Metropolis
    // test needs only Z, so they are deferred to refresh_backward (called on acceptance).
    void evaluate( const mvgb::EigenValuesBlocks& V_diag, RealType& Z, RealType& logZ,
                   std::vector<Operator>& fwd, std::vector<Operator>& shortstep,
                   std::vector<RealType>& Z_i, std::vector<std::vector<Operator>>& fwd_unc,
                   std::vector<std::vector<Operator>>& shortstep_unc ) const;
    // build the cached backward propagators from the cached (accepted) short-step propagators
    void refresh_backward();

    // borrowed references (lifetime owned by main)
    const ps::ParameterSpace& m_pspace;
    const FrequencyCovarianceCluster& m_cov;
    const mvgb::EigenValuesBlocks& m_eig;
    const mvgb::OrthogonalTransformationBlocks& m_ortho;
    SiteFields m_mean_fields;
    std::mt19937& m_engine;

    // sampling configuration
    bool m_uncoupled{};
    RealType m_pcn_beta{};
    RealType m_pcn_root{}; // sqrt(1 - pcn_beta^2)
    std::vector<std::vector<std::normal_distribution<RealType>>> m_prior_dists; // per flat block, per component
    std::uniform_real_distribution<RealType> m_uniform01{ RealType{0.}, RealType{1.} };

    // chain state (frequency space + cached imaginary-time quantities)
    mvgb::EigenValuesBlocks m_V_diag;
    RealType m_Z{};
    RealType m_logZ{};
    std::vector<Operator> m_forward, m_backward, m_shortstep;                              // coupled path
    std::vector<RealType> m_Z_i_list;                                                      // uncoupled path
    std::vector<std::vector<Operator>> m_forward_unc, m_backward_unc, m_shortstep_unc;     // uncoupled path
};
}
