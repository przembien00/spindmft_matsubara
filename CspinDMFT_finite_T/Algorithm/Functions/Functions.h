#pragma once

#include<string>
#include<fstream>
#include<Globals/Types.h>
#include"../matrices.h"

#include<Multivariate_Gaussian/Multivariate_Gaussian_Blocks.h>
#include<Multivariate_Gaussian/Symmetry_Schemes.h>
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
using TimeTrajectories = std::vector<TimeTrajectory>;
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
    TimeTrajectories sample_time_trajectories( const mvgb::EigenValuesBlocks& eig, const mvgb::OrthogonalTransformationBlocks& ortho, std::mt19937& engine, const size_t num_samples ) const;

private:
    char m_symmetry_type{};
    size_t m_num_spins{};
    size_t m_num_frequencies{};
    std::vector<size_t> m_block_sizes{};
    std::vector<Multivariate_Gaussian::SymmetricMatrix> m_covariances{};

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

void compute_spin_observables( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<Operator>& forward_propagators, const std::vector<Operator>& backward_propagators, const RealType Z, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace );

std::tuple<std::vector<RealType>, std::vector<std::vector<Operator>>, std::vector<std::vector<Operator>>> compute_uncoupled_propagators( const TimeTrajectory& Vs_of_t, const SiteFields& mean_fields, const ps::ParameterSpace& pspace );

void accumulate_uncoupled_Z( std::vector<RealType>& uncoupled_partition_functions, rtd::RunTimeData& rtdata, const std::vector<RealType>& Z_i_list, const ps::ParameterSpace& pspace );

void compute_uncoupled_spin_observables( CluCorrTen& spin_CCT_Re, CluCorrTen& spin_CCT_Im, SiteFields& spin_expvals, const std::vector<std::vector<Operator>>& forward_propagators, const std::vector<std::vector<Operator>>& backward_propagators, const std::vector<RealType>& Z_i_list, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace );


void MPI_share_results( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, SiteFields& spin_expvals, rtd::RunTimeData& rtdata, RealType& partition_function );
void MPI_share_uncoupled_results( std::vector<RealType>& partition_functions, rtd::RunTimeData& rtdata );

void normalize( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, const RealType partition_function );
void normalize_uncoupled( CluCorrTen& spin_corr_Re, CluCorrTen& spin_corr_Im, SiteFields& spin_expvals, const std::vector<RealType>& partition_functions, const size_t num_Samples );
}
