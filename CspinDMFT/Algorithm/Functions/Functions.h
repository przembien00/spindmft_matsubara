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
using VecVec = std::vector<std::vector<RealType>>;
using Corr = corr::CorrelationVector;
using CorrTen = ten::CorrelationTensor<Corr>;
using CluCorrTen = clu::CorrelationCluster<CorrTen>;

namespace Initialization
{
    void write_general_spin_matrices( const std::vector<RealType>& spin_float_list, const size_t num_HilbertSpaceDimension );

    void write_cluster_Hamiltonian( const size_t num_Spins, const size_t num_HilbertSpaceDimension, const SymmMatrix& J, const ph::SpinModel& spinspin_cmodel, const ph::ChemicalShift& chemical_shift, const ph::LocalExtraInteraction& local_extra_interaction );

    CluCorrTen generate_initial_environment_spin_correlations( const ps::ParameterSpace& pspace );
}

void propagate( Operator& TimeEvolutionOperator, const SparseObservable& old_Hamiltonian, const SparseObservable& new_Hamiltonian, const ps::ParameterSpace& pspace );

void compute_spin_correlations( CluCorrTen& spin_CCT, const Operator& TimeEvolutionOperator, const size_t time, rtd::RunTimeData& rtdata, const ps::ParameterSpace& pspace );

void MPI_share_results( CluCorrTen& spin_corr, rtd::RunTimeData& rtdata );

void normalize_and_fill_time_zero ( CluCorrTen& spin_corr, const size_t num_Samples, const std::vector<RealType>& spin_float_list );

}
