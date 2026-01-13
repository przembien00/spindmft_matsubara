#pragma once

#include<string>
#include<vector>
#include<iostream>
#include<functional>
#include<memory>
#include<Globals/Types.h>
#include<Physics/Physics.h>
#include<Observables/Clusters.h>
#include"../matrices.h"

// Forward Declaration of Model
namespace Mean_Field_Models
{
  class MeanFieldModel;
}

namespace DMFT_parameter_space
{

namespace ph = Physics;
namespace mfm = Mean_Field_Models;

using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;
using CluList = Observables::Clusters::CorrelationCluster<std::vector<RealType>>;

// ============================================================================
// ========================== PARAMETER SPACE CLASS ===========================
// ============================================================================
class ParameterSpace
{ 
public:
  // CONSTRUCTORS
  ParameterSpace() = default;
  ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank );

  // PUBLIC METHODS
  size_t get_num_Cores() const{ return static_cast<size_t>(world_size); }
  void read_initial_correlations_from_file();
  std::string create_essentials_string() const;

private:
  // PRIVATE MEMBERS
  int my_rank{};
  int world_size{};

public:
  // PUBLIC MEMBERS
  int Run_ID{};

  // ========== physical parameters ==========
  // ...concerning the mean-field model
  std::string mf_model_name{};
  std::shared_ptr< mfm::MeanFieldModel > mf_model{nullptr};

  // ...concerning the cluster
  std::string config_file{};
  size_t num_Spins{};
  size_t num_HilbertSpaceDimension{};
  std::vector<std::string> spin_list{}; // list of spin values as string
  std::vector<RealType> spin_float_list{}; // list of spin values as float
  SymmMatrix J{};                   // local spin-spin couplings
  ph::SpinModel spinspin_cmodel{};  // spin-spin coupling model
  ph::SpinModel spinmf_cmodel{};    // spin-mean-field coupling model
  ph::ChemicalShift chemical_shift{}; // chemical shift (always in z-direction)
  ph::LocalExtraInteraction local_extra_interaction{}; // extra interaction, e.g., quadrupolar 
  RealType rescale_meanfield{1.0};  // rescales spin-mean-field couplings by a factor
  RealType rescale_spinspin{1.0};   // rescales spin-spin couplings by a factor

  // ...concerning the simulation output
  IndexPairList correlation_categories{}; // correlation categories that will be computed during the simulation

  
  // ========== numerical parameters ==========
  // ...concerning the symmetry
  char symmetry_type{};

  // ...concerning the time discretization
  size_t num_TimeSteps{};
  size_t num_TimePoints{};
  RealType delta_t{};

  // ...concerning the statistics
  size_t num_SamplesPerCore{};
  size_t num_SamplesPerSet{};
  size_t num_Samples{};
  std::string seed{};
  bool adaptive_sample_size{ false }; 
  RealType statistical_error_tolerance{}; 

  // ...concerning the initial correlations 
  std::string init_corr_mode{};
  // ...if they are imported
  std::string imported_correlations_src_file{};
  std::string import_filename{};
  bool extrapolate_imported_spin_correlations{ false };
  size_t old_num_TimePoints{};
  size_t old_num_Spins{1};
  RealType old_delta_t{};
  CluList imported_correlations_linearized{}; // hierarchy: category -> direction & time (not linearized at the cluster level)
  // ...if they are generated
  ph::DiagonalSpinCorrelation init_diag_corr{};
  ph::NonDiagonalSpinCorrelation init_nondiag_corr{};

  // ...concerning the iteration
  std::string iteration_error_mode{};
  RealType absolute_iteration_error_tolerance{};
  RealType relative_iteration_error_tolerance{};
  size_t Iteration_Limit{};

    // ...concerning the eigenvalues
  RealType eigenvalue_ratio_tolerance{};
  std::string truncation_scheme_negative_eigenvalues{};

  // ...and more
  std::string matrix_exponential_computation{};


  // ========== parameters for storing and naming ==========
  std::string information_text{};
  std::string project_name{};
  std::string dst_project_name{};
  std::string filename_extension{};
  size_t num_PrintDigits{5};
};

}
