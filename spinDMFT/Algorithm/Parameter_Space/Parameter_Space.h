#pragma once

#include<string>
#include<vector>
#include<iostream>
#include<functional>
#include<Globals/Types.h>
#include<Physics/Physics.h>
#include"../matrices.h"

namespace spinDMFT::Parameter_Space
{

// ======================= NAMESPACES ============================
namespace ph = Physics;


// ===============================================================
// ============== HEADER FOR PARAMETER SPACE CLASS ===============
// ===============================================================
class ParameterSpace
{
 private:
  // PRIVATE MEMBERS
  int my_rank{};
  int world_size{};

 public:
  // CONSTRUCTORS
  ParameterSpace() = default;
  ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank );

  // PUBLIC METHODS
  size_t get_num_Cores() const{ return static_cast<size_t>(world_size); }
  void read_initial_correlations_from_file();
  std::string create_essentials_string() const;

  // PUBLIC MEMBERS
  // ========== model and physical parameters ==========
  bool self_consistency{ true };
  ph::SpinModel spin_model{};
  std::string spin{};
  RealType spin_float{};
  size_t num_HilbertSpaceDimension{};
  RealType JQ{};
  RealType JL{};
  RealType beta{};
  ph::MagneticField B{};
  ph::Noise noise{};
  ph::ExtraInteraction extra_interaction{};

  // ========== general numerical parameters ==========
  // ...concerning the symmetry
  char correlation_symmetry_type{};

  // ...concerning the time discretization
  size_t num_TimeSteps{};
  size_t num_TimePoints{};
  RealType delta_t{};

  // ...concerning the statistics
  size_t num_SamplesPerCore{};
  size_t num_SamplesPerSet{};
  size_t num_SetsPerCore{};
  size_t num_Samples{};
  std::string seed{};

  // ...concerning the iteration 
  RealType absolute_iteration_error_threshold{};
  size_t Iteration_Limit{};
  
  // ...concerning the initial correlations 
  ph::DiagonalSpinCorrelation init_diag_corr{};
  ph::NonDiagonalSpinCorrelation init_nondiag_corr{};
  bool load_initial_spin_correlations{ false };
  bool extrapolate_initial_spin_correlations{ false };
  size_t old_num_TimePoints{};
  RealType old_delta_t{};
  std::string initial_correlations_src_file{};
  std::string initial_correlations_src_directory{};
  std::vector<RealType> initial_correlations_linearized{};

  // ...concerning potential negative eigenvalues of the covariance matrix
  RealType critical_eigenvalue_ratio{};
  std::string truncation_scheme_negative_eigenvalues{};
  
  // ========== storing and naming ==========
  std::string information_text{};
  std::string project_name{};
  std::string filename_extension{};
  size_t num_PrintDigits{};
};

}
