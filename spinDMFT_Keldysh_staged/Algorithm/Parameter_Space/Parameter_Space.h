#pragma once

#include<string>
#include<array>
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
  std::string bath{"selfconsistent"};
  RealType bath_frequency{RealType{1.}};
  RealType bath_coupling{RealType{1.}};
  char bath_component{'x'};
  bool uses_harmonic_bath() const { return bath=="harmonic"; }

  // ========== general numerical parameters ==========
  // ...concerning the symmetry
  char correlation_symmetry_type{};

  // ...concerning the time discretization
  size_t num_TimeSteps{};
  size_t num_TimePoints{};
  RealType delta_t{};
  size_t num_RealTimeSteps{};
  size_t real_time_substeps{1};
  std::string real_time_evolution{"staggered-strang"};
  bool uses_staggered_strang() const
  { return num_RealTimeSteps>0; }
  size_t num_RealTimePoints{};
  RealType delta_real_t{};
  RealType Tmax{};

  // ...concerning the statistics
  size_t num_SamplesPerCore{};
  size_t num_Samples{};
  std::string seed{};

  // ...concerning Monte-Carlo sampling and correlated-chain statistics
  std::string sampling_strategy{"independent"};
  RealType mh_step_size{RealType{0.3}};
  size_t mh_burn_in{size_t{100}};
  size_t num_blocks{};
  size_t real_time_steps_per_interval() const
  { return real_time_substeps; }

  // ...concerning the iteration 
  RealType iteration_error_sigma_threshold{RealType{5.}};
  size_t Iteration_Limit{};
  size_t num_PreparationSamplesPerCore{}, num_PreparationSamples{};
  size_t Preparation_Iteration_Limit{};
  bool constant_magnetization_time{false};
  RealType covariance_tolerance{RealType{1e-10}};
  RealType takagi_tolerance{RealType{1e-8}};
  RealType minimum_phase_magnitude{RealType{1e-6}};
  
  // ...concerning the initial correlations
  ph::DiagonalSpinCorrelation init_diag_corr{};
  ph::NonDiagonalSpinCorrelation init_nondiag_corr{};
  bool load_initial_spin_correlations{ false };
  bool extrapolate_initial_spin_correlations{ false };
  size_t old_num_TimePoints{};
  RealType old_delta_t{};
  std::string initial_correlations_src_file{};
  std::vector<RealType> initial_correlations_linearized{};
  std::vector<RealType> initial_correlations_imag_linearized{};
  std::vector<RealType> initial_magnetization_linearized{};

  // ========== storing and naming ==========
  std::string information_text{};
  std::string project_name{};
  std::string filename_extension{};
  size_t num_PrintDigits{};
};

}
