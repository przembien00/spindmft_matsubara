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
  size_t num_RealTimePoints{};
  RealType delta_real_t{};
  RealType Tmax{};

  // ...concerning the statistics
  size_t num_SamplesPerCore{};
  size_t num_SamplesPerSet{};
  size_t num_SetsPerCore{};
  size_t num_Samples{};
  std::string seed{};

  // ...concerning Monte-Carlo sampling and correlated-chain statistics
  std::string sampling_strategy{"pcn"};
  bool antithetic_pairs{false};
  RealType mh_step_size{RealType{0.3}};
  size_t mh_burn_in{size_t{100}};
  RealType partition_imaginary_tolerance{RealType{1e-8}};
  size_t num_blocks{};
  std::string gaussian_factorization{"dense"};
  RealType fft_cross_frequency_cutoff{RealType{3.}};
  std::string spin_insertion_strategy{"closed-contour"};
  std::string correlation_normalization{"partition-function"};

  // ...concerning the iteration 
  RealType iteration_error_sigma_threshold{RealType{5.}};
  size_t Iteration_Limit{};
  RealType mixing_alpha{RealType{1.}};
  bool constant_magnetization_time{false};
  RealType covariance_tolerance{RealType{1e-10}};
  RealType branch_identity_tolerance{RealType{0.1}};
  RealType takagi_tolerance{RealType{1e-8}};
  RealType minimum_phase_magnitude{RealType{1e-6}};
  RealType denominator_constancy_tolerance{RealType{0.1}};
  RealType imaginary_magnetization_sigma{RealType{5.}};
  
  // ...concerning the initial correlations
  FieldVector initial_spin_expval{};
  ph::DiagonalSpinCorrelation init_diag_corr{};
  ph::NonDiagonalSpinCorrelation init_nondiag_corr{};
  bool load_initial_spin_correlations{ false };
  bool extrapolate_initial_spin_correlations{ false };
  size_t old_num_TimePoints{};
  RealType old_delta_t{};
  std::string initial_correlations_src_file{};
  std::string initial_correlations_src_directory{};
  std::vector<RealType> initial_correlations_linearized{};
  std::vector<RealType> initial_correlations_imag_linearized{};

  // ========== storing and naming ==========
  std::string information_text{};
  std::string project_name{};
  std::string filename_extension{};
  size_t num_PrintDigits{};
};

}
