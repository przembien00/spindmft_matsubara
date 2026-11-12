#pragma once
#include <functional>
#include <optional>
#include <string>
#include <vector>
#include <Physics/Physics.h>
#include "../Contour/Contour_Types.h"
namespace CspinDMFT::Mean_Field_Models { struct Model; }
namespace CspinDMFT::Parameter_Space {
class ParameterSpace {
public:
    ParameterSpace() = default;
    ParameterSpace(int argc,char* const argv[],int world_size,int rank);
    size_t get_num_Cores() const { return num_Cores; }
    void validate() const;
    std::string create_essentials_string(const Mean_Field_Models::Model&) const;
    std::string create_filename_string() const;
    std::string configuration_path() const;
    std::string output_directory() const;
    // Cluster physics follows CspinDMFT.
    std::string config_file, mf_model_name="CorrelationReplica";
    Physics::SpinModel spinspin_cmodel{"ISO"},spinmf_cmodel{"ISO"};
    Physics::Noise noise{"none",0.};
    std::string chemical_shift="none",extra_interaction="none";
    std::optional<double> uniform_field,staggered_field;
    double interaction_strength=0.,rescale_spinspin=1.,rescale_meanfield=1.,beta=1.;
    std::string bath="selfconsistent";
    double bath_frequency=1.,bath_coupling=1.;
    char bath_component='x';
    bool uses_harmonic_bath() const { return bath=="harmonic"; }
    bool self_consistency=true;
    size_t num_Spins=0;
    Contour::Layout preparation_layout{0,0,0},real_time_layout{0,0,0};
    // Contour and Monte Carlo settings follow spinDMFT_Keldysh_staged.
    char correlation_symmetry_type='D';
    bool time_translation_invariant=true,constant_magnetization_time=false;
    size_t num_TimeSteps=16,num_RealTimeSteps=10,real_time_substeps=1;
    double Tmax=1.,delta_real_t=.1;
    std::string real_time_evolution="staggered-strang";
    size_t num_SamplesPerCore=1000,num_PreparationSamplesPerCore=1000,num_blocks=32,num_SamplesPerSet=1;
    size_t Preparation_Iteration_Limit=20,Iteration_Limit=20;
    std::string seed="random",sampling_strategy="independent";
    unsigned long long generated_seed=0;
    double mh_step_size=.3;
    size_t mh_burn_in=100;
    std::string iteration_error_mode="relative",absolute_iteration_error_option="conservative";
    double absolute_iteration_error_tolerance=0.,iteration_error_sigma_threshold=5.;
    bool check_gaussian_reconstruction=false;
    double covariance_tolerance=1e-10,takagi_tolerance=1e-8,minimum_phase_magnitude=1e-6;
    std::string init_corr_mode="generate",initial_correlations_src_file,init_diag_corr="imagtime",init_nondiag_corr="zero";
    bool extrapolate_initial_spin_correlations=false;
    double correlation_periods=1.;
    std::string information_text,project_name,dst_project_name,filename_extension;
    size_t num_PrintDigits=4;
    int Run_ID=0;
private:
    size_t num_Cores=1;
};
}
namespace CspinDMFT { using ParameterSpace=Parameter_Space::ParameterSpace; }
