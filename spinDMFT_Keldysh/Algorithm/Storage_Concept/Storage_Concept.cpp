#include"Storage_Concept.h"
#include<fstream>
#include<sys/stat.h>

#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;

#include<File_Management/File_Management.h>
namespace fm = File_Management;

#include"STOC_Error_Handling.h"
namespace error = spinDMFT::Storage_Concept::Error_Handling;

namespace hdf5r = HDF5_Routines;

namespace spinDMFT::Storage_Concept
{
// ===================================================
// =============== HDF5 STORAGE CLASS ================
// ===================================================
// constructor : create folder branch and data file
HDF5_Storage::HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace, const std::string& termination ):
    m_storing_permission( my_rank == 0 ),
    m_fname_max_length( 200 ),
    m_num_TriesToBuildFile( 5 )
{
    create_folder_branch( pspace, termination );
    create_file( pspace );
}

// create the folder branch in which the data will be stored
void HDF5_Storage::create_folder_branch( const ps::ParameterSpace& pspace, const std::string& termination  )
{
    if( !m_storing_permission ){ return; } // permission request

    // determine folder branch list
    std::vector<std::string> folder_branch_list{};
    folder_branch_list.push_back( "Data" );
    if( pspace.uses_harmonic_bath() )
    {
        folder_branch_list.push_back( "HARMONIC_BATH" );
    }
    else if( !pspace.self_consistency )
    {
        folder_branch_list.push_back( "noselfcons" );
    }
    if( pspace.project_name != "" )
    {
        folder_branch_list.push_back( pspace.project_name );
    }
    if( termination == "by iteration limit" ) // additional folder, if data are not converged regularly
    {
        folder_branch_list.push_back( "NOT_CONVERGED" );
    }

    // create folder branch:
    m_filename = fm::create_folder_tree( folder_branch_list, m_fname_max_length );
}

// create the file in which the data will be stored
void HDF5_Storage::create_file( const ps::ParameterSpace& pspace )
{
    if( !m_storing_permission ){ return; } // permission request

    // create filename:
    std::string filename = pspace.spin_model.compact_info( pspace.num_PrintDigits ); // spin model info
    if( pspace.spin != "1/2" ) // add spin value to filename if unequal to 1/2
    {
        std::string spin_without_slash{pspace.spin};
        auto pos = spin_without_slash.find("/");
        if( pos != std::string::npos ) // interpret 1/2, 3/2, ... because no / in allowed in filename
        {
            spin_without_slash.replace(pos, 1, "o");
        }
        filename += "__spin=" + spin_without_slash;
    }
    if( pspace.JQ != RealType{1.} )
    {
        filename += "__JQ=" + print::remove_zeros(print::round_value_to_string(pspace.JQ,pspace.num_PrintDigits)); // rescaling
    }
    if( pspace.JL != RealType{0.} )
    {
        filename += "__JL=" + print::remove_zeros(print::round_value_to_string(pspace.JL,pspace.num_PrintDigits)); // rescaling
    }
    if( pspace.beta != RealType{0.} )
    {
         filename += "__beta=" + print::remove_zeros(print::round_value_to_string(pspace.beta,pspace.num_PrintDigits));
    }
    if( pspace.uses_harmonic_bath() )
    {
        filename += "__bath=harmonic_omega="
            +print::remove_zeros(print::round_value_to_string(
                pspace.bath_frequency,pspace.num_PrintDigits))
            +"_g="+print::remove_zeros(print::round_value_to_string(
                pspace.bath_coupling,pspace.num_PrintDigits))
            +"_axis="+std::string(1,pspace.bath_component);
    }
    if( pspace.B.m_name != "none" )
    {
        filename += "__" + pspace.B.compact_info( pspace.num_PrintDigits ); // external field info
    }
    if( pspace.noise.m_name != "none" )
    {
        filename += "__" + pspace.noise.compact_info( pspace.num_PrintDigits ); // noise info
    }
    if( pspace.extra_interaction.m_name != "none" )
    {
        filename += "__" + pspace.extra_interaction.compact_info( pspace.num_PrintDigits ); // extra_interaction info
    }
    if( pspace.constant_magnetization_time )
    {
        filename += "__mag=constant";
    }
    if( pspace.correlation_normalization=="closed-contour" )
    {
        filename += "__corrnorm=D";
    }
    if( pspace.filename_extension != "" )
    {
        filename += "_" + pspace.filename_extension;
    }

    // create the file:
    print::cut_if_too_large( filename, m_fname_max_length );
    m_filename += filename;
    size_t count = 0;
    do{
        std::string tmp = m_filename + ".hdf5";
        std::ifstream f( tmp.c_str() );
        if( !f.good() ) // then the file doesn't exist yet
        {
            m_file_id = H5Fcreate( tmp.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            if( count > 0 ) // inform about the difficulties -> creating file didn't work in the first try
            {
                std::cout << "\033[1;31mdata file already exists!\nsuccessfully created " << m_filename << " instead\033[0m\n";
            }
            break;
        }
        else // file exists already
        {
            m_filename += "X"; // adds an X to the filename and retries
        }
        f.close();
        count++;
    }while( count < m_num_TriesToBuildFile );
    m_filename += ".hdf5";
    if( count == m_num_TriesToBuildFile )
    {
        error::CREATE_FILE( m_filename, __PRETTY_FUNCTION__ );
    }
}

// store main data
void HDF5_Storage::store_main( const ps::ParameterSpace& pspace,
                              const rtd::RunTimeData& rtdata,
                              const CorrelationSet& correlations,
                              const CorrelationSet& standard_errors )
{
    if( !m_storing_permission ){ return; } // permission request

    // ++++++++++ store parameters ++++++++++
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // ========== model and physical parameters ==========
    hdf5r::store_string( ps_group_id, "self_consistency",           hdf5r::bool_to_string(pspace.self_consistency) );
    hdf5r::store_string( ps_group_id, "spin_model",                 pspace.spin_model.m_name );
    for( const auto& param : pspace.spin_model.m_parameters ) // add all parameters of the spin model
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "spin",                       pspace.spin );
    hdf5r::store_scalar( ps_group_id, "num_HilbertSpaceDimension",  pspace.num_HilbertSpaceDimension );
    hdf5r::store_scalar( ps_group_id, "JQ",                         pspace.JQ );
    hdf5r::store_scalar( ps_group_id, "JL",                         pspace.JL );
    hdf5r::store_scalar( ps_group_id, "beta",                       pspace.beta );
    hdf5r::store_string( ps_group_id, "bath",                       pspace.bath );
    hdf5r::store_scalar( ps_group_id, "bath_omega",                 pspace.bath_frequency );
    hdf5r::store_scalar( ps_group_id, "bath_coupling",              pspace.bath_coupling );
    hdf5r::store_string( ps_group_id, "bath_component",             std::string(1,pspace.bath_component) );
    hdf5r::store_string( ps_group_id, "field_source",
        pspace.uses_harmonic_bath()
        ?"prescribed thermal harmonic oscillator; no spinDMFT feedback"
        :"spinDMFT first- and second-moment self-consistency" );
    hdf5r::store_string( ps_group_id, "bath_correlation_definition",
        pspace.uses_harmonic_bath()
        ?"X(t,tau)=g^2(n_B+1)[exp(-omega*tau+i*omega*t)+exp(-omega*(beta-tau)-i*omega*t)]; tau=0 is lesser and tau=beta is greater"
        :"not applicable" );
    hdf5r::store_string( ps_group_id, "B",                          pspace.B.m_name );
    for( const auto& param : pspace.B.m_parameters ) // add all parameters of the magnetic field
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "noise",                      pspace.noise.m_name );
    for( const auto& param : pspace.noise.m_parameters ) // add all parameters of the noise
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "extra_interaction",          pspace.extra_interaction.m_name );
    for( const auto& param : pspace.extra_interaction.m_parameters ) // add all parameters of the extra interaction
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }

    // ========== general numerical parameters ==========
    // ...concerning symmetry
    hdf5r::store_scalar( ps_group_id, "correlation_symmetry_type",  pspace.correlation_symmetry_type );

    // ...concerning time discretization
    hdf5r::store_scalar( ps_group_id, "num_ImagTimeSteps",          pspace.num_TimeSteps );
    hdf5r::store_scalar( ps_group_id, "num_ImagTimePoints",         pspace.num_TimePoints );
    hdf5r::store_scalar( ps_group_id, "delta_imag_t",               pspace.delta_t );
    hdf5r::store_scalar( ps_group_id, "num_RealTimeSteps",          pspace.num_RealTimeSteps );
    hdf5r::store_scalar( ps_group_id, "num_RealTimePoints",         pspace.num_RealTimePoints );
    hdf5r::store_scalar( ps_group_id, "delta_real_t",               pspace.delta_real_t );
    hdf5r::store_scalar( ps_group_id, "Tmax",                       pspace.Tmax );
    hdf5r::store_string( ps_group_id, "equal_time_prescription",
                         "symmetric_theta_half" );
    hdf5r::store_string( ps_group_id, "equal_time_prescription_definition",
        "At coincident same-branch contour points theta(0)=1/2 is used: the canonical flat-index entry is (G greater + G lesser)/2 and its transpose reuses the identical value so E[V V^T] is complex symmetric." );

    // ...concerning statistics
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerCore",         pspace.num_SamplesPerCore );
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerSet",          pspace.num_SamplesPerSet );
    hdf5r::store_scalar( ps_group_id, "num_Samples",                pspace.num_Samples );
    hdf5r::store_scalar( ps_group_id, "num_Cores",                  pspace.get_num_Cores() );

    hdf5r::store_string( ps_group_id, "seed",                       pspace.seed );

    // ...concerning complex-field sampling
    hdf5r::store_string( ps_group_id, "sampling_strategy",
                         pspace.sampling_strategy );
    hdf5r::store_string( ps_group_id, "pcn_sampling_weight",
        pspace.sampling_strategy!="pcn"
        ?"not applicable"
        :pspace.correlation_normalization=="closed-contour"
         ?"Re D(T), where D(T)=Tr[U_+(T,0) rho_M B_-(T,0)]"
         :"Re Z_M" );
    hdf5r::store_string( ps_group_id, "antithetic_pairs",
                         hdf5r::bool_to_string(pspace.antithetic_pairs) );
    hdf5r::store_string( ps_group_id, "antithetic_pair_definition",
        !pspace.antithetic_pairs
        ?"disabled"
        :pspace.sampling_strategy=="pcn"
         ?pspace.correlation_normalization=="closed-contour"
          ?"Sign-symmetrized pCN: each Markov state evaluates r and -r, targets p0(r) Re[D_r(T)+D_-r(T)], and reweights the summed observables by that real likelihood. num_Samples counts pair states."
          :"Sign-symmetrized pCN: each Markov state evaluates r and -r, targets p0(r) Re[Z_M(r)+Z_M(-r)], and reweights the summed observables by that real likelihood. num_Samples counts pair states."
         :"Each independent real Gaussian latent draw r is evaluated as two trajectories with fluctuations Lr and -Lr. num_Samples counts both trajectories; numerator and Z_M remain accumulated into the global complex ratio." );
    if( pspace.antithetic_pairs )
        hdf5r::store_scalar( ps_group_id, "num_antithetic_pairs",
            pspace.sampling_strategy=="pcn"
            ?pspace.num_Samples:pspace.num_Samples/2 );
    hdf5r::store_string( ps_group_id, "sampling_count_unit",
        pspace.sampling_strategy=="pcn"&&pspace.antithetic_pairs
        ?"sign-symmetrized pCN pair states"
        :pspace.sampling_strategy=="pcn"
         ?"pCN Markov states"
         :"contour trajectories" );
    hdf5r::store_scalar( ps_group_id, "mh_step_size",
                         pspace.mh_step_size );
    hdf5r::store_scalar( ps_group_id, "mh_burn_in",
                         pspace.mh_burn_in );
    hdf5r::store_scalar( ps_group_id, "partition_imaginary_tolerance",
                         pspace.partition_imaginary_tolerance );
    hdf5r::store_string( ps_group_id, "gaussian_factorization",
                         pspace.gaussian_factorization );
    hdf5r::store_scalar( ps_group_id, "fft_cross_frequency_cutoff",
                         pspace.fft_cross_frequency_cutoff );
    hdf5r::store_string( ps_group_id, "gaussian_factorization_options",
        "dense: physical-grid real-lift Autonne--Takagi; svd: physical-grid canonical complex-SVD Takagi; fft: doubled-real FFT with one Matsubara plus low-real-frequency block, discarded high-frequency Matsubara-real covariance, independently sampled high {omega,-omega} Takagi blocks without a global dense factor, inverse FFT, and physical-grid restriction" );
    hdf5r::store_string( ps_group_id, "correlation_normalization",
                         pspace.correlation_normalization );
    hdf5r::store_string( ps_group_id, "magnetization_normalization",
                         pspace.correlation_normalization );
    hdf5r::store_string( ps_group_id, "normalization_strategy",
        pspace.correlation_normalization=="closed-contour"
        ?pspace.sampling_strategy=="pcn"
         ?pspace.antithetic_pairs
          ?"pCN target pi_pair(r) proportional to p0(r) Re[D_r(T)+D_-r(T)]; fixed-final-contour ratio sum[(A(r)+A(-r))/Re(D_r(T)+D_-r(T))]/sum[(D_r(T)+D_-r(T))/Re(D_r(T)+D_-r(T))] for correlations and magnetization at every time"
          :"pCN target pi(r) proportional to p0(r) Re D_r(T); fixed-final-contour ratio sum[A/Re D(T)]/sum[D(T)/Re D(T)] for correlations and magnetization at every time"
         :"bare-prior fixed-final-contour ratio: (sum numerator)/(sum D(T)) for correlations and magnetization at every time"
        :pspace.sampling_strategy=="pcn"
        ?pspace.antithetic_pairs
         ?"sign-symmetrized pCN target pi_pair(r) proportional to p0(r) Re[Z_M(r)+Z_M(-r)]; arithmetic mean of [N(r)+N(-r)]/Re[Z_M(r)+Z_M(-r)]"
         :"pCN target pi(r) proportional to p0(r) Re Z_M(r), Re Z_M>0; arithmetic mean of numerator/Re Z_M"
        :"bare-prior exact complex ratio: (sum numerator)/(sum Z_M)" );
    hdf5r::store_string( ps_group_id, "spin_insertion_strategy",
        pspace.spin_insertion_strategy=="prefix"
        ?"prefix: U_+(0,t) S B_-(t,0) for correlations and magnetization"
        :"closed contour: B_-(T,0) U_+(t,T) S U_+(t,0) for correlations and magnetization; U_+(t,T) is the forward continuation t -> T, not an inverse" );
    hdf5r::store_string( ps_group_id, "closed_contour_diagnostic_definition",
        pspace.sampling_strategy=="pcn"
        ?pspace.antithetic_pairs
         ?pspace.correlation_normalization=="closed-contour"
          ?"prefix diagnostic D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored value is the sign-symmetrized pCN mean of [D_r(t)+D_-r(t)]/Re[D_r(T)+D_-r(T)]"
          :"prefix diagnostic D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored value is the sign-symmetrized pCN mean of [D_r(t)+D_-r(t)]/Re[Z_M(r)+Z_M(-r)]"
         :pspace.correlation_normalization=="closed-contour"
          ?"prefix diagnostic D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored value is the pCN mean of D(t)/Re D(T)"
          :"prefix diagnostic D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored value is the pCN mean of D(t)/Re Z_M"
        :pspace.correlation_normalization=="closed-contour"
         ?"prefix D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored diagnostic ratio is (sum D(t))/(sum Z_M), and only its final value sum D(T) is the selected correlation and magnetization denominator at every time"
         :"prefix diagnostic D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]; stored ratio is (sum D(t))/(sum Z_M), not the closed-contour observable denominator" );
    hdf5r::store_string( ps_group_id, "connected_correlation_strategy",
        pspace.uses_harmonic_bath()
        ?"not used for field construction in prescribed harmonic-bath mode"
        :pspace.constant_magnetization_time
        ?"X_conn^{ab}(t,tau)=X^{ab}(t,tau)-m_a(0)m_b(0), complex unconjugated equilibrium projection; connect the mixed primitive before contour-block reconstruction"
        :"X_conn^{ab}(t,tau)=X^{ab}(t,tau)-m_a(t)m_b(0), complex unconjugated; connect the mixed primitive before contour-block reconstruction" );
    hdf5r::store_string( ps_group_id, "mean_field_strategy",
        pspace.uses_harmonic_bath()
        ?"zero prescribed bath mean; external magnetic field and local interactions remain active"
        :pspace.constant_magnetization_time
        ?"mean_M=mean_+(t)=mean_-(t)=JL D Re m(0), equilibrium projection"
        :"mean_M=JL D Re m(0); mean_+(t)=mean_-(t)=JL D Re m(t)" );
    hdf5r::store_string( ps_group_id, "constant_magnetization_time",
        hdf5r::bool_to_string(pspace.constant_magnetization_time) );
    hdf5r::store_scalar( ps_group_id, "num_blocks",pspace.num_blocks );
    if( pspace.sampling_strategy=="independent" )
        hdf5r::store_scalar( ps_group_id, "num_ratio_jackknife_blocks",
                             pspace.num_blocks );
    hdf5r::store_string( ps_group_id, "field_covariance",
        pspace.uses_harmonic_bath()
        ?"complex symmetric pseudo-covariance E[V V^T] from the exact contour correlation of X=g(a+a^dagger), on M,+,- physical field values"
        :"complex symmetric pseudo-covariance E[V V^T] from the connected spin correlation, on M,+,- physical field values" );
    hdf5r::store_string( ps_group_id, "contour_convention",
        "0 -> T (+), T -> 0 (-), 0 -> -i beta (M); q_M=-1, q_+=-i, q_-=+i" );
    hdf5r::store_string( ps_group_id, "field_layout",
        "flat(A,n,a)=3*(offset_A+n)+a; offsets M=0, +=N_tau+1, -=N_tau+1+N_R" );
    hdf5r::store_string( ps_group_id, "imaginary_field_grid",
        "one-sided edge grid; tau=k*delta_imag_t for k=0,...,N_tau with distinct 0+ and beta- fields" );

    // ...concerning the iteration
    hdf5r::store_scalar( ps_group_id, "iteration_error_sigma_threshold",
                         pspace.iteration_error_sigma_threshold );
    hdf5r::store_scalar( ps_group_id, "Iteration_Limit",            pspace.Iteration_Limit );
    hdf5r::store_scalar( ps_group_id, "fixed_point_mixing_alpha",   pspace.mixing_alpha );
    hdf5r::store_scalar( ps_group_id, "covariance_tolerance",       pspace.covariance_tolerance );
    hdf5r::store_scalar( ps_group_id, "branch_identity_tolerance",  pspace.branch_identity_tolerance );
    hdf5r::store_scalar( ps_group_id, "takagi_tolerance",           pspace.takagi_tolerance );
    hdf5r::store_scalar( ps_group_id, "minimum_phase_magnitude",    pspace.minimum_phase_magnitude );
    hdf5r::store_scalar( ps_group_id, "denominator_constancy_tolerance", pspace.denominator_constancy_tolerance );
    hdf5r::store_scalar( ps_group_id, "imaginary_magnetization_sigma", pspace.imaginary_magnetization_sigma );

    // ...concerning the initially inserted correlations
    hdf5r::store_string( ps_group_id, "init_diag_corr",             hdf5r::none_if_empty(pspace.init_diag_corr.m_name) );
    for( const auto& param : pspace.init_diag_corr.m_parameters ) // add all parameters of the initial diagonal correlations
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "init_nondiag_corr",          hdf5r::none_if_empty(pspace.init_nondiag_corr.m_name) );
    for( const auto& param : pspace.init_nondiag_corr.m_parameters ) // add all parameters of the initial non-diagonal correlations
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "load_initial_spin_correlations",         hdf5r::bool_to_string(pspace.load_initial_spin_correlations) );
    hdf5r::store_string( ps_group_id, "extrapolate_initial_spin_correlations",  hdf5r::bool_to_string(pspace.extrapolate_initial_spin_correlations) );
    hdf5r::store_string( ps_group_id, "initial_correlations_src_file",          hdf5r::none_if_empty(pspace.initial_correlations_src_file) );
    hdf5r::store_string( ps_group_id, "initial_correlations_src_directory",     hdf5r::none_if_empty(pspace.initial_correlations_src_directory) );

    // ========== storing and naming ==========
    hdf5r::store_string( ps_group_id, "information_text",           hdf5r::none_if_empty(pspace.information_text) );
    std::string RealType_str{};
    #ifdef USE_DOUBLE
    RealType_str = "DOUBLE";
    #elif defined USE_FLOAT
    RealType_str = "FLOAT";
    #endif
    hdf5r::store_string( ps_group_id, "RealType",                   RealType_str);
    hdf5r::store_string( ps_group_id, "system_name",                fm::get_system_name());
    hdf5r::store_string( ps_group_id, "(original) project_name",                hdf5r::none_if_empty(pspace.project_name) );
    hdf5r::store_string( ps_group_id, "(original) filename_extension",          hdf5r::none_if_empty(pspace.filename_extension) );
    H5Gclose( ps_group_id );


    // ++++++++++ store run time data ++++++++++
    auto rtd_group_id = H5Gcreate( m_file_id, "runtimedata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hdf5r::store_scalar( rtd_group_id, "generated_seed",            rtdata.generated_seed );
    hdf5r::store_list( rtd_group_id, "covariance_symmetry_errors",  rtdata.covariance_symmetry_errors );
    hdf5r::store_list( rtd_group_id, "branch_identity_errors",      rtdata.branch_identity_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_factor_reconstruction_errors",rtdata.gaussian_factor_reconstruction_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_covariance_approximation_errors",rtdata.gaussian_covariance_approximation_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_factor_latent_dimensions",      rtdata.gaussian_factor_latent_dimensions );
    hdf5r::store_list( rtd_group_id, "gaussian_largest_factorization_dimensions",rtdata.gaussian_largest_factorization_dimensions );
    hdf5r::store_list( rtd_group_id, "average_phase_magnitudes",    rtdata.average_phase_magnitudes );
    hdf5r::store_list( rtd_group_id,
        pspace.sampling_strategy=="pcn"
        ?"pcn_minimum_effective_sample_sizes"
        :"complex_weight_effective_sample_sizes",
        rtdata.effective_sample_sizes );
    hdf5r::store_list( rtd_group_id, "partition_sum_Re",            rtdata.partition_sum_Re );
    hdf5r::store_list( rtd_group_id, "partition_sum_Im",            rtdata.partition_sum_Im );
    hdf5r::store_list( rtd_group_id, "partition_abs_sum",           rtdata.partition_abs_sums );
    hdf5r::store_list( rtd_group_id, "denominator_constancy_residuals", rtdata.denominator_constancy_residuals );
    if( pspace.sampling_strategy=="pcn" )
    {
        hdf5r::store_list( rtd_group_id, "mh_acceptance_rates",       rtdata.mh_acceptance_rates );
        hdf5r::store_list( rtd_group_id, "mh_nonpositive_rejection_rates",rtdata.mh_nonpositive_rejection_rates );
        hdf5r::store_list( rtd_group_id, "maximum_relative_imaginary_sampling_weights",rtdata.maximum_relative_imaginary_sampling_weights );
        hdf5r::store_list( rtd_group_id, "blocking_curve_block_length",rtdata.blocking_curve_block_lengths );
        hdf5r::store_list( rtd_group_id, "blocking_curve_mean_error", rtdata.blocking_curve_mean_errors );
        hdf5r::store_list( rtd_group_id, "blocking_curve_max_error",  rtdata.blocking_curve_max_errors );
        hdf5r::store_list( rtd_group_id, "maximum_correlation_tau_int",rtdata.maximum_tau_int );
        hdf5r::store_list( rtd_group_id, "base_block_length_to_max_tau_ratio",rtdata.block_length_to_tau_ratios );
    }
    hdf5r::store_list( rtd_group_id, "closed_contour_ratio_Re",     rtdata.closed_contour_ratio_Re );
    hdf5r::store_list( rtd_group_id, "closed_contour_ratio_Im",     rtdata.closed_contour_ratio_Im );

    // ...concerning the statistics
    hdf5r::store_scalar( rtd_group_id, "actual_num_blocks_per_core",rtdata.num_blocks() );
    hdf5r::store_scalar( rtd_group_id, "samples_per_block",rtdata.samples_per_block() );
    if( pspace.sampling_strategy=="independent" )
    {
        hdf5r::store_scalar( rtd_group_id, "actual_num_ratio_jackknife_blocks_per_core", rtdata.num_blocks() );
        hdf5r::store_scalar( rtd_group_id, "samples_per_ratio_jackknife_block", rtdata.samples_per_block() );
    }
    const std::string contour_err_remark = pspace.sampling_strategy=="pcn"
        ?pspace.antithetic_pairs
         ?" Standard error from contiguous sign-symmetrized pCN pair-state batch means, taking the largest resolved value along the power-of-two blocking curve; blocks are pooled over independent MPI chains."
         :" Standard error from contiguous pCN batch means, taking the largest resolved value along the power-of-two blocking curve; blocks are pooled over independent MPI chains."
        :pspace.antithetic_pairs
         ?" Delete-one-block jackknife standard error of the complex ratio (sum N)/(sum Z_M), pooling equal-sized blocks of complete r,-r antithetic pairs over MPI ranks."
         :" Delete-one-block jackknife standard error of the paired complex ratio (sum N)/(sum Z_M), pooling equal-sized independent-sample blocks over MPI ranks.";
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_Re_sample_stds",
        rtdata.closed_contour_residual_Re_sample_stds );
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_Im_sample_stds",
        rtdata.closed_contour_residual_Im_sample_stds );
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_abs_sample_stds",
        rtdata.closed_contour_residual_abs_sample_stds );
    store_correlation( standard_errors.Re, rtd_group_id,
        "Re_correlation_sample_stds",
        "Standard errors of edge-grid Re X; axes t,direction_pair,tau_edge." + contour_err_remark );
    store_correlation( standard_errors.Im, rtd_group_id,
        "Im_correlation_sample_stds",
        "Standard errors of edge-grid Im X; axes t,direction_pair,tau_edge." + contour_err_remark );
    if( pspace.sampling_strategy=="pcn" )
    {
        store_correlation( rtdata.contour_tau_int.Re, rtd_group_id,
            "Re_correlation_tau_int",
            pspace.antithetic_pairs
            ?"Integrated autocorrelation time in sign-symmetrized pCN pair steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge."
            :"Integrated autocorrelation time in pCN steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge." );
        store_correlation( rtdata.contour_tau_int.Im, rtd_group_id,
            "Im_correlation_tau_int",
            pspace.antithetic_pairs
            ?"Integrated autocorrelation time in sign-symmetrized pCN pair steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge."
            :"Integrated autocorrelation time in pCN steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge." );
    }

    // Components forbidden by the symmetry are exact zeros.
    hdf5r::store_2D_tensor<RealType>( rtd_group_id, "Re_magnetization_sample_stds",
        H5_REAL_TYPE, rtdata.magnetization_time_Re_stds,
        "Axes are real_time, spin_component(x,y,z)." );
    hdf5r::store_2D_tensor<RealType>( rtd_group_id, "Im_magnetization_sample_stds",
        H5_REAL_TYPE, rtdata.magnetization_time_Im_stds,
        "Axes are real_time, spin_component(x,y,z)." );


    // ...concerning the iteration
    hdf5r::store_scalar( rtd_group_id, "num_Iterations",            rtdata.num_Iterations );
    hdf5r::store_string( rtd_group_id, "termination",               rtdata.termination );
    hdf5r::store_list( rtd_group_id, "absolute_iteration_errors",   rtdata.absolute_iteration_errors );
    hdf5r::store_list( rtd_group_id, "standardized_iteration_errors",
                       rtdata.standardized_iteration_errors );

    H5Gclose( rtd_group_id );


    // ++++++++++ store results ++++++++++
    auto results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    std::vector<std::string> direction_labels;
    for( const auto& direction : correlations.Re.front().get_direction_pairs() )
    {
        const std::string axes = "xyz";
        direction_labels.emplace_back(
            std::string(1,axes[direction[0]]) + std::string(1,axes[direction[1]]) );
    }
    hdf5r::store_string_list( results_group_id, "correlation_direction_labels",
                              direction_labels );
    std::vector<RealType> imaginary_edges(pspace.num_TimePoints),
                          real_times(pspace.num_RealTimePoints);
    for( size_t k=0;k<imaginary_edges.size();++k ) imaginary_edges[k]=k*pspace.delta_t;
    for( size_t t=0;t<real_times.size();++t ) real_times[t]=t*pspace.delta_real_t;
    hdf5r::store_list( results_group_id, "imaginary_time_edges", imaginary_edges );
    hdf5r::store_list( results_group_id, "real_times", real_times );
    store_correlation( correlations.Re, results_group_id,
        "Re_correlation",
        "Re X^{ab}(t,tau_edge)=Re <S_b(-i tau_edge) S_a(t)>; axes t,direction_pair(a,b),tau_edge." );
    store_correlation( correlations.Im, results_group_id,
        "Im_correlation",
        "Im X^{ab}(t,tau_edge); axes t,direction_pair(a,b),tau_edge." );
    hdf5r::store_2D_tensor<RealType>( results_group_id, "Re_magnetization",
        H5_REAL_TYPE, rtdata.magnetization_time_Re,
        "Re m_a(t); axes are real_time, spin_component(x,y,z)." );
    hdf5r::store_2D_tensor<RealType>( results_group_id, "Im_magnetization",
        H5_REAL_TYPE, rtdata.magnetization_time_Im,
        "Im m_a(t); axes are real_time, spin_component(x,y,z)." );
    H5Gclose( results_group_id );
}

// store all measured times
void HDF5_Storage::store_time( const tmm::DerivedTimeMeasure& tmeasure )
{
    if( !m_storing_permission ){ return; } // permission request

    // ++++++++++ store times ++++++++++
    auto td_group_id = H5Gcreate( m_file_id, "timedata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    hdf5r::store_string( td_group_id, "info", "all durations are measured in seconds" );

    // global measures
    std::for_each( tmeasure.m_global_durations.cbegin(), tmeasure.m_global_durations.cend(),
    [&td_group_id]( const tmm::DurationQuantity& q )
    {
        hdf5r::store_scalar( td_group_id, q.m_name, q.m_duration );
    } );

    // tmp measures
    std::for_each( tmeasure.m_tmp_measures.cbegin(), tmeasure.m_tmp_measures.cend(),
    [&td_group_id]( const tmm::IterationDurationQuantity& q )
    {
        hdf5r::store_list( td_group_id, q.m_name, q.m_durations );
    } );

    // tmp measure iteration averages
    std::for_each( tmeasure.m_tmp_measures.cbegin(), tmeasure.m_tmp_measures.cend(),
    [&td_group_id]( const tmm::IterationDurationQuantity& q )
    {
        hdf5r::store_scalar( td_group_id, q.m_name + " (av)", q.average() );
    } );

    // total duration and start date and time
    hdf5r::store_scalar( td_group_id, "total_duration", tmeasure.m_total_duration );
    hdf5r::store_string( td_group_id, "start_date_and_time", std::string{std::ctime(&tmeasure.m_start_date_and_time)} );

    // close group
    H5Gclose( td_group_id );
}

void HDF5_Storage::store_correlation( const ContourCorrelation& correlations,
                                     const hid_t group_id,
                                     const std::string& dataset_name,
                                     const std::string& dataset_info )
{
    if( !m_storing_permission ){ return; }
    hdf5r::store_3D_tensor<RealType>( group_id, dataset_name, H5_REAL_TYPE,
                                     correlations, dataset_info );
}

// finalize, i.e., close file
void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}

};
