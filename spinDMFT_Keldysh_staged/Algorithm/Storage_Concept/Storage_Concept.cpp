#include"Storage_Concept.h"
#include<fstream>
#include<iomanip>
#include<limits>
#include<sstream>
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
HDF5_Storage::HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace, const std::string& termination,
                            const std::string& filename_suffix ):
    m_storing_permission( my_rank == 0 ),
    m_fname_max_length( 200 ),
    m_num_TriesToBuildFile( 5 )
{
    create_folder_branch( pspace, termination );
    create_file( pspace,filename_suffix );
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
void HDF5_Storage::create_file( const ps::ParameterSpace& pspace, const std::string& filename_suffix )
{
    if( !m_storing_permission ){ return; } // permission request

    // create filename:
    std::string prefix{"staged__"};
    if(pspace.sampling_strategy=="pcn")prefix+="jointMpcn_absZ__";
    std::string filename = prefix+pspace.spin_model.compact_info( pspace.num_PrintDigits ); // spin model info
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
    // Include the requested evolution even if preparation failed before the
    // real-time stage, so those outputs also retain their run identity.
    if(pspace.real_time_evolution=="staggered-strang")
    {
        filename += "__prop=staggered-strang";
        if(pspace.real_time_substeps>1)
            filename += "__substeps="+std::to_string(pspace.real_time_substeps);
    }
    const std::string extension=pspace.filename_extension+filename_suffix;
    if( !extension.empty() )
    {
        filename += "_" + extension;
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
                              const Magnetization& magnetization, bool preparation )
{
    if( !m_storing_permission ){ return; } // permission request

    const auto& statistics=preparation?rtdata.preparation:rtdata.real_time;

    // ++++++++++ store parameters ++++++++++
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    hdf5r::store_string(ps_group_id,"joint_matsubara_pcn",hdf5r::bool_to_string(pspace.sampling_strategy=="pcn"));
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

    // ...concerning statistics
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerCore",         pspace.num_SamplesPerCore );
    hdf5r::store_scalar( ps_group_id, "num_Samples",                pspace.num_Samples );
    hdf5r::store_scalar(ps_group_id,"num_PreparationSamplesPerCore",pspace.num_PreparationSamplesPerCore);
    hdf5r::store_scalar(ps_group_id,"num_PreparationSamples",pspace.num_PreparationSamples);
    hdf5r::store_scalar(ps_group_id,"Preparation_Iteration_Limit",pspace.Preparation_Iteration_Limit);
    hdf5r::store_scalar( ps_group_id, "num_Cores",                  pspace.get_num_Cores() );

    hdf5r::store_string( ps_group_id, "seed",                       pspace.seed );

    // ...concerning complex-field sampling
    hdf5r::store_string( ps_group_id, "sampling_strategy",
                         pspace.sampling_strategy );
    hdf5r::store_string( ps_group_id, "pcn_sampling_weight",
        pspace.sampling_strategy!="pcn"
        ?"not applicable"
        :"|Z_M|; positive modulus of the Matsubara partition function" );
    hdf5r::store_string(ps_group_id,"pcn_updated_coordinates",
        pspace.sampling_strategy=="pcn"?"Matsubara-active real row space of [Re L_M; Im L_M] from the current joint Takagi factor":"not applicable");
    hdf5r::store_string(ps_group_id,"real_residual_sampling",
        pspace.sampling_strategy=="pcn"?"fresh independent standard normals projected onto the joint Matsubara null space after every production pCN step, including rejection; x=Ua+(I-UU^T)xi":"direct Gaussian draws");
    hdf5r::store_scalar( ps_group_id, "mh_step_size",
                         pspace.mh_step_size );
    hdf5r::store_scalar( ps_group_id, "mh_burn_in",
                         pspace.mh_burn_in );

    if(pspace.sampling_strategy=="pcn")
        hdf5r::store_string(ps_group_id,"pcn_row_space_rank_rule",
            "singular value > epsilon*max(2*nM,q)*largest singular value; relative reconstruction checked; rebuilt and chain burned in each iteration");
    hdf5r::store_scalar(ps_group_id,"real_time_substeps",pspace.real_time_substeps);
    hdf5r::store_scalar(ps_group_id,"delta_real_propagation_t",
        pspace.delta_real_t/static_cast<RealType>(pspace.real_time_steps_per_interval()));
    {
        // HDF5's FALSE/TRUE enum is read as a boolean by h5py.
        const hid_t boolean_type=H5Tenum_create(H5T_NATIVE_UCHAR);
        const unsigned char false_value=0,true_value=1;
        H5Tenum_insert(boolean_type,"FALSE",&false_value);
        H5Tenum_insert(boolean_type,"TRUE",&true_value);
        const hid_t space=H5Screate(H5S_SCALAR);
        const hid_t attribute=H5Acreate(ps_group_id,"constant_magnetization_time",
            boolean_type,space,H5P_DEFAULT,H5P_DEFAULT);
        const unsigned char value=pspace.constant_magnetization_time?true_value:false_value;
        H5Awrite(attribute,boolean_type,&value);
        H5Aclose(attribute);
        H5Sclose(space);
        H5Tclose(boolean_type);
    }
    hdf5r::store_scalar( ps_group_id, "num_blocks",pspace.num_blocks );
    if( pspace.sampling_strategy=="independent" )
        hdf5r::store_scalar( ps_group_id, "num_ratio_jackknife_blocks",
                             pspace.num_blocks );

    // ...concerning the iteration
    hdf5r::store_scalar( ps_group_id, "iteration_error_sigma_threshold",
                         pspace.iteration_error_sigma_threshold );
    hdf5r::store_scalar( ps_group_id, "Iteration_Limit",            pspace.Iteration_Limit );
    hdf5r::store_scalar(ps_group_id,"real_time_absolute_residual_ceiling",
        std::max(pspace.spin_float,pspace.spin_float*pspace.spin_float));
    hdf5r::store_scalar( ps_group_id, "covariance_tolerance",       pspace.covariance_tolerance );
    hdf5r::store_scalar( ps_group_id, "takagi_tolerance",           pspace.takagi_tolerance );
    hdf5r::store_scalar( ps_group_id, "minimum_phase_magnitude",    pspace.minimum_phase_magnitude );

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
    hdf5r::store_list( rtd_group_id, "covariance_symmetry_errors",  statistics.covariance_symmetry_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_factor_reconstruction_errors",statistics.gaussian_factor_reconstruction_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_covariance_approximation_errors",statistics.gaussian_covariance_approximation_errors );
    hdf5r::store_list( rtd_group_id, "gaussian_factor_latent_dimensions",      statistics.gaussian_factor_latent_dimensions );
    hdf5r::store_list( rtd_group_id, "gaussian_largest_factorization_dimensions",statistics.gaussian_largest_factorization_dimensions );
    hdf5r::store_list( rtd_group_id, "average_phase_magnitudes",    statistics.average_phase_magnitudes );
    hdf5r::store_list( rtd_group_id,
        pspace.sampling_strategy=="pcn"
        ?"pcn_minimum_effective_sample_sizes"
        :"complex_weight_effective_sample_sizes",
        statistics.effective_sample_sizes );
    hdf5r::store_list( rtd_group_id, "partition_sum_Re",            statistics.partition_sum_Re );
    hdf5r::store_list( rtd_group_id, "partition_sum_Im",            statistics.partition_sum_Im );
    hdf5r::store_list( rtd_group_id, "partition_abs_sum",           statistics.partition_abs_sums );
    hdf5r::store_list( rtd_group_id, "denominator_constancy_residuals", statistics.denominator_constancy_residuals );
    if( pspace.sampling_strategy=="pcn" )
    {
        hdf5r::store_list( rtd_group_id, "mh_acceptance_rates",       statistics.mh_acceptance_rates );
        hdf5r::store_list( rtd_group_id, "mh_zero_weight_rejection_rates",statistics.mh_nonpositive_rejection_rates );
        hdf5r::store_list( rtd_group_id, "blocking_curve_block_length",statistics.blocking_curve_block_lengths );
        hdf5r::store_list( rtd_group_id, "blocking_curve_mean_error", statistics.blocking_curve_mean_errors );
        hdf5r::store_list( rtd_group_id, "blocking_curve_max_error",  statistics.blocking_curve_max_errors );
        hdf5r::store_list( rtd_group_id, "maximum_correlation_tau_int",statistics.maximum_tau_int );
        hdf5r::store_list( rtd_group_id, "base_block_length_to_max_tau_ratio",statistics.block_length_to_tau_ratios );
    }
    hdf5r::store_list( rtd_group_id, "closed_contour_ratio_Re",     statistics.closed_contour_ratio_Re );
    hdf5r::store_list( rtd_group_id, "closed_contour_ratio_Im",     statistics.closed_contour_ratio_Im );

    // ...concerning the statistics
    hdf5r::store_scalar( rtd_group_id, "actual_num_blocks_per_core",statistics.num_blocks );
    hdf5r::store_scalar( rtd_group_id, "samples_per_block",statistics.samples_per_block );
    if( pspace.sampling_strategy=="independent" )
    {
        hdf5r::store_scalar( rtd_group_id, "actual_num_ratio_jackknife_blocks_per_core", statistics.num_blocks );
        hdf5r::store_scalar( rtd_group_id, "samples_per_ratio_jackknife_block", statistics.samples_per_block );
    }
    const std::string contour_err_remark = pspace.sampling_strategy=="pcn"
        ?" Standard error from contiguous pCN batch means, taking the largest resolved value along the power-of-two blocking curve; blocks are pooled over independent MPI chains."
        :" Delete-one-block jackknife standard error of the paired complex ratio (sum N)/(sum Z_M), pooling equal-sized independent-sample blocks over MPI ranks.";
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_Re_sample_stds",
        statistics.closed_contour_residual_Re_sample_stds );
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_Im_sample_stds",
        statistics.closed_contour_residual_Im_sample_stds );
    hdf5r::store_list( rtd_group_id, "closed_contour_residual_abs_sample_stds",
        statistics.closed_contour_residual_abs_sample_stds );
    store_correlation( statistics.contour_sample_stds.Re, rtd_group_id,
        "Re_correlation_sample_stds",
        "Standard errors of edge-grid Re X; axes t,direction_pair,tau_edge." + contour_err_remark );
    store_correlation( statistics.contour_sample_stds.Im, rtd_group_id,
        "Im_correlation_sample_stds",
        "Standard errors of edge-grid Im X; axes t,direction_pair,tau_edge." + contour_err_remark );
    if( pspace.sampling_strategy=="pcn" )
    {
        store_correlation( statistics.contour_tau_int.Re, rtd_group_id,
            "Re_correlation_tau_int",
            "Integrated autocorrelation time in pCN steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge." );
        store_correlation( statistics.contour_tau_int.Im, rtd_group_id,
            "Im_correlation_tau_int",
            "Integrated autocorrelation time in pCN steps inferred from batch-mean variance inflation; axes t,direction_pair,tau_edge." );
    }

    // Components forbidden by the symmetry are exact zeros.
    hdf5r::store_2D_tensor<RealType>( rtd_group_id, "Re_magnetization_sample_stds",
        H5_REAL_TYPE, statistics.magnetization_time_stds.Re.expand(),
        "Axes are real_time, spin_component(x,y,z)." );
    hdf5r::store_2D_tensor<RealType>( rtd_group_id, "Im_magnetization_sample_stds",
        H5_REAL_TYPE, statistics.magnetization_time_stds.Im.expand(),
        "Axes are real_time, spin_component(x,y,z)." );


    // ...concerning the iteration
    hdf5r::store_scalar( rtd_group_id, "num_Iterations",            statistics.num_Iterations );
    hdf5r::store_string( rtd_group_id, "termination",               statistics.termination );
    hdf5r::store_list( rtd_group_id, "absolute_iteration_errors",   statistics.absolute_iteration_errors );
    hdf5r::store_list( rtd_group_id, "standardized_iteration_errors",
                       statistics.standardized_iteration_errors );

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
                          real_times(correlations.Re.size());
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
        H5_REAL_TYPE, magnetization.Re.expand(),
        "Re m_a(t); axes are real_time, spin_component(x,y,z)." );
    hdf5r::store_2D_tensor<RealType>( results_group_id, "Im_magnetization",
        H5_REAL_TYPE, magnetization.Im.expand(),
        "Im m_a(t); axes are real_time, spin_component(x,y,z)." );
    H5Gclose( results_group_id );
}

void HDF5_Storage::store_preparation(const ps::ParameterSpace& p,
    const rtd::RunTimeData& data,const CorrelationSet& frozen,const Magnetization& magnetization,
    bool real_started)
{
    if(!m_storing_permission)return;
    const auto& statistics=data.preparation;
    const auto store_matsubara_correlation=[&](const ContourCorrelation& values,
        hid_t group,const std::string& name,const std::string& info)
    {
        if(values.size()!=1)
            throw std::logic_error("preparation correlations require one real-time slice");
        hdf5r::store_2D_tensor<RealType>(group,name,H5_REAL_TYPE,values.front(),
            info+"; axes direction,tau");
    };
    const hid_t group=H5Gcreate(m_file_id,"preparation",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    hdf5r::store_string(group,"termination",statistics.termination);
    hdf5r::store_string(group,"real_time_started",hdf5r::bool_to_string(real_started));
    hdf5r::store_scalar(group,"num_Iterations",statistics.num_Iterations);
    hdf5r::store_scalar(group,"absolute_iteration_error_threshold",data.matsubara_iteration_error_threshold());
    hdf5r::store_scalar(group,"samples_per_core",p.num_PreparationSamplesPerCore);
    hdf5r::store_scalar(group,"iteration_limit",p.Preparation_Iteration_Limit);
    hdf5r::store_scalar(group,"generated_seed",data.generated_seed);
    hdf5r::store_string(group,"sampling_strategy",p.sampling_strategy);
    if(p.sampling_strategy=="pcn")
    {
        hdf5r::store_string(group,"pcn_sampling_weight","|Z_M|");
        hdf5r::store_string(group,"normalization_strategy","sum[N/|Z_M|]/sum[Z_M/|Z_M|]");
        hdf5r::store_scalar(group,"mh_step_size",p.mh_step_size);
        hdf5r::store_scalar(group,"mh_burn_in",p.mh_burn_in);
        hdf5r::store_scalar(group,"actual_num_blocks_per_core",statistics.num_blocks);
        hdf5r::store_scalar(group,"samples_per_block",statistics.samples_per_block);
        hdf5r::store_list(group,"mh_acceptance_rates",statistics.mh_acceptance_rates);
        hdf5r::store_list(group,"mh_zero_weight_rejection_rates",statistics.mh_nonpositive_rejection_rates);
        hdf5r::store_list(group,"blocking_curve_block_length",statistics.blocking_curve_block_lengths);
        hdf5r::store_list(group,"blocking_curve_mean_error",statistics.blocking_curve_mean_errors);
        hdf5r::store_list(group,"blocking_curve_max_error",statistics.blocking_curve_max_errors);
        hdf5r::store_list(group,"maximum_correlation_tau_int",statistics.maximum_tau_int);
        hdf5r::store_list(group,"base_block_length_to_max_tau_ratio",statistics.block_length_to_tau_ratios);
        store_matsubara_correlation(statistics.contour_tau_int.Re,group,"Re_correlation_tau_int","pCN batch-means variance inflation");
        store_matsubara_correlation(statistics.contour_tau_int.Im,group,"Im_correlation_tau_int","pCN batch-means variance inflation");
    }
    store_matsubara_correlation(frozen.Re,group,"Re_correlation","frozen X(0,tau)");
    store_matsubara_correlation(frozen.Im,group,"Im_correlation","frozen X(0,tau)");
    store_matsubara_correlation(statistics.contour_sample_stds.Re,group,"Re_correlation_sample_stds","final preparation estimate standard errors");
    store_matsubara_correlation(statistics.contour_sample_stds.Im,group,"Im_correlation_sample_stds","final preparation estimate standard errors");
    hdf5r::store_2D_tensor<RealType>(group,"Re_magnetization",H5_REAL_TYPE,magnetization.Re.expand(),"frozen origin magnetization");
    hdf5r::store_2D_tensor<RealType>(group,"Im_magnetization",H5_REAL_TYPE,magnetization.Im.expand(),"frozen origin magnetization");
    hdf5r::store_2D_tensor<RealType>(group,"Re_magnetization_sample_stds",H5_REAL_TYPE,statistics.magnetization_time_stds.Re.expand(),"last measured origin standard errors");
    hdf5r::store_2D_tensor<RealType>(group,"Im_magnetization_sample_stds",H5_REAL_TYPE,statistics.magnetization_time_stds.Im.expand(),"last measured origin standard errors");
    hdf5r::store_list(group,"absolute_iteration_errors",statistics.absolute_iteration_errors);
    hdf5r::store_list(group,"standardized_iteration_errors",statistics.standardized_iteration_errors);
    hdf5r::store_list(group,"gaussian_factor_reconstruction_errors",statistics.gaussian_factor_reconstruction_errors);
    hdf5r::store_list(group,"average_phase_magnitudes",statistics.average_phase_magnitudes);
    hdf5r::store_list(group,"effective_sample_sizes",statistics.effective_sample_sizes);
    hdf5r::store_list(group,"partition_sum_Re",statistics.partition_sum_Re);
    hdf5r::store_list(group,"partition_sum_Im",statistics.partition_sum_Im);
    if(!data.frozen_errors.empty())hdf5r::store_list(group,"frozen_MM_residuals",data.frozen_errors);
    H5Gclose(group);
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
