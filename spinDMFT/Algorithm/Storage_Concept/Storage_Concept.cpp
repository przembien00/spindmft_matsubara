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
    if( !pspace.self_consistency )
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
void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata, const CorrTen& corr_R, const CorrTen& corr_I, const FieldVector& spin_expval )
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
    hdf5r::store_scalar( ps_group_id, "num_TimeSteps",              pspace.num_TimeSteps ); 
    hdf5r::store_scalar( ps_group_id, "num_TimePoints",             pspace.num_TimePoints ); 
    hdf5r::store_scalar( ps_group_id, "delta_t",                    pspace.delta_t );

    // ...concerning statistics
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerCore",         pspace.num_SamplesPerCore ); 
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerSet",          pspace.num_SamplesPerSet ); 
    hdf5r::store_scalar( ps_group_id, "num_Samples",                pspace.num_Samples ); 
    hdf5r::store_scalar( ps_group_id, "num_Cores",                  pspace.get_num_Cores() ); 

    hdf5r::store_string( ps_group_id, "seed",                       pspace.seed );

    // ...concerning the iteration 
    hdf5r::store_scalar( ps_group_id, "absolute_iteration_error_threshold",     pspace.absolute_iteration_error_threshold );
    hdf5r::store_scalar( ps_group_id, "Iteration_Limit",            pspace.Iteration_Limit ); 

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

    // ...concerning negative eigenvalues of the covariance matrix
    hdf5r::store_string( ps_group_id, "truncation_scheme_negative_eigenvalues", pspace.truncation_scheme_negative_eigenvalues );
    hdf5r::store_scalar( ps_group_id, "critical_eigenvalue_ratio",              pspace.critical_eigenvalue_ratio );

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

    // ...concerning the eigenvalues
    hdf5r::store_scalar( rtd_group_id, "generated_seed",            rtdata.generated_seed );
    hdf5r::store_list( rtd_group_id, "negative_eigenvalue_ratios",  rtdata.negative_eigenvalue_ratios );

    // ...concerning the statistics
    store_correlation_tensor( rtdata.sample_stds_Re, rtd_group_id, "Re_correlation_sample_stds", "Standard deviations of the correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );
    store_correlation_tensor( rtdata.sample_stds_Im, rtd_group_id, "Im_correlation_sample_stds", "Standard deviations of the correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );


    // ...concerning the iteration
    hdf5r::store_scalar( rtd_group_id, "num_Iterations",            rtdata.num_Iterations );
    hdf5r::store_string( rtd_group_id, "termination",               rtdata.termination );
    hdf5r::store_list( rtd_group_id, "absolute_iteration_errors",   rtdata.absolute_iteration_errors );
    H5Gclose( rtd_group_id );


    // ++++++++++ store results ++++++++++
    auto results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    store_correlation_tensor( corr_R, results_group_id, "Re_correlation", "Correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );
    store_correlation_tensor( corr_I, results_group_id, "Im_correlation", "Correlations <S^alpha(t)S^beta(0)>, stored according to the hierarchy alpha-beta, t" );
    hdf5r::store_scalar( results_group_id, "S_x", spin_expval[0] );
    hdf5r::store_scalar( results_group_id, "S_y", spin_expval[1] );
    hdf5r::store_scalar( results_group_id, "S_z", spin_expval[2] );
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

// store correlation tensor
void HDF5_Storage::store_correlation_tensor( const CorrTen& CT, const hid_t group_id, const std::string dataset_name, const std::string dataset_info )
{   
    if( !m_storing_permission ){ return; } // permission request

    hdf5r::store_2D_tensor<RealType>(group_id, dataset_name, H5_REAL_TYPE, CT, dataset_info);
}

// finalize, i.e., close file
void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}

};
