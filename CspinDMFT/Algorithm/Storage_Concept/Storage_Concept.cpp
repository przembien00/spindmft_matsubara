#include"Storage_Concept.h"
#include<iostream>
#include<iomanip>
#include<string>
#include<sys/types.h>
#include<sys/stat.h>
#include<fstream>

#include<Standard_Algorithms/Standard_Algorithms.h>
namespace stda = Standard_Algorithms;

#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;

#include<File_Management/File_Management.h>
namespace fm = File_Management;

#include"Error_Handling.h"
namespace error = Storage_Concept::Error_Handling;

namespace hdf5r = HDF5_Routines;

namespace Storage_Concept
{
// ==================================================================
// ======================= HDF5 STORAGE CLASS =======================
// ==================================================================
// constructor : create folder tree and data file
HDF5_Storage::HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata ):
    m_storing_permission( my_rank == 0 ),
    m_fname_max_length( 200 ),
    m_num_TriesToBuildFile( 5 )
{
    create_folder_branch( pspace, rtdata.termination );
    create_file( pspace, rtdata.eigenvalue_threshold_violated_list.back() );
}

// create the folder branch in which the data will be stored
void HDF5_Storage::create_folder_branch( const ps::ParameterSpace& pspace, const std::string& termination )
{
    if( !m_storing_permission ){ return; } // permission request

    // create folder branch list:
    std::vector<std::string> folder_branch_list{};
    folder_branch_list.push_back( "Data" );
    if( pspace.dst_project_name != "" )
    {
        folder_branch_list.push_back( pspace.dst_project_name );
    }
    if( termination == "by hand" || termination == "by iteration limit" )
    {
        folder_branch_list.push_back( "NOT_CONVERGED" );
    }
    // fbranch.push_back( mf_model_name ); CorrelationReplica is the only implemented algorithm at the moment (!)

    // create folder branch:
    m_filename = fm::create_folder_tree( folder_branch_list, m_fname_max_length );
}

// create the file in which the data will be stored
void HDF5_Storage::create_file( const ps::ParameterSpace& pspace, const bool eigenvalue_threshold_violated )
{
    if( !m_storing_permission ){ return; } // permission request

    // create filename:
    std::string filename{};
    filename += pspace.spinspin_cmodel.compact_info(pspace.num_PrintDigits);
    if( pspace.spinmf_cmodel != pspace.spinspin_cmodel )
    {
        filename += "__mf" + pspace.spinmf_cmodel.compact_info(pspace.num_PrintDigits);
    }
    filename += "__config=" + pspace.config_file;
    if( pspace.chemical_shift.m_name != "none" )
    {
        filename += "__" + pspace.chemical_shift.compact_info(pspace.num_PrintDigits);
    }
    if( pspace.local_extra_interaction.m_name != "none" )
    {
        filename += "__" + pspace.local_extra_interaction.compact_info( pspace.num_PrintDigits ); // extra_interaction info
    }

    // put rescale into the filename if the factors are NOT very close to 1
    if( abs( pspace.rescale_meanfield - RealType{1.0} ) > RealType{pow(10,-6)} ) 
    {
        filename += "__rescmf=" + print::remove_zeros(print::round_value_to_string(pspace.rescale_meanfield,pspace.num_PrintDigits));
    }
    if( abs( pspace.rescale_spinspin - RealType{1.0} ) > RealType{pow(10,-6)} ) 
    {
        filename += "__rescss=" + print::remove_zeros(print::round_value_to_string(pspace.rescale_spinspin,pspace.num_PrintDigits));
    }

    if( pspace.filename_extension != "" )
    {
        filename += "_" + pspace.filename_extension;
    }

    // create file:
    print::cut_if_too_large( filename, m_fname_max_length );
    m_filename += filename;
    if( eigenvalue_threshold_violated )
    {
        m_filename += "_ETV";
    }
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

// main storing function
void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, rtd::RunTimeData& rtdata, CluCorrTen& CCT )
{
    if( !m_storing_permission ){ return; } // permission request

    // ++++++++++ store parameters ++++++++++
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // ========== physical parameters ==========
    // ...concerning the mean-field model
    hdf5r::store_string( ps_group_id, "mf_model",                   pspace.mf_model_name );

    // ...concerning the cluster
    hdf5r::store_string( ps_group_id, "config_file",                pspace.config_file );
    hdf5r::store_scalar( ps_group_id, "num_Spins",                  pspace.num_Spins ); 
    hdf5r::store_string( ps_group_id, "spinspin_cmodel",            pspace.spinspin_cmodel.m_name );
    for( const auto& param : pspace.spinspin_cmodel.m_parameters ) // add all parameters of the spin model
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "spinmf_cmodel",              pspace.spinmf_cmodel.m_name );
    for( const auto& param : pspace.spinmf_cmodel.m_parameters ) // add all parameters of the spin model
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "chemical_shift",             pspace.chemical_shift.m_name );
    hdf5r::store_string( ps_group_id, "local_extra_interaction",    pspace.local_extra_interaction.m_name );
    for( const auto& param : pspace.local_extra_interaction.m_parameters ) // add all parameters of the local extra interaction
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_scalar( ps_group_id, "rescale_meanfield",          pspace.rescale_meanfield );
    hdf5r::store_scalar( ps_group_id, "rescale_spinspin",           pspace.rescale_spinspin );

    // ...concerning the self consistency
    hdf5r::store_2D_tensor<uint>( ps_group_id, "correlation_categories, stored according to the hierarchy n, i-j", H5T_NATIVE_INT, pspace.correlation_categories );

    // ========== numerical parameters ==========
    // ...concerning the symmetry
    hdf5r::store_scalar( ps_group_id, "correlation_symmetry_type",  pspace.symmetry_type );

    // ...concerning the time discretization
    hdf5r::store_scalar( ps_group_id, "num_TimeSteps",              pspace.num_TimeSteps );
    hdf5r::store_scalar( ps_group_id, "num_TimePoints",             pspace.num_TimePoints ); 
    hdf5r::store_scalar( ps_group_id, "delta_t",                    pspace.delta_t );

    // ...concerning the statistics
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerCore",         pspace.num_SamplesPerCore ); 
    hdf5r::store_scalar( ps_group_id, "num_SamplesPerSet",          pspace.num_SamplesPerSet );
    hdf5r::store_scalar( ps_group_id, "num_Samples",                pspace.num_Samples ); 
    hdf5r::store_scalar( ps_group_id, "num_Cores",                  pspace.get_num_Cores() ); 
    hdf5r::store_string( ps_group_id, "adaptive_sample_size",       hdf5r::bool_to_string(pspace.adaptive_sample_size) );
    hdf5r::store_scalar( ps_group_id, "statistical_error_tolerance",pspace.statistical_error_tolerance );
    hdf5r::store_string( ps_group_id, "seed type",                  pspace.seed );

    // ...concerning the initial correlations
    hdf5r::store_string( ps_group_id, "init_corr_mode",             pspace.init_corr_mode );
    // ...if they are imported
    hdf5r::store_string( ps_group_id, "extrapolate_imported_spin_correlations", hdf5r::bool_to_string(pspace.extrapolate_imported_spin_correlations) );
    hdf5r::store_string( ps_group_id, "imported_correlations_src_file",         hdf5r::none_if_empty(pspace.imported_correlations_src_file) );
    // ...if they are generated
    hdf5r::store_string( ps_group_id, "init_diag_corr",             hdf5r::none_if_empty(pspace.init_diag_corr.m_name) );
    for( const auto& param : pspace.init_diag_corr.m_parameters )
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }
    hdf5r::store_string( ps_group_id, "init_nondiag_corr",          hdf5r::none_if_empty(pspace.init_nondiag_corr.m_name) );
    for( const auto& param : pspace.init_nondiag_corr.m_parameters )
    {
        hdf5r::store_scalar( ps_group_id, param.m_name,             param.m_value );
    }

    // ...concerning the iteration
    hdf5r::store_string( ps_group_id, "iteration_error_mode",       pspace.iteration_error_mode );
    hdf5r::store_scalar( ps_group_id, "relative_iteration_error_tolerance" + rtdata.regarded("relative"), pspace.relative_iteration_error_tolerance );
    hdf5r::store_scalar( ps_group_id, "absolute_iteration_error_tolerance" + rtdata.regarded("absolute"), pspace.absolute_iteration_error_tolerance );
    hdf5r::store_scalar( ps_group_id, "Iteration_Limit",            pspace.Iteration_Limit );

    // ...concerning the eigenvalues
    hdf5r::store_string( ps_group_id, "truncation_scheme_negative_eigenvalues", pspace.truncation_scheme_negative_eigenvalues );
    hdf5r::store_scalar( ps_group_id, "eigenvalue_ratio_tolerance (manual)",    pspace.eigenvalue_ratio_tolerance );

    // ...and more
    hdf5r::store_string( ps_group_id, "matrix_exponential_computation",         pspace.matrix_exponential_computation );

    // ========== parameters for storing and naming ==========
    std::string RealType_str{};
    #ifdef USE_DOUBLE
    RealType_str = "DOUBLE";
    #elif defined USE_FLOAT 
    RealType_str = "FLOAT";
    #endif
    hdf5r::store_string( ps_group_id, "RealType",                   RealType_str);
    hdf5r::store_string( ps_group_id, "system_name",                fm::get_system_name());
    hdf5r::store_string( ps_group_id, "information_text",           hdf5r::none_if_empty(pspace.information_text) );
    hdf5r::store_string( ps_group_id, "(original) project_name",    hdf5r::none_if_empty(pspace.project_name) );
    hdf5r::store_string( ps_group_id, "(original) destination project_name",    hdf5r::none_if_empty(pspace.dst_project_name) );
    hdf5r::store_string( ps_group_id, "(original) filename_extension",          hdf5r::none_if_empty(pspace.filename_extension) );
    hdf5r::store_string( ps_group_id, "(original) filename",        m_filename );
    H5Gclose( ps_group_id );

    // ++++++++++ store runtime data ++++++++++
    auto rtd_group_id = H5Gcreate( m_file_id, "runtimedata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // ...concerning the eigenvalues
    hdf5r::store_list( rtd_group_id, "positive_eigen_values_sum_list",          rtdata.positive_eigen_values_sum_list );
    hdf5r::store_list( rtd_group_id, "negative_eigen_values_sum_list",          rtdata.negative_eigen_values_sum_list );
    hdf5r::store_list( rtd_group_id, "largest_positive_eigen_value_list",       rtdata.largest_positive_eigen_value_list );
    hdf5r::store_list( rtd_group_id, "largest_negative_eigen_value_list",       rtdata.largest_negative_eigen_value_list );
    hdf5r::store_list( rtd_group_id, "negative_eigenvalue_ratio_list",          rtdata.negative_eigenvalue_ratio_list );
    std::vector<std::string> ETV_list{};
    std::transform( rtdata.eigenvalue_threshold_violated_list.cbegin(),         rtdata.eigenvalue_threshold_violated_list.cend(), std::back_inserter(ETV_list), hdf5r::bool_to_string );
    hdf5r::store_string_list( rtd_group_id, "eigenvalue_threshold_violated_list", ETV_list );

    // ...concerning the statistics
    hdf5r::store_scalar( rtd_group_id, "generated_seed",                        rtdata.generated_seed );
    store_correlation_tensor_cluster( rtdata.sample_stds, rtd_group_id, "correlation_stds", "Standard deviations of the correlations <S^alpha_i(t)S^beta_j(0)>, stored according to the hierarchy i-j, alpha-beta, t" );
    if( pspace.adaptive_sample_size )
    {
        hdf5r::store_list( rtd_group_id, "adaptive sample size per core (including hypothetical next iteration step)", rtdata.adaptive_num_SamplesPerCore );
    }

    // ...concerning the iterations
    hdf5r::store_scalar( rtd_group_id, "num_Iterations",            rtdata.num_Iterations );
    hdf5r::store_string( rtd_group_id, "termination",               rtdata.termination );
    hdf5r::store_list( rtd_group_id, "relative_iteration_error_list" + rtdata.regarded("relative"), rtdata.relative_iteration_error_list );
    hdf5r::store_list( rtd_group_id, "absolute_iteration_error_list" + rtdata.regarded("absolute"), rtdata.absolute_iteration_error_list );
    H5Gclose( rtd_group_id );

    // ++++++++++ store results ++++++++++
    auto results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    store_correlation_tensor_cluster( CCT, results_group_id, "correlation", "Correlations <S^alpha_i(t)S^beta_j(0)>, stored according to the hierarchy i-j, alpha-beta, t" );
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
    hdf5r::store_scalar( td_group_id, "total_duration",         tmeasure.m_total_duration );
    hdf5r::store_string( td_group_id, "start_date_and_time",    std::string{std::ctime(&tmeasure.m_start_date_and_time)} );
    
    // close group
    H5Gclose( td_group_id );
}

// store correlation tensor cluster
void HDF5_Storage::store_correlation_tensor_cluster( CluCorrTen& data, const hid_t data_group_id, const std::string dataset_name, const std::string dataset_info )
{
    if( !m_storing_permission ){ return; } // permission request

    // create the data space
    auto ct_data_group_id = H5Gcreate( data_group_id, dataset_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // save the data linearized into the dataset
    data.iterate( [&]( auto& CT, auto& sites )
    {
        hdf5r::store_2D_tensor<RealType>(ct_data_group_id, std::to_string(sites[0]+1)+"-"+std::to_string(sites[1]+1), H5_REAL_TYPE, CT);
    } );
    std::string info = dataset_info;
    hdf5r::store_string( ct_data_group_id, "info", info );

    // finalize
    H5Gclose( ct_data_group_id );
}

// close the file
void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}

};