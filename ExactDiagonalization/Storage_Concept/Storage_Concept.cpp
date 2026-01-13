#include"Storage_Concept.h"

#include<fstream>
#include<Standard_Algorithms/Print_Routines.h>
//#include<sys/types.h>
#include<sys/stat.h>
#include"STOC_Error_Handling.h"

namespace SpinED::Storage_Concept
{

namespace print = Print_Routines;
namespace error = Error_Handling;
namespace hdf5r = HDF5_Routines;

// ================= CLASS IMPLEMENTATIONS =================
// HDF5 STORAGE CLASS
// constructor : create folder tree and data file
HDF5_Storage::HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace ):
    m_storing_permission( my_rank == 0 ),
    m_fname_max_length( 200 ),
    m_num_TriesToBuildFile( 5 )
{
    auto folder_branch_list = pspace.create_folder_branch_list();
    create_folder_tree( folder_branch_list );
    create_file( folder_branch_list.back() );
}


void HDF5_Storage::store_main( const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata, const mod::MeanTen& mean, const mod::CorrTen& corr )
{
    if( !m_storing_permission ){ return; } // permission request

    // ===== store parameters =====
    auto ps_group_id = H5Gcreate( m_file_id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    hdf5r::add_uint_to_group( ps_group_id, "num_Spins",                    pspace.num_Spins ); 
    hdf5r::add_uint_to_group( ps_group_id, "num_HilbertSpaceDimension",    pspace.num_HilbertSpaceDimension ); 
    std::string tmp{ pspace.couplings_filename + ".hdf5" };
    hdf5r::add_string_to_group( ps_group_id, "src_file",                   tmp );
    hdf5r::add_string_to_group( ps_group_id, "spin_model",                 pspace.spin_model->m_name );
    for( const auto& param : pspace.spin_model->return_params() ) // add all parameters of the spin model
    {
        hdf5r::add_RealType_to_group( ps_group_id, param.m_name,           param.m_value );
    }
    hdf5r::add_RealType_to_group( ps_group_id, "rescale", pspace.rescale );

    //hdf5r::add_char_to_group( ps_group_id, "mean_symmetry_type",           pspace.spin_model.m_means._symmetry_type );
    //hdf5r::add_char_to_group( ps_group_id, "correlation_symmetry_type",    pspace.spin_model.m_correlations._symmetry_type );

    hdf5r::add_uint_to_group( ps_group_id, "num_TimeSteps",                pspace.num_TimeSteps ); 
    hdf5r::add_uint_to_group( ps_group_id, "num_TimePoints",               pspace.num_TimePoints ); 
    hdf5r::add_RealType_to_group( ps_group_id, "delta_t",                  pspace.delta_t );

    hdf5r::add_string_to_group( ps_group_id, "information_text",            pspace.information_text );
    hdf5r::add_string_to_group( ps_group_id, "original project_name",       pspace.project_name );
    hdf5r::add_string_to_group( ps_group_id, "original filename_extension", pspace.filename_extension );

    H5Gclose( ps_group_id );

    // ===== store run time data =====
    auto rtd_group_id = H5Gcreate( m_file_id, "runtimedata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Gclose( rtd_group_id );

    // ===== store results =====
    m_results_group_id = H5Gcreate( m_file_id, "results", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    
    store_mean( mean );
    store_corr( corr );

    H5Gclose( m_results_group_id );
}


// store all measured times
void HDF5_Storage::store_time( const tmm::TimeMeasure& tmeasure )
{
    if( !m_storing_permission ){ return; } // permission request

    // ===== store times =====
    auto td_group_id = H5Gcreate( m_file_id, "timedata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // global measures
    std::for_each( tmeasure.m_global_durations.cbegin(), tmeasure.m_global_durations.cend(), 
    [&td_group_id]( const tmm::DurationQuantity& q )
    {
        hdf5r::add_double_to_group( td_group_id, q.m_name, q.m_duration );
    } );
    
    // total duration and start date and time
    hdf5r::add_double_to_group( td_group_id, "total_duration", tmeasure.m_total_duration );
    hdf5r::add_string_to_group( td_group_id, "start_date_and_time", std::string{std::ctime(&tmeasure.m_start_date_and_time)} );

    // close group
    H5Gclose( td_group_id );
}

void HDF5_Storage::finalize()
{
    if( !m_storing_permission ){ return; } // permission request
    H5Fclose( m_file_id );
}


// create the whole branch of folders and correspondingly extends the filename
void HDF5_Storage::create_folder_tree( const std::vector<std::string>& folder_branch_list )
{
    if( !m_storing_permission ){ return; } // permission request

    m_filename = "";
    std::for_each( folder_branch_list.cbegin(), folder_branch_list.cend()-1, 
    [this]( std::string s )
    {
        print::cut_if_too_large( s, m_fname_max_length );
        m_filename += s;
        char const *c = m_filename.data();
        if( mkdir( c, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) != -1 ) // then the folder was built (otherwise the folder was already there)
        {
            std::cout << "created folder : " << m_filename << "\n";
        }
        m_filename += "/";
    } );
}

// create the file
void HDF5_Storage::create_file( std::string filename )
{
    if( !m_storing_permission ){ return; } // permission request

    print::cut_if_too_large( filename, m_fname_max_length );
    m_filename += filename;
    uint count = 0;
    do{
        std::string tmp = m_filename + ".hdf5";
        std::ifstream f( tmp.c_str() );
        if( !f.good() ) // then the file doesn't exist yet
        {          
            m_file_id= H5Fcreate( tmp.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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


// store mean tensor
void HDF5_Storage::store_mean( const mod::MeanTen data )
{   
    if( !m_storing_permission ){ return; } // permission request

    // linearize
    std::vector<RealType> data_linearized{};
    std::for_each( data.cbegin(), data.cend(), [&data_linearized]( const auto& mean_component )
    {
        std::for_each( mean_component.cbegin(), mean_component.cend(), [&data_linearized]( const RealType& mean )
        {
            data_linearized.emplace_back( mean );
        } );
    } );

    // store
    std::vector<hsize_t> tensor_sizes{ data.size(), data.cbegin()->size()  };
    size_t tensor_rank = 2;
    auto dataspace = H5Screate_simple( tensor_rank, tensor_sizes.data(), NULL );
    auto dataset = H5Dcreate( m_results_group_id, "mean", H5_REAL_TYPE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( dataset, H5_REAL_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_linearized.data() );
    
    // add info
    std::string info = "Mean Tensor <S_1,alpha(t)>, stored according to the hierarchy alpha,t";
    hdf5r::add_string_to_group( dataset, "info", info );

    // finalize
    H5Dclose( dataset );
    H5Sclose( dataspace );
}

// store correlation tensor
void HDF5_Storage::store_corr( const mod::CorrTen data )
{   
    if( !m_storing_permission ){ return; } // permission request

    // linearize
    std::vector<RealType> data_linearized{};
    std::for_each( data.cbegin(), data.cend(), [&data_linearized]( const auto& corr_component )
    {
        std::for_each( corr_component.cbegin(), corr_component.cend(), [&data_linearized]( const RealType& corr )
        {
            data_linearized.emplace_back( corr );
        } );
    } );

    // store
    std::vector<hsize_t> tensor_sizes{ data.size(), data.cbegin()->size() };
    size_t tensor_rank = 2;
    auto dataspace = H5Screate_simple( tensor_rank, tensor_sizes.data(), NULL );
    auto dataset = H5Dcreate( m_results_group_id, "correlation", H5_REAL_TYPE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( dataset, H5_REAL_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_linearized.data() );

    // add info
    std::string info = "Correlation Tensor <S_1,alpha(t_1)S_1,beta(t_2)>, stored according to the hierarchy alphabeta,t1,t2";
    hdf5r::add_string_to_group( dataset, "info", info );

    // finalize
    H5Dclose( dataset );
    H5Sclose( dataspace );
}


};
