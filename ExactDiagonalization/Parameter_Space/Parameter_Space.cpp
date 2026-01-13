#include"Parameter_Space.h"

#include<iomanip>
#include<memory>
#include<sstream>
#include<fstream>
#include<cmath>
#include<boost/program_options.hpp>
#include<Standard_Algorithms/Print_Routines.h>
#include"../Global/MPI_Global.h"
#include"PS_Error_Handling.h"
#include"../Models/Models.h"

namespace SpinED::Parameter_Space
{

namespace bpo = boost::program_options;
namespace print = Print_Routines;
namespace error = Error_Handling;
namespace mod = SpinED::Models;

// ================= CLASS IMPLEMENTATIONS =================
// PARAMETER SPACE CLASS
// constructor : build parameter space from command line arguments
ParameterSpace::ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank ):
    my_rank( my_rank ),
    world_size( world_size )
{
    // 1a.) Define Options
    bpo::options_description description("Allowed options:");
    description.add_options()
    (
    "help", "produce help message"
    )(

    // ========== model and physical parameters ==========
    "srcfile", bpo::value<std::string>()->default_value("-"),
    "specify the file from which the couplings should be read; leave out file-ending hdf5"
    )(    
    "spinmodel", bpo::value<std::string>()->default_value("ISO_Polarized"),
    "set the spin model || options are combinations of: \
    ISO = Isotropic Heisenberg Model, \
    XXZ = Easy-Axis anisotropic Heisenberg Model \
    and \
    Disordered (infinite temperature), \
    Polarized \
    and \
    Blockwise"
    )(
    "rescale", bpo::value<RealType>()->default_value(RealType{1.0}),
    "the factor rescale is multiplied to the couplings (for example for them to be in units of JQ)"
    )(
    "lambda", bpo::value<RealType>()->default_value(RealType{0.5}),
    "set the anisotropy factor for the XXZ-model; the factor will be multiplied to z-couplings"
    )(
    "h", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the value of the polarization field if required"
    )(

    // ========== general numerical parameters ==========
    // ...concerning time discretization
    "numTimeSteps", bpo::value<uint>()->default_value(uint{5}),
    "set the number of time steps for the equidistant time discretization"
    )(
    "dt", bpo::value<RealType>()->default_value(RealType{0.1}),
    "set the step width for the equidistant time discretization"
    )(

    // ========== storing and naming ==========
    "info", bpo::value<std::string>()->default_value("-"),
    "provide some text that informs about the procedure and purpose of the computation; \
    to be able to use spaces, put the text in quotation marks"
    )(
    "project", bpo::value<std::string>()->default_value(""),
    "sort the data into a project for better order"
    )(
    "fileext", bpo::value<std::string>()->default_value(""),
    "Define an extension to the filename; it will be appended according to : filename_fileext"
    )(
    "numPrintDigits", bpo::value<uint>()->default_value(uint{4}),
    "set the value precision for printing to the terminal"
    );
    // "firstn", bpo::value<std::string>()->default_value("all"),
    // "number of spins of which the correlations and pair correlations are saved (1, 2, ..., all)"
    // )(

    // 1b.) Store command line arguments in bp options
    bpo::variables_map vm;
    bpo::store( bpo::parse_command_line(argC, argV, description), vm );
    bpo::notify( vm );

    // 2.) Output of help description and termination of the program
    if(vm.count("help"))
    {
        if(my_rank==0){
            std::cout << description << "\n";
        }
        MPI_Finalize();
        exit(0);
    }
    // 3.) Store bp options in parameter space 
    // ========== physical parameters ==========
    couplings_filename      = vm["srcfile"].as<std::string>();
    rescale = vm["rescale"].as<RealType>();
    read_SpinSystem(); // couplings and num_Spins are set
    num_HilbertSpaceDimension = std::pow( 2, num_Spins );
    model_params.lambda = vm["lambda"].as<RealType>();
    model_params.h_z = vm["h"].as<RealType>();

    // ========== model ==========
    // set up model map:
    std::map< std::string, std::shared_ptr< mod::Model > > model_map // maps a model name to a shared ptr pointing to a model
    {
        { "ISO_Disordered", std::make_shared<mod::ISO_Disordered_Model>(*this) },
        { "ISO_Disordered_Blockwise", std::make_shared<mod::ISO_Disordered_Blockwise_Model>(*this) },
        { "ISO_Polarized", std::make_shared<mod::ISO_Polarized_Model>(*this) },
        { "XXZ_Disordered", std::make_shared<mod::XXZ_Disordered_Model>(*this) },
        { "DRF_Disordered", std::make_shared<mod::DRF_Disordered_Model>(*this) },
        //{ "DRF_Disordered_Blockwise", std::make_shared<mod::DRF_Disordered_Blockwise_Model>(*this) },
    }; 
    spin_model = model_map.at( vm["spinmodel"].as<std::string>() ); // determine model from bpo

    // ========== general numerical parameters ==========
    // ...concerning time discretization
    num_TimeSteps  = vm["numTimeSteps"].as<uint>();
    num_TimePoints = num_TimeSteps + 1;
    delta_t        = vm["dt"].as<RealType>();
    // RealType Tmax  = delta_t * static_cast<RealType>(num_TimeSteps);

    // ========== saving and naming ==========
    information_text        = vm["info"].as<std::string>();
    project_name            = vm["project"].as<std::string>();
    filename_extension      = vm["fileext"].as<std::string>();
    num_PrintDigits         = vm["numPrintDigits"].as<uint>();
    // save_first_n_spins      = vm["firstn"].as<std::string>();
}


// printing method : return the essential parameters string
std::string ParameterSpace::create_essentials_string() const
{   
    size_t pre_colon_space = 35;
    std::stringstream ss{};
    ss
    << print::quantity_to_output_line( pre_colon_space, "num_Spins", std::to_string( num_Spins ) )
    << print::quantity_to_output_line( pre_colon_space, "couplings_filename", couplings_filename )
    << print::quantity_to_output_line( pre_colon_space, "num_TimeSteps" , std::to_string(num_TimeSteps) ) 
    << print::quantity_to_output_line( pre_colon_space, "delta_t"       , print::round_value_to_string(delta_t,num_PrintDigits) ) 
    << print::quantity_to_output_line( pre_colon_space, "rescale"       , print::round_value_to_string(rescale,num_PrintDigits) ) 
    << print::quantity_to_output_line( pre_colon_space, "information_text" , information_text );
    return ss.str();
}


// save data method : return the folder branch including the filename where the data will be saved
std::vector<std::string> ParameterSpace::create_folder_branch_list() const
{
    std::vector<std::string> fbranch{};
    fbranch.push_back( "Data" );
    if( project_name != "" )
    {
        fbranch.push_back( project_name );
    }
    fbranch.push_back( this->create_file_name( true ) );
    return fbranch;
}

// save data method : return the name of the file where the data will be saved
// idea : file name is created from the physical parameters
std::string ParameterSpace::create_file_name( const bool with_extension ) const
{
    std::string file_name{};
    file_name += couplings_filename + "__" + spin_model->compact_info( num_PrintDigits ) ; // spin model info
    if( rescale != RealType{1.} )
    {
        file_name += "__rescale=" + print::remove_zeros(print::round_value_to_string(rescale,num_PrintDigits)); // rescaling
    }
    if( filename_extension != "" && with_extension )
    {
        file_name += "_" + filename_extension;
    }
    return file_name;
}

/* read spin system data from hdf5 file : couplings and number of spins
this function should only be used, after the couplings_filename has been set */
void ParameterSpace::read_SpinSystem()
{
    // 1) open file and group:
    std::string total_filename = "Couplings/" + couplings_filename + ".hdf5";
    hid_t file_id = H5Fopen( total_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    hid_t group_id = H5Gopen( file_id, "/all", H5P_DEFAULT );

    // 2) read the number of spins:
    // a) open attribute and read type
    hid_t attr_id = H5Aopen( group_id, "num_Spins", H5P_DEFAULT );
    hid_t datatype_id = H5Aget_type( attr_id ); // datatype of the dataset -> uint

    // b) read attribute data
    H5Aread( attr_id, datatype_id, &num_Spins );

    // c) close attribute
    H5Tclose( datatype_id );
    H5Aclose( attr_id );

    // 3) read the coupling data:
    // a) open dataset and read type
    hid_t dataset = H5Dopen2( group_id, "/all/J_ij", H5P_DEFAULT );
    hid_t datatype = H5Dget_type( dataset ); // datatype of the dataset -> double

    // b) read linearized data from the dataset
    std::vector<double> linearized_data( num_Spins * num_Spins );
    H5Dread( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, linearized_data.data() );

    // c) close resources
    H5Tclose( datatype );
    H5Dclose( dataset );

    // 4) close resources
    H5Gclose( group_id );
    H5Fclose( file_id );

    // 5) tensorize the data:
    J.resize( num_Spins );
    for( size_t i = 0; i < num_Spins; i++ ) 
    {
        for( size_t j = 0; j < num_Spins; j++ ) 
        {
            J(i,j) = rescale * linearized_data[i*num_Spins + j];
        }
    }
}

};

