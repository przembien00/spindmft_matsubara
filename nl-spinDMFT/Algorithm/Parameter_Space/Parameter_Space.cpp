#include"Parameter_Space.h"
#include<map>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<numeric>

#include<Globals/MPI_Types.h>
#include"boost/program_options.hpp"
namespace bpo = boost::program_options;

#include"Error_Handling.h"
namespace error = DMFT_parameter_space::Error_Handling;

#include"Standard_Algorithms/Print_Routines.h"
namespace print = Print_Routines;

#include"../Mean_Field_Models/Mean_Field_Models.h"
namespace mfm = Mean_Field_Models;

#include<HDF5/HDF5_Routines.h>
namespace hdf5r = HDF5_Routines;

namespace DMFT_parameter_space
{

// Forward declaration of helper functions
IndexPairList create_whole_categories_list( const size_t& num_Spins );
std::string categories_string( const IndexPairList& ipairs );

// ============================================================================
// ========================== PARAMETER SPACE CLASS ===========================
// ============================================================================
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

    // ========== physical parameters ==========
    // ...concerning the mean-field model
    "mfmodel", bpo::value<std::string>()->default_value("CorrelationReplica"),
    "set the mean-field model || only CorrelationReplica available"
    )(

    // ...concerning the cluster
    "config", bpo::value<std::string>()->default_value(""),
    "choose the file containing the configuration data; mostly of the form <System>_csize=<Clustersize>_<Rest>"
    )(
    "spinspinmodel", bpo::value<std::string>()->default_value("ISO"),
    "set the spin coupling model || options are: ISO, Ising, XXZ, DRF, XY" 
    )(
    "lambda", bpo::value<RealType>()->default_value(RealType{2.0}),
    "set the anisotropy factor for the XXZ-model or XY-model; the factor will be multiplied to the z-couplings of a Heisenbergmodel"
    )(
    "rho", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the second anisotropy factor for the XY-model; the factor will be multiplied to the xy-couplings"
    )(
    "spinmfmodel", bpo::value<std::string>()->default_value("auto"),
    "set the model for the spin-mean-field couplings || options are: ISO, Ising, XXZ, DRF, XY" 
    )(
    "smflambda", bpo::value<RealType>()->default_value(RealType{2.0}),
    "set the anisotropy factor for the XXZ-model or XY-model; the factor will be multiplied to the z-couplings of a Heisenbergmodel"
    )(
    "smfrho", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the second anisotropy factor for the XY-model; the factor will be multiplied to the xy-couplings"
    )(
    "chemshift", bpo::value<std::string>()->default_value("none"),
    "set the chemical shift (magnetic field in z-direction) || none or comma-separated values for all spins" 
    )(
    "extraint", bpo::value<std::string>()->default_value("none"),
    "include a local extra interaction term acting at each spin site || options are: none, quadrupolar"
    )(
    "intstrength", bpo::value<RealType>()->default_value(RealType{0.0}),
    "set the strength of the extra interaction term"
    )(
    "rescmf", bpo::value<RealType>()->default_value(1.0),
    "set the rescaling factor for spin-mean-field couplings; the factor is multiplied to them" 
    )(
    "rescss", bpo::value<RealType>()->default_value(1.0),
    "set the rescaling factor for spin-spin couplings; the factor is multiplied to them" 
    )(

    // ...concerning the simulation output
    "computeonly", bpo::value<std::string>()->default_value(""),
    "provide a choice of correlations that shall be computed as comma-separated list (spin numbering starts at 1) \
    for example: 1-1,2-3,4-5 \
    if empty, all correlations will be computed"
    )(
    
    // ========== numerical parameters ==========
    // ...concerning the symmetry
    "stype", bpo::value<char>()->default_value('B'),
    "set the symmetry type || options are: \
    A = (gab=0, gxx=gyy=gzz), \
    B = (gab=0, gxx=gyy), \
    C = (gaz=0, gxx=gyy), \
    D = (no constraints)"
    )(

    // ...concerning the time discretization
    "numTimeSteps", bpo::value<size_t>()->default_value(size_t{5}),
    "set the number of time steps for the equidistant time discretization"
    )(
    "dt", bpo::value<RealType>()->default_value(RealType{0.1}),
    "set the step width for the equidistant time discretization"
    )(

    // ...concerning the statistics
    "numSamplesPerCore", bpo::value<size_t>()->default_value(size_t{100}),
    "set the number of samples per core (!) for the Monte Carlo simulation"
    )(
    "numSamplesPerSet", bpo::value<size_t>()->default_value(size_t{1}),
    "sampling the noise in sets increases the efficiency, numSamplesPerSet has to be smaller than numSamplesPerCore"
    )(
    "seed", bpo::value<std::string>()->default_value("random"),
    "set the seed for the random generator (random : seed is determined by clock)"
    )(

    // ...concerning the initial correlations 
    "initmode", bpo::value<std::string>()->default_value("import"),
    "import = import the initial spin correlations from generated spin correlations \
    generate = generate the initial spin correlations from defined analytic functions"
    )(
    // ...if they are imported
    "extrapolate", "extrapolate the imported initial spin correlations linearly to the new discretization"
    )(
    "impcorrfile", bpo::value<std::string>()->default_value(""),
    "provide the filename with path in the Data folder from which the initial correlations should be taken"
    )(
    "impcorrsrc", bpo::value<std::string>()->default_value("spinDMFT"),
    "provide the root folder of the file from which the initial correlations should be taken (spinDMFT, CspinDMFT, p-spinDMFT)"
    )(
    // ...if they are generated
    "initdcorr", bpo::value<std::string>()->default_value("exponential"),
    "set the initial diagonal correlations (gxx=gyy=gzz)"
    )(
    "initndcorr", bpo::value<std::string>()->default_value("zero"),
    "set the initial nondiagonal correlations (gab with a!=b)"
    )(
    //"corrtscale", bpo::value<RealType>()->default_value(RealType{1.0}),
    //"set the time scale of initial correlations"
    //)(
    "corrperiods", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the number of periods in case of periodic initial correlations"
    )(

    // ...concerning eigenvalues
    "eigtolerance", bpo::value<RealType>()->default_value(RealType{0.001}),
    "manually set to which ratio negative eigenvalues are tolerated ( the automatic threshold is determined from the sample size ); \
    within the algorithm always the softer one of these two thresholds is considered"
    )(
    "truncneg", bpo::value<std::string>()->default_value("set_zero"),
    "decide for the truncation scheme for negative eigenvalues || options are: set_zero, abs"
    )(

    // ...and more
    "matrixcomp", bpo::value<std::string>()->default_value("diagonalization"),
    "decide on how to compute matrix exponentials || options are: taylor (very inaccurate), diagonalization"
    )(

    // ========== parameters for storing and naming ==========
    "info", bpo::value<std::string>()->default_value(""),
    "provide some text that informs about the procedure and purpose of the computation; \
    to be able to use spaces, put the text in quotation marks"
    )(
    "project", bpo::value<std::string>()->default_value(""),
    "data are sorted in a project for better order"
    )(
    "dstproject", bpo::value<std::string>()->default_value(""),
    "sort the resulting data into a project for better order, dst project is set to project if empty"
    )(
    "fileext", bpo::value<std::string>()->default_value(""),
    "Define an extension to the filename; it will be appended according to : filename_fileext"
    );

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

    // 3.) Store bp options in parameter space (note that several parameters are set later)
    // ========== physical parameters ==========
    // ...concerning the mean-field model
    mf_model_name              = vm["mfmodel"].as<std::string>();
    if( mf_model_name != "CorrelationReplica" ){ exit(0); } // first implement new model and change the folder branch construction

    // ...concerning the cluster
    config_file             = vm["config"].as<std::string>();
    if( config_file == "" ){ error::CONFIG_NOT_SET( __PRETTY_FUNCTION__ ); } // config file must be set!
    spinspin_cmodel = ph::SpinModel{ vm["spinspinmodel"].as<std::string>(), vm["lambda"].as<RealType>(), vm["rho"].as<RealType>() };
    std::string smf = vm["spinmfmodel"].as<std::string>();
    spinmf_cmodel = (smf == "auto") ? spinspin_cmodel : ph::SpinModel{ smf, vm["smflambda"].as<RealType>(), vm["smfrho"].as<RealType>() };
    chemical_shift = ph::ChemicalShift{ vm["chemshift"].as<std::string>() };
    rescale_meanfield = vm["rescmf"].as<RealType>(); // the rescaling happens later in the models
    rescale_spinspin  = vm["rescss"].as<RealType>(); // the rescaling happens later in the models

    // ========== numerical parameters ==========
    // ...concerning symmetry
    symmetry_type           = vm["stype"].as<char>();

    // ...concerning time discretization
    delta_t                 = vm["dt"].as<RealType>();
    num_TimeSteps           = vm["numTimeSteps"].as<size_t>();
    num_TimePoints          = num_TimeSteps+1;

    // ...concerning statistics
    num_SamplesPerCore      = vm["numSamplesPerCore"].as<size_t>();
    num_SamplesPerSet       = vm["numSamplesPerSet"].as<size_t>();
    num_SetsPerCore         = (size_t) std::ceil( static_cast<double>(num_SamplesPerCore) / static_cast<double>(num_SamplesPerSet) ); 
    num_Samples             = world_size * num_SamplesPerCore;
    seed                    = vm["seed"].as<std::string>();

    // ...concerning the initial correlations 
    init_corr_mode = vm["initmode"].as<std::string>();
    // ...if they are generated
    RealType Tmax           = delta_t * static_cast<RealType>(num_TimeSteps);
    init_diag_corr    = ph::DiagonalSpinCorrelation{ vm["initdcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    init_nondiag_corr = ph::NonDiagonalSpinCorrelation{ vm["initndcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    // rest is set below
   
    // ...concerning eigenvalues
    eigenvalue_ratio_tolerance = vm["eigtolerance"].as<RealType>();
    truncation_scheme_negative_eigenvalues = vm["truncneg"].as<std::string>();

    // ...and more
    matrix_exponential_computation = vm["matrixcomp"].as<std::string>();

    // ========== parameters for storing and naming ==========
    information_text        = vm["info"].as<std::string>();
    project_name            = vm["project"].as<std::string>();
    std::string dst_project_name_tmp = vm["dstproject"].as<std::string>();
    dst_project_name = (dst_project_name_tmp == "") ? project_name : dst_project_name_tmp;
    filename_extension      = vm["fileext"].as<std::string>();

    // ========== parameters that could not be set before ==========
    // determine initial correlations src folder if necessary
    if( init_corr_mode == "import" )
    { 
        if( vm.count("extrapolate") )
        {
            extrapolate_imported_spin_correlations = true;
        }
        imported_correlations_src_file = vm["impcorrfile"].as<std::string>(); 
        imported_correlations_src_directory = vm["impcorrsrc"].as<std::string>(); 
        init_diag_corr.m_name = "";
        init_nondiag_corr.m_name = "";
    }

    // construct and initialize the mean field model and hand back general physical parameters from it
    std::map< std::string, std::shared_ptr<mfm::MeanFieldModel> > model_map // maps a model name to a shared ptr pointing to a model
    { 
        { "CorrelationReplica", std::make_shared<mfm::CorrelationReplicaModel>() },
    };
    mf_model = model_map.at( mf_model_name );
    mf_model->constructor_function( *this );

    // interpret computeonly
    if( compute_only_categories.size()==0 ) // otherwise, this parameter has been set by config file
    {
        std::string compute_only = vm["computeonly"].as<std::string>();
        if( compute_only == "" ) // compute all correlations
        {
            compute_only_categories = create_whole_categories_list(num_Spins);
        }
    }

    // determine local extra interaction
    local_extra_interaction = ph::LocalExtraInteraction{ vm["extraint"].as<std::string>(), vm["intstrength"].as<RealType>(), num_Spins };

    // import initial correlations according to given scheme
    if( init_corr_mode == "import" )
    {   
        read_initial_correlations_from_file();
    }
}

// return the essential parameters string
std::string ParameterSpace::create_essentials_string() const
{
    size_t pre_colon_space = 20;
    std::stringstream ss{};
    ss << "I compute the following correlations (numbering starts at 1): " << categories_string(compute_only_categories) << "\n";
    if( init_corr_mode == "import" )
    {
        ss << "I import the spin correlations " << categories_string(import_only_categories) << " from \'" << total_import_filename << "\'\n";
    }
    ss
    << print::quantity_to_output_line( pre_colon_space, "mf_model"      , mf_model_name )
    << print::quantity_to_output_line( pre_colon_space, "symmetry_type" , std::string(1, symmetry_type) )
    << print::quantity_to_output_line( pre_colon_space, "num_TimeSteps" , std::to_string(num_TimeSteps) ) 
    << print::quantity_to_output_line( pre_colon_space, "delta_t"       , print::round_value_to_string(delta_t,num_PrintDigits) ) 
    << print::quantity_to_output_line( pre_colon_space, "num_Samples"   , std::to_string(num_Samples) )
    << print::quantity_to_output_line( pre_colon_space, "config_file" , config_file )
    << print::quantity_to_output_line( pre_colon_space, "information_text", information_text ) 
    << print::quantity_to_output_line( pre_colon_space, "num_Spins"     , std::to_string(num_Spins) )
    << print::quantity_to_output_line( pre_colon_space, "spin_list"     , print::concatenate_string_with_delimiter(spin_list, ',') );
    if( num_Spins < 5 )
    {
        ss << print::quantity_to_output_line( pre_colon_space, "J"      , std::string{std::to_string(num_Spins) + "x" + std::to_string(num_Spins) + "-Matrix"} );
        ss << J;
    }
    else 
    {
        ss << print::quantity_to_output_line( pre_colon_space, "J(0,1)" , print::round_value_to_string(J(0,1),num_PrintDigits) );
    }
    ss 
    << print::quantity_to_output_line( pre_colon_space, "spinspin_cmodel info"   , spinspin_cmodel.compact_info(num_PrintDigits) )
    << print::quantity_to_output_line( pre_colon_space, "spinmf_cmodel info" , spinmf_cmodel.compact_info(num_PrintDigits) )
    << print::quantity_to_output_line( pre_colon_space, "chemical_shift info" , chemical_shift.compact_info(num_PrintDigits) )
    << print::quantity_to_output_line( pre_colon_space, "rescale_meanfield", print::round_value_to_string(rescale_meanfield,num_PrintDigits) )
    << print::quantity_to_output_line( pre_colon_space, "rescale_spinspin",  print::round_value_to_string(rescale_spinspin,num_PrintDigits) );
    return ss.str();
}

// read the initial correlations in semi-linearized form from file
void ParameterSpace::read_initial_correlations_from_file()
{ 
    // 1) open file and group:
    total_import_filename = "../" + imported_correlations_src_directory + "/Data/" + imported_correlations_src_file + ".hdf5";
    hid_t file_id = H5Fopen( total_import_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( file_id < 0 )
    {
        error::IMPORT_FILE_NOT_FOUND( __PRETTY_FUNCTION__, total_import_filename );
    }
    hid_t group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );

    // 2) import parameters of loaded data and check if they match with the new ones:
    char old_correlation_symmetry_type{};
    std::string old_RealType{};
    hdf5r::import_scalar( group_id, "num_TimePoints", old_num_TimePoints );
    hdf5r::import_scalar( group_id, "delta_t", old_delta_t );
    hdf5r::import_scalar( group_id, "correlation_symmetry_type", old_correlation_symmetry_type );
    hdf5r::import_string( group_id, "RealType", old_RealType );

    // check discretization (may only mismatch if extrapolation is used)
    if( !extrapolate_imported_spin_correlations )
    {
        if( old_num_TimePoints != num_TimePoints || old_delta_t != delta_t )
        {
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        }
    }

    // check symmetry and type (has to match in any case)
    if( old_correlation_symmetry_type != symmetry_type )
    {
        error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
    }
    else // check RealType (also has to match in any case)
    {
        #ifdef USE_DOUBLE
        if( old_RealType != "DOUBLE" )
        {
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        }
        #endif
        #ifdef USE_FLOAT 
        if( old_RealType != "FLOAT" )
        {
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        }
        #endif
    }

    // 3) read the correlation data:
    if( imported_correlations_src_directory == "spinDMFT" )
    {
        group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );
        imported_correlations_linearized.resize(1);
        hdf5r::import_ND_tensor_linearized( group_id, "/results/correlation", imported_correlations_linearized[0] );
        std::for_each( imported_correlations_linearized[0].cbegin(), imported_correlations_linearized[0].cend(), []( const RealType& v )
        {
            if( std::isnan( v ) ) // check data for nans
            {
                error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ );
            }
        } );
        H5Gclose( group_id );
    }
    else if( imported_correlations_src_directory == "CspinDMFT" )
    {
        // a) read num_Spins
        group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );
        hdf5r::import_scalar( group_id, "num_Spins", old_num_Spins );
        H5Gclose( group_id );

        // b) open group and dataset and read type
        group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );
        hid_t sub_group_id = H5Gopen( file_id, "/results/correlation", H5P_DEFAULT );

        // c) read linearized data for each (desired) correlation from the dataset
        // loop over spins of the old data set
        std::for_each( import_only_categories.cbegin(), import_only_categories.cend(), [&](const IndexPair& ipair)
        {
            std::string dataset_name = "/results/correlation/" + std::to_string(ipair[0]) + "-" + std::to_string(ipair[1]);
            std::vector<RealType> linearized_correlation_tensor{}; // linearized correlation tensor for single pair of spins (i,j)
            hdf5r::import_ND_tensor_linearized( sub_group_id, dataset_name, linearized_correlation_tensor );
            std::for_each( linearized_correlation_tensor.cbegin(), linearized_correlation_tensor.cend(), []( const RealType& v )
            {
                if( std::isnan( v ) ) // check data for nans
                {
                    error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ );
                }
            } );
            imported_correlations_linearized.emplace_back( std::move(linearized_correlation_tensor) );
        } );

        // d) close resources
        H5Gclose( sub_group_id );
        H5Gclose( group_id );
    }
    else if( imported_correlations_src_directory == "p-spinDMFT" )
    {
        // a) read num_Sites
        group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );
        hdf5r::import_scalar( group_id, "num_Sites", old_num_Spins );
        H5Gclose( group_id );

        // b) open group and dataset and read type
        group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );
        hid_t sub_group_id = H5Gopen( file_id, "/results/correlation", H5P_DEFAULT );

        // c) read linearized data for each (desired) correlation from the dataset
        // loop over spins of the old data set, note that i has to be equal to j for p-spinDMFT
        std::for_each( import_only_categories.cbegin(), import_only_categories.cend(), [&](const IndexPair& ipair)
        {
            std::string dataset_name = "/results/correlation/" + std::to_string(ipair[0]) + "-" + std::to_string(ipair[0]);
            std::vector<RealType> linearized_correlation_tensor{}; // linearized correlation tensor for single autocorrelation (i,i)
            hdf5r::import_ND_tensor_linearized( sub_group_id, dataset_name, linearized_correlation_tensor );
            std::for_each( linearized_correlation_tensor.cbegin(), linearized_correlation_tensor.cend(), []( const RealType& v )
            {
                if( std::isnan( v ) ) // check data for nans
                {
                    error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ );
                }
            } );
            imported_correlations_linearized.emplace_back( std::move(linearized_correlation_tensor) );
        } );

        // d) close resources
        H5Gclose( sub_group_id );
        H5Gclose( group_id );
    }
    else
    {
        error::INIT_CORRELATIONS_SRC_FOLDER( __PRETTY_FUNCTION__, imported_correlations_src_directory );
    }

    // 4) close resources
    H5Fclose( file_id );
}

// ============================================================================
// ============================= HELPER FUNCTIONS =============================
// ============================================================================
// create a list of all possible correlations based on num_Spins
IndexPairList create_whole_categories_list( const size_t& num_Spins )
{
    IndexPairList ipairs{};
    for( size_t i=0; i<num_Spins; ++i )
    {
        for( size_t j=i; j<num_Spins; ++j )
        {
            ipairs.emplace_back( IndexPair{i,j} );
        }
    }
    return ipairs;
}

// transform IndexPairList to string
std::string categories_string( const IndexPairList& ipairs )
{
    std::string str{};
    for( auto ipair : ipairs ){ str += std::to_string(ipair[0]+1) + "-" + std::to_string(ipair[1]+1) + ","; }
    str.pop_back(); // remove last comma
    return str;
}

}