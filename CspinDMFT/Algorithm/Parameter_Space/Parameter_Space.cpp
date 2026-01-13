#include"Parameter_Space.h"
#include<map>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<numeric>
#include<random>
#include<algorithm>

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
    "runID", bpo::value<std::string>()->default_value("random"),
    "via the run ID parallel jobs can be distinguished from each other \
    simulations can be terminated after their current iteration step by creating the file terminate_<runID> in the main directory \
    the code does NOT check whether an inserted runID aready exists (user's responsibility)"
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
    "include a local extra interaction term acting at each spin sites || options are: none, quadrupolar"
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
    "adaptive", "adaptively increase the sample size in each iteration step to fulfill staterrtolerance"
    )(
    "staterrtolerance", bpo::value<RealType>()->default_value(RealType{0.01}),
    "set the absolute error tolerance for the statistical error, which is employed in case of an adaptive sample size"
    )(

    // ...concerning the initial correlations 
    "initmode", bpo::value<std::string>()->default_value("generate"),
    "import = import the initial spin correlations from generated spin correlations \
    generate = generate the initial spin correlations from defined analytic functions"
    )(
    // ...if they are imported
    "extrapolate", "extrapolate the imported initial spin correlations linearly to the new discretization"
    )(
    "impcorrfile", bpo::value<std::string>()->default_value(""),
    "provide the filename with path in the Data folder from which the initial correlations should be taken"
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

    // ...concerning the iteration
    "itermode", bpo::value<std::string>()->default_value("absolute"),
    "decide which iteration error should be regarded for the termination condition, options are : absolute, relative"
    )(
    "absiterror", bpo::value<std::string>()->default_value("conservative"),
    "set the absolute tolerance for the iteration error (conservative = 1/4sqrt(3*sample_size))"
    )(
    "reliterror", bpo::value<RealType>()->default_value(RealType{2.0}),
    "set the tolerance for the iteration error relative to the sample std devs"
    )(
    "iterlimit", bpo::value<size_t>()->default_value(20),
    "maximum number of iteration steps (iteration is automatically terminated once the limit is reached)"
    )(

    // ...concerning the eigenvalues
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
    // ========== determine the ID of this run ==========
    std::string Run_ID_str = vm["runID"].as<std::string>();
    if( Run_ID_str == "random" )
    {
        if( my_rank == 0 )
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> distribution(1, 10000);
            Run_ID = distribution(gen); // draw the Run_ID randomly
            std::cout << "\033[1;32mID of this Run: " << std::to_string(Run_ID) << "\033[0m\n";
        }
        MPI_Bcast( &Run_ID, 1, MPI_INT, 0, MPI_COMM_WORLD ); // broadcast Run_ID to all the other ranks
    }
    else
    {
        Run_ID = std::stoi( Run_ID_str );
    }

    // ========== physical parameters ==========
    // ...concerning the mean-field model
    mf_model_name              = vm["mfmodel"].as<std::string>();
    if( mf_model_name != "CorrelationReplica" ){ exit(0); } // before considering a new model: first implement new model and change the folder branch construction

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
    // ...concerning the symmetry
    symmetry_type           = vm["stype"].as<char>();

    // ...concerning the time discretization
    delta_t                 = vm["dt"].as<RealType>();
    num_TimeSteps           = vm["numTimeSteps"].as<size_t>();
    num_TimePoints          = num_TimeSteps+1;

    // ...concerning the statistics
    num_SamplesPerCore      = vm["numSamplesPerCore"].as<size_t>();
    num_SamplesPerSet       = vm["numSamplesPerSet"].as<size_t>();
    num_Samples             = world_size * num_SamplesPerCore;
    seed                    = vm["seed"].as<std::string>();
    if(vm.count("adaptive")){ adaptive_sample_size = true; }
    statistical_error_tolerance = vm["staterrtolerance"].as<RealType>();

    // ...concerning the initial correlations 
    init_corr_mode = vm["initmode"].as<std::string>();
    // ...if they are generated
    RealType Tmax           = delta_t * static_cast<RealType>(num_TimeSteps);
    init_diag_corr    = ph::DiagonalSpinCorrelation{ vm["initdcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    init_nondiag_corr = ph::NonDiagonalSpinCorrelation{ vm["initndcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    // ...if they are imported
    if( init_corr_mode == "import" ) 
    { 
        if( vm.count("extrapolate") )
        {
            extrapolate_imported_spin_correlations = true;
        }
        imported_correlations_src_file = vm["impcorrfile"].as<std::string>(); // determine src folder if necessary
        init_diag_corr.m_name = "";
        init_nondiag_corr.m_name = "";
    }

    // ...concerning the iteration
    iteration_error_mode    = vm["itermode"].as<std::string>();
    relative_iteration_error_tolerance = vm["reliterror"].as<RealType>();
    std::string abs_str = vm["absiterror"].as<std::string>();
    if( abs_str == "conservative" )
    {
        absolute_iteration_error_tolerance = RealType{0.25} / ( std::sqrt( RealType{3.0} * static_cast<RealType>(num_Samples) ) );
    }
    else
    {
        absolute_iteration_error_tolerance = std::stof( abs_str );
    }
    Iteration_Limit         = vm["iterlimit"].as<size_t>();
   
    // ...concerning the eigenvalues
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
    // construct and initialize the mean field model and hand back general physical parameters from it
    std::map< std::string, std::shared_ptr<mfm::MeanFieldModel> > model_map // maps a model name to a shared ptr pointing to a model
    { 
        { "CorrelationReplica", std::make_shared<mfm::CorrelationReplicaModel>() },
    };
    mf_model = model_map.at( mf_model_name );
    mf_model->constructor_function( *this );

    // determine local extra interaction
    local_extra_interaction = ph::LocalExtraInteraction{ vm["extraint"].as<std::string>(), vm["intstrength"].as<RealType>(), num_Spins };

    // import initial correlations according to given scheme
    if( init_corr_mode == "import" )
    {   
        read_initial_correlations_from_file();
    }
}

// read the initial correlations in semi-linearized form from file
void ParameterSpace::read_initial_correlations_from_file()
{
    // 1) open file and group:
    import_filename = "Data/" + imported_correlations_src_file + ".hdf5";
    hid_t file_id = H5Fopen( import_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( file_id < 0 )
    {
        error::IMPORT_FILE_NOT_FOUND( __PRETTY_FUNCTION__, import_filename );
    }
    hid_t group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );

    // 2) import parameters of loaded data and check if they match with the new ones:
    char old_correlation_symmetry_type{};
    std::string old_RealType{};
    IndexPairList old_correlation_categories{};
    hdf5r::import_scalar( group_id, "num_TimePoints", old_num_TimePoints );
    hdf5r::import_scalar( group_id, "delta_t", old_delta_t );
    hdf5r::import_scalar( group_id, "correlation_symmetry_type", old_correlation_symmetry_type );
    hdf5r::import_string( group_id, "RealType", old_RealType );
    hdf5r::import_list( group_id, "correlation_categories", old_correlation_categories );

    if( !extrapolate_imported_spin_correlations ) // check discretization (may only mismatch if extrapolation is used)
    {
        if( old_num_TimePoints != num_TimePoints || old_delta_t != delta_t )
        {
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        }
    }
    if( old_correlation_symmetry_type != symmetry_type ) // check symmetry and type (has to match in any case)
    {
        error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
    }
    #ifdef USE_DOUBLE // check RealType (also has to match in any case)
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
    if( !std::equal(old_correlation_categories.cbegin(),old_correlation_categories.cend(),correlation_categories.cbegin()) ) // check cc
    {
        error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
    }

    // 3) read the correlation data:
    // open group and dataset
    group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );
    hid_t sub_group_id = H5Gopen( file_id, "/results/correlation", H5P_DEFAULT );

    // read linearized data for each (desired) correlation from the dataset
    imported_correlations_linearized = CluList{correlation_categories};
    imported_correlations_linearized.iterate( [&]( std::vector<RealType>& lin_CT, const IndexPair& ij )
    {
        std::string dataset_name = "/results/correlation/" + std::to_string(ij[0]) + "-" + std::to_string(ij[1]);
        hdf5r::import_ND_tensor_linearized( sub_group_id, dataset_name, lin_CT );
        std::for_each( lin_CT.cbegin(), lin_CT.cend(), []( const RealType& x )
        {
            if( std::isnan( x ) ) // check data for nan's
            {
                error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ );
            }
        } );
    } );

    // 4) close resources
    H5Gclose( sub_group_id );
    H5Gclose( group_id );
    H5Fclose( file_id );
}

// return the essential parameters string
std::string ParameterSpace::create_essentials_string() const
{
    size_t pre_colon_space = 20;
    std::stringstream ss{};
    ss << "The following correlations enter the self consistency (numbering starts at 1): " << categories_string(correlation_categories) << "\n";
    if( init_corr_mode == "import" )
    {
        ss << "The initial environment spin correlations are imported from from \'" << import_filename << "\'\n";
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


// ============================================================================
// ============================= HELPER FUNCTIONS =============================
// ============================================================================
// transform IndexPairList to string
std::string categories_string( const IndexPairList& ipairs )
{
    std::string str{};
    for( auto ipair : ipairs ){ str += std::to_string(ipair[0]+1) + "-" + std::to_string(ipair[1]+1) + ","; }
    str.pop_back(); // remove last comma
    return str;
}

}