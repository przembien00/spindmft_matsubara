#include"Parameter_Space.h"

#include<iomanip>
#include<memory>
#include<sstream>
#include<fstream>
#include<cmath>
#include<hdf5.h>
#include<boost/program_options.hpp>
namespace bpo = boost::program_options;

#include<Globals/MPI_Types.h>

#include<HDF5/HDF5_Routines.h>
namespace hdf5r = HDF5_Routines;

#include<Standard_Algorithms/Print_Routines.h>
namespace print = Print_Routines;

#include<Observables/Tensors.h>
namespace ten = Observables::Tensors;

#include<Physics/Spin.h>
namespace sp = Physics::Spin;

#include"PS_Error_Handling.h"
namespace error = spinDMFT::Parameter_Space::Error_Handling;

namespace spinDMFT::Parameter_Space
{

// ===============================================================
// ==================== PARAMETER SPACE CLASS ====================
// ===============================================================
// constructor : build parameter space from command line arguments
ParameterSpace::ParameterSpace( const int argC, char* const argV[], const int world_size, const int my_rank ):
    my_rank( my_rank ),
    world_size( world_size )
{
    // 1a.) Define Options
    bpo::options_description description("All Options:");
    bpo::options_description description_help("Help options:");
    bpo::options_description description_physics("General options concerning physics:");
    bpo::options_description description_numerics("General options concerning numerics:");
    bpo::options_description description_storing("General options concerning storing and naming:");

    description_help.add_options()
    (
    "help", "show all options"
    )(
    "helpNum", "show allowed options concerning numerics"
    );

    description_physics.add_options()
    (
    "noselfcons", "perform only one iteration, and ignore the self-consistency (however, the sc-equations are applied to the initial spin-correlations)"
    )(
    // ========== model and physical parameters ==========
    "spinmodel", bpo::value<std::string>()->default_value("ISO"),
    "set the spin model || options are : \
    ISO = Isotropic Heisenberg Model, \
    XXZ = Easy-Axis anisotropic Heisenberg Model "
    )(
    "spin", bpo::value<std::string>()->default_value("1/2"),
    "set the spin value as fraction or float (1/2, 1, 1.5, ...)"
    )(
    "lambda", bpo::value<RealType>()->default_value(RealType{2.0}),
    "set the anisotropy factor for the XXZ-model or XY-model; the factor will be multiplied to the z-couplings of a Heisenbergmodel"
    )(
    "rho", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the second anisotropy factor for the XY-model; the factor will be multiplied to the xy-couplings"
    )(
    "JQ", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the value of JQ (quadratic coupling factor)"
    )(
    "JL", bpo::value<RealType>()->default_value(RealType{0.}),
    "set the value of JL (linear coupling factor)"
    )
    (
    "Bname", bpo::value<std::string>()->default_value("none"),
    "set the external field type (direction): none, x, y, z, arb" 
    )(
    "Babs", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the absolute value of the external field B (Larmor frequency in units of JQ)"
    )(
    "Bphi", bpo::value<RealType>()->default_value(RealType{0.0}),
    "set the polar angle of the external field B" 
    )(    
    "Btheta", bpo::value<RealType>()->default_value(RealType{0.0}),
    "set the azimuth angle of the external field B"
    )( 
    "dontrescaleB", "rescaling factor is not multiplied to the field B"
    )( 
    "Cname", bpo::value<std::string>()->default_value("none"),
    "set the (static) noise field type: none, z" // non-static not implemented
    )(
    "C", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the (static) noise field variance (units of JQ^2)"
    )(
    "extraint", bpo::value<std::string>()->default_value("none"),
    "include a local extra interaction term || options are: none, quadrupolar"
    )(
    "intstrength", bpo::value<RealType>()->default_value(RealType{0.0}),
    "set the strength of the extra interaction term"
    )(
    "beta", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the value of the inverse temperature"
    );   

    description_numerics.add_options()
    (
    // ========== general numerical parameters ==========
    // ...concerning the symmetry
    "cstype", bpo::value<char>()->default_value('D'),
    "set the correlation symmetry type || options are: \
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
    "set the step width for the equidistant time discretization (time in units of hbar/JQ)"
    )(
    // ...concerning the statistics
    "numSamplesPerCore", bpo::value<size_t>()->default_value(size_t{100}),
    "set the number of samples per core (!) for the Monte Carlo simulation"
    )(
    "numSamplesPerSet", bpo::value<size_t>()->default_value(size_t{100}),
    "set the number of samples per set for the Monte Carlo simulation"
    )(
    "seed", bpo::value<std::string>()->default_value("random"),
    "set the seed for the random generator (random : seed is determined by clock)"
    )(
    // ...concerning the iteration 
    "reliterror", bpo::value<RealType>()->default_value(RealType{5.0}),
    "set the the threshold for the iteration error relative to the estimated Monte-Carlo std devs"
    )(
    "iterlimit", bpo::value<size_t>()->default_value(size_t{20}),
    "set the the maximum number of iterations of the self-consistency problem"
    )(
    // ...concerning the initially inserted correlations
    "initdcorr", bpo::value<std::string>()->default_value("imagtime"),
    "set the initial diagonal correlations (gxx=gyy=gzz)"
    )(
    "initndcorr", bpo::value<std::string>()->default_value("zero"),
    "set the initial nondiagonal correlations (gxx=gyy=gzz)"
    )(
    "corrperiods", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the number of periods in case of periodical initial correlations"
    )(
    "loadinit", "import the initial spin correlations from existing spin correlations in the data folder"
    )(
    "extrapolate", "extrapolate the imported initial spin correlations linearly to a new discretization"
    )(
    "initcorrfile", bpo::value<std::string>()->default_value("[auto]"),
    "provide the filename from which the initial correlations should be taken"
    )(
    "initcorrsrc", bpo::value<std::string>()->default_value("spinDMFT"),
    "provide the root folder of the file from which the initial correlations should be taken (spinDMFT, CspinDMFT)"
    )(
    
    // ...concerning potential negative eigenvalues of the covariance matrix
    "critneg", bpo::value<RealType>()->default_value(RealType{0.001}),
    "set the prefactor for the critical ratio between the positive and negative eigenvalue sum; in case of a violation the computation terminates"
    )(
    "truncneg", bpo::value<std::string>()->default_value("set_zero"),
    "decide for the truncation scheme for negative eigenvalues || options are: set_zero, abs"
    );

    description_storing.add_options()
    (
    // ========== storing and naming ==========
    "info", bpo::value<std::string>()->default_value(""),
    "provide some text that informs about the procedure and purpose of the computation; \
    to be able to use spaces, put the text in quotation marks"
    )(
    "project", bpo::value<std::string>()->default_value(""),
    "sort the data into a project for better order"
    )(
    "fileext", bpo::value<std::string>()->default_value(""),
    "Define an extension to the filename; it will be appended according to : filename_fileext"
    )(
    "numPrintDigits", bpo::value<size_t>()->default_value(size_t{4}),
    "set the value precision for printing to the terminal"
    );

    // 1b.) Store command line arguments in bp options
    description.add( description_help );
    description.add( description_physics );
    description.add( description_numerics );
    description.add( description_storing );
    bpo::variables_map vm;
    bpo::store( bpo::parse_command_line(argC, argV, description), vm );
    bpo::notify( vm );

    // 2.) Output of help descriptions and termination of the program
    if( vm.count("help") || vm.count("helpNum") )
    {
        if( vm.count("help") && my_rank==0 ){
            std::cout << description << "\n";
        }
        else if( vm.count("helpNum") && my_rank==0 ){
            std::cout << description_numerics << "\n";
        }
        MPI_Finalize();
        exit(0);
    }

    // 3.) Store bp options in parameter space 
    // ========== model and physical parameters ==========
    if( vm.count("noselfcons") )
    {
        self_consistency = false;
    }
    spin_model  = ph::SpinModel{ vm["spinmodel"].as<std::string>(), vm["lambda"].as<RealType>(), vm["rho"].as<RealType>() };
    spin        = vm["spin"].as<std::string>();
    spin_float  = sp::convert_spin_string_to_float( spin );
    num_HilbertSpaceDimension = std::lround(2*spin_float + 1);
    JQ          = vm["JQ"].as<RealType>();
    JL          = vm["JL"].as<RealType>();
    RealType Babs = vm["Babs"].as<RealType>();
    if( !vm.count("dontrescaleB") )
    {
        Babs = JQ * Babs; // rescale the external field
    }
    B           = ph::MagneticField{ vm["Bname"].as<std::string>(), Babs, vm["Bphi"].as<RealType>(), vm["Btheta"].as<RealType>() };
    noise       = ph::Noise{ vm["Cname"].as<std::string>(), vm["C"].as<RealType>() };
    extra_interaction = ph::ExtraInteraction{ vm["extraint"].as<std::string>(), vm["intstrength"].as<RealType>() };
    beta = vm["beta"].as<RealType>();

    // ========== general numerical parameters ==========
    // ...concerning the symmetry
    correlation_symmetry_type = vm["cstype"].as<char>();

    // ...concerning the time discretization
    num_TimeSteps  = vm["numTimeSteps"].as<size_t>();
    num_TimePoints = num_TimeSteps + 1;
    delta_t        = beta / static_cast<RealType>(num_TimeSteps);
    RealType Tmax  = beta;

    // ...concerning the statistics
    num_SamplesPerCore      = vm["numSamplesPerCore"].as<size_t>();
    num_SamplesPerSet       = vm["numSamplesPerSet"].as<size_t>();
    num_SetsPerCore         = (size_t) std::ceil( static_cast<RealType>(num_SamplesPerCore) / static_cast<RealType>(num_SamplesPerSet) ); 
    num_Samples             = world_size * num_SamplesPerCore;
    seed                    = vm["seed"].as<std::string>();

    // ...concerning the initially inserted correlations and the iteration 
    absolute_iteration_error_threshold = vm["reliterror"].as<RealType>() / std::sqrt( RealType{3.0} * static_cast<RealType>(num_Samples) ) * spin_float * (spin_float+1.) / 3.;
    Iteration_Limit         = vm["iterlimit"].as<size_t>();
    init_diag_corr    = ph::DiagonalSpinCorrelation{ vm["initdcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    init_nondiag_corr = ph::NonDiagonalSpinCorrelation{ vm["initndcorr"].as<std::string>(), Tmax, vm["corrperiods"].as<RealType>() };
    // initial correlations are perhaps read from file below

    // ...concerning potential negative eigenvalues of the covariance matrix
    critical_eigenvalue_ratio               = vm["critneg"].as<RealType>(); // / std::sqrt( static_cast<RealType>(48*num_Samples) );
    truncation_scheme_negative_eigenvalues  = vm["truncneg"].as<std::string>();

    // ========== saving and naming ==========
    information_text        = vm["info"].as<std::string>();
    project_name            = vm["project"].as<std::string>();
    filename_extension      = vm["fileext"].as<std::string>();
    num_PrintDigits         = vm["numPrintDigits"].as<size_t>();

    initial_spin_expval = FieldVector{0., 0., 0.};
    if( vm.count("loadinit") )
    {
        load_initial_spin_correlations = true;
        if( vm.count("extrapolate") )
        {
            extrapolate_initial_spin_correlations = true;
        }
        initial_correlations_src_file = vm["initcorrfile"].as<std::string>();
        initial_correlations_src_directory = vm["initcorrsrc"].as<std::string>();
        read_initial_correlations_from_file();
    }

}

// method for importing the initial correlations
void ParameterSpace::read_initial_correlations_from_file()
{ 
    // 1) open file and group:
    std::string total_filename = "../" + initial_correlations_src_directory + "/Data/" + initial_correlations_src_file + ".hdf5";
    hid_t file_id = H5Fopen( total_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( file_id < 0 ){ error::IMPORT_FILE_NOT_FOUND( __PRETTY_FUNCTION__, total_filename ); }
    hid_t group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );
    // 2) import parameters of loaded data and check if they match with the new ones:
    char old_correlation_symmetry_type{}; 
    std::string old_RealType{};
    hdf5r::import_scalar( group_id, "num_TimePoints", old_num_TimePoints );
    hdf5r::import_scalar( group_id, "delta_t", old_delta_t );
    hdf5r::import_scalar( group_id, "correlation_symmetry_type", old_correlation_symmetry_type );
    hdf5r::import_string( group_id, "RealType", old_RealType );
    // check discretization (may only mismatch if extrapolation is used)
    if( !extrapolate_initial_spin_correlations )
    {
        if( old_num_TimePoints != num_TimePoints )
        {
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        }
    }
    // check symmetry and type (has to match in any case)
    if( old_correlation_symmetry_type != correlation_symmetry_type )
    {
        error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
    }
//    else // check RealType (also has to match in any case)
//    {
//        #ifdef USE_DOUBLE
//        if( old_RealType != "DOUBLE" )
//        {
//            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
//        }
//        #endif
//        #ifdef USE_FLOAT 
//        if( old_RealType != "FLOAT" )
//        {
//            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
//        }
//        #endif
//    }
    // 3) import the correlation data:
    group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );    
    if( initial_correlations_src_directory == "spinDMFT" )
    {
        hdf5r::import_ND_tensor_linearized( group_id, "/results/Re_correlation", initial_correlations_linearized );
        std::for_each( initial_correlations_linearized.cbegin(), initial_correlations_linearized.cend(), []( const RealType& v ) 
        {
            if( std::isnan( v ) )
            {
                error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ ); // check data for nan's
            }
        } );
    }
    else
    {
        error::SRC_DIRECTORY( initial_correlations_src_directory, __PRETTY_FUNCTION__ );
    }
    // 4) import spin expectation values:
    hdf5r::import_scalar( group_id, "S_x", initial_spin_expval[0] );
    hdf5r::import_scalar( group_id, "S_y", initial_spin_expval[1] );
    hdf5r::import_scalar( group_id, "S_z", initial_spin_expval[2] );

    // 5) close resources:
    H5Gclose( group_id );
    H5Fclose( file_id );
}

// printing method : return the essential parameters string
std::string ParameterSpace::create_essentials_string() const
{   
    size_t pre_colon_space = 35;
    std::stringstream ss{};

    if( !self_consistency )
    {
        ss << "I don't perform self-consistent iterations in this simulation!\n";
    }
    ss
    << print::quantity_to_output_line( pre_colon_space, "spin_model_string"     , spin_model.compact_info( num_PrintDigits ) )
    << print::quantity_to_output_line( pre_colon_space, "spin"          , spin )
    << print::quantity_to_output_line( pre_colon_space, "JQ"            , print::remove_zeros(print::round_value_to_string(JQ,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "JL"            , print::remove_zeros(print::round_value_to_string(JL,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "beta"          , print::remove_zeros(print::round_value_to_string(beta,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "B_string"      , B.compact_info( num_PrintDigits ) )
    << print::quantity_to_output_line( pre_colon_space, "noise_string"  , noise.compact_info( num_PrintDigits ) )
    << print::quantity_to_output_line( pre_colon_space, "extra_interaction_string", extra_interaction.compact_info( num_PrintDigits ) )
    << print::quantity_to_output_line( pre_colon_space, "correlation_symmetry_type" , std::string(1, correlation_symmetry_type) )
    << print::quantity_to_output_line( pre_colon_space, "num_TimeSteps" , std::to_string(num_TimeSteps) ) 
    << print::quantity_to_output_line( pre_colon_space, "delta_t"       , print::remove_zeros(print::round_value_to_string(delta_t,num_PrintDigits)) ) 
    << print::quantity_to_output_line( pre_colon_space, "num_Samples"   , std::to_string(num_Samples) )
    << print::quantity_to_output_line( pre_colon_space, "information_text" , information_text );
    return ss.str();
}

};

