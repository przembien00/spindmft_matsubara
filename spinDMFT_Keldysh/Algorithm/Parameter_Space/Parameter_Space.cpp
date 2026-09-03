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
    )(
    "bath", bpo::value<std::string>()->default_value("selfconsistent"),
    "field covariance source: selfconsistent or harmonic"
    )(
    "bathOmega", bpo::value<RealType>()->default_value(RealType{1.0}),
    "frequency omega_0 > 0 of the prescribed harmonic bath"
    )(
    "bathCoupling", bpo::value<RealType>()->default_value(RealType{1.0}),
    "coupling g >= 0 in the prescribed bath operator X=g(a+a^dagger)"
    )(
    "bathComponent", bpo::value<char>()->default_value('x'),
    "spin component coupled to the prescribed harmonic bath: x, y, or z; harmonic mode requires cstype=D"
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
    "numImagTimeSteps", bpo::value<size_t>()->default_value(size_t{5}),
    "set the number of imaginary-time steps on [0,beta]"
    )(
    "dt", bpo::value<RealType>()->default_value(RealType{0.1}),
    "set the step width for the equidistant time discretization (time in units of hbar/JQ)"
    )(
    "numRealTimeSteps", bpo::value<size_t>()->default_value(size_t{10}),
    "set the number of real-time steps on [0,Tmax]"
    )(
    "Tmax", bpo::value<RealType>()->default_value(RealType{1.0}),
    "set the maximum real time"
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
    "samplingStrategy", bpo::value<std::string>()->default_value("pcn"),
    "Monte-Carlo strategy: pcn targets p0(r) Re Z_M(r); independent retains the bare-Gaussian complex-ratio estimator"
    )(
    "antitheticPairs",
    "evaluate Gaussian latent states as r and -r pairs; independent mode counts total trajectories, while pCN mode counts sign-symmetrized pair states"
    )(
    "mhStepSize", bpo::value<RealType>()->default_value(RealType{0.3}),
    "pCN proposal parameter beta in (0,1]: r'=sqrt(1-beta^2)r+beta xi"
    )(
    "mhBurnIn", bpo::value<size_t>()->default_value(size_t{100}),
    "pCN burn-in steps per MPI chain after every bath/self-consistency update"
    )(
    "partitionImagTolerance", bpo::value<RealType>()->default_value(RealType{1e-8}),
    "maximum observed |Im Z_M|/|Re Z_M| allowed by the real-positive pCN target assumption"
    )(
    "numBlocks", bpo::value<size_t>()->default_value(size_t{32}),
    "number of contiguous blocks per core; pCN errors use a multi-scale batch-means plateau"
    )(
    "gaussianFactorization", bpo::value<std::string>()->default_value("dense"),
    "complex pseudo-covariance factorization: dense, svd, or fft"
    )(
    "fftCrossFrequencyCutoff", bpo::value<RealType>()->default_value(RealType{3.}),
    "for gaussianFactorization=fft, retain dense Matsubara-real coupling for |omega| up to this cutoff; a negative value disables truncation"
    )(
    "spinInsertionStrategy", bpo::value<std::string>()->default_value("full-contour"),
    "real-time spin insertion: full-contour uses B_-(T,0) U_+(t,T) S U_+(t,0); prefix uses U_+(0,t) S B_-(t,0)"
    )(
    // ...concerning the iteration
    "reliterror", bpo::value<RealType>()->default_value(RealType{5.0}),
    "maximum pointwise |F(x)-x| in units of its sampling-strategy-aware standard error"
    )(
    "iterlimit", bpo::value<size_t>()->default_value(size_t{20}),
    "set the the maximum number of iterations of the self-consistency problem"
    )(
    "mixingAlpha", bpo::value<RealType>()->default_value(RealType{1.}),
    "mix the physical magnetization and edge-grid contour primitive with alpha in (0,1]"
    )(
    "constantMagnetization",
    "project m_a(t) to m_a(0) in first- and second-moment self-consistency; retain the measured m_a(t) for output and diagnostics"
    )(
    "covarianceTolerance", bpo::value<RealType>()->default_value(RealType{1e-10}),
    "maximum final covariance transpose-symmetry residual"
    )(
    "branchIdentityTolerance", bpo::value<RealType>()->default_value(RealType{0.1}),
    "maximum pre-mirroring branch transpose-identity residual"
    )(
    "takagiTolerance", bpo::value<RealType>()->default_value(RealType{1e-8}),
    "maximum Autonne--Takagi reconstruction residual"
    )(
    "minimumPhaseMagnitude", bpo::value<RealType>()->default_value(RealType{1e-6}),
    "minimum resolved magnitude |sum Z_M|/sum |Z_M|"
    )(
    "denominatorConstancyTolerance", bpo::value<RealType>()->default_value(RealType{0.1}),
    "maximum max_t |sum D(t)/sum Z_M-1| convergence diagnostic"
    )(
    "imagMagSigma", bpo::value<RealType>()->default_value(RealType{5.}),
    "allowed imaginary magnetization in units of its jackknife standard error"
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
    "initcorrsrc", bpo::value<std::string>()->default_value("spinDMFT_Keldysh"),
    "deprecated compatibility option; initial correlations are read from this solver's Data folder"
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
    bath=vm["bath"].as<std::string>();
    bath_frequency=vm["bathOmega"].as<RealType>();
    bath_coupling=vm["bathCoupling"].as<RealType>();
    bath_component=vm["bathComponent"].as<char>();
    if( bath!="selfconsistent" && bath!="harmonic" )
        throw std::invalid_argument("bath must be selfconsistent or harmonic");
    if( uses_harmonic_bath() ) self_consistency=false;
    RealType Babs = vm["Babs"].as<RealType>();
    if( !vm.count("dontrescaleB") && !uses_harmonic_bath() )
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
    if( uses_harmonic_bath() )
    {
        if( beta<=RealType{0.} )
            throw std::invalid_argument("the prescribed harmonic bath requires beta > 0");
        if( bath_frequency<=RealType{0.} )
            throw std::invalid_argument("bathOmega must be positive");
        if( bath_coupling<RealType{0.} )
            throw std::invalid_argument("bathCoupling must be non-negative");
        if( bath_component!='x' && bath_component!='y' && bath_component!='z' )
            throw std::invalid_argument("bathComponent must be x, y, or z");
        if( correlation_symmetry_type!='D' )
            throw std::invalid_argument(
                "the single-component harmonic bath requires cstype=D");
    }

    // ...concerning the time discretization
    num_TimeSteps  = vm["numImagTimeSteps"].as<size_t>();
    if( num_TimeSteps == 0 ) throw std::invalid_argument("numImagTimeSteps must be positive");
    num_TimePoints = num_TimeSteps + 1;
    delta_t        = beta / static_cast<RealType>(num_TimeSteps);
    num_RealTimeSteps  = vm["numRealTimeSteps"].as<size_t>();
    if( num_RealTimeSteps == 0 ) throw std::invalid_argument("numRealTimeSteps must be positive");
    num_RealTimePoints = num_RealTimeSteps + 1;
    Tmax                = vm["Tmax"].as<RealType>();
    if( Tmax < RealType{0.} ) throw std::invalid_argument("Tmax must be non-negative");
    delta_real_t        = Tmax / static_cast<RealType>(num_RealTimeSteps);
    // ...concerning the statistics
    num_SamplesPerCore      = vm["numSamplesPerCore"].as<size_t>();
    num_SamplesPerSet       = vm["numSamplesPerSet"].as<size_t>();
    num_SetsPerCore         = (size_t) std::ceil( static_cast<RealType>(num_SamplesPerCore) / static_cast<RealType>(num_SamplesPerSet) ); 
    num_Samples             = world_size * num_SamplesPerCore;
    seed                    = vm["seed"].as<std::string>();
    sampling_strategy       = vm["samplingStrategy"].as<std::string>();
    antithetic_pairs        = vm.count("antitheticPairs")>0;
    mh_step_size            = vm["mhStepSize"].as<RealType>();
    mh_burn_in              = vm["mhBurnIn"].as<size_t>();
    partition_imaginary_tolerance=vm["partitionImagTolerance"].as<RealType>();
    num_blocks              = vm["numBlocks"].as<size_t>();
    if( sampling_strategy!="pcn"&&sampling_strategy!="independent" )
        throw std::invalid_argument("samplingStrategy must be pcn or independent");
    if( antithetic_pairs&&sampling_strategy=="independent"
        &&num_SamplesPerCore%2!=0 )
        throw std::invalid_argument(
            "independent antitheticPairs requires an even numSamplesPerCore on every MPI rank");
    if( mh_step_size<=RealType{0.}||mh_step_size>RealType{1.} )
        throw std::invalid_argument("mhStepSize must lie in (0,1]");
    if( partition_imaginary_tolerance<RealType{0.} )
        throw std::invalid_argument("partitionImagTolerance must be non-negative");
    gaussian_factorization  = vm["gaussianFactorization"].as<std::string>();
    if( gaussian_factorization!="dense" && gaussian_factorization!="svd"
        && gaussian_factorization!="fft" )
        throw std::invalid_argument(
            "gaussianFactorization must be dense, svd, or fft" );
    fft_cross_frequency_cutoff=vm["fftCrossFrequencyCutoff"].as<RealType>();
    spin_insertion_strategy=vm["spinInsertionStrategy"].as<std::string>();
    if( spin_insertion_strategy!="full-contour"
        && spin_insertion_strategy!="prefix" )
        throw std::invalid_argument(
            "spinInsertionStrategy must be full-contour or prefix" );

    // ...concerning the initially inserted correlations and the iteration 
    iteration_error_sigma_threshold=vm["reliterror"].as<RealType>();
    Iteration_Limit         = vm["iterlimit"].as<size_t>();
    mixing_alpha            = vm["mixingAlpha"].as<RealType>();
    constant_magnetization_time=vm.count("constantMagnetization")>0;
    if( mixing_alpha<=RealType{0.} || mixing_alpha>RealType{1.} )
        throw std::invalid_argument("mixingAlpha must lie in (0,1]");
    covariance_tolerance    = vm["covarianceTolerance"].as<RealType>();
    branch_identity_tolerance=vm["branchIdentityTolerance"].as<RealType>();
    takagi_tolerance        = vm["takagiTolerance"].as<RealType>();
    minimum_phase_magnitude = vm["minimumPhaseMagnitude"].as<RealType>();
    denominator_constancy_tolerance=vm["denominatorConstancyTolerance"].as<RealType>();
    imaginary_magnetization_sigma=vm["imagMagSigma"].as<RealType>();
    if( iteration_error_sigma_threshold<RealType{0.}
        ||covariance_tolerance<RealType{0.} || branch_identity_tolerance<RealType{0.}
        ||takagi_tolerance<RealType{0.} ||minimum_phase_magnitude<RealType{0.}
        ||minimum_phase_magnitude>RealType{1.}
        ||denominator_constancy_tolerance<RealType{0.}
        ||imaginary_magnetization_sigma<RealType{0.} )
        throw std::invalid_argument("Keldysh convergence and diagnostic tolerances must be non-negative and minimumPhaseMagnitude <= 1");
    init_diag_corr    = ph::DiagonalSpinCorrelation{ vm["initdcorr"].as<std::string>(), beta, vm["corrperiods"].as<RealType>() };
    init_nondiag_corr = ph::NonDiagonalSpinCorrelation{ vm["initndcorr"].as<std::string>(), beta, vm["corrperiods"].as<RealType>() };
    // initial correlations are perhaps read from file below

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
        initial_correlations_src_directory = "spinDMFT_Keldysh";
        read_initial_correlations_from_file();
    }

}

// method for importing the initial correlations
void ParameterSpace::read_initial_correlations_from_file()
{ 
    // 1) open file and group:
    std::string total_filename = "Data/" + initial_correlations_src_file + ".hdf5";
    hid_t file_id = H5Fopen( total_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( file_id < 0 ){ error::IMPORT_FILE_NOT_FOUND( __PRETTY_FUNCTION__, total_filename ); }
    hid_t group_id = H5Gopen( file_id, "/parameters", H5P_DEFAULT );
    // 2) import parameters of loaded data and check if they match with the new ones:
    char old_correlation_symmetry_type{}; 
    std::string old_RealType{};
    if( H5Aexists(group_id,"num_ImagTimePoints")>0 )
        hdf5r::import_scalar( group_id, "num_ImagTimePoints", old_num_TimePoints );
    else
        hdf5r::import_scalar( group_id, "num_TimePoints", old_num_TimePoints );
    if( H5Aexists(group_id,"delta_imag_t")>0 )
        hdf5r::import_scalar( group_id, "delta_imag_t", old_delta_t );
    else
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
    H5Gclose( group_id );
    group_id = H5Gopen( file_id, "/results", H5P_DEFAULT );
    if( initial_correlations_src_directory == "spinDMFT_Keldysh" )
    {
        const std::string re_name=H5Lexists(group_id,"Re_correlation_edge",H5P_DEFAULT)>0
            ?"/results/Re_correlation_edge":"/results/Re_correlation";
        hdf5r::import_ND_tensor_linearized( group_id, re_name,
                                            initial_correlations_linearized );
        if( H5Lexists(group_id,"Im_correlation_edge",H5P_DEFAULT)>0 )
            hdf5r::import_ND_tensor_linearized(
                group_id,"/results/Im_correlation_edge",
                initial_correlations_imag_linearized );
        else if( H5Lexists(group_id,"Im_correlation",H5P_DEFAULT)>0 )
            hdf5r::import_ND_tensor_linearized(
                group_id,"/results/Im_correlation",
                initial_correlations_imag_linearized );
        else
            initial_correlations_imag_linearized.assign(
                initial_correlations_linearized.size(),RealType{0.});
        auto reject_nan=[]( const RealType& v )
        {
            if( std::isnan( v ) )
                error::INIT_CORRELATIONS_CONTAIN_NANS( __PRETTY_FUNCTION__ );
        };
        std::for_each(initial_correlations_linearized.cbegin(),
                      initial_correlations_linearized.cend(),reject_nan);
        std::for_each(initial_correlations_imag_linearized.cbegin(),
                      initial_correlations_imag_linearized.cend(),reject_nan);
    }
    else
    {
        error::SRC_DIRECTORY( initial_correlations_src_directory, __PRETTY_FUNCTION__ );
    }
    // 4) import m(0). New files store magnetization only as a trajectory;
    // retain support for the legacy scalar checkpoint schema.
    if( H5Lexists(group_id,"Re_magnetization",H5P_DEFAULT)>0 )
    {
        std::vector<RealType> magnetization_linearized;
        hdf5r::import_ND_tensor_linearized(
            group_id,"/results/Re_magnetization",magnetization_linearized );
        if( magnetization_linearized.size()<initial_spin_expval.size() )
            error::INIT_CORRELATIONS_MISMATCH( __PRETTY_FUNCTION__ );
        std::copy_n(magnetization_linearized.cbegin(),initial_spin_expval.size(),
                    initial_spin_expval.begin());
    }
    else
    {
        hdf5r::import_scalar( group_id, "S_x", initial_spin_expval[0] );
        hdf5r::import_scalar( group_id, "S_y", initial_spin_expval[1] );
        hdf5r::import_scalar( group_id, "S_z", initial_spin_expval[2] );
    }

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
        ss << "I don't perform self-consistent iterations in this simulation!\n";
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
    << print::quantity_to_output_line( pre_colon_space, "num_ImagTimeSteps", std::to_string(num_TimeSteps) )
    << print::quantity_to_output_line( pre_colon_space, "delta_t"       , print::remove_zeros(print::round_value_to_string(delta_t,num_PrintDigits)) ) 
    << print::quantity_to_output_line( pre_colon_space, "num_RealTimeSteps", std::to_string(num_RealTimeSteps) )
    << print::quantity_to_output_line( pre_colon_space, "Tmax", print::remove_zeros(print::round_value_to_string(Tmax,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "delta_real_t", print::remove_zeros(print::round_value_to_string(delta_real_t,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "equal_time_prescription", "symmetric_theta_half" )
    << print::quantity_to_output_line( pre_colon_space, "num_Samples"   , std::to_string(num_Samples) )
    << print::quantity_to_output_line( pre_colon_space, "sampling_strategy", sampling_strategy )
    << print::quantity_to_output_line( pre_colon_space, "antithetic_pairs", print::bool_to_string(antithetic_pairs) )
    << print::quantity_to_output_line( pre_colon_space, "mh_step_size", print::remove_zeros(print::round_value_to_string(mh_step_size,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "mh_burn_in", std::to_string(mh_burn_in) )
    << print::quantity_to_output_line( pre_colon_space, "gaussian_factorization", gaussian_factorization )
    << print::quantity_to_output_line( pre_colon_space, "fft_cross_frequency_cutoff", std::to_string(fft_cross_frequency_cutoff) )
    << print::quantity_to_output_line( pre_colon_space, "iteration_error_sigma_threshold", print::remove_zeros(print::round_value_to_string(iteration_error_sigma_threshold,num_PrintDigits)) )
    << print::quantity_to_output_line( pre_colon_space, "constant_magnetization_time", print::bool_to_string(constant_magnetization_time) )
    << print::quantity_to_output_line( pre_colon_space, "information_text" , information_text );
    return ss.str();
}

};
