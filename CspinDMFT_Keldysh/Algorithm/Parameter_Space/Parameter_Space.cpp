#include "Parameter_Space.h"
#include "../Mean_Field_Models/Mean_Field_Models.h"
#include <Physics/Spin.h>
#include "../Storage_Concept/HDF5_Utilities.h"
#include <boost/program_options.hpp>
#include <mpi.h>
#include <filesystem>
#include <iostream>
#include <random>
#include <limits>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <algorithm>
#include <queue>

namespace CspinDMFT::Parameter_Space {
namespace {
bool one_of(const std::string& s,std::initializer_list<const char*> values) {
    for(auto v:values) if(s==v) return true;
    return false;
}
void relative_path(const std::string& value,const char* option) {
    auto path=std::filesystem::path(value);
    if(path.is_absolute()) throw std::invalid_argument(std::string(option)+" must be relative");
    for(const auto& part:path) if(part=="..") throw std::invalid_argument(std::string(option)+" cannot contain ..");
}
}
ParameterSpace::ParameterSpace(int argc,char* const argv[],int world_size,int rank):num_Cores(world_size) {
    namespace po=boost::program_options;
    std::string spinspin="ISO",spinmf="auto",noise_name="none",run_id="random";
    double lambda=2.,rho=1.,smflambda=2.,smfrho=1.,noise_variance=1.;
    bool joint=false;
    size_t preparation_samples=0,preparation_limit=0;
    po::options_description all("All options"),help("Help"),physics("Cluster physics"),numerics("Numerics"),storing("Storing and naming");
    help.add_options()("help","show all options")("helpNum","show numerical options");
    physics.add_options()
      ("mfmodel",po::value(&mf_model_name)->default_value(mf_model_name),"CorrelationReplica")
      ("config",po::value(&config_file),"configuration basename in Configuration_Data/<project>, or existing HDF5 path")
      ("spinspinmodel",po::value(&spinspin)->default_value(spinspin),"ISO, Ising, XXZ, DRF, XY, het")
      ("spinmodel",po::value<std::string>(),"Keldysh alias for spinspinmodel")
      ("lambda",po::value(&lambda)->default_value(lambda),"local z anisotropy")
      ("rho",po::value(&rho)->default_value(rho),"local xy mixing")
      ("spinmfmodel",po::value(&spinmf)->default_value(spinmf),"auto follows spinspinmodel; otherwise independent bath model")
      ("smflambda",po::value(&smflambda)->default_value(smflambda),"bath z anisotropy")
      ("smfrho",po::value(&smfrho)->default_value(smfrho),"bath xy mixing")
      ("chemshift",po::value(&chemical_shift)->default_value(chemical_shift),"none or comma-separated z fields, one per site")
      ("unifield",po::value<double>(),"uniform z field applied equally to all sites; alternative to chemshift")
      ("stagfield",po::value<double>(),"staggered z field (+h at site 1); requires a connected bipartite coupling graph")
      ("extraint",po::value(&extra_interaction)->default_value(extra_interaction),"none or quadrupolar")
      ("intstrength",po::value(&interaction_strength)->default_value(interaction_strength),"quadrupolar coefficient: H_extra=3*intstrength*sum_i (S_i^z)^2")
      ("rescmf",po::value(&rescale_meanfield)->default_value(rescale_meanfield),"scale bath coupling")
      ("rescss",po::value(&rescale_spinspin)->default_value(rescale_spinspin),"scale local spin-spin coupling")
      ("Cname",po::value(&noise_name)->default_value(noise_name),"independent static site noise: none or z")
      ("C",po::value(&noise_variance)->default_value(noise_variance),"static noise variance")
      ("beta",po::value(&beta)->default_value(beta),"inverse temperature")
      ("noselfcons","sample the initial prescribed field without updating feedback")
      ("bath",po::value(&bath)->default_value(bath),"selfconsistent or harmonic (independent identical baths per site)")
      ("bathOmega",po::value(&bath_frequency)->default_value(bath_frequency),"harmonic bath frequency")
      ("bathCoupling",po::value(&bath_coupling)->default_value(bath_coupling),"harmonic bath coupling")
      ("bathComponent",po::value(&bath_component)->default_value(bath_component),"x, y or z");
    numerics.add_options()
      ("cstype",po::value(&correlation_symmetry_type)->default_value(correlation_symmetry_type),"A, B, C or D; D retains all directions")
      ("stype",po::value<char>(),"CspinDMFT alias for cstype")
      ("timeTranslationInvariant",po::value(&time_translation_invariant)->default_value(true),"equilibrium reconstruction; must be true")
      ("numImagTimeSteps",po::value(&num_TimeSteps)->default_value(num_TimeSteps),"imaginary intervals on [0,beta]")
      ("numRealTimeSteps",po::value(&num_RealTimeSteps)->default_value(num_RealTimeSteps),"measurement intervals on [0,Tmax]")
      ("Tmax",po::value(&Tmax)->default_value(Tmax),"maximum real time")
      ("realTimeSubsteps",po::value(&real_time_substeps)->default_value(real_time_substeps),"positive Strang microsteps per interval")
      ("realTimeEvolution",po::value(&real_time_evolution)->default_value(real_time_evolution),"staggered-strang")
      ("numSamplesPerCore",po::value(&num_SamplesPerCore)->default_value(num_SamplesPerCore),"production samples per rank and iteration")
      ("preparationSamplesPerCore",po::value(&preparation_samples)->default_value(0),"preparation samples; 0 uses numSamplesPerCore")
      ("numSamplesPerSet",po::value(&num_SamplesPerSet)->default_value(1),"independent field batch size")
      ("numBlocks",po::value(&num_blocks)->default_value(num_blocks),"requested contiguous blocks per rank")
      ("samplingStrategy",po::value(&sampling_strategy)->default_value(sampling_strategy),"independent or pcn")
      ("jointMatsubaraPCN",po::value(&joint)->default_value(false)->implicit_value(true),"pCN alias on current joint Matsubara row space")
      ("mhStepSize",po::value(&mh_step_size)->default_value(mh_step_size),"pCN proposal parameter in (0,1]")
      ("mhBurnIn",po::value(&mh_burn_in)->default_value(mh_burn_in),"pCN burn-in after each covariance update")
      ("seed",po::value(&seed)->default_value(seed),"random or unsigned integer seed")
      ("reliterror",po::value(&iteration_error_sigma_threshold)->default_value(iteration_error_sigma_threshold),"iteration residual threshold in standard-error units; sets conservative preparation scale")
      ("itermode",po::value(&iteration_error_mode)->default_value(iteration_error_mode),"relative (Keldysh) or absolute")
      ("absiterror",po::value(&absolute_iteration_error_option)->default_value(absolute_iteration_error_option),"conservative or numeric absolute threshold")
      ("iterlimit",po::value(&Iteration_Limit)->default_value(Iteration_Limit),"real-time iteration limit")
      ("preparationIterlimit",po::value(&preparation_limit)->default_value(0),"preparation limit; 0 uses iterlimit")
      ("constantMagnetization","use prepared m(0) in real-time feedback; retain measured m(t)")
      ("covarianceTolerance",po::value(&covariance_tolerance)->default_value(covariance_tolerance),"transpose symmetry tolerance")
      ("takagiTolerance",po::value(&takagi_tolerance)->default_value(takagi_tolerance),"Takagi reconstruction tolerance")
      ("minimumPhaseMagnitude",po::value(&minimum_phase_magnitude)->default_value(minimum_phase_magnitude),"minimum resolved |sum Z|/sum |Z|")
      ("initmode",po::value(&init_corr_mode)->default_value(init_corr_mode),"generate or import")
      ("initdcorr",po::value(&init_diag_corr)->default_value(init_diag_corr),"imagtime, const, exponential, Gaussian, cos, isolated or zero")
      ("initndcorr",po::value(&init_nondiag_corr)->default_value(init_nondiag_corr),"zero, sin or isolated")
      ("corrperiods",po::value(&correlation_periods)->default_value(correlation_periods),"number of initial cosine/sine periods on [0,beta]")
      ("loadinit","import initial preparation correlations")
      ("initcorrfile",po::value(&initial_correlations_src_file),"initial correlation file; [auto] uses the generated destination name")
      ("impcorrfile",po::value<std::string>(),"CspinDMFT alias for initcorrfile")
      ("extrapolate","interpolate imported imaginary correlations to this beta grid");
    storing.add_options()
      ("project",po::value(&project_name)->default_value(""),"configuration source and default data project")
      ("dstproject",po::value(&dst_project_name),"destination project; defaults to project")
      ("fileext",po::value(&filename_extension)->default_value(""),"filename suffix")
      ("info",po::value(&information_text)->default_value(""),"description stored in parameters")
      ("numPrintDigits",po::value(&num_PrintDigits)->default_value(num_PrintDigits),"numeric filename precision")
      ("runID",po::value(&run_id)->default_value(run_id),"stop after an iteration by creating terminate_<runID>");
    all.add(help).add(physics).add(numerics).add(storing);
    po::variables_map vm; po::store(po::parse_command_line(argc,argv,all),vm); po::notify(vm);
    if(vm.count("help") || vm.count("helpNum")) {
        std::ostringstream message;
        message<<(vm.count("help")?all:numerics)<<'\n';
        Print_Routines::print_R0(rank,message.str());
        MPI_Finalize();
        std::exit(EXIT_SUCCESS);
    }
    if(vm.count("spinmodel")) {
        if(!vm["spinspinmodel"].defaulted()) throw std::invalid_argument("choose spinmodel or spinspinmodel");
        spinspin=vm["spinmodel"].as<std::string>();
    }
    if(vm.count("stype")) {
        if(!vm["cstype"].defaulted()) throw std::invalid_argument("choose stype or cstype");
        correlation_symmetry_type=vm["stype"].as<char>();
    }
    for(const auto& s:{spinspin,spinmf=="auto"?spinspin:spinmf})
        if(!one_of(s,{"ISO","Ising","XXZ","DRF","XY","het"})) throw std::invalid_argument("unknown spin model: "+s);
    if(vm.count("unifield") && vm.count("stagfield")) throw std::invalid_argument("choose unifield or stagfield");
    if((vm.count("unifield") || vm.count("stagfield")) && chemical_shift!="none")
        throw std::invalid_argument("choose unifield, stagfield or chemshift");
    if(vm.count("unifield")) uniform_field=vm["unifield"].as<double>();
    if(vm.count("stagfield")) staggered_field=vm["stagfield"].as<double>();
    if((uniform_field && !std::isfinite(*uniform_field)) || (staggered_field && !std::isfinite(*staggered_field)))
        throw std::invalid_argument("unifield and stagfield must be finite");
    if(!one_of(noise_name,{"none","z"}) || noise_variance<0.)
        throw std::invalid_argument("invalid static noise");
    for(double x:{lambda,rho,smflambda,smfrho,noise_variance})
        if(!std::isfinite(x)) throw std::invalid_argument("nonfinite physics parameter");
    spinspin_cmodel=Physics::SpinModel(spinspin,lambda,rho);
    spinmf_cmodel=spinmf=="auto"?spinspin_cmodel:Physics::SpinModel(spinmf,smflambda,smfrho);
    noise=Physics::Noise(noise_name,noise_variance);
    self_consistency=!vm.count("noselfcons") && !uses_harmonic_bath();
    constant_magnetization_time=vm.count("constantMagnetization");
    if(joint) {
        if(!vm["samplingStrategy"].defaulted() && sampling_strategy!="pcn") throw std::invalid_argument("conflicting pCN options");
        sampling_strategy="pcn";
    } else if(!vm["jointMatsubaraPCN"].defaulted() && sampling_strategy=="pcn") throw std::invalid_argument("conflicting pCN options");
    num_PreparationSamplesPerCore=preparation_samples?preparation_samples:num_SamplesPerCore;
    Preparation_Iteration_Limit=preparation_limit?preparation_limit:Iteration_Limit;
    if(!vm.count("dstproject")) dst_project_name=project_name;
    if(vm.count("impcorrfile")) {
        if(vm.count("initcorrfile")) throw std::invalid_argument("choose impcorrfile or initcorrfile");
        initial_correlations_src_file=vm["impcorrfile"].as<std::string>();
    }
    if(vm.count("loadinit")) init_corr_mode="import";
    if(init_diag_corr=="constant") init_diag_corr="const";
    if(init_diag_corr=="gaussian") init_diag_corr="Gaussian";
    extrapolate_initial_spin_correlations=vm.count("extrapolate");
    if(seed=="random") {
        if(rank==0) { std::random_device rd; generated_seed=(static_cast<unsigned long long>(rd())<<32)^rd(); }
    } else {
        if(seed.empty() || seed.find_first_not_of("0123456789")!=std::string::npos) throw std::invalid_argument("seed must be random or an unsigned integer");
        generated_seed=std::stoull(seed);
    }
    MPI_Bcast(&generated_seed,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
    if(rank==0) Run_ID=run_id=="random"?int(generated_seed%10000+1):std::stoi(run_id);
    MPI_Bcast(&Run_ID,1,MPI_INT,0,MPI_COMM_WORLD);
    if(absolute_iteration_error_option!="conservative") {
        size_t end=0; absolute_iteration_error_tolerance=std::stod(absolute_iteration_error_option,&end);
        if(end!=absolute_iteration_error_option.size() || !(absolute_iteration_error_tolerance>0)) throw std::invalid_argument("invalid absiterror");
    }
    validate(); delta_real_t=Tmax/double(num_RealTimeSteps);
    {
        using namespace HDF5_Utilities;
        Handle configuration(H5Fopen(configuration_path().c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
        const double sites=read(configuration,"num_Spins",{},true)[0];
        if(sites<1 || sites>12 || sites!=std::round(sites))
            throw std::invalid_argument("num_Spins must be an integer in [1,12]");
        num_Spins=static_cast<size_t>(sites);
        if(uniform_field || staggered_field) {
            std::vector<int> color(num_Spins,0);
            if(staggered_field) {
                const auto couplings=read(configuration,"spin-spin couplings",{num_Spins,num_Spins});
                std::fill(color.begin(),color.end(),-1);
                std::queue<size_t> pending;
                color[0]=0; pending.push(0);
                while(!pending.empty()) {
                    const size_t i=pending.front(); pending.pop();
                    for(size_t j=0;j<num_Spins;++j) if(std::abs(couplings[i*num_Spins+j])>1e-12) {
                        if(color[j]<0) { color[j]=1-color[i]; pending.push(j); }
                        else if(color[j]==color[i]) throw std::invalid_argument("stagfield requires a bipartite coupling graph");
                    }
                }
                for(int sublattice:color) if(sublattice<0)
                    throw std::invalid_argument("stagfield requires a connected coupling graph");
            }
            const double h=uniform_field?*uniform_field:*staggered_field;
            std::ostringstream fields;
            fields<<std::setprecision(std::numeric_limits<double>::max_digits10);
            for(size_t i=0;i<num_Spins;++i) {
                if(i) fields<<',';
                fields<<(color[i]==0?h:-h);
            }
            chemical_shift=fields.str();
        }
    }
    preparation_layout={num_TimeSteps+1,0,3*num_Spins};
    real_time_layout={num_TimeSteps+1,num_RealTimeSteps*real_time_substeps,3*num_Spins};
    preparation_layout.dimension(); real_time_layout.dimension();
    if(init_corr_mode=="import" && (initial_correlations_src_file.empty() || initial_correlations_src_file=="[auto]"))
        initial_correlations_src_file=(std::filesystem::path(output_directory())/(create_filename_string()+".hdf5")).string();
}
void ParameterSpace::validate() const {
    if(config_file.empty() || mf_model_name!="CorrelationReplica") throw std::invalid_argument("config is required; mfmodel must be CorrelationReplica");
    if(num_TimeSteps<3 || !num_RealTimeSteps || !real_time_substeps || !(Tmax>0) || !(beta>0)
       || num_SamplesPerCore<4 || num_PreparationSamplesPerCore<4 || num_blocks<2 || !Iteration_Limit || !Preparation_Iteration_Limit || !num_SamplesPerSet)
        throw std::invalid_argument("invalid grids, samples or iteration limits (N_tau>=3, samples>=4, blocks>=2)");
    if(!std::isfinite(beta+Tmax+interaction_strength+rescale_meanfield+rescale_spinspin+mh_step_size+
                      iteration_error_sigma_threshold+covariance_tolerance+takagi_tolerance+minimum_phase_magnitude+absolute_iteration_error_tolerance))
        throw std::invalid_argument("nonfinite numerical parameter");
    if(num_SamplesPerCore>std::numeric_limits<size_t>::max()/num_Cores || num_PreparationSamplesPerCore>std::numeric_limits<size_t>::max()/num_Cores)
        throw std::invalid_argument("total sample count overflows size_t");
    if(num_RealTimeSteps>(std::numeric_limits<size_t>::max()/6-1)/real_time_substeps || num_TimeSteps>std::numeric_limits<size_t>::max()/6-1)
        throw std::invalid_argument("contour grid overflow");
    if(!one_of(sampling_strategy,{"independent","pcn"}) || !(mh_step_size>0 && mh_step_size<=1)) throw std::invalid_argument("invalid sampling strategy or mhStepSize");
    if(!one_of(iteration_error_mode,{"absolute","relative"}) || iteration_error_sigma_threshold<0
       || covariance_tolerance<0 || !(takagi_tolerance>0) || minimum_phase_magnitude<0 || minimum_phase_magnitude>1)
        throw std::invalid_argument("invalid convergence or diagnostic tolerance");
    if(!time_translation_invariant || real_time_evolution!="staggered-strang") throw std::invalid_argument("this cluster solver supports equilibrium staggered-strang only");
    if(std::string("ABCD").find(correlation_symmetry_type)==std::string::npos) throw std::invalid_argument("cstype must be A, B, C or D");
    if(!one_of(extra_interaction,{"none","quadrupolar"}) || !one_of(bath,{"selfconsistent","harmonic"}) || !(bath_frequency>0) || bath_coupling<0
       || !std::isfinite(bath_frequency+bath_coupling) || std::string("xyz").find(bath_component)==std::string::npos)
        throw std::invalid_argument("invalid local interaction or bath");
    if(!one_of(init_corr_mode,{"generate","import"}) || !one_of(init_diag_corr,{"imagtime","zero","const","exponential","Gaussian","cos","isolated"}) || !one_of(init_nondiag_corr,{"zero","sin","isolated"}))
        throw std::invalid_argument("invalid initialization option");
    if(!(correlation_periods>0) || !std::isfinite(correlation_periods)) throw std::invalid_argument("corrperiods must be positive and finite");
    if(num_PrintDigits<1 || num_PrintDigits>12) throw std::invalid_argument("numPrintDigits must be in [1,12]");
    relative_path(project_name,"project"); relative_path(dst_project_name,"dstproject");
    if(filename_extension.find_first_of("/\\")!=std::string::npos) throw std::invalid_argument("fileext must not contain a path separator");
}
std::string ParameterSpace::configuration_path() const {
    namespace fs=std::filesystem;
    fs::path direct(config_file);
    if(fs::is_regular_file(direct)) return direct.string();
    if(direct.is_absolute()) return direct.string();
    fs::path path=fs::path("Configuration_Data")/project_name/direct;
    if(path.extension()!=".hdf5" && path.extension()!=".h5") path += ".hdf5";
    return path.string();
}
std::string ParameterSpace::output_directory() const { return (std::filesystem::path("Data")/dst_project_name).string(); }
std::string ParameterSpace::create_essentials_string(const Mean_Field_Models::Model& model) const {
    constexpr size_t width=35;
    const auto number=[&](double value) {
        return Print_Routines::remove_zeros(Print_Routines::round_value_to_string(value,num_PrintDigits));
    };
    std::ostringstream output;
    const auto line=[&](const std::string& label,const std::string& value) {
        output<<Print_Routines::quantity_to_output_line(width,label,value);
    };
    output<<"\033[1;32mID of this Run: "<<Run_ID<<"\033[0m\n";
    output<<"The following correlations enter the self consistency (numbering starts at 1): ";
    for(size_t n=0;n<model.correlation_categories.size();++n) {
        const auto ij=model.correlation_categories[n];
        if(n) output<<',';
        output<<ij[0]+1<<'-'<<ij[1]+1;
    }
    output<<'\n';
    if(!self_consistency) output<<"I don't perform self-consistent iterations in this simulation!\n";
    if(init_corr_mode=="import")
        output<<"The initial correlations and magnetizations are imported from '"<<initial_correlations_src_file<<"'\n";
    else {
        line("init_diag_corr",init_diag_corr);
        line("init_nondiag_corr",init_nondiag_corr);
    }
    line("mf_model",mf_model_name);
    line("config_file",config_file);
    line("num_Spins",std::to_string(model.sites));
    line("num_HilbertSpaceDimension",std::to_string(model.identity.rows()));
    line("uncoupled_spins",rescale_spinspin==0. || std::all_of(model.J.begin(),model.J.end(),[](double j) { return j==0.; })?"true":"false");
    std::vector<std::string> spins;
    for(double spin:model.spins) spins.push_back(Physics::Spin::convert_spin_float_to_string(spin));
    line("spin_list",Print_Routines::concatenate_string_with_delimiter(spins,','));
    if(model.sites<5) {
        line("J",std::to_string(model.sites)+"x"+std::to_string(model.sites)+"-Matrix");
        blaze::DynamicMatrix<double,blaze::rowMajor> J(model.sites,model.sites);
        for(size_t i=0;i<model.sites;++i) for(size_t j=0;j<model.sites;++j) J(i,j)=model.J[i*model.sites+j];
        output<<std::setprecision(num_PrintDigits)<<J;
    } else line("J(0,1)",number(model.J[1]));
    line("spinspin_cmodel info",spinspin_cmodel.compact_info(num_PrintDigits));
    line("spinmf_cmodel info",spinmf_cmodel.compact_info(num_PrintDigits));
    line("chemical_shift info",chemical_shift);
    line("rescale_meanfield",number(rescale_meanfield));
    line("rescale_spinspin",number(rescale_spinspin));
    line("beta",number(beta));
    line("noise_string",noise.compact_info(num_PrintDigits));
    line("bath",bath);
    if(uses_harmonic_bath()) {
        line("bath_frequency",number(bath_frequency));
        line("bath_coupling",number(bath_coupling));
        line("bath_component",std::string(1,bath_component));
    }
    line("extra_interaction_string",extra_interaction);
    if(extra_interaction!="none") line("intstrength",number(interaction_strength));
    line("correlation_symmetry_type",std::string(1,correlation_symmetry_type));
    line("time_translation_invariant",time_translation_invariant?"true":"false");
    line("num_ImagTimeSteps",std::to_string(num_TimeSteps));
    line("delta_imag_t",number(beta/double(num_TimeSteps)));
    line("num_RealTimeSteps",std::to_string(num_RealTimeSteps));
    line("real_time_substeps",std::to_string(real_time_substeps));
    line("delta_real_propagation_t",std::to_string(delta_real_t/double(real_time_substeps)));
    line("Tmax",number(Tmax));
    line("delta_real_t",number(delta_real_t));
    line("num_Cores",std::to_string(num_Cores));
    line("num_Samples",std::to_string(num_SamplesPerCore*num_Cores));
    line("num_SamplesPerCore",std::to_string(num_SamplesPerCore));
    line("num_PreparationSamples",std::to_string(num_PreparationSamplesPerCore*num_Cores));
    line("num_PreparationSamplesPerCore",std::to_string(num_PreparationSamplesPerCore));
    line("sampling_strategy",sampling_strategy);
    line("staged_field_sampling","full-contour Takagi");
    if(sampling_strategy=="pcn") {
        line("pcn_target","p0(a)|Z_M(L_M U a)|; fresh joint null-space draws");
        line("mh_step_size",number(mh_step_size));
        line("mh_burn_in",std::to_string(mh_burn_in));
    }
    line("iteration_error_mode",iteration_error_mode);
    line("iteration_error_sigma_threshold",number(iteration_error_sigma_threshold));
    line("Preparation_Iteration_Limit",std::to_string(Preparation_Iteration_Limit));
    line("Iteration_Limit",std::to_string(Iteration_Limit));
    line("information_text",information_text);
    return output.str();
}
std::string ParameterSpace::create_filename_string() const {
    const auto number=[&](double v) { std::ostringstream os; os<<std::setprecision(num_PrintDigits)<<v; return os.str(); };
    std::string name=spinspin_cmodel.compact_info(num_PrintDigits);
    if(spinspin_cmodel!=spinmf_cmodel) name+="__mf"+spinmf_cmodel.compact_info(num_PrintDigits);
    name+="__config="+std::filesystem::path(config_file).stem().string()+"__beta="+number(beta);
    if(uniform_field) name+="__unifield="+number(*uniform_field);
    else if(staggered_field) name+="__stagfield="+number(*staggered_field);
    else if(chemical_shift!="none") name+="__chemshift="+chemical_shift;
    if(extra_interaction!="none") name+="__extraint="+extra_interaction+"_strength="+number(interaction_strength);
    if(rescale_meanfield!=1.) name+="__rescmf="+number(rescale_meanfield);
    if(rescale_spinspin!=1.) name+="__rescss="+number(rescale_spinspin);
    if(noise.m_name!="none") name+="__"+noise.compact_info(num_PrintDigits);
    if(uses_harmonic_bath()) name+="__bath=harmonic_omega="+number(bath_frequency)+"_g="+number(bath_coupling)+"_component="+bath_component;
    if(constant_magnetization_time) name+="__mag=constant";
    if(real_time_substeps!=1) name+="__substeps="+std::to_string(real_time_substeps);
    if(!filename_extension.empty()) name+="_"+filename_extension;
    if(name.size()>220) name.resize(220);
    return name;
}
}
