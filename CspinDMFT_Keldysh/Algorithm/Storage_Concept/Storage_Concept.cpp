#include "../main_header.h"
#include "HDF5_Utilities.h"
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace CspinDMFT::Storage_Concept {
namespace {
using namespace HDF5_Utilities;
void string(hid_t g,const char* name,const std::string& value) {
    Handle type(H5Tcopy(H5T_C_S1),H5Tclose); check(H5Tset_size(type,value.size()+1));
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle attribute(H5Acreate2(g,name,type,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(attribute,type,value.c_str()));
}
void scalar(hid_t g,const char* name,double value) {
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle a(H5Acreate2(g,name,H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(a,H5T_NATIVE_DOUBLE,&value));
}
void integer(hid_t g,const char* name,unsigned long long value) {
    Handle space(H5Screate(H5S_SCALAR),H5Sclose);
    Handle a(H5Acreate2(g,name,H5T_STD_U64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    check(H5Awrite(a,H5T_NATIVE_ULLONG,&value));
}
void list(hid_t g,const char* name,const std::vector<double>& values) {
    hsize_t size=values.size();
    Handle space(H5Screate_simple(1,&size,nullptr),H5Sclose);
    Handle a(H5Acreate2(g,name,H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
    if(!values.empty()) check(H5Awrite(a,H5T_NATIVE_DOUBLE,values.data()));
}
void dataset(hid_t g,const std::string& name,const std::vector<hsize_t>& shape,const std::vector<double>& values) {
    size_t n=1; for(auto d:shape) n*=d;
    if(n!=values.size()) throw std::logic_error("output shape mismatch");
    Handle space(H5Screate_simple(int(shape.size()),shape.data(),nullptr),H5Sclose);
    Handle data(H5Dcreate2(g,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    if(n) check(H5Dwrite(data,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()));
}
Handle group(hid_t file,const char* name) { return Handle(H5Gcreate2(file,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Gclose); }
template<class Correlations>
void correlations(hid_t g,const char* name,const Correlations& state,char symmetry,bool imaginary,bool preparation) {
    auto output=group(g,name); const auto dirs=direction_pairs(symmetry);
    string(output,"info",preparation?"Axes direction,tau; sites are one-based i-j.":"Axes real_time,direction,tau; sites are one-based i-j. X=<S_j^b(-i tau) S_i^a(t)>.");
    for(auto ij:state.categories) {
        const auto i=ij[0],j=ij[1];
        std::vector<double> values;
        for(size_t t=0;t<time_points(state);++t) for(auto ab:dirs) for(size_t k=0;k<state.imaginary_points;++k) {
            C v=correlation(state,t,3*i+ab[0],3*j+ab[1],k); values.push_back(imaginary?v.imag():v.real());
        }
        std::vector<hsize_t> shape{dirs.size(),state.imaginary_points};
        if(!preparation) shape.insert(shape.begin(),time_points(state));
        dataset(output,std::to_string(i+1)+"-"+std::to_string(j+1),shape,values);
    }
}
void category_dataset(hid_t g,const char* name,const std::vector<IndexPair>& categories) {
    std::vector<unsigned int> values;
    for(auto ij:categories) { values.push_back(static_cast<unsigned int>(ij[0])); values.push_back(static_cast<unsigned int>(ij[1])); }
    hsize_t shape[2]{categories.size(),2};
    Handle space(H5Screate_simple(2,shape,nullptr),H5Sclose);
    Handle data(H5Dcreate2(g,name,H5T_STD_U32LE,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    if(!values.empty()) check(H5Dwrite(data,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()));
}
template<class Magnetization>
void store_magnetization(hid_t g,const char* name,const Magnetization& moments,bool imaginary) {
    std::vector<double> values;
    for(size_t t=0;t<time_points(moments);++t) for(size_t a=0;a<3*moments.sites;++a) {
        C value=magnetization(moments,t,a); values.push_back(imaginary?value.imag():value.real());
    }
    dataset(g,name,{time_points(moments),moments.sites,3},values);
}
template<class Correlations>
void errors(hid_t group_id,const BasicStageData<Correlations>& e,char symmetry,bool prep) {
    correlations(group_id,"Re_correlation_sample_stds",e.correlation_error_re,symmetry,false,prep);
    correlations(group_id,"Im_correlation_sample_stds",e.correlation_error_im,symmetry,false,prep);
    store_magnetization(group_id,"Re_magnetization_sample_stds",e.magnetization_error_re,false);
    store_magnetization(group_id,"Im_magnetization_sample_stds",e.magnetization_error_im,false);
}
void history(hid_t g,const Run_Time_Data::StageData& stage) {
    string(g,"termination",stage.termination); integer(g,"num_Iterations",stage.num_Iterations);
    list(g,"absolute_iteration_errors",stage.absolute_iteration_errors); list(g,"standardized_iteration_errors",stage.standardized_iteration_errors);
    list(g,"average_phase_magnitudes",stage.average_phase_magnitudes); list(g,"gaussian_factor_reconstruction_errors",stage.gaussian_factor_reconstruction_errors);
    list(g,"covariance_symmetry_errors",stage.covariance_symmetry_errors); list(g,"partition_sum_Re",stage.partition_sum_Re); list(g,"partition_sum_Im",stage.partition_sum_Im);
    list(g,"effective_sample_sizes",stage.effective_sample_sizes); list(g,"mh_acceptance_rates",stage.mh_acceptance_rates);
}
void labels(hid_t g,char symmetry) {
    Handle type(H5Tcopy(H5T_C_S1),H5Tclose); check(H5Tset_size(type,3));
    auto dirs=direction_pairs(symmetry); hsize_t size=dirs.size();
    Handle space(H5Screate_simple(1,&size,nullptr),H5Sclose);
    Handle d(H5Dcreate2(g,"correlation_direction_labels",type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    std::vector<char> data(3*size,0);
    for(size_t i=0;i<size;++i) { data[3*i]="xyz"[dirs[i][0]]; data[3*i+1]="xyz"[dirs[i][1]]; }
    check(H5Dwrite(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data.data()));
}
}
HDF5_Storage::HDF5_Storage(int rank,const ParameterSpace& p,const std::string& termination,const std::string& suffix) {
    if(rank!=0) return;
    namespace fs=std::filesystem;
    fs::path folder=p.output_directory();
    if(termination=="by iteration limit" || termination=="by hand") folder/="NOT_CONVERGED";
    fs::create_directories(folder);
    std::string base=p.create_filename_string()+suffix;
    for(size_t tries=0;tries<100;++tries) {
        fs::path candidate=folder/(base+".hdf5");
        if(fs::exists(candidate)) { base+="X"; continue; }
        hid_t file;
        H5E_BEGIN_TRY { file=H5Fcreate(candidate.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT); } H5E_END_TRY;
        if(file>=0) {
            m_file=file; m_filename=candidate.string();
            if(tries>0) print::print_R0(rank,"\033[1;31mdata file already exists!\nsuccessfully created "+m_filename+" instead\033[0m\n");
            break;
        }
        if(!fs::exists(candidate)) throw std::runtime_error("could not create HDF5 file: "+candidate.string());
        base+="X";
    }
    if(m_file<0) throw std::runtime_error("could not allocate unique output filename");
}
HDF5_Storage::~HDF5_Storage() { if(m_file>=0) H5Fclose(m_file); }
void HDF5_Storage::finalize() { if(m_file>=0) { HDF5_Utilities::check(H5Fclose(m_file)); m_file=-1; } }
template<class Correlations>
void HDF5_Storage::store_main(const ParameterSpace& p,const Model& m,const Run_Time_Data::RunTimeData& runtime,const Correlations& values,const MagnetizationFor<Correlations>& moments) {
    if(m_file<0) return;
    const auto& e=runtime.stage<Correlations>();
    auto ps=group(m_file,"parameters");
    string(ps,"mf_model",p.mf_model_name); string(ps,"config_file",p.config_file);
    string(ps,"configuration_src_file",std::filesystem::path(p.configuration_path()).filename().string()); integer(ps,"num_Spins",m.sites); integer(ps,"num_HilbertSpaceDimension",m.identity.rows());
    string(ps,"spinspin_cmodel",p.spinspin_cmodel.m_name); string(ps,"spinmf_cmodel",p.spinmf_cmodel.m_name);
    for(auto param:p.spinspin_cmodel.m_parameters) scalar(ps,param.m_name.c_str(),param.m_value);
    for(auto param:p.spinmf_cmodel.m_parameters) scalar(ps,("smf_"+param.m_name).c_str(),param.m_value);
    string(ps,"chemical_shift",p.chemical_shift); string(ps,"local_extra_interaction",p.extra_interaction); scalar(ps,"intstrength",p.interaction_strength);
    if(p.uniform_field) scalar(ps,"unifield",*p.uniform_field);
    if(p.staggered_field) scalar(ps,"stagfield",*p.staggered_field);
    string(ps,"noise",p.noise.m_name); for(auto param:p.noise.m_parameters) scalar(ps,param.m_name.c_str(),param.m_value);
    scalar(ps,"beta",p.beta); scalar(ps,"rescale_meanfield",p.rescale_meanfield); scalar(ps,"rescale_spinspin",p.rescale_spinspin);
    string(ps,"bath",p.bath); scalar(ps,"bath_frequency",p.bath_frequency); scalar(ps,"bath_coupling",p.bath_coupling); string(ps,"bath_component",std::string(1,p.bath_component));
    integer(ps,"num_ImagTimeSteps",p.num_TimeSteps); integer(ps,"num_ImagTimePoints",p.num_TimeSteps+1); scalar(ps,"delta_imag_t",p.beta/double(p.num_TimeSteps));
    integer(ps,"num_RealTimeSteps",p.num_RealTimeSteps); integer(ps,"num_RealTimePoints",p.num_RealTimeSteps+1); scalar(ps,"delta_real_t",p.delta_real_t); scalar(ps,"Tmax",p.Tmax);
    string(ps,"correlation_symmetry_type",std::string(1,p.correlation_symmetry_type)); string(ps,"time_translation_invariant","yes");
    string(ps,"constant_magnetization_time",p.constant_magnetization_time?"yes":"no"); string(ps,"self_consistency",p.self_consistency?"yes":"no");
    integer(ps,"real_time_substeps",p.real_time_substeps); scalar(ps,"delta_real_propagation_t",p.delta_real_t/double(p.real_time_substeps));
    integer(ps,"num_SamplesPerCore",p.num_SamplesPerCore); integer(ps,"num_Samples",p.num_SamplesPerCore*p.get_num_Cores());
    integer(ps,"num_PreparationSamplesPerCore",p.num_PreparationSamplesPerCore); integer(ps,"num_PreparationSamples",p.num_PreparationSamplesPerCore*p.get_num_Cores());
    integer(ps,"num_Cores",p.get_num_Cores()); integer(ps,"num_SamplesPerSet",p.num_SamplesPerSet); integer(ps,"num_blocks",p.num_blocks);
    string(ps,"seed",p.seed); string(ps,"sampling_strategy",p.sampling_strategy); scalar(ps,"mh_step_size",p.mh_step_size); integer(ps,"mh_burn_in",p.mh_burn_in);
    if(p.sampling_strategy=="pcn") string(ps,"pcn_updated_coordinates","current joint Takagi Matsubara row space; fresh null-space real latents each production step");
    string(ps,"adaptive_sample_size","no");
    integer(ps,"Iteration_Limit",p.Iteration_Limit); integer(ps,"Preparation_Iteration_Limit",p.Preparation_Iteration_Limit);
    string(ps,"iteration_error_mode",p.iteration_error_mode); scalar(ps,"iteration_error_sigma_threshold",p.iteration_error_sigma_threshold);
    scalar(ps,"absolute_iteration_error_threshold",runtime.preparation_threshold(p,m));
    scalar(ps,"covariance_tolerance",p.covariance_tolerance); scalar(ps,"takagi_tolerance",p.takagi_tolerance); scalar(ps,"minimum_phase_magnitude",p.minimum_phase_magnitude);
    scalar(ps,"corrperiods",p.correlation_periods); string(ps,"init_corr_mode",p.init_corr_mode); string(ps,"init_diag_corr",p.init_diag_corr); string(ps,"init_nondiag_corr",p.init_nondiag_corr);
    string(ps,"initial_correlations_src_file",p.initial_correlations_src_file); string(ps,"extrapolate_initial_spin_correlations",p.extrapolate_initial_spin_correlations?"yes":"no");
    string(ps,"information_text",p.information_text); string(ps,"project_name",p.project_name); string(ps,"destination project_name",p.dst_project_name);
    string(ps,"filename_extension",p.filename_extension); string(ps,"filename",std::filesystem::path(m_filename).filename().string()); string(ps,"RealType","double"); integer(ps,"Run_ID",p.Run_ID);
    list(ps,"spin lengths",m.spins);
    category_dataset(ps,"correlation_categories, stored according to the hierarchy n, i-j",m.correlation_categories);
    dataset(ps,"spin-spin couplings",{m.sites,m.sites},m.J); dataset(ps,"expectation weights",{m.sites,m.sites},m.expectation_weights);
    dataset(ps,"correlation weights",{m.sites,m.sites,m.sites,m.sites},m.correlation_weights);
    auto rtd=group(m_file,"runtimedata");
    history(rtd,runtime.real_time); list(rtd,"num_samples_per_core_history",runtime.real_time.sample_counts); integer(rtd,"generated_seed",runtime.generated_seed);
    integer(rtd,"actual_num_blocks_per_core",e.actual_num_blocks); integer(rtd,"samples_per_block",e.samples_per_block);
    list(rtd,"blocking_curve_block_length",e.blocking_lengths); list(rtd,"blocking_curve_max_error",e.blocking_errors);
    std::vector<double> re,im;
    for(C z:e.closure) { re.push_back(z.real()); im.push_back(z.imag()); }
    list(rtd,"closed_contour_ratio_Re",re); list(rtd,"closed_contour_ratio_Im",im);
    errors(rtd,e,p.correlation_symmetry_type,false);
    auto results=group(m_file,"results"); labels(results,p.correlation_symmetry_type);
    correlations(results,"Re_correlation",values,p.correlation_symmetry_type,false,false);
    correlations(results,"Im_correlation",values,p.correlation_symmetry_type,true,false);
    store_magnetization(results,"Re_magnetization",moments,false); store_magnetization(results,"Im_magnetization",moments,true);
}
void HDF5_Storage::store_preparation(const ParameterSpace& p,const Model& m,const Run_Time_Data::RunTimeData& rt,const MatsubaraCorrelations& values,const ClusterMagnetization& moments) {
    if(m_file<0) return;
    const auto& e=rt.preparation;
    auto prep=group(m_file,"preparation"); history(prep,rt.preparation);
    scalar(prep,"absolute_iteration_error_threshold",rt.preparation_threshold(p,m));
    integer(prep,"samples_per_core",e.samples_per_core); integer(prep,"iteration_limit",p.Preparation_Iteration_Limit);
    string(prep,"sampling_strategy",p.sampling_strategy); integer(prep,"actual_num_blocks_per_core",e.actual_num_blocks); integer(prep,"samples_per_block",e.samples_per_block);
    correlations(prep,"Re_correlation",values,p.correlation_symmetry_type,false,true);
    correlations(prep,"Im_correlation",values,p.correlation_symmetry_type,true,true);
    store_magnetization(prep,"Re_magnetization",moments,false); store_magnetization(prep,"Im_magnetization",moments,true); errors(prep,e,p.correlation_symmetry_type,true);
}
void HDF5_Storage::store_time(const Time_Measure::DerivedTimeMeasure& time) {
    if(m_file<0) return;
    auto g=group(m_file,"timedata"); string(g,"info","all durations are measured in seconds"); scalar(g,"total_duration",time.total_duration);
    list(g,"preparation iteration",time.preparation_iterations); list(g,"real-time iteration",time.real_time_iterations);
    string(g,"start_date_and_time",std::string(std::ctime(&time.start_date_and_time)));
}
}

namespace CspinDMFT::Storage_Concept {
template void HDF5_Storage::store_main(const ParameterSpace&,const Model&,const Run_Time_Data::RunTimeData&,const MatsubaraCorrelations&,const ClusterMagnetization&);
template void HDF5_Storage::store_main(const ParameterSpace&,const Model&,const Run_Time_Data::RunTimeData&,const ContourCorrelations&,const TimeMagnetization&);
}
