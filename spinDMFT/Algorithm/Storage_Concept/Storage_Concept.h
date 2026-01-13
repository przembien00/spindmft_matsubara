#pragma once

#include<HDF5/HDF5_Routines.h>
#ifdef USE_DOUBLE
static hid_t H5_REAL_TYPE = H5T_IEEE_F64LE;
#endif 
#ifdef USE_FLOAT 
static hid_t H5_REAL_TYPE = H5T_IEEE_F32LE;
#endif
#include<Globals/Types.h>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include"../Parameter_Space/Parameter_Space.h"
#include"../Time_Measure/Time_Measure.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace spinDMFT::Storage_Concept
{

namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace tmm = Time_Measure;
namespace ps = spinDMFT::Parameter_Space;
namespace rtd = spinDMFT::Run_Time_Data;

using CorrTen = ten::CorrelationTensor<corr::CorrelationVector>;

// =================================================================
// ================= HEADER FOR HDF5 STORAGE CLASS =================
// =================================================================
class HDF5_Storage
{
 public:
    HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace, const std::string& termination );
    void store_main( const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata, const CorrTen& corr_R, const CorrTen& corr_I, const FieldVector& spin_expval );
    void store_time( const tmm::DerivedTimeMeasure& tmeasure );
    void finalize();

 private:
    const bool m_storing_permission{true}; // permission for the class instance to create folders/files and store data
    const size_t m_fname_max_length{};
    const size_t m_num_TriesToBuildFile{};
    hid_t m_file_id;
    std::string m_filename{};

    void create_folder_branch( const ps::ParameterSpace& pspace, const std::string& termination );
    void create_file( const ps::ParameterSpace& pspace );
    void store_correlation_tensor( const CorrTen& CT, const hid_t group_id, const std::string dataset_name, const std::string dataset_info );
};

};