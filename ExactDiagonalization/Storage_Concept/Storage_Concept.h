#pragma once

#include"../Global/RealType_Global.h"
#include"HDF5/HDF5_Routines.h"
#ifdef USE_DOUBLE
static hid_t H5_REAL_TYPE = H5T_IEEE_F64LE;
#endif 
#ifdef USE_FLOAT 
static hid_t H5_REAL_TYPE = H5T_IEEE_F32LE;
#endif
#include"HDF5/HDF5_Routines_RealType.h"
#include<spinDMFT_Tensors/Tensors.h>
#include"../Time_Measure/Time_Measure.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"
#include"../Models/Models.h"


namespace SpinED::Storage_Concept
{

namespace tmm = Time_Measure;
namespace ps = SpinED::Parameter_Space;
namespace rtd = SpinED::Run_Time_Data;
namespace mod = SpinED::Models;

// ================= CLASS DEFINITIONS =================
// HDF5 STORAGE CLASS
class HDF5_Storage
{
 public:
    HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace );

    void store_main( const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata, const mod::MeanTen& mean, const mod::CorrTen& corr );
    void store_time( const tmm::TimeMeasure& tmeasure );

    void finalize();

 private:
    const bool m_storing_permission{true}; // permission for the class instance to create folders/files and store data
    const size_t m_fname_max_length{};
    const size_t m_num_TriesToBuildFile{};
    hid_t m_file_id;
    hid_t m_results_group_id;
    std::string m_filename{};

    void create_folder_tree( const std::vector<std::string>& folder_branch_list );
    void create_file( std::string filename );
    void store_mean( const mod::MeanTen data );
    void store_corr( const mod::CorrTen data );
};


};