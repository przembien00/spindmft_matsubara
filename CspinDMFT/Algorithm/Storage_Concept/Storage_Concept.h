#pragma once

#include<Globals/Types.h>
#include<HDF5/HDF5_Routines.h>
#ifdef USE_DOUBLE
static hid_t H5_REAL_TYPE = H5T_IEEE_F64LE;
#endif 
#ifdef USE_FLOAT 
static hid_t H5_REAL_TYPE = H5T_IEEE_F32LE;
#endif
#include<memory>
#include<vector>
#include<string>
#include<Observables/Correlations.h>
#include<Observables/Tensors.h>
#include<Observables/Clusters.h>
#include"../Time_Measure/Time_Measure.h"
#include"../Parameter_Space/Parameter_Space.h"
#include"../Run_Time_Data/Run_Time_Data.h"

namespace Storage_Concept
{

namespace ps = DMFT_parameter_space;
namespace tmm = Time_Measure;
namespace rtd = Run_Time_Data;
namespace corr = Observables::Correlations;
namespace ten = Observables::Tensors;
namespace clu = Observables::Clusters;
using CorrTen = ten::CorrelationTensor<corr::CorrelationVector>;
using CluCorrTen = clu::CorrelationCluster<CorrTen>;

// ==================================================================
// ======================= HDF5 STORAGE CLASS =======================
// ==================================================================
/* stores all the different data and parameters produced by the CspinDMFT algorithm */
class HDF5_Storage
{ 
 public:
    HDF5_Storage( const int my_rank, const ps::ParameterSpace& pspace, const rtd::RunTimeData& rtdata );

    void store_main( const ps::ParameterSpace& pspace, rtd::RunTimeData& rtdata, CluCorrTen& CCT );
    void store_time( const tmm::DerivedTimeMeasure& tmeasure );
    void finalize();

 private:
    const bool m_storing_permission{true}; // permission for the class instance to create folders/files and store data
    const size_t m_fname_max_length{};
    const size_t m_num_TriesToBuildFile{};
    hid_t m_file_id;
    std::string m_filename{};

    void create_folder_branch( const ps::ParameterSpace& pspace, const std::string& termination );
    void create_file( const ps::ParameterSpace& pspace, const bool eigenvalue_threshold_violated );
    void store_correlation_tensor_cluster( CluCorrTen& data, const hid_t data_group_id, const std::string dataset_name, const std::string dataset_info );
};

};



