#pragma once
#include "../Run_Time_Data/Run_Time_Data.h"
#include "../Time_Measure/Time_Measure.h"
namespace CspinDMFT::Storage_Concept {
class HDF5_Storage {
public:
    HDF5_Storage(int rank,const ParameterSpace&,const std::string& termination,const std::string& suffix="");
    template<class Correlations>
    void store_main(const ParameterSpace&,const Model&,const Run_Time_Data::RunTimeData&,const Correlations&,const MagnetizationFor<Correlations>&);
    void store_preparation(const ParameterSpace&,const Model&,const Run_Time_Data::RunTimeData&,const MatsubaraCorrelations&,const ClusterMagnetization&);
    void store_time(const Time_Measure::DerivedTimeMeasure&);
    void finalize();
    ~HDF5_Storage();
    const std::string& filename() const { return m_filename; }
    HDF5_Storage(const HDF5_Storage&)=delete;
private:
    long long m_file=-1;
    std::string m_filename;
};
}
