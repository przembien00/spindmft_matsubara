#include"Time_Measure/Time_Measure.h"

#pragma once

namespace Time_Measure
{

// ================= CLASS DEFINITIONS =================
// SPIN ED TIME MEASURE CLASS
/* derive a time measure class from the time measure base class and define an internal measure */
class SpinEDTimeMeasure : public TimeMeasure
{
 public: 
    // CONSTRUCTORS
    SpinEDTimeMeasure( const uint my_rank_, const uint world_size_ );

    // PUBLIC METHODS
    void finalize();

    // PUBLIC MEMBERS
    std::vector<DurationQuantity> m_tmp_measures{};

 private:
    DurationType internal_measure( const std::string what, const bool print ) override;
};

}
