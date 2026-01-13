#include"Time_Measure/Time_Measure.h"

#pragma once

namespace Time_Measure
{
// ===============================================================
// ============ HEADER FOR DERIVED TIME MEASURE CLASS ============
// ===============================================================
// derive a time measure class from the time measure base class and define an internal measure
class DerivedTimeMeasure : public TimeMeasure
{
 public: 
    // CONSTRUCTORS
    DerivedTimeMeasure( const int my_rank, const int world_size );

 private:
    DurationType internal_measure( const std::string what, const bool print ) override;
};

}
