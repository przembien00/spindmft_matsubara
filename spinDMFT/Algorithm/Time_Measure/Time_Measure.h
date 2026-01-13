#include"Time_Measure/Time_Measure.h"

#pragma once

namespace Time_Measure
{

// ===============================================================
// ======== HEADER FOR ITERATION DURATION QUANTITY CLASS =========
// ===============================================================
// stores the duration of a certain process/quantity for different iteration steps
struct IterationDurationQuantity
{
  IterationDurationQuantity( const std::string& quantity_name ) : m_name( quantity_name ) {};
  DurationType average() const
  { 
    return std::accumulate(m_durations.cbegin(), m_durations.cend(), DurationType{0.}) / static_cast<DurationType>(m_durations.size()); 
  }
  std::string m_name{};                     // name of the process/quantity
  std::vector<DurationType> m_durations{};  // vector of the measured durations during the iterations
};

// ===============================================================
// ============ HEADER FOR DERIVED TIME MEASURE CLASS ============
// ===============================================================
// derive a time measure class from the time measure base class and define an internal measure
class DerivedTimeMeasure : public TimeMeasure
{
 public: 
    DerivedTimeMeasure( const int my_rank, const int world_size );
    void new_iteration();

    std::vector<IterationDurationQuantity> m_tmp_measures{};

 private:
    DurationType internal_measure( const std::string what, const bool print ) override;
};

}
