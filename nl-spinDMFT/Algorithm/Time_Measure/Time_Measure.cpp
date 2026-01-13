#include"Time_Measure.h" 
#include<iomanip>
#include<algorithm>
#include<Globals/MPI_Types.h>

namespace Time_Measure
{
// ===============================================================
// ================== DERIVED TIME MEASURE CLASS =================
// ===============================================================
// constructor : start time measurement and define durations
DerivedTimeMeasure::DerivedTimeMeasure( const int my_rank, const int world_size ): 
    TimeMeasure( my_rank, world_size ) 
{
    // add global measures
    m_global_durations.emplace_back( DurationQuantity{"initialization"} );
    m_global_durations.emplace_back( DurationQuantity{"mean-field distribution"} );
    m_global_durations.emplace_back( DurationQuantity{"mean-field sampling"} );
    m_global_durations.emplace_back( DurationQuantity{"time evolution"} );
    m_global_durations.emplace_back( DurationQuantity{"Monte-Carlo simulation"} );
    m_global_durations.emplace_back( DurationQuantity{"MPI communication"} );
    m_global_durations.emplace_back( DurationQuantity{"finalization"} );
    m_global_durations.emplace_back( DurationQuantity{"storing"} );
}

// measure the duration for 'what'
DurationType DerivedTimeMeasure::internal_measure( const std::string what, const bool print )
{
    TimeType time = MPI_Wtime();
    DurationType duration{};
    if( m_loop_level == 0 ) // global duration
    {
        duration = time - m_global_measure_time;
        auto iterator = find_measure_in( what, m_global_durations );
        iterator->m_duration = duration; // save the measured duration
        m_global_measure_time = time; // update the global measure time
    }
    else // higher loop levels -> tmp duration
    {
        duration = time - m_tmp_measure_times[m_loop_level-1];
        auto iterator = find_measure_in( what, m_global_durations );
        iterator->m_duration += duration; // add up to the measured duration
        m_tmp_measure_times[m_loop_level-1] = time; // update the tmp measure time
    }
    if( print ) // print the time measurement
    {
        print_measure( what, duration );
    }
    return duration;
}

};