#include"Time_Measure.h" 

#include<iomanip>
#include<algorithm>
#include"../Global/MPI_Global.h"

namespace Time_Measure
{

// MPI_Datatype MPI_TIMETYPE = MPI_DOUBLE;
namespace error = Error_Handling;
namespace print = Print_Routines;

// ================= CLASS IMPLEMENTATIONS =================
// CLASS SPIN ED TIME MEASURE
// constructor : start time measurement
SpinEDTimeMeasure::SpinEDTimeMeasure( const uint my_rank_, const uint world_size_ ): 
    TimeMeasure( my_rank_, world_size_ ) 
{
    // add global measures
    m_global_durations.emplace_back( DurationQuantity{"construction"} );
    m_global_durations.emplace_back( DurationQuantity{"diagonalization"} );
    m_global_durations.emplace_back( DurationQuantity{"transformation"} );
    m_global_durations.emplace_back( DurationQuantity{"autocorrelations"} );
    m_global_durations.emplace_back( DurationQuantity{"saving"} );

    // generate temporaries (for duration measurements in deeper loops)
    m_tmp_measures.emplace_back( DurationQuantity{"quantities_at"} );
}

// measure the duration for 'what' internally
DurationType SpinEDTimeMeasure::internal_measure( const std::string what, const bool print )
{
    TimeType time = MPI_Wtime();
    DurationType duration{};
    if( m_loop_level == 0 ) // global duration
    {
        duration = time - m_global_measure_time;
        auto iterator = find_measure_in( what, m_global_durations );
        iterator->m_duration = duration;
        m_global_measure_time = time; // update the global measure time
    }
    else // time loop level -> tmp duration
    {
        duration = time - m_tmp_measure_times[0];
        auto iterator = find_measure_in( what, m_tmp_measures );
        iterator->m_duration += duration; // add up to the corresponding duration
        m_tmp_measure_times[0] = time; // update the tmp measure time
    }
    if( print ) // print the time measurement
    {
        print_measure( what, duration );
    }
    return duration;
}

// some finalization
void SpinEDTimeMeasure::finalize()
{   
    // do smt, e.g. core averaging
}


};