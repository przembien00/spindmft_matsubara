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
    m_global_durations.emplace_back( DurationQuantity{"self consistency"} );
    m_global_durations.emplace_back( DurationQuantity{"storing"} );

    // generate temporaries (for duration measurements in deeper loops)
    m_tmp_measures.emplace_back( IterationDurationQuantity{"mean-field distribution"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"mean-field sampling"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"time propagation"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"expectation values"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"Monte-Carlo simulation"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"MPI communication"} );
    m_tmp_measures.emplace_back( IterationDurationQuantity{"iteration-step finalization"} );  
}

// measure the duration for 'what' internally
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
        auto iterator = find_measure_in( what, m_tmp_measures );
        *(--(iterator->m_durations.end())) += duration; // add up to the iteration duration
        m_tmp_measure_times[m_loop_level-1] = time; // update the tmp measure time
    }
    if( print ) // print the time measurement
    {
        print_measure( what, duration );
    }
    return duration;
}

// proceed to the next iteration step
void DerivedTimeMeasure::new_iteration()
{
    std::for_each( m_tmp_measures.begin(), m_tmp_measures.end(), []( auto& q ){ q.m_durations.emplace_back(0.); } );
}

};