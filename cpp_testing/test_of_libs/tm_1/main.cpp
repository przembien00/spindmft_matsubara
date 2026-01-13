using RealType = double;

#include<iostream>
#include"Time_Measure.h"

namespace Time_Measure 
{

/* derive a time measure class from the time measure base class and define internal measure
what should happen in which loop level? */
class TestTimeMeasure : public TimeMeasure
{
 public: 
    TestTimeMeasure( const uint my_rank_, const uint world_size_ ): TimeMeasure( my_rank_, world_size_ ) 
    {
        // add global measures
        m_global_durations.emplace_back( DurationQuantity{"preparation"} );
        m_global_durations.emplace_back( DurationQuantity{"simulation"} );
        m_global_durations.emplace_back( DurationQuantity{"finalizing"} );
    }

 private:
    DurationType internal_measure( const std::string what, const bool print ) override
    {
        TimeType time = MPI_Wtime();
        DurationType duration{};
        if( m_loop_counter == 0 ) // global duration
        {
            duration = time - m_global_measure_time;
            auto iterator = find_measure_in( what, m_global_durations );
            iterator->m_duration = duration;
            m_global_measure_time = MPI_Wtime(); // update the global measure time
        }
        else if( m_loop_counter == 1 ) // first loop level
        {
            duration = time - m_tmp_measure_times[0];
            m_tmp_measure_times[0] = MPI_Wtime(); // update the tmp measure time
        }
        if( print ) // print the time measurement
        {
            print_measure( what, duration );
        }
        return duration;
    }

};

};

namespace tmm = Time_Measure;

int main() 
{
    const int loop_range = 10;
    tmm::TestTimeMeasure my_clock( 0, 0 );
    Time_Measure::SimpleDurationEstimator my_estimator( 0, "firstloop", loop_range, true );

    // some preparation
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    my_clock.measure( "preparation", true );

    // some simulation
    my_clock.enter_loop();
    for( uint i = 0; i < loop_range; ++i )
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        my_estimator.obtain( my_clock.measure( "firstloop" ) );
    }
    my_clock.leave_loop();
    my_clock.measure( "simulation", true );

    // some finalization
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    my_clock.measure( "finalizing", true );
    my_clock.stop();
}
