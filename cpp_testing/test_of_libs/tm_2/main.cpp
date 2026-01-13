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

        // add first loop measures 
        m_first_loop_durations.emplace_back( DurationQuantity{"1st"} );
    }

    std::vector<DurationQuantity> m_first_loop_durations{};

 private:
    DurationType internal_measure( const std::string what, const bool print ) override
    {
        TimeType time = MPI_Wtime();
        DurationType duration{};
        if( m_loop_level == 0 ) // global duration
        {
            duration = time - m_global_measure_time;
            auto iterator = find_measure_in( what, m_global_durations );
            iterator->m_duration = duration;
            m_global_measure_time = time; // update the global reference time
        }
        else if( m_loop_level == 1 ) // first loop duration
        {
            duration = time - m_tmp_measure_times[m_loop_level-1];
            auto iterator = find_measure_in( what, m_first_loop_durations );
            iterator->m_duration += duration;
            m_tmp_measure_times[m_loop_level-1] = time; // update the first loop reference time
        }
        else
        {
            duration = time - m_tmp_measure_times[m_loop_level-1];
            m_tmp_measure_times[m_loop_level-1] = time;
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
    const uint first_loop_range = 10;
    const uint second_loop_range = 50;
    tmm::TestTimeMeasure my_clock( 0, 0 );
    std::vector<uint> loop_sizes{ first_loop_range, first_loop_range*second_loop_range };
    std::vector<std::string> loop_names{ "1st", "2nd" };
    Time_Measure::DurationEstimator my_estimator( 0, "doubleloop", loop_sizes, loop_names, true );

    // some preparation
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    my_clock.measure( "preparation", true );

    // some simulation
    my_clock.enter_loop();
    for( uint i = 0; i < first_loop_range; ++i )
    {
        my_clock.update_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        my_estimator.obtain( my_clock.measure( "1st" ) );
        my_clock.enter_loop();
        for( uint i = 0; i < second_loop_range; ++i )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(30));
            my_estimator.obtain( my_clock.measure( "2nd" ) );
        }
        my_clock.leave_loop();
    }
    my_clock.leave_loop();
    my_clock.measure( "simulation", true );

    // some finalization
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    my_clock.measure( "finalizing", true );
    my_clock.stop();

    // return duration of 1st
    std::cout << "1st duration: " << my_clock.m_first_loop_durations[0].m_duration << "\n";
}
