#include<vector>
#include<string>
#include<chrono>
#include<thread>
#include<mpi.h>
#include<numeric>
#include"../Standard_Algorithms/Print_Routines.h"
#include"../Standard_Algorithms/Standard_Algorithms.h"
#include"TM_Error_Handling.h"

#pragma once

namespace Time_Measure
{

// ================= NAMESPACES ========================
namespace stda = Standard_Algorithms;
namespace print = Print_Routines;
namespace error = Error_Handling;

// ================= USING STATEMENTS ==================
using TimeType = double;
using DurationType = double;


// =====================================================
// ============== DURATION QUANTITY CLASS ==============
// =====================================================
// class that stores the duration of a certain process/quantity
struct DurationQuantity
{
  explicit DurationQuantity( const std::string& quantity_name ) : m_name( quantity_name ) {};
  DurationQuantity( const std::string& quantity_name, const DurationType value ) : m_name( quantity_name ), m_duration( value ) {};
  std::string m_name{};       // name of the process/quantity
  DurationType m_duration{};  // value of the measured duration
};


// =====================================================
// ========== HEADER OF TIME MEASURE CLASS =============
// =====================================================
/* base class for measuring the times as well as the durations of different computations
- for-loops are entered by enter() and left by leave()
- a duration measurement is carried out with measure() (in or outside of loops), 
- one might want to use update_time() which entails that the subsequent duration measurement is done with respect to the time of this function call; otherwise duration measurements are always done with regard to the previous measurement
- in the end the time measurement is m_stopped with stop()
- two different measures are distinguished: m_global_measure_time and m_tmp_measure_times
global measures are measure times of some specific computation with a name, whose durations are always stored
temporary measures are local times within some nested for loops; their elements are constantly overwritten during the corresponding loop and deleted if the loop is left; hence, storing or using specific measurements has to be done within the loops */
class TimeMeasure
{ 
 protected:
  // PROTECTED MEMBERS
  const int my_rank;
  const int world_size;
  const uint m_pre_colon_space;
  uint m_loop_level{0};
  bool m_stopped{false};

 public:
  // CONSTRUCTORS
  TimeMeasure( const int my_rank_, const int world_size_ );

  // PUBLIC METHODS
  void enter_loop();
  void update_time();
  void leave_loop();
  DurationQuantity measure( const std::string what, const bool print = false );
  void stop(); 

  // PUBLIC MEMBERS
  const TimeType m_start_time;              // time of the constructor call 
  const std::time_t m_start_date_and_time;  // date and time of the constructor call
  TimeType m_end_time;                        // time where stop has been called
  DurationType m_total_duration;              // total duration of the job
  
  TimeType m_global_measure_time{};              // time of the previous global measurement (global reference time)
  std::vector<TimeType> m_tmp_measure_times{};   // times of the previous measurement in the respective loop (local reference time)

  std::vector< DurationQuantity > m_global_durations{}; // global durations
  
 protected:
  // PROTECTED METHODS 
  virtual DurationType internal_measure( const std::string what, const bool print ) = 0;
  void print_measure( const std::string what, const DurationType duration ) const; 
};


// =====================================================
// ===== HEADER OF SIMPLE DURATION ESTIMATOR CLASS =====
// =====================================================
/* class that estimates the duration of a for loop and prints a progress bar if desired 
- the estimator obtains the duration of a single loop run and simply accumulates it with the loop size to estimate the total duration
- the estimated duration is not updated as the loop is carried out */
class SimpleDurationEstimator
{ 
 private:
  // PRIVATE MEMBERS
  const uint my_rank;
  const uint m_pre_colon_space;
  const std::string m_name;
  const uint m_steps;
  const bool m_show_progress;  
  const uint m_bar_width;
  bool m_estimation_done{false};
  uint m_current_step{};
  DurationType m_single_step_duration{};
  DurationType m_estimated_total{};

 public:
  // CONSTRUCTORS
  SimpleDurationEstimator( const uint my_rank, const std::string name, 
    const uint multiplicity, const bool show_progress = false );

  // PUBLIC METHODS
  void obtain( const DurationQuantity& single_step_duration );
  void reset();

 private:
  // PRIVATE METHODS 
  void estimate(); 
  void show_progress() const;
};


// =====================================================
// ======== HEADER OF DURATION ESTIMATOR CLASS =========
// =====================================================
/* class that estimates the duration of nested for-loops and prints a progress bar if desired 
- this more complex estimator obtains the durations of a couple of computations that occur in different (nested) loops as well as the corresponding loop sizes and accumulates everything to estimate the total duration
- the estimated duration is not updated as the loop is carried out */
class DurationEstimator
{ 
 private:
  // PRIVATE MEMBERS
  const uint my_rank;
  const uint m_pre_colon_space;
  const std::string m_name;
  const std::vector<size_t> m_loop_sizes;
  std::vector<uint> m_loop_counters{};
  const bool m_show_progress;  
  const uint m_bar_width;
  bool m_estimation_done{false};
  std::vector<DurationQuantity> m_duration_quantities{};
  std::vector<double> m_fractions{};
  DurationType m_estimated_total{};

 public:
  // CONSTRUCTORS
  DurationEstimator( const uint my_rank, const std::string name, 
    const std::vector<size_t>& loop_sizes, const std::vector<std::string>& loop_names, const bool show_progress = false );

  // PUBLIC METHODS
  void obtain( const DurationQuantity& single_step_duration );
  void reset();

 private:
  // PRIVATE METHODS 
  void estimate(); 
  void show_progress() const;
};


// ================= FUNCTIONS =================
template<typename T> // requirement : m_name needs to be member of T
typename std::vector<T>::iterator find_measure_in( const std::string what, std::vector<T>& list_of_measures )
{
    auto iterator = std::find_if( list_of_measures.begin(), list_of_measures.end(), 
    [&what]( auto& q )
    {
        return (what == q.m_name);
    } );
    if( iterator != list_of_measures.end() ) // what has been found
    {
        return iterator;
    }
    else // what has not been found
    {
        error::TIME_MEASURE( what, std::string(__PRETTY_FUNCTION__) );
        return list_of_measures.end();
    }
}


// =====================================================
// ========== IMPLEMENTATION OF TIME MEASURE ===========
// =====================================================
// base constructor
inline TimeMeasure::TimeMeasure( const int my_rank_, const int world_size_ ):
    my_rank( my_rank_ ),
    world_size( world_size_ ),
    m_pre_colon_space( 35 ),
    m_start_time( MPI_Wtime() ),
    m_start_date_and_time( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) )
{
    m_global_measure_time = m_start_time;
}

// enter a loop : add new temporary time measurement and increase counter
inline void TimeMeasure::enter_loop()
{
    m_tmp_measure_times.emplace_back(MPI_Wtime());
    ++m_loop_level;
}

// update the outermost time; this may be done in the beginning of loops
inline void TimeMeasure::update_time()
{
    if( m_loop_level == 0 )
    {
        m_global_measure_time = MPI_Wtime();
    }
    else
    {
        *(--m_tmp_measure_times.end()) = MPI_Wtime();
    }
}

// leave a loop : delete temporary time measurement and decrease counter
inline void TimeMeasure::leave_loop()
{
    m_tmp_measure_times.pop_back();
    --m_loop_level;
}

// measure the duration for 'what' (calling private measure method)
inline DurationQuantity TimeMeasure::measure( const std::string what, const bool print )
{
    if( !m_stopped )
    {
        return DurationQuantity{ what, internal_measure( what, print ) };
    }
    else
    {
        error::TIME_STOPPED( __PRETTY_FUNCTION__ ) ;
        return DurationQuantity{ what, DurationType{0.} }; // suppress warnings
    }
}

// stop the whole measurement : compute the total time, print the measurement results at core 0 and forbid the usage of measure
inline void TimeMeasure::stop()
{   
    // check whether all entered loops have been left
    if( m_loop_level != 0 )
    {
        error::TIME_STOPPED_IN_LOOP( m_loop_level, __PRETTY_FUNCTION__ );
    }

    // measure the end time and total duration
    m_end_time = MPI_Wtime();
    m_total_duration = m_end_time - m_start_time;

    // print the total simulation time
    if( my_rank == 0 )
    {
        std::string what = "everything";
        std::stringstream ss{};
        ss << "in total it took me " << m_total_duration << " seconds";
        std::cout << print::quantity_to_output_line( m_pre_colon_space, "\033[1;36m" + what + " done (0)\033[0m", ss.str() );
    }

    // set m_stopped to true
    m_stopped = true;
}

// printing method : print the measurement results at core 0
inline void TimeMeasure::print_measure( const std::string what, const DurationType duration ) const
{
    if( my_rank == 0 )
    {
        std::stringstream ss{};
        ss << "it took me " << duration << " seconds";
        std::cout << print::quantity_to_output_line( m_pre_colon_space, "\033[1;33m" + what + " done (0)\033[0m", ss.str() );
    }
}


// =====================================================
// === IMPLEMENTATION OF SIMPLE DURATION ESTIMATOR =====
// =====================================================
// constructor
inline SimpleDurationEstimator::SimpleDurationEstimator( const uint my_rank, const std::string name, const uint steps, const bool show_progress ):
    my_rank( my_rank ),
    m_pre_colon_space( 35 ),
    m_name( name ),
    m_steps( steps ),
    m_show_progress( show_progress ),
    m_bar_width( 50 )
{}

// obtain the duration of a process (name is not checked)
inline void SimpleDurationEstimator::obtain( const DurationQuantity& single_step_duration )
{
    ++m_current_step;
    if( !m_estimation_done ) // if the time of the loop has not yet been estimated
    {
        m_single_step_duration = single_step_duration.m_duration;
        estimate();
        m_estimation_done = true;
    }
    if( m_show_progress ) // if demanded show progress bar
    {
        show_progress();
    }
}

// reset the time estimator so that it can be reused
inline void SimpleDurationEstimator::reset()
{
    m_estimation_done = false;
    m_current_step = 0;
    m_single_step_duration = DurationType{0.};
    m_estimated_total = DurationType{0.};
}

// compute and print the duration estimate
inline void SimpleDurationEstimator::estimate()
{
    // estimate:
    m_estimated_total = m_single_step_duration * static_cast<DurationType>(m_steps);

    // print:
    if( my_rank == 0 )
    {
        std::stringstream ss{};
        ss << m_estimated_total << " seconds";
        std::cout << print::quantity_to_output_line( m_pre_colon_space, std::string{"\033[1;33m"} 
            + std::string{"estimated duration of "} + m_name + std::string{"\033[0m"}, ss.str() );
    }
}

// print the progress bar
inline void SimpleDurationEstimator::show_progress() const
{
    if( my_rank == 0 )
    {
        double fraction = static_cast<double>( m_current_step ) / static_cast<double>( m_steps );
        uint filled = static_cast<uint>( fraction * m_bar_width );

        std::cout << "\r[";
        for ( uint i = 0; i < m_bar_width; ++i ) {
            if (i < filled) {
                std::cout << "=";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << (uint)( fraction * double{100.0} ) << "%";
        std::cout.flush();
        if( m_current_step == m_steps )
        {
            std::cout << "\n";
        }
    }

}


// =====================================================
// ======= IMPLEMENTATION OF DURATION ESTIMATOR ========
// =====================================================
// constructor
inline DurationEstimator::DurationEstimator( const uint my_rank, const std::string name, 
    const std::vector<size_t>& loop_sizes, const std::vector<std::string>& loop_names, const bool show_progress ):
    my_rank( my_rank ),
    m_pre_colon_space( 35 ),
    m_name( name ),
    m_loop_sizes( loop_sizes ),
    m_show_progress( show_progress ),
    m_bar_width( 50 )
{
    // resize and fill counters with zeros:
    std::transform( m_loop_sizes.cbegin(), m_loop_sizes.cend(), std::back_inserter( m_loop_counters ), []( auto q ){ return 0; } );
    // resize and fill duration quantities:
    std::transform( m_loop_sizes.cbegin(), m_loop_sizes.cend(), loop_names.cbegin(), std::back_inserter( m_duration_quantities ), 
        []( const auto q, const std::string& name )
    { 
        return DurationQuantity{ name, 0. }; 
    } );
}

// obtain the duration of a process
inline void DurationEstimator::obtain( const DurationQuantity& single_step_duration )
{
    // determine the obtained quantity and increase the corresponding loop counter:
    auto iterator = find_measure_in( single_step_duration.m_name, m_duration_quantities );
    uint q_index = std::distance( m_duration_quantities.begin(), iterator );
    ++m_loop_counters[ q_index ];

    if( !m_estimation_done ) // if the time of the loop has not yet been estimated...
    {
        m_duration_quantities[ q_index ].m_duration = single_step_duration.m_duration;
        if( m_loop_counters.back() == 1 ) // ...and the last loop level has just been reached
        {

            estimate();
            m_estimation_done = true;
        }
    }
    if( m_show_progress && q_index==0 ) // if demanded show progress bar
    {
        show_progress();
    }
}

// reset the time estimator so that it can be reused
inline void DurationEstimator::reset()
{
    m_estimation_done = false;
    std::for_each( m_duration_quantities.begin(), m_duration_quantities.end(), []( auto& q ){ q = DurationQuantity{ q.m_name, 0. }; } );
    std::fill( m_loop_counters.begin(), m_loop_counters.end(), 0 );
    m_estimated_total = DurationType{0.};
}

// compute and print the duration estimate
inline void DurationEstimator::estimate()
{
    // estimate:
    stda::for_2each( m_duration_quantities.cbegin(), m_duration_quantities.cend(), m_loop_sizes.cbegin(), 
        [this]( auto& q, auto& num )
    {
        m_estimated_total += static_cast<DurationType>(num) * q.m_duration;
    } );

    // print:
    if( my_rank == 0 )
    {
        std::stringstream ss{};
        ss << m_estimated_total << " seconds";
        std::cout << print::quantity_to_output_line( m_pre_colon_space, std::string{"\033[1;33m"} 
            + std::string{"estimated duration of "} + m_name + std::string{"\033[0m"}, ss.str() );
    }
}

// print the progress bar (with regard to the outermost loop)
inline void DurationEstimator::show_progress() const
{
    if( my_rank == 0 )
    {
        // progress bar should not be shown, if the single step durations are yet unknown
        if( std::find( m_loop_counters.cbegin(), m_loop_counters.cend(), uint{0} ) != m_loop_counters.cend() )
        {
            return;
        }

        double fraction = static_cast<double>( m_loop_counters[0] ) / static_cast<double>( m_loop_sizes[0] );
        uint filled = static_cast<uint>(fraction * m_bar_width);

        std::cout << "\r[";
        for ( uint i = 0; i < m_bar_width; ++i ) {
            if (i < filled) {
                std::cout << "=";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << (uint)( fraction * double{100.0} ) << "%";
        std::cout.flush();
        if( m_loop_counters[0] == m_loop_sizes[0] )
        {
            std::cout << "\n";
        }
    }
}

};
