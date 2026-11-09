#include "Time_Measure.h"
#include <Standard_Algorithms/Print_Routines.h>
#include <numeric>
#include <sstream>
namespace CspinDMFT::Time_Measure {
void DerivedTimeMeasure::finish_iteration(bool preparation) {
    auto elapsed=std::chrono::duration<double>(Clock::now()-last).count();
    (preparation?preparation_iterations:real_time_iterations).push_back(elapsed);
}
void DerivedTimeMeasure::print_measure(const std::string& what,double duration) const {
    std::ostringstream value; value<<"it took me "<<duration<<" seconds";
    Print_Routines::print_R0(rank,Print_Routines::quantity_to_output_line(35,
        "\033[1;33m"+what+" done (0)\033[0m",value.str()));
}
void DerivedTimeMeasure::measure(const std::string& what) {
    const auto now=Clock::now();
    print_measure(what,std::chrono::duration<double>(now-last_measure).count());
    last_measure=now;
}
void DerivedTimeMeasure::finish_stage(bool preparation) const {
    const auto& iterations=preparation?preparation_iterations:real_time_iterations;
    print_measure(preparation?"Matsubara self consistency":"real-time self consistency",
        std::accumulate(iterations.begin(),iterations.end(),0.));
}
void DerivedTimeMeasure::stop() {
    total_duration=std::chrono::duration<double>(Clock::now()-start).count();
    std::ostringstream value; value<<"in total it took me "<<total_duration<<" seconds";
    Print_Routines::print_R0(rank,Print_Routines::quantity_to_output_line(35,
        "\033[1;36meverything done (0)\033[0m",value.str()));
}
}
