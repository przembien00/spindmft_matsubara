#pragma once
#include <chrono>
#include <ctime>
#include <vector>
#include <string>
namespace CspinDMFT::Time_Measure {
class DerivedTimeMeasure {
    using Clock=std::chrono::steady_clock;
    Clock::time_point start=Clock::now(),last=start,last_measure=start;
    int rank;
    void print_measure(const std::string&,double) const;
public:
    explicit DerivedTimeMeasure(int my_rank=0):rank(my_rank) {}
    std::time_t start_date_and_time=std::time(nullptr);
    std::vector<double> preparation_iterations,real_time_iterations;
    double total_duration=0.;
    void new_iteration() { last=last_measure=Clock::now(); }
    void measure(const std::string&);
    void finish_stage(bool preparation) const;
    void finish_iteration(bool preparation);
    void stop();
};
}
