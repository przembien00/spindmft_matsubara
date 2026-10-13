#include"../Run_Time_Data/Ratio_Statistics.h"

#include<cmath>
#include<iostream>

namespace rtd = spinDMFT::Run_Time_Data;
namespace
{
int require( bool condition,const char* message )
{
    if( condition ) return 0;
    std::cerr<<"FAILED: "<<message<<'\n'; return 1;
}
}

int main()
{
    int failures{};
    const std::vector<ComplexType> Z{{1.,1.},{2.,-0.5},{0.5,0.25},{1.5,-0.2}};
    const std::vector<ComplexType> N{{2.,-1.},{-0.5,0.4},{1.2,0.7},{0.1,-0.8}};
    ComplexType sum_Z{},sum_N{}; RealType sum_abs{};
    for( size_t i=0;i<Z.size();++i )
    { sum_Z+=Z[i]; sum_N+=N[i]; sum_abs+=std::abs(Z[i]); }
    const auto exact=rtd::exact_complex_ratio(sum_N,sum_Z,sum_abs);
    failures+=require(std::abs(exact-sum_N/sum_Z)<RealType{1e-14},
                      "bare estimator uses the complex Z_M sum");
    failures+=require(std::abs(exact-sum_N/sum_abs)>RealType{1e-3},
                      "phase-quenched absolute-Z ratio is distinct");
    std::vector<RealType> block_abs;
    for( const auto z:Z ) block_abs.push_back(std::abs(z));
    const auto replicates=rtd::delete_one_block_ratios(
        sum_N,sum_Z,sum_abs,N,Z,block_abs);
    failures+=require(replicates.size()==Z.size(),"delete-one-block replicate count");
    failures+=require(rtd::jackknife_component_error(replicates,false)>RealType{0.},
                      "real ratio jackknife error");
    failures+=require(rtd::jackknife_component_error(replicates,true)>RealType{0.},
                      "imaginary ratio jackknife error");
    failures+=require(rtd::jackknife_absolute_error(replicates)>RealType{0.},
                      "absolute ratio jackknife error");
    std::vector<ComplexType> paired_differences(replicates.size());
    for( size_t b=0;b<replicates.size();++b )
        paired_differences[b]=replicates[b]-replicates[b];
    failures+=require(rtd::jackknife_component_error(
                          paired_differences,false)==RealType{0.},
                      "paired difference preserves exact common-mode cancellation");
    failures+=require(rtd::jackknife_absolute_error(
                          paired_differences)==RealType{0.},
                      "paired absolute difference preserves exact common-mode cancellation");
    failures+=require(!rtd::denominator_resolved(ComplexType{},sum_abs),
                      "zero complex denominator is rejected");
    return failures==0?0:1;
}
