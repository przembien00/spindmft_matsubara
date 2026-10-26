#include"Ratio_Statistics.h"

#include<cmath>
#include<limits>
#include<stdexcept>

namespace spinDMFT::Run_Time_Data
{

bool denominator_resolved( ComplexType denominator, RealType absolute_weight_sum )
{
    if( !std::isfinite(std::real(denominator))
        || !std::isfinite(std::imag(denominator))
        || !std::isfinite(absolute_weight_sum) ) return false;
    const RealType scale=std::max(RealType{1.},absolute_weight_sum);
    const RealType tolerance=RealType{100.}*std::numeric_limits<RealType>::epsilon()*scale;
    return std::abs(denominator)>tolerance;
}

ComplexType exact_complex_ratio( ComplexType numerator, ComplexType denominator,
                                 RealType absolute_weight_sum )
{
    if( !denominator_resolved(denominator,absolute_weight_sum) )
        throw std::runtime_error("the complex partition-function denominator is unresolved");
    return numerator/denominator;
}

}
