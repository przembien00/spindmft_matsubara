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

std::vector<ComplexType> delete_one_block_ratios(
    ComplexType total_numerator, ComplexType total_denominator,
    RealType total_absolute_weight,
    const std::vector<ComplexType>& block_numerators,
    const std::vector<ComplexType>& block_denominators,
    const std::vector<RealType>& block_absolute_weights )
{
    const size_t blocks=block_numerators.size();
    if( block_denominators.size()!=blocks || block_absolute_weights.size()!=blocks )
        throw std::invalid_argument("inconsistent ratio-jackknife block arrays");
    std::vector<ComplexType> result(blocks);
    for( size_t b=0;b<blocks;++b )
        result[b]=exact_complex_ratio(total_numerator-block_numerators[b],
                                      total_denominator-block_denominators[b],
                                      total_absolute_weight-block_absolute_weights[b]);
    return result;
}

RealType jackknife_component_error( const std::vector<ComplexType>& values,
                                    bool imaginary )
{
    if( values.size()<2 ) return RealType{0.};
    RealType mean{};
    for( const auto value:values ) mean+=imaginary?std::imag(value):std::real(value);
    mean/=static_cast<RealType>(values.size());
    RealType square_sum{};
    for( const auto value:values )
    {
        const RealType x=imaginary?std::imag(value):std::real(value);
        square_sum+=(x-mean)*(x-mean);
    }
    const RealType B=static_cast<RealType>(values.size());
    return std::sqrt((B-RealType{1.})/B*square_sum);
}

RealType jackknife_absolute_error( const std::vector<ComplexType>& values )
{
    if( values.size()<2 ) return RealType{0.};
    RealType mean{};
    for( const auto value:values ) mean+=std::abs(value);
    mean/=static_cast<RealType>(values.size());
    RealType square_sum{};
    for( const auto value:values )
    {
        const RealType magnitude=std::abs(value);
        square_sum+=(magnitude-mean)*(magnitude-mean);
    }
    const RealType B=static_cast<RealType>(values.size());
    return std::sqrt((B-RealType{1.})/B*square_sum);
}

}
