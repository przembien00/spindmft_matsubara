#pragma once
#include"../../Run_Time_Data/Ratio_Statistics.h"
#include<vector>

namespace spinDMFT::Run_Time_Data
{
std::vector<ComplexType> delete_one_block_ratios(
    ComplexType total_numerator, ComplexType total_denominator,
    RealType total_absolute_weight,
    const std::vector<ComplexType>& block_numerators,
    const std::vector<ComplexType>& block_denominators,
    const std::vector<RealType>& block_absolute_weights );

RealType jackknife_component_error( const std::vector<ComplexType>& replicates,
                                    bool imaginary );

RealType jackknife_absolute_error( const std::vector<ComplexType>& replicates );

}
