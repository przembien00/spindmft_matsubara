#pragma once

#include<vector>

#include<Globals/Types.h>

namespace spinDMFT::Run_Time_Data
{

bool denominator_resolved( ComplexType denominator, RealType absolute_weight_sum );

ComplexType exact_complex_ratio( ComplexType numerator, ComplexType denominator,
                                 RealType absolute_weight_sum );

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
