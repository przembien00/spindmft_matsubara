#pragma once

#include<Globals/Types.h>

namespace spinDMFT::Run_Time_Data
{

bool denominator_resolved( ComplexType denominator, RealType absolute_weight_sum );

ComplexType exact_complex_ratio( ComplexType numerator, ComplexType denominator,
                                 RealType absolute_weight_sum );

}
