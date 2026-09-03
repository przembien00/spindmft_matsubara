#pragma once

#include<Globals/Types.h>

#include"Contour_Types.h"

namespace spinDMFT::Contour
{

ComplexType edge_value( const CorrelationSet& correlations, size_t real_time,
                        size_t first_spin, size_t second_spin,
                        size_t imaginary_edge );

ComplexType greater_value( const CorrelationSet& correlations,
                           size_t first_spin, size_t second_spin,
                           std::ptrdiff_t time_difference );

ComplexType lesser_value( const CorrelationSet& correlations,
                          size_t first_spin, size_t second_spin,
                          std::ptrdiff_t time_difference );

// Raw contour-ordered spin correlation C^{AB}_{ab}. This function reconstructs
// every branch block from the single edge-grid primitive X(t,tau).
ComplexType branch_correlation( const CorrelationSet& correlations,
                                const ContourLayout& layout,
                                const ContourIndex& first,
                                const ContourIndex& second );

}
