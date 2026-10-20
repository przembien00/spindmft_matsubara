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

// Smooth primitive interpolation in units of the measurement-grid spacing.
// Uses at most four neighboring knots (linear/quadratic on short grids).
ComplexType interpolated_edge_value(const CorrelationSet& correlations,
    RealType real_coordinate,size_t a,size_t b,size_t imaginary_edge);

// Direct covariance in the staggered (M,eta,nu) basis. Only smooth kernel
// values are interpolated; retarded support is applied at the actual nodes.
ComplexType staggered_correlation(const CorrelationSet& correlations,
    const StaggeredIndex& first,const StaggeredIndex& second,size_t substeps);

}
