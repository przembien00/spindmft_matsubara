#pragma once
#include "Contour_Types.h"
namespace CspinDMFT::Contour {
C kernel(const MatsubaraCorrelations&,Node,Node,size_t substeps);
C kernel(const ContourCorrelations&,Node,Node,size_t substeps);
template<class Correlations>
Matrix covariance(const Correlations&,const Layout&,size_t substeps);
}
namespace CspinDMFT { using Contour::kernel; using Contour::covariance; }
