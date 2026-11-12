#pragma once
#include <Complex_Gaussian/Complex_Gaussian.h>
#include "../Contour/Contour_Types.h"
#include "../Parameter_Space/Parameter_Space.h"
namespace CspinDMFT::Mean_Field_Models {
struct Model {
    size_t sites=0;
    std::vector<double> spins,J,expectation_weights,correlation_weights;
    std::vector<Matrix> operators;
    std::vector<IndexPair> correlation_categories;
    Matrix hamiltonian,identity;
    Physics::FieldMatrix coupling{ {1.,0.,0.},{0.,1.,0.},{0.,0.,1.} };
    void initialize(const ParameterSpace&);
    double weight(size_t,size_t,size_t,size_t) const;
};
Model read_model(const ParameterSpace&);
template<class Correlations>
Correlations connected(const Correlations&,const MagnetizationFor<Correlations>&);
template<class Correlations>
Correlations bath_primitive(const Model&,const Correlations&,const MagnetizationFor<Correlations>&);
template<class Correlations>
Correlations prescribed_bath(const Model&,const ParameterSpace&,const Correlations& grid);
template<class Correlations>
Matrix self_consistent_equations(const Model&,const ParameterSpace&,const Correlations&,
    const MagnetizationFor<Correlations>&,const Contour::Layout&);
template<class Correlations>
spinDMFT::Functions::SymmetricCovarianceSource covariance_source(const Model&,const ParameterSpace&,
    const Correlations&,const MagnetizationFor<Correlations>&,const Contour::Layout&);
template<class Magnetization>
std::vector<C> mean_field(const Model&,const Magnetization&);

}

namespace CspinDMFT {
using Mean_Field_Models::Model;
using Mean_Field_Models::read_model;
using Mean_Field_Models::connected;
using Mean_Field_Models::bath_primitive;
using Mean_Field_Models::prescribed_bath;
using Mean_Field_Models::mean_field;
using Mean_Field_Models::self_consistent_equations;
}
