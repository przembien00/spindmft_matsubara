#include "Mean_Field_Models/Mean_Field_Models.h"
#include <cmath>
#include <stdexcept>

using namespace Mean_Field_Models;
struct TestModel : CorrelationReplicaModel {
    TestModel() {
        num_Spins=2;
        correlation_categories={{0,0},{0,1}};
        mf_expectation_weights=Matrix(2,2,0.);
        spinmf_cmodel.coupling_matrix=blaze::IdentityMatrix<RealType>(3);
    }
};
int main() {
    TestModel model;
    std::vector<double> weights(16,0.);
    weights[0]=3.; // W[0,0,0,0]: diagonal, counted once.
    weights[1]=0.5; // W[0,0,0,1]: representative, counted twice.
    weights[2]=0.5; // Symmetric reverse entry must not be counted a second time.
    model.interpret_model_specific_parameters(weights);
    CluCorrTen source{IndexPairList{{0,0},{0,1}}, 'D', size_t{3}};
    source(0,0)[0]=Corr{3,0.25};
    source(0,1)[0]=Corr{3,0.125};
    SiteFields moments{{0.1,0.,0.},{0.2,0.,0.}};
    const auto [mean,bath]=model.self_consistency(source,moments);
    const double expected=3.*(0.25-0.1*0.1)+(0.125-0.1*0.2);
    for(size_t t=0;t<3;++t)
        if(std::abs(bath(0,0)(0,0)[t]-expected)>1e-13)
            throw std::runtime_error("incorrect diagonal/off-diagonal connected feedback multiplicity");
}
