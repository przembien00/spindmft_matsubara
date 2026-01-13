using RealType = double;

#include<iostream>
#include<Correlation_Tensor.h>
#include<Multivariate_Gaussian_Blocks.h>
#include<Covariance_Filling_Schemes.h>

namespace ct = spinDMFT::Correlation_Tensor;
namespace mvgb = Multivariate_Gaussian::Blocks;
namespace cfs = mvgb::Covariance_Filling_Schemes;

int main()
{   
    // 1.) create correlation tensor
    char symmetry_type{'D'};
    ct::CorrelationTensor<std::vector<RealType>> my_corr( symmetry_type );
    size_t p = 1;
    std::for_each( my_corr.begin(), my_corr.end(), [&p]( auto& x )
    { 
        x = std::vector<RealType>{ 10.0 - p, 9.1 - p, 8.2 - p };
        ++p;
    } );

    // 2.) create filling scheme
    cfs::CorrelationVectorTensorFillingScheme<ct::CorrelationTensor<std::vector<RealType>>> my_scheme{ my_corr, symmetry_type };

    // 3.) create and fill Covariance Matrix Blocks
    mvgb::CovarianceMatrixBlocks my_cov{};
    my_cov.fill_from_scheme( my_scheme );
    my_cov.print();
}