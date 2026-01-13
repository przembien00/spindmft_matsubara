using RealType = double;

#include<iostream>
#include<Multivariate_Gaussian.h>

namespace mvg = Multivariate_Gaussian;

int main()
{
    // 1.) ...inserting a correlation function to a sub triangular matrix
    std::vector<RealType> my_corr{ 1.0, 0.7, 0.45, 0.1 };
    mvg::CovarianceMatrix Cov{};
    Cov.resize( 6 );
    Cov.fill_correlationvector_to_triangle( my_corr, 1, 2 );
    Cov.print();
}