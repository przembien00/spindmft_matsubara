using RealType = double;

#include<iostream>
#include<Multivariate_Gaussian.h>

namespace mvg = Multivariate_Gaussian;

int main()
{
    // Random Number Generator
    std::random_device my_seed{};
    std::mt19937 engine( my_seed() );

    // Sampling Noise from a multivariate Gaussian distribution...
    // 1.) ...inserting a correlation function
    std::vector<RealType> my_corr{ 1.0, 0.7, 0.45, 0.25, 0.1, 0.0 };
    size_t num_noises = 5;
    auto my_gaussnoise = mvg::sample_multivariateGaussian_noise( my_corr, engine, num_noises );
    std::cout << "Created " << num_noises << " noises from a correlation function: \n";
    my_gaussnoise.print();

    // 2.) ...inserting a symmetric matrix
    mvg::SymmetricMatrix my_mat{ {1.0, 0.99}, {0.99, 1.0} };
    std::vector<RealType> mean_values{ 17.0, 7.0 };
    num_noises = 4;
    my_gaussnoise = mvg::sample_multivariateGaussian_noise( my_mat, engine, num_noises );
    my_gaussnoise.add_mean( mean_values );
    std::cout << "Created " << num_noises << " noises from a symmetric matrix: \n";
    my_gaussnoise.print();

    // 3.) ...inserting a covariance matrix
    mvg::CovarianceMatrix my_cov{ my_mat };
    RealType mean_value{5.0};
    num_noises = 3;
    my_gaussnoise = mvg::sample_multivariateGaussian_noise( my_cov, engine, num_noises );
    my_gaussnoise.add_const_mean( mean_value );
    std::cout << "Created " << num_noises << " noises from a covariance matrix: \n";
    my_gaussnoise.print();
}