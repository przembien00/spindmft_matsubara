using RealType = double;

#include<iostream>
#include<Multivariate_Gaussian_Blocks.h>

namespace mgb = Multivariate_Gaussian::Blocks;
namespace mgb = Multivariate_Gaussian::Blocks;

using CorrelationVectorTensor = std::vector<std::vector<RealType>>;

int main()
{
    // Random Number Generator
    std::random_device my_seed{};
    std::mt19937 engine( my_seed() );

    // 1.) create Symmetric Matrix Blocks
    char symmetry_type{'B'};
    mgb::SymmetricMatrix block_xx{ { 1.0, 0.5, 0.25 }, { 0.5, 1.0, 0.5 }, { 0.25, 0.5, 1.0 } };
    mgb::SymmetricMatrix block_zz{ { 1.0, 0.99, 0.98 }, { 0.99, 1.0, 0.99 }, { 0.98, 0.99, 1.0 } };
    mgb::SymmetricMatrixBlocks my_blocks{ block_zz, block_xx };

    // 2.) copy fill it to Covariance Matrix Blocks
    mgb::CovarianceMatrixBlocks my_cov_1{ my_blocks };
    // my_cov_1.print();

    // 3.) move fill it to Covariance Matrix Blocks
    mgb::CovarianceMatrixBlocks my_cov_2{ std::move(my_blocks) };
    my_cov_2.print();
    // std::cout << my_blocks[0](0,0) << std::endl; -> gives memory error as it should

    // 4.) diagonalize
    mgb::EigenValuesBlocks my_eig{};
    mgb::OrthogonalTransformationBlocks my_ortho{};
    my_cov_2.diagonalize( my_eig, my_ortho );
    
    // 5.) write into diagonal normal distributions
    mgb::DiagonalBasisNormalDistributionsBlocks my_ndist{ my_eig };
    my_ndist.print();

    // 6.) sample Gaussian Noise
    std::vector<size_t> num_noises_per_block{ 2, 1 }; // x,y || z
    mgb::GaussianNoiseVectorsBlocks my_noise{ my_ndist, num_noises_per_block, engine };
    my_noise.basis_transform( my_ortho );
    my_noise.print();

}