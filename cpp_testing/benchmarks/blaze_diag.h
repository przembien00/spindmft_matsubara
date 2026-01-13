#include<complex>
#include<blaze/Math.h>
#include<string>
#include<random>

std::string benchmark_name = "blaze_diag";

typedef double RealType;
typedef std::complex<RealType> ComplexType;

// Quantum Operators
typedef blaze::HermitianMatrix<blaze::DynamicMatrix<ComplexType,blaze::rowMajor>> Observable;
typedef blaze::SymmetricMatrix<blaze::DynamicMatrix<RealType,blaze::rowMajor>> RealObservable;
typedef blaze::DynamicMatrix<ComplexType,blaze::rowMajor> Operator;
typedef blaze::DynamicMatrix<RealType,blaze::rowMajor> RealOperator;

typedef blaze::DynamicVector<RealType,blaze::columnVector> EigenValues;

/*
// Sparse Quantum Operators
typedef blaze::DiagonalMatrix<Operator> DiagonalOperator;
typedef blaze::DiagonalMatrix<RealOperator> RealDiagonalOperator;
typedef blaze::HermitianMatrix<blaze::CompressedMatrix<ComplexType,blaze::rowMajor>> SparseObservable;
typedef blaze::SymmetricMatrix<blaze::CompressedMatrix<RealType,blaze::rowMajor>> SparseRealObservable;
typedef blaze::CompressedMatrix<ComplexType,blaze::rowMajor> SparseOperator;
typedef blaze::CompressedMatrix<RealType,blaze::rowMajor> SparseRealOperator;
*/

// Function to benchmark
void operation() 
{
    uint size = 1000;
    RealObservable M(size);

    // fill with random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1,1);

    for(uint i = 0; i < size; ++i) 
    {
        for(uint j = i; j < size; ++j) 
        {
            M(i,j) = dist(gen);
        }
    }

    EigenValues lambda;
    RealOperator U;
    blaze::eigen(M,lambda,U);
}
