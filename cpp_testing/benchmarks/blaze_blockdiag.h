#include<complex>
#include<blaze/Math.h>
#include<string>
#include<random>

std::string benchmark_name = "blaze_blockdiag";

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
    RealObservable M1(size/2);
    RealObservable M2(size/2);

    // fill with random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1,1);
    for(uint i = 0; i < size/2; ++i) 
    {
        for(uint j = i; j < size/2; ++j) 
        {
            M1(i,j) = dist(gen);
        }
    }
    for(uint i = 0; i < size/2; ++i) 
    {
        for(uint j = i; j < size/2; ++j) 
        {
            M2(i,j) = dist(gen);
        }
    }

    EigenValues lambda1,lambda2;
    RealOperator U1,U2;
    blaze::eigen(M1,lambda1,U1);
    blaze::eigen(M2,lambda2,U2);
}
