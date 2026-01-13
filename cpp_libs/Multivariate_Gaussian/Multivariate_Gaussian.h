#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"MVG_Error_Handling.h"
#include"../Matrices/Diagonalization.h"

namespace Multivariate_Gaussian
{

namespace error = Multivariate_Gaussian::Error_Handling;
namespace stda = Standard_Algorithms;
namespace diag = Matrices::Diagonalization;

// ================== USING STATEMENTS =================
using EigenValues = blaze::DynamicVector<RealType, blaze::columnVector>;
using SymmetricMatrix = blaze::SymmetricMatrix<blaze::DynamicMatrix<RealType, blaze::rowMajor>>;
using NormalDistributions = std::vector<std::normal_distribution<RealType>>;
using OrthogonalTransformation = blaze::DynamicMatrix<RealType, blaze::rowMajor>;
using NoiseMatrix = blaze::DynamicMatrix<RealType, blaze::columnMajor>;


// =====================================================
// ======== HEADER FOR COVARIANCE MATRIX CLASS =========
// =====================================================
/* builds a covariance matrix that is enforced to be symmetric (wrapper around SymmetricMatrix)
the matrix is not enforced to be positive (user responsibility);
it can be filled from a matrix or from correlation vectors and it can be diagonalized 
vectors (v1, v2, v3, ...) are inserted according to the following scheme:
( v1 v2 v3 ... )
( v2 v1 v2 ... )
( v3 v2 v1 ... )
( .  .  .  ... )
matrices and vectors may also be inserted into a subblock of the matrix providing the start row and column */
class CovarianceMatrix
{
public:
    // CONSTRUCTORS
    CovarianceMatrix() = default;
    explicit CovarianceMatrix( const size_t& size ) { resize( size ); };
    template<typename CorrelationVector> // needs to be an iterable object
    explicit CovarianceMatrix( const CorrelationVector& corr ) { this->fill_correlationvector( corr ); };
    explicit CovarianceMatrix( const SymmetricMatrix& mat ) { m_cov = mat; };
    explicit CovarianceMatrix( SymmetricMatrix&& mat ) { m_cov = std::move( mat ); };

    // PUBLIC METHODS
    template<typename CorrelationVector> // needs to be an iterable object
    void fill_correlationvector( const CorrelationVector& corr );
    template<typename CorrelationVector> // needs to be an iterable object
    void fill_correlationvector_to_subtriangle( const CorrelationVector& corr, const size_t srow, const size_t scol );
    template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
    void fill_correlationmatrix( const CorrelationMatrix& corr );
    template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
    void fill_correlationmatrix_to_subsquare( const CorrelationMatrix& corr, const size_t srow, const size_t scol );
    template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
    void fill_symmetriccorrelationmatrix_to_diagonal_subblock( const CorrelationMatrix& corr, const size_t srow );

    void diagonalize( EigenValues& eig, OrthogonalTransformation& ortho ) const;
    void resize( const size_t size ) { m_cov.resize( size ); };
    const size_t size() const { return m_cov.rows(); };

    void print() const { std::cout << m_cov; };

private:
    // PRIVATE MEMBERS
    SymmetricMatrix m_cov{};
};


// =====================================================
// ======== HEADER FOR NORMAL DISTRIBUTION CLASS =======
// =====================================================
/* builds a vector of normal distributions (wrapper around NormalDistributions) setting their mean values to zero and their variances from a vector of eigenvalues; throws an error if some eigenvalues are negative */
class DiagonalBasisNormalDistributions
{
public:
    // CONSTRUCTORS
    DiagonalBasisNormalDistributions();
    explicit DiagonalBasisNormalDistributions( const EigenValues& eig ) { this->fill( eig ); };

    // PUBLIC METHODS
    void fill( const EigenValues& eig );
    void resize( const size_t size ){ m_ndist.resize( size ); }
    const size_t size() const { return m_ndist.size(); }
    void print() const;

    // ITERATORS
    auto begin() { return m_ndist.begin(); };
    auto end() { return m_ndist.end(); };
    auto cbegin() const { return m_ndist.cbegin(); };
    auto cend() const { return m_ndist.cend(); };

private:
    // PRIVATE MEMBERS
    NormalDistributions m_ndist{};
};


// =====================================================
// ====== HEADER FOR GAUSSIAN NOISE VECTORS CLASS ======
// =====================================================
/* builds noise vectors from a given vector of distributions and writes them in the columns of a matrix (wrapper around NoiseMatrix);
The matrix has m_noise.rows() number of rows (noise degrees of freedom) and m_noise.cols() number of columns (number of drawn noise vectors) */
class GaussianNoiseVectors 
{
public:
    // CONSTRUCTORS
    GaussianNoiseVectors();
    GaussianNoiseVectors( DiagonalBasisNormalDistributions& ndist, std::mt19937& engine ) { this->fill( ndist, engine ); };
    GaussianNoiseVectors( DiagonalBasisNormalDistributions& ndist, std::mt19937& engine, const size_t num_noises ) { this->fill( ndist, engine, num_noises );};

    // PUBLIC METHODS
    void fill( DiagonalBasisNormalDistributions& ndist, std::mt19937& engine ) { this->fill( ndist, engine, 1 ); };
    void fill( DiagonalBasisNormalDistributions& ndist, std::mt19937& engine, const size_t num_noises );
    void basis_transform( const OrthogonalTransformation& ortho ) { m_noise = ortho * m_noise; };
    void print() const { std::cout << m_noise; };

    // OPERATORS
    RealType& operator()( const size_t at, const size_t noise_vector_index ) { return m_noise( at, noise_vector_index ); };
    const RealType& operator()( const size_t at, const size_t noise_vector_index ) const { return m_noise( at, noise_vector_index ); };

    // ITERATORS
    auto begin( const size_t row ) { return m_noise.begin( row ); };
    auto end( const size_t row ) { return m_noise.end( row ); };
    auto cbegin( const size_t row ) const { return m_noise.cbegin( row ); };
    auto cend( const size_t row ) const { return m_noise.cend( row ); };

private:
    // PRIVATE MEMBERS
    NoiseMatrix m_noise{};
};


// =======================================================
// ====== IMPLEMENTATION OF COVARIANCE MATRIX CLASS ======
// =======================================================
/*  fills the covariance matrix from a correlationvector 
also automatically resizes corresponding to the correlationvector size
the lower triangular matrix is explicitly filled and due to symmetry also the rest */
template<typename CorrelationVector> // needs to be an iterable object
void CovarianceMatrix::fill_correlationvector( const CorrelationVector& corr )
{   
    this->resize( corr.size() );
    for( size_t row = 0; row < corr.size(); ++row ) // fill from top to bottom (small to large)
    {
        std::copy( corr.crbegin() + (corr.size()-row-1), corr.crend(), m_cov.begin( row ) );
    }
}
/* fills a part of the covariance matrix from a correlationvector 
an upper sub-triangular matrix is explicitly filled (and due to symmetry also the corresponding lower sub-triangular matrix)
the point (srow, scol) defines the upper left corner of the triangle */
template<typename CorrelationVector> // needs to be an iterable object
void CovarianceMatrix::fill_correlationvector_to_subtriangle( const CorrelationVector& corr, const size_t srow, const size_t scol )
{   
    if( (int)corr.size() > (int)m_cov.rows() - (int)srow || (int)corr.size() > (int)m_cov.rows() - (int)scol ) // check container size
    {
        error::INVALID_SIZE( corr.size(), __PRETTY_FUNCTION__ );
    }
    for( size_t row = 0; row < corr.size(); ++row ) // fill from top to bottom (small to large)
    {
        std::copy( corr.crbegin() + (corr.size()-row-1), corr.crend(), m_cov.begin( srow + row ) + scol );
    }
}

/* fills the covariance matrix from a correlationmatrix 
also automatically resizes corresponding to the correlationmatrix size
symmetry of the inserted correlationmatrix is assumed and not checked (users responsibility) */
template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
void CovarianceMatrix::fill_correlationmatrix( const CorrelationMatrix& corr )
{
    this->resize( corr.rows() );
    for( size_t row = 0; row < corr.rows(); ++row ) 
    {
        std::copy( corr.cbegin(row) + row, corr.cend(row), m_cov.begin(row) + row ); // copy only the upper triangle (rest is autofilled)
    }
}
/* fills a part of the covariance matrix from a CorrelationMatrix 
a subblock is explicitly filled (and due to symmetry of the covariance matrix potentially also mirrored indices)
if the subblock has overlap with diagonal elements of the covariance matrix, the inserted CorrelationMatrix must be symmetric around a certain diagonal; this is not explicitly checked in this method (users responsibility) */
template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
void CovarianceMatrix::fill_correlationmatrix_to_subsquare( const CorrelationMatrix& corr, const size_t srow, const size_t scol )
{
    if( (int)corr.rows() > (int)m_cov.rows() - (int)srow || (int)corr.rows() > (int)m_cov.rows() - (int)scol ) // check container size
    {
        error::INVALID_SIZE( corr.rows(), __PRETTY_FUNCTION__ );
    }
    for( size_t row = 0; row < corr.rows(); ++row ) 
    {
        std::copy( corr.cbegin( row ), corr.cend( row ), m_cov.begin( srow + row ) + scol ); // copy the whole matrix to the subblock (in addition, some mirrored indices in the covariance matrix may also be filled due to symmetry)
    }
}
/* fills a diagonal subblock of the covariance matrix from a symmetric CorrelationMatrix 
the upper sub-triangular matrix is explicitly filled (and due to symmetry also the mirrored indices)
this function is a more efficient alternative to fill_correlationmatrix_to_subsquare for diagonal subblocks
for non-diagonal subblocks it cannot be used 
not also that the required symmetry of CorrelationMatrix is not checked (users responsibility) */
template<typename CorrelationMatrix> // needs to be a rowwise iterable square matrix
void CovarianceMatrix::fill_symmetriccorrelationmatrix_to_diagonal_subblock( const CorrelationMatrix& corr, const size_t srow )
{
    size_t scol = srow; 
    if( (int)corr.rows() > (int)m_cov.rows() - (int)srow || (int)corr.rows() > (int)m_cov.rows() - (int)scol ) // check container size
    {
        error::INVALID_SIZE( corr.rows(), __PRETTY_FUNCTION__ );
    }
    for( size_t row = 0; row < corr.rows(); ++row ) 
    {
        std::copy( corr.cbegin( row ) + row, corr.cend( row ), m_cov.begin( srow + row ) + (scol + row) );
    }
}

// diagonalizes the covariance matrix and fills Eigenvalues and Eigenvectors into eig and ortho
inline void CovarianceMatrix::diagonalize( EigenValues& eig, OrthogonalTransformation& ortho ) const
{
    diag::diagonalize_real( m_cov, eig, ortho );
    ortho = blaze::trans( ortho ); // D = O^T M O
}


// =======================================================
// === IMPLEMENTATION OF DIAG BASIS NORMAL DIST CLASS ====
// =======================================================
// fill the normal distributions setting the averages to zero and the standard deviations to the square roots of the eigenvalues handed over
inline void DiagonalBasisNormalDistributions::fill( const EigenValues& eig )
{
    m_ndist.clear();
    std::transform( eig.cbegin(), eig.cend(), std::back_inserter( m_ndist ),
    []( const RealType& eig ) -> std::normal_distribution<RealType>
    {
        if( eig < RealType{0.} )
        {
            error::NEGATIVE_VARIANCE( eig, std::string(__PRETTY_FUNCTION__) );
        }
        return std::normal_distribution<RealType>( RealType{0.}, std::sqrt(eig) );
    } );
}

// print
inline void DiagonalBasisNormalDistributions::print() const
{
    std::cout << "Containing " << m_ndist.size() << " normal distributions with zero mean and stddevs\n";
    size_t index{1};
    std::for_each( m_ndist.cbegin(), m_ndist.cend(),
    [&index]( const auto& dist ){
        std::cout << "sigma_" << index++ << " = " << dist.stddev() << "\n";
    } );
}


// =======================================================
// ==== IMPLEMENTATION OF GAUSSIAN NOISE VECTORS CLASS ===
// =======================================================
// fills each noise vector (column of m_noise) from a set of normal distributions and Mersenne Twister random generator
inline void GaussianNoiseVectors::fill( DiagonalBasisNormalDistributions& ndist, std::mt19937& engine, const size_t num_noises )
{
    m_noise.resize( ndist.size(), num_noises );
    for( size_t column = 0; column < num_noises; ++column )
    {
        stda::for_2each( ndist.begin(), ndist.end(), m_noise.begin( column ),
        [&engine]( auto& dist, auto& n )
        {
            n = dist(engine);
        } );
    }
}


};