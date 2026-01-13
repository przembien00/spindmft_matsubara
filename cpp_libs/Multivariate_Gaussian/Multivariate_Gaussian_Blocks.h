#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"../Standard_Algorithms/Standard_Algorithms.h"
#include"Multivariate_Gaussian.h"
#include"MVG_Error_Handling.h"

namespace Multivariate_Gaussian::Blocks
{

// ================== NAMESPACES =======================
namespace error = Multivariate_Gaussian::Error_Handling;
namespace stda = Standard_Algorithms;

// ================== USING STATEMENTS =================
using EigenValuesBlocks = std::vector<EigenValues>;
using OrthogonalTransformationBlocks = std::vector<OrthogonalTransformation>;
using SymmetricMatrixBlocks = std::vector<SymmetricMatrix>;


// =====================================================
// ===== HEADER FOR COVARIANCE MATRIX BLOCKS CLASS =====
// =====================================================
// This is simply a wrapper around a std::vector of CovarianceMatrix
class CovarianceMatrixBlocks
{
public:
    // CONSTRUCTORS
    CovarianceMatrixBlocks() = default;
    explicit CovarianceMatrixBlocks( const SymmetricMatrixBlocks& matb );
    explicit CovarianceMatrixBlocks( SymmetricMatrixBlocks&& matb );

    // PUBLIC METHODS
    template<typename FillingScheme>
    void fill_from_scheme( FillingScheme& scheme ) { scheme.fill( m_cov ); };
    void diagonalize( EigenValuesBlocks& eig, OrthogonalTransformationBlocks& ortho ) const;
    void resize( const std::vector<size_t>&& block_sizes );
    const size_t size() const { return m_cov.size(); };
    void print() const;

    // ITERATORS
    auto begin() { return m_cov.begin(); };
    auto end() { return m_cov.end(); };
    auto cbegin() const { return m_cov.cbegin(); };
    auto cend() const { return m_cov.cend(); };

private:
    // PRIVATE MEMBERS
    std::vector<CovarianceMatrix> m_cov{};
};


// =====================================================
// ==== HEADER FOR NORMAL DISTRIBUTION BLOCKS CLASS ====
// =====================================================
// this is a wrapper around std::vector<DiagonalBasisNormalDistributions> with zero mean value
class DiagonalBasisNormalDistributionsBlocks
{
public:
    // CONSTRUCTORS
    DiagonalBasisNormalDistributionsBlocks() = default; 
    explicit DiagonalBasisNormalDistributionsBlocks( const EigenValuesBlocks& eig ){ this->fill( eig ); };

    // PUBLIC METHODS
    void fill( const EigenValuesBlocks& eig );
    void resize( const std::vector<size_t>& block_sizes );
    const size_t size() const { return m_ndist.size(); };
    void print() const;

    // ITERATORS
    auto begin() { return m_ndist.begin(); };
    auto end() { return m_ndist.end(); };
    auto cbegin() const { return m_ndist.cbegin(); };
    auto cend() const { return m_ndist.cend(); };

private:
    // PRIVATE MEMBERS
    std::vector<DiagonalBasisNormalDistributions> m_ndist{};
};


// =====================================================
// === HEADER FOR GAUSSIAN NOISE VECTORS BLOCKS CLASS ==
// =====================================================
/* builds blocks of noise vectors from given blocks of vectors of distributions 
(wrapper around std::vector<GaussianNoiseVectors>)*/
class GaussianNoiseVectorsBlocks
{
public:
    // CONSTRUCTORS
    GaussianNoiseVectorsBlocks() = default;
    GaussianNoiseVectorsBlocks( const std::vector<size_t>& num_noises_per_block,  DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine, const OrthogonalTransformationBlocks& ortho );
    GaussianNoiseVectorsBlocks( const std::vector<size_t>&& num_noises_per_block, DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine, const OrthogonalTransformationBlocks& ortho );

    // PUBLIC METHODS
    void fill( const std::vector<size_t>& num_noises_per_block, DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine );
    void basis_transform( const OrthogonalTransformationBlocks& ortho );
    void print() const;

    // OPERATORS
    GaussianNoiseVectors& operator[]( const size_t index ) { return m_noise[index]; };
    const GaussianNoiseVectors& operator[]( const size_t index ) const { return m_noise[index]; };

    // ITERATORS
    auto begin() { return m_noise.begin(); };
    auto end() { return m_noise.end(); };
    auto cbegin() const { return m_noise.cbegin(); };
    auto cend() const { return m_noise.cend(); };

private:
    // PRIVATE MEMBERS
    std::vector<GaussianNoiseVectors> m_noise{};
};


// =====================================================
// == IMPLEMENTATION OF COVARIANCE MATRIX BLOCKS CLASS =
// =====================================================
// copy constructor for symmetric matrix blocks
inline CovarianceMatrixBlocks::CovarianceMatrixBlocks( const SymmetricMatrixBlocks& matb )
{
    std::transform( matb.cbegin(), matb.cend(), std::back_inserter(m_cov), [this]( const SymmetricMatrix& symm_block )
    {
        return CovarianceMatrix{ symm_block };
    } );
}

// move constructor for symmetric matrix blocks
inline CovarianceMatrixBlocks::CovarianceMatrixBlocks( SymmetricMatrixBlocks&& matb )
{
    std::transform( matb.begin(), matb.end(), std::back_inserter(m_cov), [this]( SymmetricMatrix& symm_block )
    {
        return CovarianceMatrix{ std::move( symm_block ) };
    } );
}

// diagonalization of each block
inline void CovarianceMatrixBlocks::diagonalize( EigenValuesBlocks& eig, OrthogonalTransformationBlocks& ortho ) const
{
    eig.resize( m_cov.size() );
    ortho.resize( m_cov.size() );
    stda::for_3each( m_cov.cbegin(), m_cov.cend(), eig.begin(), ortho.begin(),
    []( const CovarianceMatrix& cov_block, EigenValues& eig_block, OrthogonalTransformation& ortho_block )
    {
        cov_block.diagonalize( eig_block, ortho_block );
    } );
}

// clear and resize
inline void CovarianceMatrixBlocks::resize( const std::vector<size_t>&& block_sizes )
{
    m_cov.clear();
    std::transform( block_sizes.cbegin(), block_sizes.cend(), std::back_inserter( m_cov ), [this]( const size_t block_size )
    {
        return CovarianceMatrix{ block_size };
    } );
}

// print
inline void CovarianceMatrixBlocks::print() const
{
    size_t block_index{};
    std::for_each( m_cov.cbegin(), m_cov.cend(),
    [&block_index]( const CovarianceMatrix& cov_block )
    {
        std::cout << "Covariance Matrix Block No. " << block_index++ << ":\n";
        cov_block.print();
        std::cout << "\n";
    } );
}


// =====================================================
// ===== IMPLEMENTATION OF NORMAL DIST BLOCKS CLASS ====
// =====================================================
// clear everything and fill eigen values blocks
inline void DiagonalBasisNormalDistributionsBlocks::fill( const EigenValuesBlocks& eig )
{
    m_ndist.clear();
    std::for_each( eig.cbegin(), eig.cend(),
    [this]( const EigenValues& eig_block )
    {
        m_ndist.emplace_back( DiagonalBasisNormalDistributions{ eig_block } );
    } );
}

// print
inline void DiagonalBasisNormalDistributionsBlocks::print() const
{
    size_t block_index{};
    std::for_each( m_ndist.cbegin(), m_ndist.cend(),
    [&block_index]( const DiagonalBasisNormalDistributions& dist_block )
    {
        std::cout << "Normal Distributions Block No. " << block_index++ << ":\n";
        dist_block.print();
        std::cout << "\n";
    } );
}


// =====================================================
// === IMPLEMENTATION OF GAUSS NOISE VEC BLOCKS CLASS ==
// =====================================================
// constructor from normal distribution blocks
inline GaussianNoiseVectorsBlocks::GaussianNoiseVectorsBlocks( const std::vector<size_t>& num_noises_per_block, DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine, const OrthogonalTransformationBlocks& ortho )
{
    this->fill( num_noises_per_block, ndist, engine );
    this->basis_transform( ortho );
}

// constructor from normal distribution blocks (moved)
inline GaussianNoiseVectorsBlocks::GaussianNoiseVectorsBlocks( const std::vector<size_t>&& num_noises_per_block, DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine, const OrthogonalTransformationBlocks& ortho )
{
    this->fill( num_noises_per_block, ndist, engine );
    this->basis_transform( ortho );
}

// fill from normal distribution blocks
inline void GaussianNoiseVectorsBlocks::fill( const std::vector<size_t>& num_noises_per_block, DiagonalBasisNormalDistributionsBlocks& ndist, std::mt19937& engine )
{    
    if( ndist.size() != num_noises_per_block.size() )
    {
        error::INVALID_SIZE( ndist.size(), __PRETTY_FUNCTION__ );
    }
    m_noise.clear();
    std::transform( ndist.begin(), ndist.end(), num_noises_per_block.cbegin(), std::back_inserter( m_noise ),
    [this,&engine]( auto& dist_block, const size_t num_noises )
    {
        return GaussianNoiseVectors{ dist_block, engine, num_noises };
    } );
}

// transform each block into a new basis
inline void GaussianNoiseVectorsBlocks::basis_transform( const OrthogonalTransformationBlocks& ortho )
{
    stda::for_2each( m_noise.begin(), m_noise.end(), ortho.cbegin(),
    []( GaussianNoiseVectors& noise_block, const OrthogonalTransformation& ortho_block )
    {
        noise_block.basis_transform( ortho_block );
    } );
}

// prints each noise blocks
inline void GaussianNoiseVectorsBlocks::print() const
{
    size_t block_index{};
    std::for_each( m_noise.cbegin(), m_noise.cend(),
    [&block_index]( const GaussianNoiseVectors& noise_block )
    {
        std::cout << "Noise Matrix Block No. " << block_index++ << ":\n";
        noise_block.print();
        std::cout << "\n";
    } );
}


};
