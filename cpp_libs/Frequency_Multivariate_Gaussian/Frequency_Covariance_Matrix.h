#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<memory>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"MVG_Error_Handling.h"
#include"../Matrices/Diagonalization.h"
#include"../Standard_Algorithms/Standard_Algorithms.h"

namespace Frequency_Multivariate_Gaussian
{

// ============================ NAMESPACES ===============================
namespace error = Frequency_Multivariate_Gaussian::Error_Handling;
namespace stda = Standard_Algorithms;
namespace diag = Matrices::Diagonalization;

// ========================= USING STATEMENTS ============================
using EigenValues = blaze::StaticVector<RealType,3UL,blaze::columnVector>;
using SymmetricMatrix = blaze::SymmetricMatrix<blaze::StaticMatrix<RealType,3UL,3UL,blaze::rowMajor>>;
using OrthogonalTransformation = blaze::StaticMatrix<RealType,3UL,3UL,blaze::rowMajor>;
using EigenValuesList = std::vector<EigenValues>;
using SymmetricMatrixList = std::vector<SymmetricMatrix>;
using OrthogonalTransformationList = std::vector<OrthogonalTransformation>;


// =======================================================================
// ============ HEADER FOR FREQUENCY COVARIANCE MATRIX CLASS  ============
// =======================================================================
class FrequencyCovarianceMatrix
{
public:
    // CONSTRUCTORS
    FrequencyCovarianceMatrix() = default;
    explicit FrequencyCovarianceMatrix( const size_t& size ) { resize( size ); };
    template<typename CorrelationTensor> // needs to be an iterable object
    explicit FrequencyCovarianceMatrix( const CorrelationTensor& corr, const char symmetry_type ) { this->fourier_fill_correlationtensor( corr, symmetry_type ); };
    explicit FrequencyCovarianceMatrix( const SymmetricMatrixList& mat ) { m_cov = mat; };

    // PUBLIC METHODS
    template<typename CorrelationTensor> // needs to be an iterable object
    void fourier_fill_correlationtensor( const CorrelationTensor& corr, const char symmetry_type );
    template<typename CorrelationVector> // needs to be an iterable object
    std::vector<RealType> fourier_transform( const CorrelationVector& corr );
    inline void diagonalize( EigenValuesList& eig, OrthogonalTransformationList& ortho ) const;
    void resize( const size_t size ) { m_cov.resize( size ); };
    const size_t size() const { return m_cov.size(); };

    // void print() const { std::cout << m_cov; };

private:
    // PRIVATE MEMBERS
    SymmetricMatrixList m_cov{};
};

// =======================================================================
// ========================== HELPER FUNCTIONS ===========================
// =======================================================================

template<typename CorrelationVector>
std::vector<RealType> FrequencyCovarianceMatrix::fourier_transform( const CorrelationVector& corr )
{
    size_t N = corr.size() - 1;
    std::vector<RealType> corr_omega( N, RealType{0.} );
    for( size_t k = 0; k < N; ++k )
    {
        for( size_t n = 0; n < N; ++n )
        {
            corr_omega[k] += corr[n] * std::cos( (2.*M_PI*k*n)/RealType(N) );
        }
        // corr_omega[k] /= RealType(N);
    }
    return corr_omega;
}

// ========================================================================
// ========== IMPLEMENTATION OF FREQUENCY COVARIANCE MATRIX CLASS =========
// ========================================================================
// resize and fill std::vector<SymmetricMatrix> from a CorrelationTensor
template<typename CorrelationTensor>
void FrequencyCovarianceMatrix::fourier_fill_correlationtensor( const CorrelationTensor& corr, const char symmetry_type )
{
    size_t bsize = corr->get_xx().size() - 1; // base size
    resize( bsize );
    switch( symmetry_type )
    {
        case 'A':
        {
            std::vector<RealType> C_omega = fourier_transform( corr->get_xx() );
            for( size_t i=0; i<bsize; ++i )
            {
                m_cov[i] = SymmetricMatrix{ {C_omega[i],0.,0.}, 
                                           {0.,C_omega[i],0.}, 
                                           {0.,0.,C_omega[i]} };
            }
            break;
        }
        case 'B':
        {
            std::vector<RealType> C_xx_omega = fourier_transform( corr->get_xx() );
            std::vector<RealType> C_zz_omega = fourier_transform( corr->get_zz() );
            for( size_t i=0; i<bsize; ++i )
            {
                m_cov[i] = SymmetricMatrix{ {C_xx_omega[i],0.,0.},
                                           {0.,C_xx_omega[i],0.},
                                           {0.,0.,C_zz_omega[i]} };
            }
            break;
        }
        case 'C':
        {
            std::vector<RealType> C_xx_omega = fourier_transform( corr->get_xx() );
            std::vector<RealType> C_xy_omega = fourier_transform( corr->get_xy() );
            std::vector<RealType> C_zz_omega = fourier_transform( corr->get_zz() );
            for( size_t i=0; i<bsize; ++i )
            {
                m_cov[i] = SymmetricMatrix{ {C_xx_omega[i], C_xy_omega[i], 0.},
                                           {C_xy_omega[i], C_xx_omega[i], 0.},
                                           {0., 0., C_zz_omega[i]} };
            }
            break;
        }
        case 'D':
        {
            std::vector<RealType> C_xx_omega = fourier_transform( corr->get_xx() );
            std::vector<RealType> C_xy_omega = fourier_transform( corr->get_xy() );
            std::vector<RealType> C_xz_omega = fourier_transform( corr->get_xz() );
            // std::vector<RealType> C_yx_omega = fourier_transform( corr->get_yx() );
            std::vector<RealType> C_yy_omega = fourier_transform( corr->get_yy() );
            std::vector<RealType> C_yz_omega = fourier_transform( corr->get_yz() );
            // std::vector<RealType> C_zx_omega = fourier_transform( corr->get_zx() );
            // std::vector<RealType> C_zy_omega = fourier_transform( corr->get_zy() );
            std::vector<RealType> C_zz_omega = fourier_transform( corr->get_zz() );
            for( size_t i=0; i<bsize; ++i )
            {
                m_cov[i] = SymmetricMatrix{ {C_xx_omega[i], C_xy_omega[i], C_xz_omega[i]},
                                           {C_xy_omega[i], C_yy_omega[i], C_yz_omega[i]},
                                           {C_xz_omega[i], C_yz_omega[i], C_zz_omega[i]} };
            }
            break;
        }
        default:
            error::SYMMETRY_TYPE( symmetry_type, __PRETTY_FUNCTION__ );
    }
}

inline void FrequencyCovarianceMatrix::diagonalize( EigenValuesList& eig, OrthogonalTransformationList& ortho ) const
{
    eig.clear();
    ortho.clear();
    std::for_each( m_cov.cbegin(), m_cov.cend(), 
    [&]( const SymmetricMatrix& mat )
    {
        EigenValues evals;
        OrthogonalTransformation otrafo;
        diag::diagonalize_real( mat, evals, otrafo );
        eig.push_back( evals );
        ortho.push_back( otrafo );
    } );
}

}