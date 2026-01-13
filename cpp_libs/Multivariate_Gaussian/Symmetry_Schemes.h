/* No classes or functions in the namespaces Multivariate_Gaussian and Multivariate_Gaussian::Blocks know
about symmetrytypes, since they are formulated in an abstract way. In order to fill the abstract 
class instances with spin or mean-field correlations, one requires filling schemes, who do know about the 
symmetry types; these schemes are implemented here */
#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<memory>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"../Standard_Algorithms/Standard_Algorithms.h"
#include"Multivariate_Gaussian.h"
#include"Multivariate_Gaussian_Blocks.h"
#include"MVG_Error_Handling.h"

namespace Multivariate_Gaussian::Blocks::Symmetry_Schemes
{

// ================== NAMESPACES =========================================
namespace error = Multivariate_Gaussian::Error_Handling;
namespace stda = Standard_Algorithms;

// ================== USING STATEMENTS ===================================
using Vec = blaze::StaticVector<RealType,3UL,blaze::columnVector>;
using VecVec = std::vector<Vec>;


// =======================================================================
// ====== HEADER FOR CORRELATIONVECTOR TENSOR FILLING SCHEME CLASS  ======
// =======================================================================
// HEADER FOR 
// This is a class for filling correlationvector tensors to covariance matrix blocks
// CorrelationVectorTensor needs to be an iterable container of iterable vectors
template<typename CorrelationVectorTensor>
class CorrelationVectorTensorFillingScheme
{
 public:
    // CONSTRUCTORS
    CorrelationVectorTensorFillingScheme( const std::shared_ptr<const CorrelationVectorTensor>& ten_ptr, const char symmetry_type );

    // PUBLIC METHODS
    void fill( std::vector<CovarianceMatrix>& cov );

 private:
    // PRIVATE MEMBERS
    const char m_symmetry_type{};
    std::shared_ptr<const CorrelationVectorTensor> m_container{nullptr};

    // PRIVATE CONSTRUCTORS 
    CorrelationVectorTensorFillingScheme();
};


// =======================================================================
// === HEADER FOR CORRELATIONVECTOR TENSOR CLUSTER FILLING SCHEME CLASS ==
// =======================================================================
// This is a class for filling correlationvector tensors to covariance matrix blocks
// CorrelationVectorTensorCluster needs to be an iterable container of iterable vectors
template<typename CorrelationVectorTensorCluster>
class CorrelationVectorTensorClusterFillingScheme
{
 public:
    // CONSTRUCTORS
    CorrelationVectorTensorClusterFillingScheme( const std::shared_ptr<const CorrelationVectorTensorCluster>& clu_ptr, const char symmetry_type, const size_t num_Spins );

    // PUBLIC METHODS
    void fill( std::vector<CovarianceMatrix>& cov );

 private:
    // PRIVATE MEMBERS
    const char m_symmetry_type{};
    const size_t m_num_Spins{};
    std::shared_ptr<const CorrelationVectorTensorCluster> m_container{nullptr};

    // PRIVATE CONSTRUCTORS 
    CorrelationVectorTensorClusterFillingScheme();
};


// =======================================================================
// ========== HEADER FOR CORRELATIONMATRIX FILLING SCHEME CLASS ==========
// =======================================================================
// This is a class for a filling scheme filling correlationmatrices to covariance matrix blocks
// CorrelationMatrixTensor needs to be an iterable container of rowwise iterable square matrices
template<typename CorrelationMatrixTensor>
class CorrelationMatrixTensorFillingScheme
{
 public:
    // CONSTRUCTORS 
    CorrelationMatrixTensorFillingScheme( const std::shared_ptr<const CorrelationMatrixTensor>& ten_ptr, const char symmetry_type );

    // PUBLIC METHODS
    void fill( std::vector<CovarianceMatrix>& cov );

 private:
    // PRIVATE MEMBERS
    const char m_symmetry_type{};
    std::shared_ptr<const CorrelationMatrixTensor> m_container{nullptr};

    // PRIVATE CONSTRUCTORS 
    CorrelationMatrixTensorFillingScheme(); // shall not be called
};


// =======================================================================
// ============ HEADER FOR NOISE TENSOR READER SCHEME CLASS ==============
// =======================================================================
// This is a class for reading out noise vectors from a GaussianNoiseVectorsBlocks instance
class NoiseTensorReaderScheme
{
public:
    // CONSTRUCTORS 
     NoiseTensorReaderScheme( const std::shared_ptr<const GaussianNoiseVectorsBlocks>& noise, const char symmetry_type, const size_t num_Samples, const size_t noise_size );

    // PUBLIC METHODS
    Vec read();

private:
    // PRIVATE MEMBERS
    const char m_symmetry_type{};
    const size_t m_num_Samples{};
    const size_t m_noise_size{};
    int m_sample_counter{};
    size_t m_component_counter{};
    std::shared_ptr<const GaussianNoiseVectorsBlocks> m_noise{nullptr};
    NoiseMatrix::ConstIterator m_ptr_x{};
    NoiseMatrix::ConstIterator m_ptr_y{};
    NoiseMatrix::ConstIterator m_ptr_z{};

    // PRIVATE METHODS
    void reset_ptrs();

    // PRIVATE CONSTRUCTORS 
    NoiseTensorReaderScheme(); // shall not be called
};


// =======================================================================
// ========= HEADER FOR NOISE TENSOR CLUSTER READER SCHEME CLASS =========
// =======================================================================
class NoiseTensorClusterReaderScheme
{
public:
    // CONSTRUCTORS 
     NoiseTensorClusterReaderScheme( const std::shared_ptr<const GaussianNoiseVectorsBlocks>& noise, const char symmetry_type, const size_t num_Samples, const size_t sub_noise_size, const size_t num_noises );

    // PUBLIC METHODS
    VecVec read();

private:
    // PRIVATE MEMBERS
    const char m_symmetry_type{};
    const size_t m_num_Samples{};
    const size_t m_sub_noise_size{};
    const size_t m_num_noises{};
    int m_sample_counter{};
    size_t m_component_counter{};
    std::shared_ptr<const GaussianNoiseVectorsBlocks> m_noise{nullptr};
    std::vector<NoiseMatrix::ConstIterator> m_ptrs_x{};
    std::vector<NoiseMatrix::ConstIterator> m_ptrs_y{};
    std::vector<NoiseMatrix::ConstIterator> m_ptrs_z{};

    // PRIVATE METHODS
    void reset_ptrs();

    // PRIVATE CONSTRUCTORS 
    NoiseTensorClusterReaderScheme(); // shall not be called
};


// =======================================================================
// =========================== HELPER FUNCTIONS ==========================
// =======================================================================
// initialize a vector of covariance matrix blocks to match a certain symmetry type
inline void initialize_covariance( std::vector<CovarianceMatrix>& cov, const char symmetry_type, size_t base_size )
{
    switch( symmetry_type )
    {
        case 'A':
        {
            cov.resize( 1 );
            cov[0].resize( base_size ); // Block 0 : xx
            break;
        }
        case 'B':
        {
            cov.resize( 2 );
            cov[0].resize( base_size ); // Block 0 : xx
            cov[1].resize( base_size ); // Block 1 : zz
            break;
        }
        case 'C':
        {
            cov.resize( 2 );
            cov[0].resize( 2*base_size ); // Block 0 : xx xy yx yy
            cov[1].resize( base_size );   // Block 1 : zz
            break;
        }
        case 'D':
        {
            cov.resize( 1 );
            cov[0].resize( 3*base_size ); // Block 0 : xx xy xz yx yy yz zx zy zz
            break;
        }
        default:
            error::SYMMETRY_TYPE( symmetry_type, __PRETTY_FUNCTION__ );
    }   
}


// determine the number of noise samples per block depending on the symmetry type
inline std::vector<size_t> return_num_noises_per_block( const char symmetry_type, const size_t num_Samples )
{
    std::vector<size_t> num_noises_per_block;
    switch( symmetry_type )    
    {
        case 'A':
        {
            num_noises_per_block = std::vector<size_t>{ 3*num_Samples }; // Vx & Vy & Vz -> block 0
            break;
        }
        case 'B':
        {
            num_noises_per_block = std::vector<size_t>{ 2*num_Samples, num_Samples }; // Vx & Vy -> block 0 || Vz -> block 1
            break;
        }
        case 'C':
        {
            num_noises_per_block = std::vector<size_t>{ num_Samples, num_Samples }; // Vxy -> block 0 || Vz -> block 1
            break;
        }
        case 'D':
        {
            num_noises_per_block = std::vector<size_t>{ num_Samples }; // Vxyz -> block 0
            break;
        }
        default:
            error::SYMMETRY_TYPE( symmetry_type, __PRETTY_FUNCTION__ );
    }
    return num_noises_per_block;
}


// =======================================================================
// === IMPLEMENTATION OF CORRELATIONVECTOR TENSOR FILLING SCHEME CLASS ===
// =======================================================================
// constructor from a CorrelationVectorTensor(Tensor) and symmetry_type
template<typename CorrelationVectorTensor>
CorrelationVectorTensorFillingScheme<CorrelationVectorTensor>::CorrelationVectorTensorFillingScheme( const std::shared_ptr<const CorrelationVectorTensor>& ten_ptr, const char symmetry_type ):
    m_symmetry_type( symmetry_type ),
    m_container( ten_ptr )
{}


// resize and fill std::vector<CovarianceMatrix> from the stored CorrelationVectorTensor
template<typename CorrelationVectorTensor>
void CorrelationVectorTensorFillingScheme<CorrelationVectorTensor>::fill( std::vector<CovarianceMatrix>& cov )
{
    size_t bsize = m_container->get_xx().size(); // base size
    initialize_covariance( cov, m_symmetry_type, bsize );
    switch( m_symmetry_type )
    {
        case 'A':
        {
            cov[0].fill_correlationvector( m_container->get_xx() ); // Block 0 : xx
            break;
        }
        case 'B':
        {
            cov[0].fill_correlationvector( m_container->get_xx() ); // Block 0 : xx
            cov[1].fill_correlationvector( m_container->get_zz() ); // Block 1 : zz
            break;
        }
        case 'C':
        {
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_xx(), 0, 0          );  
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_xy(), 0, bsize      ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_yx(), bsize, 0      ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_yy(), bsize, bsize  );
            cov[1].fill_correlationvector( m_container->get_zz() ); // Block 1 : zz
            break;
        }
        case 'D':
        {
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_xx(), 0, 0            ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_xy(), 0, bsize        );
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_xz(), 0, 2*bsize      ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_yx(), bsize, 0        ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_yy(), bsize, bsize    );
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_yz(), bsize, 2*bsize  ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_zx(), 2*bsize, 0      ); 
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_zy(), 2*bsize, bsize  );
            cov[0].fill_correlationvector_to_subtriangle( m_container->get_zz(), 2*bsize, 2*bsize);
            break;
        }
        default:
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    }
}


// =======================================================================
// ======== CORRELATIONVECTOR TENSOR CLUSTER FILLING SCHEME CLASS ========
// =======================================================================
// constructor from a CorrelationVectorTensor(Tensor) and symmetry_type
template<typename CorrelationVectorTensorCluster>
CorrelationVectorTensorClusterFillingScheme<CorrelationVectorTensorCluster>::CorrelationVectorTensorClusterFillingScheme( const std::shared_ptr<const CorrelationVectorTensorCluster>& clu_ptr, const char symmetry_type, const size_t num_Spins ):
    m_symmetry_type( symmetry_type ),
    m_num_Spins( num_Spins ),
    m_container( clu_ptr )
{}

// resize and fill std::vector<CovarianceMatrix> from the stored CorrelationVectorTensor
template<typename CorrelationVectorTensorCluster>
void CorrelationVectorTensorClusterFillingScheme<CorrelationVectorTensorCluster>::fill( std::vector<CovarianceMatrix>& cov )
{
    size_t vsize = m_container->operator[](0).get_xx().size(); // number of time points
    size_t bsize = vsize * m_num_Spins; // base size
    initialize_covariance( cov, m_symmetry_type, bsize );
    for( size_t i=0; i<m_num_Spins; ++i )
    {
        for( size_t j=0; j<m_num_Spins; ++j )
        {
            switch( m_symmetry_type )
            {
                case 'A':
                {
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xx(), i*vsize, j*vsize ); // Block 0 : xx
                    break;
                }
                case 'B':
                {
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xx(), i*vsize, j*vsize ); // Block 0 : xx
                    cov[1].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_zz(), i*vsize, j*vsize ); // Block 1 : zz
                    break;
                }
                case 'C':
                {
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xx(), i*vsize, j*vsize ); // Block 0 ...
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xy(), i*vsize, j*vsize + bsize ); 
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_yx(), i*vsize + bsize, j*vsize ); 
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_yy(), i*vsize + bsize, j*vsize + bsize );
                    cov[1].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_zz(), i*vsize, j*vsize ); // Block 1 : zz
                    break;
                }
                case 'D':
                {
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xx(), i*vsize, j*vsize         );
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xy(), i*vsize, j*vsize +bsize  );
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_xz(), i*vsize, j*vsize +2*bsize);
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_yx(), i*vsize +bsize, j*vsize         );
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_yy(), i*vsize +bsize, j*vsize +bsize  );
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_yz(), i*vsize +bsize, j*vsize +2*bsize);
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_zx(), i*vsize +2*bsize, j*vsize       ); 
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_zy(), i*vsize +2*bsize, j*vsize+bsize  );
                    cov[0].fill_correlationvector_to_subtriangle( m_container->operator()(i,j).get_zz(), i*vsize +2*bsize, j*vsize+2*bsize);
                    break;
                }
                default:
                    error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
            }
        }
    }
}


// =======================================================================
// ============ CORRELATIONMATRIX TENSOR FILLING SCHEME CLASS ============
// =======================================================================
// constructor from a CorrelationMatrixTensor(Tensor) and symmetry_type
template<typename CorrelationMatrixTensor>
 CorrelationMatrixTensorFillingScheme<CorrelationMatrixTensor>:: CorrelationMatrixTensorFillingScheme( const std::shared_ptr<const CorrelationMatrixTensor>& ten_ptr, const char symmetry_type ):
    m_symmetry_type( symmetry_type ),
    m_container( ten_ptr )
{}

// resize and fill std::vector<CovarianceMatrix> from the stored CorrelationMatrixTensor
template<typename CorrelationMatrixTensor> 
void CorrelationMatrixTensorFillingScheme<CorrelationMatrixTensor>::fill( std::vector<CovarianceMatrix>& cov )
{
    size_t bsize = m_container->get_xx().rows(); // sub block size
    initialize_covariance( cov, m_symmetry_type, bsize );
    switch( m_symmetry_type )
    {
        case 'A':
        {
            cov[0].fill_correlationmatrix( m_container->get_xx() ); // Block 0 : xx
            break;
        }
        case 'B':
        {
            cov[0].fill_correlationmatrix( m_container->get_xx() ); // Block 0 : xx
            cov[1].fill_correlationmatrix( m_container->get_zz() ); // Block 1 : zz
            break;
        }
        case 'C':
        {
            cov[0].fill_symmetriccorrelationmatrix_to_diagonal_subblock(  m_container->get_xx(), 0            );  
            cov[0].fill_correlationmatrix_to_subsquare(             m_container->get_xy(), 0, bsize           ); 
            cov[0].fill_symmetriccorrelationmatrix_to_diagonal_subblock(  m_container->get_yy(), bsize, bsize );
            cov[1].fill_correlationmatrix(                          m_container->get_zz() ); // Block 1 : zz
            break;
        }
        case 'D':
        {
            cov[0].fill_symmetriccorrelationmatrix_to_diagonal_subblock(  m_container->get_xx(), 0        ); 
            cov[0].fill_correlationmatrix_to_subsquare(             m_container->get_xy(), 0, bsize       );
            cov[0].fill_correlationmatrix_to_subsquare(             m_container->get_xz(), 0, 2*bsize     ); 
            cov[0].fill_symmetriccorrelationmatrix_to_diagonal_subblock(  m_container->get_yy(), bsize    );
            cov[0].fill_correlationmatrix_to_subsquare(             m_container->get_yz(), bsize, 2*bsize ); 
            cov[0].fill_symmetriccorrelationmatrix_to_diagonal_subblock(  m_container->get_zz(), 2*bsize  );
            break;
        }
        default:
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    }
}


// =======================================================================
// ================== NOISE TENSOR READER SCHEME CLASS ===================
// =======================================================================
// constructor
inline NoiseTensorReaderScheme::NoiseTensorReaderScheme( const std::shared_ptr<const GaussianNoiseVectorsBlocks>& noise, const char symmetry_type, const size_t num_Samples, const size_t noise_size ):
    m_symmetry_type( symmetry_type ),
    m_num_Samples( num_Samples ),
    m_noise_size( noise_size ),
    m_noise( noise )
{
    m_sample_counter = -1;
    reset_ptrs();
}

// reset iterators
inline void NoiseTensorReaderScheme::reset_ptrs()
{
    m_component_counter = 0;
    m_sample_counter++;
    switch( m_symmetry_type )
    {
        case 'A':
        {
            m_ptr_x = m_noise->operator[](0).cbegin(m_sample_counter);
            m_ptr_y = m_noise->operator[](0).cbegin(m_sample_counter+m_num_Samples);
            m_ptr_z = m_noise->operator[](0).cbegin(m_sample_counter+2*m_num_Samples);
            break;
        }
        case 'B':
        {
            m_ptr_x = m_noise->operator[](0).cbegin(m_sample_counter);
            m_ptr_y = m_noise->operator[](0).cbegin(m_sample_counter+m_num_Samples);
            m_ptr_z = m_noise->operator[](1).cbegin(m_sample_counter);
            break;
        }
        case 'C':
        {
            m_ptr_x = m_noise->operator[](0).cbegin(m_sample_counter);
            m_ptr_y = m_noise->operator[](0).cbegin(m_sample_counter)+m_noise_size;
            m_ptr_z = m_noise->operator[](1).cbegin(m_sample_counter);
            break;
        }
        case 'D':
        {
            m_ptr_x = m_noise->operator[](0).cbegin(m_sample_counter);
            m_ptr_y = m_noise->operator[](0).cbegin(m_sample_counter)+m_noise_size;
            m_ptr_z = m_noise->operator[](0).cbegin(m_sample_counter)+2*m_noise_size;
            break;
        }
        default: // type invalid
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    };
}

// reader methods
inline Vec NoiseTensorReaderScheme::read()
{
    if( m_component_counter != m_noise_size )
    {
        m_component_counter++;
        return Vec{ *m_ptr_x++, *m_ptr_y++, *m_ptr_z++ };
    }
    else // return and go to new sample
    {
        if( static_cast<size_t>(m_sample_counter) < m_num_Samples-1 ) // not last sample
        {
            reset_ptrs();
            return read();
        }
        else
        {
            return Vec{ *m_ptr_x, *m_ptr_y, *m_ptr_z }; // no increment
        }
    }
}


// =======================================================================
// =============== NOISE TENSOR CLUSTER READER SCHEME CLASS ==============
// =======================================================================
// constructor
inline NoiseTensorClusterReaderScheme::NoiseTensorClusterReaderScheme( const std::shared_ptr<const GaussianNoiseVectorsBlocks>& noise, const char symmetry_type, const size_t num_Samples, const size_t sub_noise_size, const size_t num_noises ):
    m_symmetry_type( symmetry_type ),
    m_num_Samples( num_Samples ),
    m_sub_noise_size( sub_noise_size ),
    m_num_noises( num_noises ),
    m_noise( noise )
{
    m_ptrs_x.resize( num_noises );
    m_ptrs_y.resize( num_noises );
    m_ptrs_z.resize( num_noises );
    m_sample_counter = -1;
    reset_ptrs();
}

// reset iterators
inline void NoiseTensorClusterReaderScheme::reset_ptrs()
{
    m_component_counter = 0;
    m_sample_counter++;
    switch( m_symmetry_type )
    {
        case 'A':
        {
            stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&,noise_index=0]( auto& ptr_x, auto& ptr_y, auto& ptr_z ) mutable
            {
                ptr_x = m_noise->operator[](0).cbegin(m_sample_counter) + noise_index*m_sub_noise_size;
                ptr_y = m_noise->operator[](0).cbegin(m_sample_counter+m_num_Samples) + noise_index*m_sub_noise_size;
                ptr_z = m_noise->operator[](0).cbegin(m_sample_counter+2*m_num_Samples) + noise_index++*m_sub_noise_size;
            } );
            break;
        }
        case 'B':
        {
            stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&,noise_index=0]( auto& ptr_x, auto& ptr_y, auto& ptr_z ) mutable
            {
                ptr_x = m_noise->operator[](0).cbegin(m_sample_counter) + noise_index*m_sub_noise_size;
                ptr_y = m_noise->operator[](0).cbegin(m_sample_counter+m_num_Samples) + noise_index*m_sub_noise_size;
                ptr_z = m_noise->operator[](1).cbegin(m_sample_counter) + noise_index++*m_sub_noise_size;
            } );
            break;
        }
        case 'C':
        {
            stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&,noise_index=0]( auto& ptr_x, auto& ptr_y, auto& ptr_z ) mutable
            {
                ptr_x = m_noise->operator[](0).cbegin(m_sample_counter) + noise_index*m_sub_noise_size;
                ptr_y = m_noise->operator[](0).cbegin(m_sample_counter) + (m_num_noises+noise_index)*m_sub_noise_size;
                ptr_z = m_noise->operator[](1).cbegin(m_sample_counter) + noise_index++*m_sub_noise_size;
            } );
            break;
        }
        case 'D':
        {
            stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&,noise_index=0]( auto& ptr_x, auto& ptr_y, auto& ptr_z ) mutable
            {
                ptr_x = m_noise->operator[](0).cbegin(m_sample_counter) + noise_index*m_sub_noise_size;
                ptr_y = m_noise->operator[](0).cbegin(m_sample_counter) + (m_num_noises+noise_index)*m_sub_noise_size;
                ptr_z = m_noise->operator[](0).cbegin(m_sample_counter) + (2*m_num_noises+noise_index++)*m_sub_noise_size;
            } );
            break;
        }
        default: // type invalid
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    };
}

// reader methods
inline VecVec NoiseTensorClusterReaderScheme::read()
{
    if( m_component_counter != m_sub_noise_size )
    {
        m_component_counter++;
        VecVec Vs{};
        stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&]( auto& ptr_x, auto& ptr_y, auto& ptr_z )
        { 
            Vs.emplace_back( Vec{ *ptr_x++, *ptr_y++, *ptr_z++ } );
        } );
        return Vs;
    }
    else // return and go to new sample
    {
        if( static_cast<size_t>(m_sample_counter) < m_num_Samples-1 ) // not last sample
        {
            reset_ptrs();
            return read();
        }
        else
        {
            VecVec Vs{};
            stda::for_3each( m_ptrs_x.begin(), m_ptrs_x.end(), m_ptrs_y.begin(), m_ptrs_z.begin(), [&]( auto& ptr_x, auto& ptr_y, auto& ptr_z )
            { 
                Vs.emplace_back( Vec{ *ptr_x, *ptr_y, *ptr_z } ); // no increment
            } );
            return Vs;
        }
    }
}

};
