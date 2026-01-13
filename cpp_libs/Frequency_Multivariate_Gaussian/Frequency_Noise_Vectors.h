#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<random>
#include<memory>
#include<blaze/Math.h>
#include"../Globals/Types.h"
#include"../Standard_Algorithms/Standard_Algorithms.h"

namespace Frequency_Multivariate_Gaussian
{

// ============================ NAMESPACES ==============================
namespace stda = Standard_Algorithms;

// ========================= USING STATEMENTS ===========================
using Vector = blaze::StaticVector<RealType,3UL,blaze::columnVector>;
using NormalDistributions = std::vector<std::normal_distribution<RealType>>;

// ======================================================================
// ================ HEADER FOR NORMAL DISTRIBUTION CLASS ================
// ======================================================================
/* builds a vector of normal distributions (wrapper around NormalDistributions) setting their mean values to zero and their variances from a vector of eigenvalues; throws an error if some eigenvalues are negative */
class NormalDistributionsList
{
public:
    // CONSTRUCTORS
    NormalDistributionsList();
    explicit NormalDistributionsList( const EigenValuesList& eig ) { this->fill( eig ); };

    // PUBLIC METHODS
    void fill( const EigenValuesList& eig );
    void resize( const size_t size ){ m_ndist.resize( size ); }
    const size_t size() const { return m_ndist.size(); }
    void print() const;

    // ITERATORS
    auto begin() { return m_ndist.begin(); };
    auto end() { return m_ndist.end(); };
    auto cbegin() const { return m_ndist.cbegin(); };
    auto cend() const { return m_ndist.cend(); };

    // OPERATORS
    NormalDistributions& operator()( const size_t at ) { return m_ndist[at]; };
    const NormalDistributions& operator()( const size_t at ) const { return m_ndist[at]; };

    private:
    // PRIVATE MEMBERS
    std::vector<NormalDistributions> m_ndist{};
};

// ======================================================================
// ================= HEADER FOR NOISE VECTORS CLASS  ====================
// ======================================================================
class FrequencyNoiseVectors
{

public:
    // CONSTRUCTORS
    FrequencyNoiseVectors() = default;
    FrequencyNoiseVectors( NormalDistributionsList& dist, OrthogonalTransformationList& ortho, std::mt19937& engine, const size_t num_noises );
    
    // PUBLIC METHODS
    inline void draw( NormalDistributionsList& dist, OrthogonalTransformationList& ortho, std::mt19937& engine, const size_t num_noises );
    inline void fourier_back_transform();
    // void print() const { std::cout << m_noise; };

    // OPERATORS
    std::vector<Vector>& operator()( const size_t noise_vector_index ) { return m_noise[noise_vector_index]; };
    const std::vector<Vector>& operator()( const size_t noise_vector_index ) const { return m_noise[noise_vector_index]; };


private:
    // PRIVATE MEMBERS
    std::vector<std::vector<Vector>> m_noise{};

};

// ======================================================================
// ============== IMPLEMENTATION OF NOISE VECTORS CLASS  ================
// ======================================================================

inline FrequencyNoiseVectors::FrequencyNoiseVectors( NormalDistributionsList& dist, OrthogonalTransformationList& ortho, std::mt19937& engine, const size_t num_noises )
{
    this->draw( dist, ortho, engine, num_noises );
    this->fourier_back_transform();
}

inline void FrequencyNoiseVectors::draw( NormalDistributionsList& dist, OrthogonalTransformationList& ortho, std::mt19937& engine, const size_t num_noises )
{
    m_noise.resize( num_noises );
    std::for_each( m_noise.begin(), m_noise.end(), [&]( auto& noise_sample )
    {
        noise_sample.resize( ortho.size() );
        for( size_t block = 0; block < ortho.size(); ++block )
        {
            Vector diag_noise;
            for(int i=0; i<3; i++)
            {
                diag_noise[i] = dist(block)[i]( engine );
            }

            noise_sample[block] = ortho[block] * diag_noise;
        }
    } );
}

inline void FrequencyNoiseVectors::fourier_back_transform()
{
    std::for_each( m_noise.begin(), m_noise.end(), [&]( auto& noise_sample )
    {
        size_t N = noise_sample.size();
            std::vector<Vector> ft_sample;
            for( size_t t = 0; t < N; t++ )
            {
                Vector result{0., 0., 0.};
                result += noise_sample[0] / std::sqrt(RealType(N)); // n = 0 term
                for( size_t n = 1; n < N/2 + 1; n++ )
                {
                    result += noise_sample[n] * std::cos( (2.*M_PI*t*n)/RealType(N) ) / std::sqrt(RealType(N)/RealType(2.));
                }
                for( size_t n = N/2 + 1; n < N; n++ )
                {
                    result += - noise_sample[n] * std::sin( (2.*M_PI*t*n)/RealType(N) ) / std::sqrt(RealType(N)/RealType(2.));
                }
                ft_sample.push_back( result );
            }

            noise_sample = ft_sample;
    } );
}

// ======================================================================
// ================= IMPLEMENTATION OF NORMAL DIST CLASS ================
// ======================================================================
// fill the normal distributions setting the averages to zero and the standard deviations to the square roots of the eigenvalues handed over
inline void NormalDistributionsList::fill( const EigenValuesList& eig )
{
    m_ndist.clear();
    std::for_each( eig.cbegin(), eig.cend(),
    [&]( const auto& evals ){
        NormalDistributions dist_vec;
        for(int i=0; i<3; i++)
        {
            dist_vec.emplace_back( 0., std::sqrt( evals[i] ) );
        }
        m_ndist.push_back( dist_vec );
    } );
}

// print the normal distributions
inline void NormalDistributionsList::print() const
{
    std::cout << "Containing " << m_ndist.size() << " normal distributions with zero mean and stddevs\n";
    size_t index{1};
    std::for_each( m_ndist.cbegin(), m_ndist.cend(),
    [&index]( const auto& dist ){
        for(int i=0; i<3; i++)
        {
            std::cout << "sigma_" << index++ << " = " << dist[i].stddev() << "\n";
        }
    } );
}

}