/* Dim has to be defined as constexpr size_t to include this header */
#pragma once
#include<iostream>
#include<vector>
#include<functional>
#include<blaze/Math.h>
#include"Standard_Algorithms/Numerics.h"

namespace Physics::Couplings
{

// ===================== NAMESPACES ===========================
namespace num = Standard_Algorithms::Numerics;

// ===================== USING STATEMENTS =====================
using Spin = blaze::StaticVector<double, Dim>;
using Distance = blaze::StaticVector<double, Dim>;
using Direction = blaze::StaticVector<double, Dim>;
using CouplingFunction = std::function<double(const Spin&, const Spin&)>;
using Cluster = std::vector<Spin>;
using Ensemble = std::vector<Spin>;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;
using Matrix = blaze::DynamicMatrix<double, blaze::rowMajor>;
using SymmMatrix = blaze::SymmetricMatrix<Matrix>;
using DistanceFunction = std::function<Distance(const Spin&,const Spin&)>;


// compute the coupling matrix for a given container of spins and coupling function
template<typename SpinContainer>
SymmMatrix compute_coupling_matrix( const SpinContainer& list, const CouplingFunction& coupling )
{
    SymmMatrix J(list.size());
    for( size_t i = 0; i < list.size() - 1; i++ )
    {
        for( size_t j = i + 1; j < list.size(); j++ )
        {
            J(i,j) = coupling(list[i],list[j]);
        }
    }
    return J;
}


// ============================================================
// ========== Compute the coupling between two spins ==========
// ============================================================

// ====================== N dimensions ========================
// ...for spins with isotropic dipolar couplings
inline double coupling_iso( const Spin& s1, const Spin& s2 )
{
    auto dist = s2 - s1;
    double r = length(dist);
    return 1.0/std::pow( r, 3 );
}

// ...for spin lattices with nearest neighbor interactions ( assuming lattice constant = 1.0 )
inline double coupling_nearest_neighbor( const Spin& s1, const Spin& s2 )
{
    double r = length( s2 - s1 );
    if( !num::is_equal(r,double{1.0}) )
    {
        return double{0.0};
    }
    else
    {
        return double{1.0};
    }
}

// ...for spin lattices with nearest neighbor interactions ( assuming lattice constant = 1.0 )
inline double coupling_nearest_neighbor_PBC( const Spin& s1, const Spin& s2, const DistanceFunction& PBC_distance )
{
    double r = length( PBC_distance(s2,s1) );
    if( !num::is_equal(r,double{1.0}) )
    {
        return double{0.0};
    }
    else
    {
        return double{1.0};
    }
}


// ====================== 2 dimensions ========================
// ...for surface spins with dipolar couplings in the doubly rotating frame 
inline double coupling_ddrf( const Spin& s1, const Spin& s2 )
{
    auto dist = s2 - s1;
    double r = length(dist);
    return 1.0/std::pow( r, 3 ) * std::cos( 2*std::acos(dist[0]/r) ); // careful: acos(rx/r) is not the actual angle but for cos(2phi) this does not matter
}

// ...for surface spins with dipolar couplings in the doubly rotating frame tilted at an angle phi_0
inline double coupling_ddrf_tilted( const double phi_0, const Spin& s1, const Spin& s2 )
{
    auto dist = s2 - s1;
    double r = length(dist);
    double phi{};
    if( dist[1] > 0 ) // this distinguishing is required because acos(rx/r) is not the actual angle
    {
        phi = std::acos(dist[0]/r);
    }
    else 
    {
        phi = -std::acos(dist[0]/r); // phi = 2pi - acos(...), but the 2pi can be omitted here
    }
    return 1.0/std::pow(r,3) * std::cos(2.0*(phi-phi_0));
}

// ...for spins with dipolar couplings in the doubly rotating frame on a square of length 1
// J = 1/r^3 * cos(2 phi)
inline double periodic_coupling_ddrf( const Spin& s1, const Spin& s2 )
{
    auto direct_dist = s2 - s1;
    Distance minimum_dist{};
    if( length(direct_dist) > 0.5 )
    {
        std::vector<Distance> dist{};
        Distance ex{ 1.0, 0.0 };
        Distance ey{ 0.0, 1.0 };
        dist.emplace_back( direct_dist - ey - ex );
        dist.emplace_back( direct_dist - ey );
        dist.emplace_back( direct_dist - ey + ex );
        dist.emplace_back( direct_dist - ex );
        dist.emplace_back( direct_dist );
        dist.emplace_back( direct_dist + ex );
        dist.emplace_back( direct_dist + ey - ex );
        dist.emplace_back( direct_dist + ey );
        dist.emplace_back( direct_dist + ey + ex );
        minimum_dist = *std::min_element( dist.begin(), dist.end(), []( const Distance& a, const Distance& b )
        {
            return length( a ) < length( b );
        } );
    }
    else{
        minimum_dist = direct_dist;
    }
    double r = length(minimum_dist);
    return 1.0/std::pow( r, 3 ) * std::cos( 2*std::acos( minimum_dist[0] / r ) );
}

// ...for spins with isotropic dipolar coupling on a square of length 1
// J = 1/r^3
inline double periodic_coupling_iso( const Spin& s1, const Spin& s2 )
{
    auto direct_dist = s2 - s1;
    Distance minimum_dist{};
    if( length(direct_dist) > 0.5 )
    {
        std::vector<Distance> dist{};
        Distance ex{ 1.0, 0.0 };
        Distance ey{ 0.0, 1.0 };
        dist.emplace_back( direct_dist - ey - ex );
        dist.emplace_back( direct_dist - ey );
        dist.emplace_back( direct_dist - ey + ex );
        dist.emplace_back( direct_dist - ex );
        dist.emplace_back( direct_dist );
        dist.emplace_back( direct_dist + ex );
        dist.emplace_back( direct_dist + ey - ex );
        dist.emplace_back( direct_dist + ey );
        dist.emplace_back( direct_dist + ey + ex );
        minimum_dist = *std::min_element( dist.begin(), dist.end(), []( const Distance& a, const Distance& b )
        {
            return length( a ) < length( b );
        } );
    }
    else{
        minimum_dist = direct_dist;
    }
    double r = length(minimum_dist);
    return 1.0/std::pow( r, 3 );
}

// periodic boundary conditions on a parallelogram spanned by a1, a2
inline Distance PBC_parallelogram_distance( const Spin& s1, const Spin& s2, const Direction& a1, const Direction& a2 )
{
    auto direct_dist = s2-s1;
    std::vector<Distance> dists{};
    for( int k=-1; k<=1; ++k )
    {
        for( int l=-1; l<=1; ++l )
        {
            dists.emplace_back( direct_dist + static_cast<double>(k)*a1 + static_cast<double>(l)*a2 );
        }
    }
    return *std::min_element( dists.begin(), dists.end(), []( const Distance& a, const Distance& b )
    {
        return length( a ) < length( b );
    } );
}


// ====================== 3 dimensions ========================
// compute the spherical coordinates normal vector
inline blaze::StaticVector<double,3> normal_vector( const double& phi, const double& theta )
{
    return blaze::StaticVector<double,3>{ std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta) };
}

// ...for spins with dipolar couplings in the rotating frame, nB is the direction of the magnetic field
// J = (1-3cos^2(theta)) * 1/r^3 | theta = nB*n(s1,s2)
inline double coupling_drf_3D( const Spin& s1, const Spin& s2, const Direction& nB )
{
    Direction dist = s2 - s1;
    double r = length(dist);
    double cos_theta_sq = std::pow( double{1.}/r * (blaze::trans(dist) * nB), 2 );
    return double{1.}/std::pow(r,3) * ( 1-3*cos_theta_sq );
}


};