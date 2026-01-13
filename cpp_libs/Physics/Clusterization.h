/* Dim has to be defined as constexpr size_t to include this header */
#pragma once
#include<iostream>
#include<vector>
#include<tuple>
#include<functional>
#include<blaze/Math.h>
#include"Standard_Algorithms/Standard_Algorithms.h"
#include"Standard_Algorithms/Numerics.h"


namespace Physics::Clusterization
{

// ===================== NAMESPACES ===========================
namespace stda = Standard_Algorithms;
namespace num = stda::Numerics;

// ===================== USING STATEMENTS =====================
using Spin = blaze::StaticVector<double, Dim>;
using Distance = blaze::StaticVector<double, Dim>;
using Direction = blaze::StaticVector<double, Dim>;
using CouplingFunction = std::function<double(const Spin&, const Spin&)>;
using Cluster = std::vector<Spin>;
using Ensemble = std::vector<Spin>;
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;
using IndexList = std::vector<size_t>;
using Matrix = blaze::DynamicMatrix<double,blaze::rowMajor>;
using SymmMatrix = blaze::SymmetricMatrix<Matrix>;


// =============================================================
// =================== MANUAL CLUSTERIZATION ===================
// =============================================================
// assign the first NGamma spins in the ensemble to the cluster and the rest to the environment
inline std::tuple<Cluster,Cluster,IndexList> auto_clusterize( const Ensemble& ensemble, const size_t NGamma )
{
    // 0.) initialize
    Cluster cluster{};
    Cluster environment{};
    IndexList cluster_indices{};

    // 1.) add the spins
    for( size_t s = 0; s < NGamma; ++s ) // the considered first spins are added to the central cluster
    {
        cluster.emplace_back( ensemble[s] );
        cluster_indices.emplace_back( s );
    }
    for( size_t s = NGamma; s < ensemble.size(); ++s ) // the other spins are added to the rest
    {
        environment.emplace_back( ensemble[s] );
    }

    // 2.) return
    return std::make_tuple( cluster, environment, cluster_indices );
}


// add single spin manually to the cluster by ensemble index
inline void manually_add_spin_to_cluster( Cluster& cluster, Cluster& environment, IndexList& cluster_indices, const Ensemble& ensemble, const size_t ensemble_index )
{
    // remove spin from environment:
    auto new_spin = ensemble[ensemble_index];
    auto it_to_new_spin = std::find_if( environment.cbegin(), environment.cend(), [&new_spin]( const auto& s )
    {
        return num::is_equal( s, new_spin );
    } );
    environment.erase( it_to_new_spin );

    // add spin to cluster:
    cluster.emplace_back( new_spin );

    // add index to cluster indices:
    cluster_indices.emplace_back( ensemble_index ); 
}



// =============================================================
// ================= CSPINBASED CLUSTERIZATION =================
// =============================================================
// central-spin-based clusterization scheme for lattice systems
inline std::tuple<Cluster,Cluster,IndexList> clusterize_ordered_cspinbased( const Ensemble& ensemble, const size_t NGamma, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 0.) initialize
    Cluster cluster{};
    Cluster environment{};
    IndexList cluster_indices{};

    // 1.) add the first num_fixed spins and determine the central spin
    for( size_t s = 0; s < num_fixed; ++s )
    {
        cluster_indices.emplace_back( s );
        cluster.emplace_back( ensemble[s] );
    }
    auto central_spin = cluster[0];


    // 2.) determine the couplings to the central spin
    std::vector<double> quad_coupling_to_central_spin(num_fixed,double{0.}); // set intra-cluster couplings and self-coupling manually to zero
    for( auto it_to_spin=ensemble.cbegin()+num_fixed; it_to_spin!=ensemble.cend(); ++it_to_spin )
    {
        quad_coupling_to_central_spin.emplace_back( pow(coupling(*it_to_spin,central_spin),2) );
    }

    // 3.) add the spins shellwise
    auto grouped_index_lists = stda::argsort_flip_extended(quad_coupling_to_central_spin,[&](double J1, double J2){return num::is_equal(J1,J2);});
    size_t shell = 0;
    while( NGamma!=cluster_indices.size() )
    {
        auto indices = grouped_index_lists[shell++];
        std::sort( indices.begin(), indices.end(), [&ensemble]( size_t index_1, size_t index_2 ) -> bool
        {
            size_t dir = Dim;
            while( --dir >= 0 )
            {
                if( !num::is_equal(ensemble[index_1][dir],ensemble[index_2][dir]) )
                {
                    return ensemble[index_1][dir] > ensemble[index_2][dir]; // compare position components, since they are not equal
                }
            }
            return false;
        } );

        for( size_t i=0; i<indices.size(); ++i )
        {
            cluster_indices.emplace_back( indices[i] );
            cluster.emplace_back( ensemble[indices[i]] );
            if( NGamma==cluster_indices.size() )
            {
                if( i!=indices.size()-1 )
                {
                    std::cout << "\033[1;36mWarning: cluster size to small, I cannot add the whole shell of equally coupled spins. Missing spins: ";
                    for( size_t p=i+1; p<indices.size(); ++p )
                    {
                        std::cout << indices[p] << " ";
                    }
                    std::cout << "\033[0m\n";
                }
                break;
            }
        }
    }
         
    // 4.) add the other spins to the rest
    for( size_t s=0; s<ensemble.size(); ++s )
    {
        if( std::find(cluster_indices.begin(),cluster_indices.end(),s) == cluster_indices.end() ) // the considered spin is not part of the central cluster
        {
            environment.emplace_back( ensemble[s] );
        }
    }

    // 5.) return
    return std::make_tuple(cluster, environment, cluster_indices);   
}


// central-spin-based clusterization scheme for inhomogeneous systems
inline std::tuple<Cluster,Cluster,IndexList> clusterize_inhomogeneous_cspinbased( const Ensemble& ensemble, const size_t NGamma, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 0.) initialize
    Cluster cluster{};
    Cluster environment{};
    IndexList cluster_indices{};

    // 1.) add the first num_fixed spins and determine the central spin
    for( size_t s=0; s<num_fixed; ++s )
    {
        cluster_indices.emplace_back( s );
        cluster.emplace_back( ensemble[s] );
    }
    auto central_spin = cluster[0];

    // 2.) determine the couplings to the central spin
    std::vector<double> quad_coupling_to_central_spin(num_fixed,double{0.}); // set intra-cluster couplings and self-coupling manually to zero
    for( auto it_to_spin=ensemble.cbegin()+num_fixed; it_to_spin!=ensemble.cend(); ++it_to_spin )
    {
        quad_coupling_to_central_spin.emplace_back( pow(coupling(*it_to_spin,central_spin),2) );
    }

    // 3.) add the strongest coupled spins to the cluster (assumption : there are no exactly equal couplings)
    for( size_t s=num_fixed; s<NGamma; ++s )
    {
        auto iterator_to_max = std::max_element(quad_coupling_to_central_spin.begin(), quad_coupling_to_central_spin.end());
        *iterator_to_max = double{0.};
        size_t index_of_strongest_coupled = std::distance( quad_coupling_to_central_spin.begin(), iterator_to_max );
        cluster_indices.emplace_back( index_of_strongest_coupled );
        cluster.emplace_back( ensemble[index_of_strongest_coupled] );
    }
        
    // 4.) add the other spins to the rest
    for( size_t s = 0; s < ensemble.size(); ++s )
    {
        if( std::find( cluster_indices.begin(), cluster_indices.end(), s) == cluster_indices.end() ) // the considered spin is not part of the central cluster
        {
            environment.emplace_back( ensemble[s] );
        }
    }

    // 5.) return
    return std::make_tuple(cluster, environment, cluster_indices);
}



// ====================================================================
// === CLUSTERBASED CLUSTERIZATION (ONLY FOR INHOMOGENEOUS SYSTEMS) ===
// ====================================================================
// cluster-based clusterization scheme
inline std::tuple<Cluster,Cluster,IndexList> clusterize_inhomogeneous_clusterbased( const Ensemble& ensemble, const size_t NGamma, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 0.) initialize
    Cluster cluster{};
    Cluster environment{};
    IndexList cluster_indices{};

    // 1.) add the first num_fixed spins and determine the central spin
    for( size_t s=0; s<num_fixed; ++s )
    {
        cluster_indices.emplace_back( s );
        cluster.emplace_back( ensemble[s] );
    }
    auto central_spin = cluster[0];

    // 2.) add the spin strongest coupled to the cluster one by one
    for( size_t s=num_fixed; s<NGamma ; s++ )
    {
        // determine the spin-cluster couplings (assumption : there are no exactly equal couplings)
        std::vector<double> Delta{}; // needs to be maximized 
        for( size_t k=0; k<ensemble.size(); ++k ) // iterate over the whole ensemble...
        {   
            if( std::find( cluster_indices.begin(), cluster_indices.end(), k ) == cluster_indices.end() ) // ...except the spins that are already part of the cluster
            {
                auto spin_k = ensemble[k];
                double Delta_k = std::accumulate( cluster.cbegin(), cluster.cend(), double{0.}, [coupling,&spin_k]( double sum, const Spin & spin_i )
                {
                    return sum + std::abs( coupling( spin_k, spin_i ) ); // compute the cluster coupling (linear sum of the mod couplings)
                } );
                Delta.emplace_back( Delta_k ); 
            }
            else
            {
                Delta.emplace_back( double{0.} ); 
            }
        }

        // search for maximum coupling and add it to the cluster
        size_t index = std::distance( Delta.begin(), std::max_element(Delta.begin(), Delta.end()) );
        cluster_indices.emplace_back( index );
        cluster.emplace_back( ensemble[index] );
    }
        
    // 3.) add the other spins to the rest
    for( size_t s = 0; s < ensemble.size(); ++s )
    {
        if( std::find( cluster_indices.begin(), cluster_indices.end(), s) == cluster_indices.end() ) // the considered spin is not part of the central cluster
        {
            environment.emplace_back( ensemble[s] );
        }
    }

    // 4.) return
    return std::make_tuple(cluster, environment, cluster_indices);
}


// cluster-based clusterization scheme with variable cluster size
inline std::tuple<Cluster,Cluster,IndexList> clusterize_inhomogeneous_clusterbased_vsize( const Ensemble& ensemble, const std::array<size_t,2> NGamma_limits, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 1.) create clusters for different NGamma
    std::vector<std::tuple<Cluster,Cluster,IndexList>> clusterizations{};
    for( size_t NGamma=NGamma_limits[0]; NGamma<=NGamma_limits[1]; ++NGamma )
    {
        clusterizations.emplace_back( clusterize_inhomogeneous_clusterbased( ensemble, NGamma, coupling, num_fixed ) );
    }

    // 2.) find most optimal clusterization
    std::vector<double> Delta_ls{};
    for( const auto& clusterization: clusterizations )
    {
        // compute all Delta_ks
        auto cluster = std::get<0>(clusterization);
        auto environment = std::get<1>(clusterization);
        std::vector<double> Delta_ks{};
        for( const auto& spin_k: environment )
        {
            Delta_ks.emplace_back( std::accumulate( cluster.cbegin(), cluster.cend(), double{0.}, [&](auto sum, auto spin_i)
            {
                return sum + std::abs(coupling(spin_i,spin_k));
            }) );
        }

        // maximum Delta_k = Delta_l
        Delta_ls.emplace_back( *std::max_element(Delta_ks.begin(), Delta_ks.end()) );
    }
    // best clusterization = minimum Delta_l
    size_t index = std::distance(Delta_ls.begin(), std::min_element(Delta_ls.begin(), Delta_ls.end()));

    // 3.) return
    return clusterizations[index];
}



// ===============================================================
// ================== HYBRID CLUSTERIZATION ======================
// ===============================================================
/* computes the quantity for the criterion in the hybrid clusterization strategy
assumes spin_k is not in cluster and that the central spin is the first spin of the cluster */
inline double Delta_hybrid( const Spin& spin_k, const Cluster& cluster, const CouplingFunction& coupling )
{
    double Delta_cluster_k = std::accumulate( cluster.cbegin(), cluster.cend(), 0., [coupling,&spin_k]( double sum, const Spin & spin_i )
    {
        return sum + std::abs( coupling( spin_k, spin_i ) ); // compute the cluster coupling (linear sum of the mod couplings)
    } );
    double Delta_cspin_k = std::abs(coupling(spin_k,cluster[0]));
    double current_csize = static_cast<double>( cluster.size() );
    double alpha = std::pow(current_csize,-3.0/static_cast<double>(Dim));
    return alpha*Delta_cluster_k + (1-alpha)*Delta_cspin_k; 
}


// hybrid clusterization scheme
inline std::tuple<Cluster,Cluster,IndexList> clusterize_hybrid( const Ensemble& ensemble, const size_t NGamma, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 0.) initialize
    Cluster cluster{}, environment{};
    IndexList cluster_indices{};

    // 1.) add the first num_fixed spins and determine the central spin
    for( size_t s=0; s<num_fixed && s<NGamma; ++s )
    {
        cluster_indices.emplace_back( s );
        cluster.emplace_back( ensemble[s] );
    }
    auto central_spin = cluster[0];

    // 2.) add the spins
    while( cluster_indices.size()<NGamma )
    {
        // a.) determine the spin-cluster couplings
        std::vector<double> Delta{}; // needs to be maximized 
        for( size_t k=0; k<ensemble.size(); ++k ) // iterate over the whole ensemble...
        {   
            if( std::find( cluster_indices.begin(), cluster_indices.end(), k ) == cluster_indices.end() ) // ...except the spins that are already part of the cluster
            {
                Delta.emplace_back( Delta_hybrid(ensemble[k],cluster,coupling) ); 
            }
            else
            {
                Delta.emplace_back( double{0.} ); 
            }
        }

        // b.) find maximally coupled shell of spins (in inhomogeneous systems shell=spin)
        auto indices = stda::argsort_flip_extended(Delta,[&](double D1, double D2){return num::is_equal(D1,D2);})[0]; // list of indices of the spins with maximum Delta_k

        // c.) sort spins of the shell according to some unique scheme
        std::sort(indices.begin(), indices.end(), [&ensemble](size_t index_1, size_t index_2) -> bool
        {
            size_t dir = Dim;
            while( --dir >= 0 )
            {
                if( !num::is_equal(ensemble[index_1][dir],ensemble[index_2][dir]) )
                {
                    return ensemble[index_1][dir] > ensemble[index_2][dir]; // compare position components, since they are not equal
                }
            }
            return false;
        } );

        // d.) add the spins of the shell
        for( size_t i=0; i<indices.size(); ++i )
        {
            cluster_indices.emplace_back( indices[i] );
            cluster.emplace_back( ensemble[indices[i]] );
            if( NGamma==cluster_indices.size() )
            {
                if( i!=indices.size()-1 )
                {
                    std::cout << "\033[1;36mWarning: cluster size to small, I cannot add the whole shell of equally coupled spins. Missing spins: ";
                    for( size_t p=i+1; p<indices.size(); ++p )
                    {
                        std::cout << indices[p] << " ";
                    }
                    std::cout << "\033[0m\n";
                }
                break;
            }
        }
    }
        
    // 3.) add the other spins to the rest
    for( size_t s = 0; s < ensemble.size(); ++s )
    {
        if(std::find( cluster_indices.begin(), cluster_indices.end(), s) == cluster_indices.end()) // the considered spin is not part of the central cluster
        {
            environment.emplace_back( ensemble[s] );
        }
    }

    // 4.) return
    return std::make_tuple(cluster, environment, cluster_indices);
}


/* hybrid clusterization scheme with variable cluster size 
mainly used for inhomogeneous systems, not very useful for lattices because many cluster sizes are implausible there */
std::tuple<Cluster,Cluster,IndexList>  clusterize_hybrid_vsize( const Ensemble& ensemble, const size_t NGamma_max, const CouplingFunction& coupling, const size_t num_fixed = 1 )
{
    // 1.) create clusters for different NGamma
    std::vector<std::tuple<Cluster,Cluster,IndexList>> clusterizations{};
    for( size_t NGamma=num_fixed; NGamma<=NGamma_max; ++NGamma )
    {
        clusterizations.emplace_back( clusterize_hybrid( ensemble, NGamma, coupling, num_fixed ) );
    }

    // 2.) find most optimal clusterization
    std::vector<double> rho_ns{};
    for( const auto& clusterization: clusterizations )
    {
        auto cluster = std::get<0>(clusterization);
        auto environment = std::get<1>(clusterization);
        std::vector<double> Delta_hybrid_ks{};
        for( const auto& spin_k: environment ) // compute all Delta_hybrid for all k
        {
            Delta_hybrid_ks.emplace_back( Delta_hybrid(spin_k,cluster,coupling) );
        }
        rho_ns.emplace_back( *std::max_element(Delta_hybrid_ks.begin(),Delta_hybrid_ks.end()) ); // maximum Delta_hybrid_k = rho_n
    }
    // best clusterization = minimum rho_n
    size_t index = std::distance( rho_ns.begin(), std::min_element(rho_ns.begin(),rho_ns.end()) );

    // 3.) return
    return clusterizations[index];
}


}