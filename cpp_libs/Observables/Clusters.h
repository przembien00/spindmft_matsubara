#pragma once

#include<vector>
#include<string>
#include<functional>
#include"Globals/Types.h"
#include"Standard_Algorithms/Standard_Algorithms.h"
#include"O_Error_Handling.h"

namespace Observables::Clusters
{

// ===================== NAMESPACES ==========================
namespace error = Observables::Error_Handling;
namespace stda = Standard_Algorithms;


// ===================== USING STATEMENTS ====================
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;


// ===========================================================
// =============== CORRELATION CLUSTER HEADER ================
// ===========================================================
/* contains a cluster of correlation entities Gij=Gji which are stored linearly in a single vector 
either, all possible correlations for a given number of spins are stored according to
G11 G12 ... G1N G22 ... G2N ... GNN
or, a specific subset of correlations is stored 
m_site_pairs stores the corresponding index pairs ij (always with i<j, but all functions accept also index pairs i>j)
*/
template<typename Correlation>
class CorrelationCluster
{
 public:
    // CONSTRUCTORS
    CorrelationCluster() = default;
    explicit CorrelationCluster( const size_t num_Spins );
    template<typename ParameterPack>
    CorrelationCluster( const size_t num_Spins, const ParameterPack& params );
    template<typename Parameter1, typename Parameter2>
    CorrelationCluster( const size_t num_Spins, const Parameter1& p1, const Parameter2& p2 );
    explicit CorrelationCluster( const IndexPairList& list );
    template<typename ParameterPack>
    CorrelationCluster( const IndexPairList& list, const ParameterPack& params );
    template<typename Parameter1, typename Parameter2>
    CorrelationCluster( const IndexPairList& list, const Parameter1& p1, const Parameter2& p2 );

    // OPERATORS
    Correlation& operator[]( const size_t linear_index ){ return m_cluster[linear_index]; }
    const Correlation& operator[]( const size_t linear_index ) const { return m_cluster[linear_index]; }
    Correlation& operator()( const size_t spin_i, const size_t spin_j ) { return m_cluster[position(IndexPair{spin_i,spin_j})]; }
    const Correlation& operator()( const size_t spin_i, const size_t spin_j ) const { return m_cluster[position(IndexPair{spin_i,spin_j})]; }
    
    // ITERATORS
    auto begin() { return m_cluster.begin(); }
    auto cbegin() const { return m_cluster.cbegin(); }
    auto end() { return m_cluster.end(); }
    auto cend() const { return m_cluster.cend(); }

    // PUBLIC METHODS
    template<typename BinaryFunction>
    void iterate( BinaryFunction f ); 
    template<typename BinaryFunction>
    void const_iterate( BinaryFunction f ) const; 
    template<typename TrinaryFunction>
    void iterate2( CorrelationCluster<Correlation>& other, TrinaryFunction f ); 
    template<typename TrinaryFunction>
    void const_iterate2( const CorrelationCluster<Correlation>& other, TrinaryFunction f ) const; 

    // GET FUNCTIONS
    size_t size() const { return m_size; }
    size_t get_linear_index( const size_t spin_i, const size_t spin_j ) const { return position(IndexPair{spin_i,spin_j}); }
    IndexPair get_site_pair( const size_t linear_index ) const { return m_site_pairs[linear_index]; }
    const IndexPairList& get_site_pairs() const { return m_site_pairs; }
    bool is_diagonal( size_t linear_index ) const { auto ip = m_site_pairs[linear_index]; return ip[0]==ip[1]; }

    void print( size_t my_rank = 0 ) const;

 private:
    // PRIVATE MEMBERS
    std::vector<Correlation> m_cluster{};
    IndexPairList m_site_pairs{};
    size_t m_size{};

    // PRIVATE METHODS
    void initialize_site_pairs_from_num_Spins( const size_t num_Spins );
    IndexPairList reduce_and_order( const IndexPairList& list ) const;
    IndexPair ascending_order( const IndexPair& ip ) const { return (ip[0]<ip[1]) ? ip : IndexPair{ip[1],ip[0]}; } // swap if i>j
    size_t position( const IndexPair& ip ) const { return std::distance( m_site_pairs.cbegin(), std::find(m_site_pairs.cbegin(),m_site_pairs.cend(),ascending_order(ip)) );}
};


// ===========================================================
// ===== IMPLEMENTATION OF CORRELATION CLUSTER CLASS =========
// ===========================================================
// initialize CorrelationCluster from number of spins, leave sub-elements uninitialized
template<typename Correlation>
CorrelationCluster<Correlation>::CorrelationCluster( const size_t num_Spins ):
    m_size((num_Spins*(num_Spins+1)) / 2)
{
    m_cluster.resize( m_size );
    initialize_site_pairs_from_num_Spins( num_Spins );
}

// initialize CorrelationCluster from number of spins, initialize correlations from ParameterPack
template<typename Correlation>
template<typename ParameterPack>
CorrelationCluster<Correlation>::CorrelationCluster( const size_t num_Spins, const ParameterPack& params ):
    m_size((num_Spins*(num_Spins+1)) / 2)
{
    m_cluster.resize( m_size, Correlation{params} );

}

// initialize CorrelationCluster from number of spins, initialize correlations from two template parameters
template<typename Correlation>
template<typename Parameter1, typename Parameter2>
CorrelationCluster<Correlation>::CorrelationCluster( const size_t num_Spins, const Parameter1& p1, const Parameter2& p2 ):
    m_size((num_Spins*(num_Spins+1)) / 2)
{
    m_cluster.resize( m_size, Correlation{p1,p2} );
    initialize_site_pairs_from_num_Spins( num_Spins );
}

// initialize CorrelationCluster from index list, leave sub-elements uninitialized
template<typename Correlation>
CorrelationCluster<Correlation>::CorrelationCluster( const IndexPairList& list):
    m_site_pairs(reduce_and_order(list)),
    m_size(list.size())
{
    m_cluster.resize( m_size );
}

// initialize CorrelationCluster from index list, initialize correlations from ParameterPack
template<typename Correlation>
template<typename ParameterPack>
CorrelationCluster<Correlation>::CorrelationCluster( const IndexPairList& list, const ParameterPack& params ):
    m_site_pairs(reduce_and_order(list)),
    m_size(list.size())
{
    m_cluster.resize( m_size, Correlation{params} );
}

// initialize CorrelationCluster from index list, initialize correlations from two template parameters
template<typename Correlation>
template<typename Parameter1, typename Parameter2>
CorrelationCluster<Correlation>::CorrelationCluster( const IndexPairList& list, const Parameter1& p1, const Parameter2& p2 ):
    m_site_pairs(reduce_and_order(list)),
    m_size(list.size())
{
    m_cluster.resize( m_size, Correlation{p1,p2} );
}

// iterate over m_cluster and m_site_pairs simultaneously using a binary function
template<typename Correlation>
template<typename BinaryFunction>
void CorrelationCluster<Correlation>::iterate( BinaryFunction f )
{
    stda::for_2each( m_cluster.begin(), m_cluster.end(), m_site_pairs.cbegin(), f );
}

// const iterate over m_cluster and m_site_pairs simultaneously using a binary function
template<typename Correlation>
template<typename BinaryFunction>
void CorrelationCluster<Correlation>::const_iterate( BinaryFunction f ) const
{
    stda::for_2each( m_cluster.cbegin(), m_cluster.cend(), m_site_pairs.cbegin(), f );
}

// iterate over m_cluster, m_cluster of another CorrelationCluster and m_site_pairs simultaneously using a trinary function
// similarity of this and the other CorrelationCluster is not checked (user responsibility)
template<typename Correlation>
template<typename TrinaryFunction>
void CorrelationCluster<Correlation>::iterate2( CorrelationCluster<Correlation>& other, TrinaryFunction f )
{
    stda::for_3each( m_cluster.begin(), m_cluster.end(), other.begin(), m_site_pairs.cbegin(), f );
}

// const iterate over m_cluster, m_cluster of another CorrelationCluster and m_site_pairs simultaneously using a trinary function
// similarity of this and the other CorrelationCluster is not checked (user responsibility)
template<typename Correlation>
template<typename TrinaryFunction>
void CorrelationCluster<Correlation>::const_iterate2( const CorrelationCluster<Correlation>& other, TrinaryFunction f ) const
{
    stda::for_3each( m_cluster.cbegin(), m_cluster.cend(), other.cbegin(), m_site_pairs.cbegin(), f );
}

template<typename Correlation>
void CorrelationCluster<Correlation>::print( size_t my_rank ) const
{
    if( my_rank == 0 )
    {
        std::string object_name = "=== Correlation Cluster ===\n";
        std::cout << object_name;
        const_iterate( [my_rank, i=0]( const auto& c, const auto& ip ) mutable
        {
            std::cout << std::to_string(i++) << ": index pair " << std::to_string(ip[0]) << "-" << std::to_string(ip[1]) << "\n";
            c.print();
            std::cout << "\n";
        } );
        std::cout << std::string(object_name.size()-1,'=') << "\n\n";
    }
}

template<typename Correlation>
void CorrelationCluster<Correlation>::initialize_site_pairs_from_num_Spins( const size_t num_Spins )
{
    for( size_t spin_i = 0; spin_i < num_Spins; spin_i++ )
    {
        for( size_t spin_j = spin_i; spin_j < num_Spins; spin_j++ )
        {
            m_site_pairs.emplace_back(IndexPair{spin_i,spin_j}); // spin_i < spin_j ensured
        }
    }
}

template<typename Correlation>
IndexPairList CorrelationCluster<Correlation>::reduce_and_order( const IndexPairList& list ) const
{
    IndexPairList new_list{};
    std::transform( list.cbegin(), list.cend(), std::back_inserter(new_list), [&]( const auto& ip ){ return ascending_order(ip); } );
    std::sort( new_list.begin(), new_list.end(), []( const auto& ip1, const auto& ip2 )
    { 
        return ip1[0] < ip2[0] || (ip1[0] == ip2[0] && ip1[1] < ip2[1]);
    } );
    new_list.erase(std::unique(new_list.begin(), new_list.end()), new_list.end());
    return new_list;
}

};