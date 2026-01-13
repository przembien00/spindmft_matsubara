/* Dim has to be defined as constexpr size_t to include this header */
#pragma once
#include<iostream>
#include<vector>
#include<functional>
#include<blaze/Math.h>
#include"Standard_Algorithms/Numerics.h"

namespace Physics::Mean_Field_Couplings
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
using Matrix = blaze::DynamicMatrix<double,blaze::rowMajor>;
using SymmMatrix = blaze::SymmetricMatrix<Matrix>;
using RVector = blaze::DynamicVector<double, blaze::columnVector>;
using MatrixOfVector = blaze::DynamicMatrix<RVector, blaze::rowMajor>;
using SymmMatrixOfVector = blaze::SymmetricMatrix<MatrixOfVector>;
using MatrixOfMatrix = blaze::DynamicMatrix<Matrix,blaze::rowMajor>;
using SymmMatrixOfMatrix = blaze::SymmetricMatrix<MatrixOfMatrix>;


// ============================================================
// ======= Coupling constants and coordination numbers ========
// ============================================================
/* compute the quadratic coupling constant JQ^2 of spin with index i */
inline double JQ_sq( const size_t i, const Ensemble& ensemble, const CouplingFunction J )
{
    size_t j{};
    auto spin_i = ensemble[i];
    return std::accumulate( ensemble.begin(), ensemble.end(), double{0.}, [&]( double sum, const Spin& spin_j )
        {
            return sum + ((i!=j++) ? std::pow(J(spin_i,spin_j),2) : double{0.});
        } );
}

/* compute the nth-order coupling constant Jn^n = sum_j J_ij^n of spin with index i */
inline double J_nth( const size_t n, const size_t i, const Ensemble& ensemble, const CouplingFunction J )
{
    size_t j{};
    auto spin_i = ensemble[i];
    return std::accumulate( ensemble.begin(), ensemble.end(), double{0.}, [&]( double sum, const Spin& spin_j )
        {
            return sum + ((i!=j++) ? std::pow(J(spin_i,spin_j),n) : double{0.});
        } );
}

/* compute the cross coupling constant JQ_ij^2 of two spins with indices i and j */
inline double JQ_cross_sq( const size_t i, const size_t j, const Ensemble& ensemble, const CouplingFunction J )
{
    size_t k{};
    Spin spin_i = ensemble[i], spin_j = ensemble[j];
    return std::accumulate( ensemble.begin(), ensemble.end(), double{0.}, [&]( double sum, const Spin& spin_k )
        {
            bool condition = (i!=k && j!=k);
            k += 1;
            return sum + (condition ? J(spin_i,spin_k)*J(spin_j,spin_k) : double{0.});
        } );
}

/* compute the effective coordination number z_eff = JQ^4 / JT^4 of spin with index i */
inline double zeff( const size_t i, const Ensemble& ensemble, const CouplingFunction J )
{
    return std::pow(JQ_sq(i,ensemble,J),2) / J_nth(4,i,ensemble,J);
}


// ============================================================
// ==================== Coupling tensors ======================
// ============================================================
/* extract a row from a blaze matrix */
inline std::vector<double> extract_row( const Matrix& C, const size_t row_index )
{
    std::vector<double> row{};
    for( uint col_index=0; col_index<C.columns(); ++col_index )
    {
        row.emplace_back(C(row_index,col_index));
    }
    return row;
}

/* Compute the mean-field-coupling tensor in correlation replica approximation (CRA) for an inhomogeneous system (required for CspinDMFT) 
inhomogeneous means that no pair of distances or couplings matches */
inline SymmMatrixOfMatrix compute_JCRA_sq( const Cluster& Gamma, const Cluster& environment, const CouplingFunction coupling )
{    
    // create ensemble (cluster + environment)
    Ensemble ensemble = Gamma;
    ensemble.insert( ensemble.end(), environment.begin(), environment.end() );

    // treat cluster sizes of 1 
    if( Gamma.size() == 1 )
    {
        SymmMatrixOfMatrix JCRA_sq(1);
        JCRA_sq(0,0) = Matrix{ {JQ_sq(0,ensemble,coupling)} };
        return JCRA_sq;
    }

    // 1.) compute S_k and C_kl for Gamma (k<=l)
    std::vector<double> S_Gamma(Gamma.size()); // autocorrelation quantity
    Matrix C_Gamma = blaze::ZeroMatrix<double,blaze::rowMajor>(Gamma.size(),Gamma.size()); // pair-correlation quantity
    std::vector<double> C_Gamma_linearized{}; // linearized version of C_Gamma (allows for more efficient iteration)
    for( uint k = 0; k < Gamma.size(); ++k )
    {
        for( uint l = k; l < Gamma.size(); ++l )
        {
            if( k==l ) // autocorrelations
            {
                S_Gamma[k] = JQ_sq(k,ensemble,coupling);
            }   
            else // pair-correlations 
            {
                C_Gamma(k,l) = std::pow(coupling(Gamma[k],Gamma[l]),2);
                C_Gamma_linearized.emplace_back(C_Gamma(k,l));
            }
        }
    }

    // 2.) compute the cutoff
    double C_cutoff = double{0.5} * (*std::min_element(C_Gamma_linearized.cbegin(),C_Gamma_linearized.cend()));

    // 3.) compute S_p and C_pq for the environment (p<=q)
    std::vector<double> S_env(environment.size()); // autocorrelation quantity
    Matrix C_env = blaze::ZeroMatrix<double,blaze::rowMajor>(environment.size(),environment.size()); // pair-correlation quantity
    for( uint p = 0; p < environment.size(); ++p )
    {
        for( uint q = p; q < environment.size(); ++q )
        {
            if( p==q ) // autocorrelations
            {
                S_env[p] = JQ_sq(p+Gamma.size(),ensemble,coupling);
            }   
            else // pair-correlations 
            {
                C_env(p,q) = std::pow(coupling(environment[p],environment[q]),2);
            }
        }
    }

    // 4.) compute JCRA_sq
    SymmMatrixOfMatrix JCRA_sq = blaze::declsym( MatrixOfMatrix(Gamma.size(), Gamma.size(), blaze::ZeroMatrix<double,blaze::rowMajor>(Gamma.size(),Gamma.size())) );
    for( uint p = 0; p < environment.size(); ++p )
    {
        for( uint q = p; q < environment.size(); ++q ) // iterate only over half of the pair correlations due to symmetry (factor 2 later)
        {
            // a.) no contribution from pair-correlations below the cutoff           
            if( p!=q ) // pair-correlations
            {
                if( C_env(p,q) < C_cutoff )
                {
                    continue;
                }
            }

            // b.) determine {kl}=f({pq}) from CRA
            size_t k, l;     
            if( p==q ) // autocorrelations 
            {
                // (i) determine the differences
                std::vector<double> sqdiff{};
                auto S_env_p = S_env[p];
                std::transform( S_Gamma.cbegin(), S_Gamma.cend(), std::back_inserter( sqdiff ), [&S_env_p]( const double& S_Gamma_k )
                {
                    return std::pow(S_Gamma_k-S_env_p,2);
                } );
                // (ii) locate the index (pair) with the minimum difference
                k = std::distance( sqdiff.begin(), std::min_element( sqdiff.begin(), sqdiff.end() ) );
                l = k;
            }
            else // pair-correlations 
            {
                // (i) determine the sqared differences
                double C_env_pq = C_env(p,q);
                auto sqdiff = blaze::map(C_Gamma,[&C_env_pq](double C_Gamma_kl){return std::pow(C_Gamma_kl-C_env_pq,2);}); // due to the cutoff the minimum is never found at the zero elements of C_env
                // (ii) locate the index pair {kl} with the minimum difference
                blaze::DynamicVector<double,blaze::columnVector> rowmins = blaze::min<blaze::rowwise>(sqdiff);
                #ifdef WORKSTATION
                k = std::argmin(rowmins);
                #else
                k = blaze::argmin(rowmins);
                #endif
                auto row = extract_row(sqdiff,k);
                l = std::distance( row.cbegin(), std::min_element(row.cbegin(),row.cend()) );
            }

            // c.) for all {ij} add the result to JCRA_sq
            for( uint i = 0; i < Gamma.size(); ++i )
            {
                for( uint j = 0; j < Gamma.size(); ++j )
                {
                    double Jsq = coupling(Gamma[i],environment[p]) * coupling(Gamma[j],environment[q]);
                    JCRA_sq(i,j)(k,l) += Jsq * ((p==q) ? double{1.} : double{2.}); // factor 2 for pair corr's
                }
            }
        }
    }

    // 5.) return
    return JCRA_sq;
}



/* Definition of Category class required for "compute_JCR_sq" below 
Categorization of an in-cluster correlation {kl} through the computation of absolute distance r_kl and coupling J_kl_sq 
This should not be used for lattices with more-atomic bases
PBC not usable here yet */
struct Category
{
    // members
    size_t k{};
    size_t l{};
    double r_kl{};
    double J_kl_sq{};
    bool is_auto{}; // true if k==l

    // constructor
    Category( size_t k_, size_t l_, const Cluster& cluster, const CouplingFunction& J ):
        k(k_),
        l(l_)
    {
        if( k==l )
        {
            is_auto = true;
        }
        else
        {
            is_auto = false;
            r_kl = blaze::norm(cluster[k]-cluster[l]);
            J_kl_sq = std::pow(J(cluster[k],cluster[l]),2);
        }
    }

    // comparison operators
    bool operator==(const Category& other) const
    {
        if( (is_auto && !other.is_auto) || (!is_auto && other.is_auto) ) // compairing auto with pair
        {
            return false;
        }
        return num::is_equal(r_kl,other.r_kl,double{1e-6}) && num::is_equal(J_kl_sq,other.J_kl_sq,double{1e-6});
    }
    bool operator!=(const Category& other) const
    {
        return !(*this == other);
    }

    // printing routine
    void print() const 
    {
        std::cout << "Category {k,l} = {" << k << "," << l << "} :\n";
        std::cout << "r_kl   = " << r_kl << "\n";
        std::cout << "J_kl^2 = " << J_kl_sq << "\n";
    }   
};

/* Compute the mean-field-coupling tensor using correlation replicas (CR) in case of a lattice system (required for CspinDMFT)
The algorithm categorizes and finds correlation replicas comparing the absolute distances and couplings between spins 
This should not be used for lattices with more-atomic bases
PBC not usable here yet */
inline std::tuple<SymmMatrixOfMatrix,IndexPairList> compute_JCR_sq( const Cluster& Gamma, const Cluster& environment, const CouplingFunction coupling )
{    
    // treat cluster sizes of 1 
    if( Gamma.size() == 1 )
    {
        Ensemble ensemble = Gamma;
        ensemble.insert( ensemble.end(), environment.begin(), environment.end() );
        SymmMatrixOfMatrix JCR_sq(1);
        JCR_sq(0,0) = Matrix{ {JQ_sq(0,ensemble,coupling)} };
        return std::make_tuple(JCR_sq, IndexPairList{IndexPair{0,0}});
    }

    // 1.) categorize the in-cluster correlations (k<=l)
    std::vector<Category> categories{};
    for( uint k = 0; k < Gamma.size(); ++k )
    {
        for( uint l = k; l < Gamma.size(); ++l )
        {
            Category new_cat(k,l,Gamma,coupling);
            if( std::find(categories.cbegin(),categories.cend(),new_cat) == categories.cend() ) // then the new category does not exist yet
            {
                categories.emplace_back(new_cat);
            }
        }
    }

    // 2.) compute JCR_sq
    SymmMatrixOfMatrix JCR_sq = blaze::declsym( MatrixOfMatrix(Gamma.size(), Gamma.size(), blaze::ZeroMatrix<double,blaze::rowMajor>(Gamma.size(),Gamma.size())) );
    for( uint p = 0; p < environment.size(); ++p )
    {
        for( uint q = p; q < environment.size(); ++q ) // iterate only over half of the pair correlations due to symmetry (factor 2 later)
        {
            // a.) categorize the correlation
            Category new_cat(p,q,environment,coupling);
            auto it_to_incluster_category = std::find(categories.cbegin(),categories.cend(),new_cat);
            if( it_to_incluster_category == categories.cend() ) // then the category does not exist
            {
                continue;
            }

            // b.) find {kl}
            size_t k = it_to_incluster_category->k;
            size_t l = it_to_incluster_category->l;

            // c.) for all {ij} add the result to JCR_sq
            for( uint i = 0; i < Gamma.size(); ++i )
            {
                for( uint j = 0; j < Gamma.size(); ++j )
                {
                    double Jsq = coupling(Gamma[i],environment[p]) * coupling(Gamma[j],environment[q]);
                    JCR_sq(i,j)(k,l) += Jsq * ((p==q) ? double{1.} : double{2.}); // factor 2 for pair corr's
                }
            }
        }
    }

    // 3.) extract double indices from categories
    IndexPairList indices{};
    for( const auto& cat: categories )
    {
        indices.emplace_back( IndexPair{cat.k,cat.l} );
    }

    // 4.) return
    return std::make_tuple(JCR_sq, indices);
}



/* Compute the mean-field-coupling tensor which couples the cluster to the bath (required for nl-spinDMFT)
assumption: only autocorrelations contribute to the mean-field dynamics and there is only a single sort of autocorrelatipns */
inline SymmMatrixOfVector compute_Jbath_sq_simple( const Cluster& Gamma, const Cluster& environment, const CouplingFunction& coupling )
{
    SymmMatrixOfVector Jbath_sq( Gamma.size());
    for( size_t i = 0; i < Gamma.size(); ++i ) // outer (symmetric) matrix : rows
    {
        for( size_t j = i; j < Gamma.size(); ++j ) // outer (symmetric) matrix : cols
        {
            Jbath_sq(i,j).resize( 1, double{0.} ); // size 1 because only the autocorrelation contributes
            for( size_t k = 0; k < environment.size(); ++k )
            {
                double J_sq = coupling(Gamma[i],environment[k]) * coupling(Gamma[j],environment[k]); // J_ik * J_jk
                Jbath_sq(i,j)[0] += J_sq; // Jbath^2 = sum_k  J_ik * J_jk
            }
        }
    }
    return Jbath_sq;
}


// todo: perhaps implement compute_JCRav_sq to computes the mf-coupling tensor in CRav for lattice systems


};