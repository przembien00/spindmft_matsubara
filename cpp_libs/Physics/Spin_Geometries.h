/* Dim has to be defined as constexpr size_t to include this header */
#pragma once
#include<iostream>
#include<vector>
#include<functional>
#include<blaze/Math.h>
#include <random>
#include"Standard_Algorithms/Standard_Algorithms.h"
#include"Standard_Algorithms/Numerics.h"

namespace Physics::Spin_Geometries
{

// =================== NAMESPACES ====================================
namespace stda = Standard_Algorithms;
namespace num = stda::Numerics;

// =================== USING STATEMENTS ==============================
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
using Constraint = std::function<bool(const Spin&)>;
using Basis_Vector = Distance;
using Group = std::vector<Spin>;
using Unit_Cell = std::vector<Group>;


// ===================================================================
// ================= INHOMOGENEOUS ENSEMBLE CLASS ====================
// ===================================================================
/* This class generates a Dim-dimensional cube of randomly positioned spins
- the cube center is at (0,0,...)^T
- the ensemble can be constructed completely random or with some constraints, e.g., a set of fixed spins and or a general constraint such as a minimum relative distance between a pair of spins
- the final ensemble contains a list of a fixed number of spins (for example the central spin and some of its neighbors) followed by the remaining spins */
class InhomogeneousEnsemble 
{
 private:
    Ensemble m_ensemble{};
    size_t seed{};
    size_t num_Spins{};
    double L{};

 public:
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++++++++++ completely random constructor +++++++++++
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
    InhomogeneousEnsemble( const size_t seed_, const size_t num_Spins_, const double density = 1.0 ):
        seed(seed_),
        num_Spins(num_Spins_),
        L(std::pow(static_cast<double>(num_Spins)/density,1.0/static_cast<double>(Dim)))
    {
        // initialize random generator:
        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> uni_dist(0.0, L);

        // draw spins:
        m_ensemble.resize(num_Spins);
        std::vector<double> distances{};
        for(auto& spin: m_ensemble)
        {
            for(auto& c: spin)
            {
                c = uni_dist(generator)-L/2; // draw component randomly
            }
            distances.emplace_back( std::accumulate(spin.cbegin(),spin.cend(),double{0.0}, 
                [](double sum, const double& ri)
            {
                return sum + std::pow(ri,2.0);
            }) );
        }

        // put central spin first:
        auto index_of_central_spin = std::distance(distances.cbegin(),std::min_element(distances.cbegin(),distances.cend()));
        auto it_to_central_spin = m_ensemble.begin() + index_of_central_spin;
        auto central_spin = *it_to_central_spin;
        m_ensemble.erase( it_to_central_spin );
        m_ensemble.insert( m_ensemble.begin(), central_spin );
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++ random constructor with some spins fixed and optional constaint +++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    InhomogeneousEnsemble( const size_t seed_, const size_t num_Spins_, const Cluster& fixed_spins, const double r_min = 0.0, const Constraint& special_constraint = [](auto){return true;}, const double density = 1.0 ):
        seed(seed_),
        num_Spins(num_Spins_),
        L(std::pow(static_cast<double>(num_Spins)/density,1.0/static_cast<double>(Dim)))
    {
        // initialize random generator:
        std::mt19937 generator(seed);
        std::uniform_real_distribution<double> uni_dist(0.0, L);

        // create minimum distance constraint:
        Constraint min_dist_constraint = [](auto){return true;};
        if( !num::is_equal(r_min,double{0.0}) )
        {
            min_dist_constraint = [this,r_min]( const Spin& new_spin )
            {
                for(const auto& spin: this->m_ensemble)
                {
                    if( length(spin-new_spin)<r_min )
                    {
                        return false;
                    }
                }
                return true;
            };
        }

        // add fixed spins:
        m_ensemble.resize(num_Spins);
        auto it_to_ensemble = m_ensemble.begin();
        for(const auto& spin: fixed_spins)
        {
            *it_to_ensemble++ = spin;
        }
        
        // draw remaining spins:
        for(; it_to_ensemble<m_ensemble.end(); ++it_to_ensemble)
        {
            Spin new_spin{};
            bool valid = false;
            while( !valid ) // draw the spin many times until it fulfills the constraints
            {
                for(auto& c: new_spin)
                {
                    c = uni_dist(generator)-L/2; // draw component randomly
                }
                if( min_dist_constraint(new_spin) && special_constraint(new_spin) )
                {
                    valid = true;
                }
            }
            *it_to_ensemble = new_spin;
        }
    }

    Ensemble get_ensemble() const
    {
        return m_ensemble;
    }

    void print() const
    {
        std::cout << Dim << "D inhomogeneous spin ensemble (" << num_Spins << " spins)\n";
        for( const auto& s : m_ensemble )
        {
            std::cout << s << "\n"; 
        }
    }
};


// ===================================================================
// =========== HELPER FUNCTIONS FOR LATTICE CLASS BELOW ==============
// ===================================================================
// check if two blaze vectors are equal
template<typename BlazeVector>
bool is_equal( const BlazeVector& v1, const BlazeVector& v2 )
{
    return num::is_zero_soft( blaze::norm( v1 - v2 ) );
}

// check if a spin is in a group
inline bool is_in( const Spin& s, const Group& g )
{
    return std::any_of( g.cbegin(), g.cend(), [s](const auto& s_p){return is_equal(s, s_p);} );
}

// print a group
inline void print_group( const Group& g )
{
    for( const auto& s : g )
    {
        std::cout << "(";
        for( size_t dir=0; dir<Dim; ++dir )
        {
            std::cout << s[dir] << ", ";
        }  
        std::cout << ")\n";
    }
}

// print a cell
inline void print_cell( const Unit_Cell& uc )
{
    size_t g_number{1};
    for( const auto& g : uc )
    {
        std::cout << "Group No. " << g_number++ << ":\n";
        print_group( g );
        std::cout << "\n";
    }
}

// operator for shifting a whole cell by some vector
inline Unit_Cell operator+( const Unit_Cell& uc, const Basis_Vector& shift )
{
    auto new_cell = uc;
    for( auto& g : new_cell )
    {
        for( auto& s : g )
        {
            s += shift;
        }
    }
    return new_cell;
}

// function that updates indices, required in Lattice class constructor
inline void update_indices( std::vector<int>& index_tuple, const int direction, const int max_index )
{
    if(direction<0){ throw std::runtime_error("dimension index exceeds limits."); }
    if(index_tuple[direction]!=max_index)
    { 
        ++index_tuple[direction]; // increase and terminate
    }
    else{ 
        index_tuple[direction]=-max_index; // reset periodically
        update_indices(index_tuple,direction-1,max_index); // update index of higher dimension
    } 
};

// ===================================================================
// ======================= LATTICE CLASS =============================
// ===================================================================
/* this class generates a lattice in Dim dimensions:
1) several groups of spins are build (they can have different sizes)
2) the groups are combined to a unit cell 
3) the unit cell and basis vectors are inserted into the Lattice class, which generates the lattice up to a desired number of cells
The final lattice contains the spin positions ordered into cells that are ordered into groups. Having two grouping objects, namely groups and unit cells, has some advantages for the clusterizations. A group (or later cluster) does not necessarily form a cell, by which a total lattice can be generated. A unit cell of combined groups, however, is able to do so. The cell and group structure is removed in the returned ensemble (central group, central spin first) */
class Lattice 
{
 private:
    std::vector<Unit_Cell> m_lattice{}; // vector of unit cells = vector of vector of groups
    size_t num_CellsPerDir{};
    size_t num_Cells{};
    size_t num_GroupsPerCell{};
    size_t num_SpinsPerCell{};
    size_t num_Spins{};
    size_t index_of_central_cell{};

 public:
    Lattice( const Unit_Cell& uc, const int num_CellsInDir, const std::array<Basis_Vector,Dim>& a ):
        num_CellsPerDir( 2*num_CellsInDir+1 ),
        num_Cells( std::pow(num_CellsPerDir,Dim) ),
        num_GroupsPerCell( uc.size() ),
        num_SpinsPerCell( std::accumulate(uc.cbegin(),uc.cend(),size_t{0},[](size_t sum, const Group& g){return sum+g.size();}) ),
        num_Spins( num_SpinsPerCell * std::pow(num_CellsPerDir,Dim) )
    {
        // construct lattice:
        std::vector<int> current_index_tuple(Dim,-num_CellsInDir);
        std::vector<int> last_index_tuple(Dim,num_CellsInDir);
        while(true)
        {
            // add cell
            Basis_Vector shift{};
            for(size_t dir=0; dir<Dim; ++dir)
            {
                shift += current_index_tuple[dir]*a[dir];
            }
            add_cell( uc + shift );
           
            // recursively update index tuple
            if( current_index_tuple == last_index_tuple )
            {
                break;
            }
            else 
            {
                update_indices(current_index_tuple,Dim-1,num_CellsInDir); 
            }
        }

        // central cell:
        index_of_central_cell = static_cast<size_t>((float)num_Cells/2.);
    }

    void add_cell( const Unit_Cell&& cell )
    {
        m_lattice.emplace_back( cell );
    }

    Ensemble get_ensemble() const
    {
        // first add the central cell
        Ensemble ensemble{};
        for(const auto& group: m_lattice[index_of_central_cell])
        {
            for(const auto& spin: group)
            {
                ensemble.emplace_back(spin);
            }
        }

        // then the other cells
        size_t cell_index{};
        for(const auto& cell: m_lattice)
        {
            if( cell_index++ != index_of_central_cell )
            {
                for(const auto& group: cell)
                {
                    for(const auto& spin: group)
                    {
                        ensemble.emplace_back(spin);
                    }
                }
            }
        }

        // return
        return ensemble;
    }

    void print() const
    {
        size_t c_number{1};
        std::cout << Dim << "D-Lattice:\n";
        for( const auto& cell : m_lattice )
        {
            std::cout << "==================\n";
            std::cout << "\nCell No. " << c_number++ << ":\n";
            print_cell( cell );
        }
    }
};


// ===================================================================
// ============== IMPLEMENTATION OF SPECIFIC LATTICES ================
// ===================================================================
// simple cubic lattice, works in arbitrary dimensions
inline Ensemble create_SimpleCubic( const size_t num_CellsInDir, const double lattice_constant )
{
    std::array<Basis_Vector,Dim> a;
    for( size_t dir=0; dir<Dim; ++dir )
    {
        Basis_Vector a_i{};
        a_i[dir] = lattice_constant;
        a[dir] = a_i;
    }
    Unit_Cell uc{Group{Spin{}}};
    Lattice l(uc,num_CellsInDir,a);
    return l.get_ensemble();
}

// triangular lattice, works only in Dim=2
inline std::array<Basis_Vector,Dim> get_Triangular_lattice_vectors( const double lattice_constant )
{
    return std::array<Basis_Vector,Dim>{ Basis_Vector{lattice_constant,0.0} , Basis_Vector{0.5*lattice_constant,lattice_constant*std::sqrt(3.0)/2.0} };
}
inline Ensemble create_Triangular( const size_t num_CellsInDir, const double lattice_constant )
{
    if( Dim != 2 )
    {
        throw std::runtime_error( "Triangular lattice has to be 2D, but dimension is " + std::to_string(Dim) );
    };
    Unit_Cell uc{Group{Spin{}}};
    Lattice l(uc,num_CellsInDir,get_Triangular_lattice_vectors(lattice_constant));
    return l.get_ensemble();
}

}