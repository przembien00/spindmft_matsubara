#pragma once

#include<vector>
#include<iostream>
#include"Globals/Types.h"
#include"Globals/Matrix_Types.h"
#include"O_Error_Handling.h"

namespace Observables::Magnetization
{
// ===================== NAMESPACES =======================
namespace error = Observables::Error_Handling;


// ===================== USING STATEMENTS =================
using SiteFields = std::vector<FieldVector>;


// ================ FUNCTION IMPLEMENTATIONS ================
/* directions (0=x, 1=y, 2=z) of the thermal first moments <S^a> permitted by the symmetry type
A: full spin isotropy, all first moments vanish
B, C: axial symmetry about z, only <S^z> survives
D: no constraints */
inline std::vector<size_t> determine_magnetization_directions( const char symmetry_type )
{
    switch( symmetry_type )
    {
        case 'A':
        {
            return {};
        }
        case 'B':
        case 'C':
        {
            return { 2 };
        }
        case 'D':
        {
            return { 0, 1, 2 };
        }
        default:
        {
            error::SYMMETRY_TYPE( symmetry_type, __PRETTY_FUNCTION__ );
            return {};
        }
    }
}


// ========================================================
// ============= MAGNETIZATION VECTOR HEADER ==============
// ========================================================
/* contains the thermal first moments <S^a> of a single site, storing only the components the
symmetry type permits (the analogue of CorrelationTensor for the second moments; a Vector
rather than a Tensor because the magnetization carries no imaginary-time dependence).
Components forbidden by the symmetry have no storage, so they are zero by construction and
cannot drift under Monte-Carlo noise. */
class MagnetizationVector
{
 public:
    // CONSTRUCTORS
    MagnetizationVector() = default;
    explicit MagnetizationVector( const char symmetry_type );
    MagnetizationVector( const char symmetry_type, const FieldVector& full ); // projects the full vector

    // PUBLIC METHODS
    FieldVector expand() const; // full 3-vector with zeros for the forbidden components
    void print( const size_t my_rank = 0 ) const;

    // GET FUNCTIONS
    size_t size() const { return m_components.size(); }
    char get_symmetry() const { return m_symmetry_type; }
    size_t get_direction( const size_t linear_index ) const { return m_directions[linear_index]; }
    const std::vector<size_t>& get_directions() const { return m_directions; }

    // ITERATORS
    auto begin() { return m_components.begin(); }
    auto end() { return m_components.end(); }
    auto cbegin() const { return m_components.cbegin(); }
    auto cend() const { return m_components.cend(); }

    // OPERATORS
    RealType& operator[]( const size_t linear_index ) { return m_components[linear_index]; }
    const RealType& operator[]( const size_t linear_index ) const { return m_components[linear_index]; }
    MagnetizationVector& operator+=( const MagnetizationVector& other );
    MagnetizationVector& operator*=( const RealType& factor );

 private:
    // PRIVATE MEMBERS
    std::vector<RealType> m_components{};
    std::vector<size_t> m_directions{};
    char m_symmetry_type{};
};


// ========================================================
// ========= MAGNETIZATION VECTOR IMPLEMENTATION ==========
// ========================================================
// constructor from symmetry type, components initialized to zero
inline MagnetizationVector::MagnetizationVector( const char symmetry_type ):
    m_directions( determine_magnetization_directions(symmetry_type) ),
    m_symmetry_type( symmetry_type )
{
    m_components.resize( m_directions.size(), RealType{0.} );
}

// constructor from symmetry type and a full 3-vector, keeping only the permitted components
inline MagnetizationVector::MagnetizationVector( const char symmetry_type, const FieldVector& full ):
    MagnetizationVector( symmetry_type )
{
    for( size_t c = 0; c < m_components.size(); ++c )
    {
        m_components[c] = full[ m_directions[c] ];
    }
}

// return the full 3-vector, with zeros for the components the symmetry forbids
inline FieldVector MagnetizationVector::expand() const
{
    FieldVector full{ 0., 0., 0. };
    for( size_t c = 0; c < m_components.size(); ++c )
    {
        full[ m_directions[c] ] = m_components[c];
    }
    return full;
}

// print
inline void MagnetizationVector::print( const size_t my_rank ) const
{
    if( my_rank == 0 )
    {
        const std::string direction_names = "xyz";
        for( size_t c = 0; c < m_components.size(); ++c )
        {
            std::cout << "<S^" << direction_names[ m_directions[c] ] << "> = " << m_components[c] << " ";
        }
        std::cout << "\n";
    }
}

// add another magnetization vector
inline MagnetizationVector& MagnetizationVector::operator+=( const MagnetizationVector& other )
{
    if( m_components.size() != other.size() )
    {
        error::SIZE_MISMATCH( __PRETTY_FUNCTION__ );
    }
    for( size_t c = 0; c < m_components.size(); ++c )
    {
        m_components[c] += other[c];
    }
    return *this;
}

// multiply assign with factor
inline MagnetizationVector& MagnetizationVector::operator*=( const RealType& factor )
{
    for( auto& component : m_components )
    {
        component *= factor;
    }
    return *this;
}


};
