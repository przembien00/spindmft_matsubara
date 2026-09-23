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
// ============= MAGNETIZATION TENSOR HEADER ==============
// ========================================================
/* Time-dependent single-site magnetization.  The outer index is the real-time
point and every entry is a symmetry-reduced MagnetizationVector. */
template<typename Magnetization>
class MagnetizationTensor
{
 public:
    // CONSTRUCTORS
    MagnetizationTensor() = default;
    MagnetizationTensor( const char symmetry_type, const size_t num_TimePoints );
    MagnetizationTensor( const char symmetry_type, const SiteFields& full );

    // PUBLIC METHODS
    SiteFields expand() const;
    void print( const size_t my_rank = 0 ) const;

    // GET FUNCTIONS
    size_t size() const { return m_tensor.size(); }
    bool empty() const { return m_tensor.empty(); }
    char get_symmetry() const { return m_symmetry_type; }
    size_t num_components() const { return m_directions.size(); }
    const std::vector<size_t>& get_directions() const { return m_directions; }

    // ITERATORS
    auto begin() { return m_tensor.begin(); }
    auto end() { return m_tensor.end(); }
    auto cbegin() const { return m_tensor.cbegin(); }
    auto cend() const { return m_tensor.cend(); }

    // OPERATORS
    Magnetization& operator[]( const size_t time ) { return m_tensor[time]; }
    const Magnetization& operator[]( const size_t time ) const { return m_tensor[time]; }
    Magnetization& front() { return m_tensor.front(); }
    const Magnetization& front() const { return m_tensor.front(); }
    Magnetization& back() { return m_tensor.back(); }
    const Magnetization& back() const { return m_tensor.back(); }
    MagnetizationTensor& operator+=( const MagnetizationTensor& other );
    MagnetizationTensor& operator*=( const RealType& factor );

 private:
    std::vector<Magnetization> m_tensor{};
    std::vector<size_t> m_directions{};
    char m_symmetry_type{};
};

template<typename Magnetization>
MagnetizationTensor<Magnetization> operator*(
    const RealType& factor, const MagnetizationTensor<Magnetization>& tensor );
template<typename Magnetization>
MagnetizationTensor<Magnetization> operator*(
    const MagnetizationTensor<Magnetization>& tensor, const RealType& factor );


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


// ========================================================
// ========= MAGNETIZATION TENSOR IMPLEMENTATION ==========
// ========================================================
template<typename Magnetization>
inline MagnetizationTensor<Magnetization>::MagnetizationTensor(
    const char symmetry_type, const size_t num_TimePoints ):
    m_tensor( num_TimePoints, Magnetization{symmetry_type} ),
    m_directions( determine_magnetization_directions(symmetry_type) ),
    m_symmetry_type( symmetry_type )
{}

template<typename Magnetization>
inline MagnetizationTensor<Magnetization>::MagnetizationTensor(
    const char symmetry_type, const SiteFields& full ):
    m_directions( determine_magnetization_directions(symmetry_type) ),
    m_symmetry_type( symmetry_type )
{
    m_tensor.reserve( full.size() );
    for( const auto& value : full )
    {
        m_tensor.emplace_back( symmetry_type, value );
    }
}

template<typename Magnetization>
inline SiteFields MagnetizationTensor<Magnetization>::expand() const
{
    SiteFields full{};
    full.reserve( m_tensor.size() );
    for( const auto& value : m_tensor )
    {
        full.emplace_back( value.expand() );
    }
    return full;
}

template<typename Magnetization>
inline void MagnetizationTensor<Magnetization>::print( const size_t my_rank ) const
{
    if( my_rank == 0 )
    {
        for( size_t t = 0; t < m_tensor.size(); ++t )
        {
            std::cout << "t[" << t << "]: ";
            m_tensor[t].print( my_rank );
        }
    }
}

template<typename Magnetization>
inline MagnetizationTensor<Magnetization>& MagnetizationTensor<Magnetization>::operator+=(
    const MagnetizationTensor<Magnetization>& other )
{
    if( m_symmetry_type != other.get_symmetry() || m_tensor.size() != other.size() )
    {
        error::SIZE_MISMATCH( __PRETTY_FUNCTION__ );
    }
    for( size_t t = 0; t < m_tensor.size(); ++t )
    {
        m_tensor[t] += other[t];
    }
    return *this;
}

template<typename Magnetization>
inline MagnetizationTensor<Magnetization>& MagnetizationTensor<Magnetization>::operator*=(
    const RealType& factor )
{
    for( auto& value : m_tensor )
    {
        value *= factor;
    }
    return *this;
}

template<typename Magnetization>
inline MagnetizationTensor<Magnetization> operator*(
    const RealType& factor, const MagnetizationTensor<Magnetization>& tensor )
{
    MagnetizationTensor<Magnetization> result{ tensor };
    result *= factor;
    return result;
}

template<typename Magnetization>
inline MagnetizationTensor<Magnetization> operator*(
    const MagnetizationTensor<Magnetization>& tensor, const RealType& factor )
{
    return factor * tensor;
}


};
