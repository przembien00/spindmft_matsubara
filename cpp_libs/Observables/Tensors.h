#pragma once

#include<vector>
#include<string>
#include<iostream>
#include<functional>
#include<type_traits>
#include"Standard_Algorithms/Standard_Algorithms.h"
#include"O_Error_Handling.h"

namespace Observables::Tensors
{

// ================== NAMESPACES ========================
namespace error = Observables::Error_Handling;
namespace stda = Standard_Algorithms;


// ================== USING STATEMENTS ==================
using IndexPair = std::array<size_t,2>;
using IndexPairList = std::vector<IndexPair>;


// ================== FORWARD DECLARATIONS ==============
template<typename Mean>
class MeanTensor;
template<typename Correlation>
class CorrelationTensor;


// ================ OPERATOR DEFINITIONS ================
// Outstream Operator for CorrelationTensor
template<typename Correlation2>
std::ostream& operator<<( std::ostream& os, const CorrelationTensor<Correlation2>& tensor );


// ======================================================
// ============== CORRELATION TENSOR HEADER =============
// ======================================================
/* Abstract class that holds correlations in different directions 
depending on the symmetry 
A: saving gxx (=gyy=gzz), rest is 0
B: saving gxx (=gyy) & gzz, rest is 0
C: saving gxx (=gyy) & gxy & gyx & gzz, rest is 0
D: saving gxx & gxy & ... & gzy & gzz (all) */
template<typename Correlation>
class CorrelationTensor
{
using It = typename std::vector<Correlation>::iterator;
using constIt = typename std::vector<Correlation>::const_iterator;

 public:
    // CONSTRUCTORS
    CorrelationTensor() = default;
    explicit CorrelationTensor( const char symmetry_type );
    template<typename ParameterPack>
    CorrelationTensor( const char symmetry_type, const ParameterPack& params );
    template<typename ParameterPack>
    CorrelationTensor( const ParameterPack& params );
    CorrelationTensor( const CorrelationTensor& CT, const std::function<Correlation(const CorrelationTensor&, const IndexPair&)>& transformation_scheme );

    // PUBLIC METHODS
    void resize( const char symmetry_type );
    template<typename CorrelationContainer>
    void initialized_fill( CorrelationContainer&& CC );
    template<typename ThreeTimesThreeElements, typename ComponentOperator> 
    void operate( const ThreeTimesThreeElements& container, const ComponentOperator& operation );
    void print( const size_t my_rank = 0 ) const;

    template<typename BinaryFunction>
    void iterate( BinaryFunction f ); 
    template<typename BinaryFunction>
    void const_iterate( BinaryFunction f ) const; 
    template<typename BinaryFunction>
    void iterate2( CorrelationTensor<Correlation>& other, BinaryFunction f ); 
    template<typename Function>
    void iterate3( CorrelationTensor<Correlation>& other1, CorrelationTensor<Correlation>& other2, Function f ); 
    template<typename BinaryFunction>
    void const_iterate2( const CorrelationTensor<Correlation>& other, BinaryFunction f ) const; 

    // GET FUNCTIONS
    Correlation& get_xx();
    Correlation& get_xy();
    Correlation& get_xz();
    Correlation& get_yx();
    Correlation& get_yy();
    Correlation& get_yz();
    Correlation& get_zx();
    Correlation& get_zy();
    Correlation& get_zz();
    Correlation get_xx() const;
    Correlation get_xy() const;
    Correlation get_xz() const;
    Correlation get_yx() const;
    Correlation get_yy() const;
    Correlation get_yz() const;
    Correlation get_zx() const;
    Correlation get_zy() const;
    Correlation get_zz() const;
    size_t size() const { return m_direction_pairs.size(); };
    char get_symmetry() const { return m_symmetry_type; }
    IndexPair get_direction_pair( const size_t linear_index ) const { return m_direction_pairs[linear_index]; }
    const IndexPairList& get_direction_pairs() const { return m_direction_pairs; }
    bool is_diagonal( const size_t linear_index ) const;

    // ITERATORS
    auto begin() { return m_tensor.begin(); }
    auto end() { return m_tensor.end(); }
    auto cbegin() const { return m_tensor.cbegin(); }
    auto cend() const { return m_tensor.cend(); }

    // OPERATORS
    Correlation& operator[]( const uint index ) { return m_tensor[index]; }
    const Correlation& operator[]( const uint index ) const { return m_tensor[index]; }
    Correlation& operator()( const uint a, const uint b );
    Correlation operator()( const uint a, const uint b ) const;
    CorrelationTensor& operator+=(const CorrelationTensor& other);

    // FRIENDS
    template<typename Correlation2>
    friend std::ostream& operator<<( std::ostream& os, const CorrelationTensor<Correlation2>& tensor );

    IndexPairList m_direction_pairs{};
 private:
    // PRIVATE MEMBERS
    std::vector<Correlation> m_tensor{};

    char m_symmetry_type{};

    // PRIVATE METHODS
    void set_direction_pairs();
    Correlation zero_correlation( const std::string name ) const;
};


// ================ FUNCTION IMPLEMENTATIONS ================
inline size_t determine_number_of_correlations( const char symmetry_type )
{
    switch( symmetry_type )
    {
        case 'A':
        {
            return 1;            
        }
        case 'B':
        {
            return 2;
        }
        case 'C':
        {
            return 4;
        }
        case 'D':
        {
            return 9;
        }
        default:
        {
            error::SYMMETRY_TYPE( symmetry_type, __PRETTY_FUNCTION__ );
            return 0;
        }
    }
}

inline IndexPairList get_all_direction_pairs()
{
    IndexPairList all_dirs{};
    for( size_t i = 0; i < 3; ++i )
    {
        for( size_t j = 0; j < 3; ++j )
        {
            all_dirs.emplace_back( IndexPair{i,j} );
        }
    }
    return all_dirs;
}


// ================ OPERATOR IMPLEMENTATIONS ================
// Outstream Operator for CorrelationTensor
template<typename Correlation2>
std::ostream& operator<<( std::ostream& os, const CorrelationTensor<Correlation2>& tensor )
{
    switch( tensor.m_symmetry_type )
    {
        case 'A':
        {
            os << "g_XX = g_YY = g_ZZ:\n" << tensor.get_zz();
            break;
        }
        case 'B':
        {
            os << "g_XX = g_YY:\n" << tensor.get_xx();
            os << "\ng_ZZ:\n" << tensor.get_zz();
            break;
        }
        case 'C':
        {
            os << "g_XX = g_YY:\n"<< tensor.get_xx();
            os << "\ng_XY:\n" << tensor.get_xy();
            os << "\ng_YX:\n" << tensor.get_yx();
            os << "\ng_ZZ:\n" << tensor.get_zz();
            break;
        }
        case 'D':
        {
            os << "g_XX:\n" << tensor.get_xx();
            os << "\ng_XY:\n" << tensor.get_xy();
            os << "\ng_XZ:\n" << tensor.get_xz();
            os << "\ng_YX:\n" << tensor.get_yx();
            os << "\ng_YY:\n" << tensor.get_yy();
            os << "\ng_YZ:\n" << tensor.get_yz();
            os << "\ng_ZX:\n" << tensor.get_zx();
            os << "\ng_ZY:\n" << tensor.get_zy();
            os << "\ng_ZZ:\n" << tensor.get_zz();
            break;
        }
        default:
            error::SYMMETRY_TYPE( tensor.m_symmetry_type, __PRETTY_FUNCTION__ );
    }
    return os;
}

// scalar multiplication (scalar * tensor)
template<typename Correlation, typename Scalar>
CorrelationTensor<Correlation> operator*(const Scalar& s, const CorrelationTensor<Correlation>& tensor)
{
    CorrelationTensor<Correlation> out = tensor; // copy
    std::transform(out.cbegin(), out.cend(), out.begin(),
                   [&s](const Correlation& c)->Correlation { return s * c; });
    return out;
}

// ==========================================================================
// =============== CORRELATION TENSOR CLASS IMPLEMENTATION ==================
// ==========================================================================
// constructor from symmetry type 
template<typename Correlation>
CorrelationTensor<Correlation>::CorrelationTensor( const char symmetry_type )
{
    resize( symmetry_type );
    set_direction_pairs();
} 

// constructor initializating all Correlations equally
template<typename Correlation>
template<typename ParameterPack>
CorrelationTensor<Correlation>::CorrelationTensor( const char symmetry_type, const ParameterPack& params )
{
    resize( symmetry_type );
    std::for_each( m_tensor.begin(), m_tensor.end(),
    [&params]( Correlation& corr )
    {
        corr = Correlation{ params };
    } );
    set_direction_pairs();

}

// constructor initializing all Correlations equally
template<typename Correlation>
template<typename ParameterPack>
CorrelationTensor<Correlation>::CorrelationTensor( const ParameterPack& params )
{
    resize( params.m_symmetry_type );
    std::for_each( m_tensor.begin(), m_tensor.end(),
    [&params]( Correlation& corr )
    {
        corr = Correlation{ params.m_sub_parameter_pack };
    } );
    set_direction_pairs();
}

// constructor initializing from other CorrelationTensor using a transformation scheme
template<typename Correlation>
CorrelationTensor<Correlation>::CorrelationTensor( const CorrelationTensor& CT, const std::function<Correlation(const CorrelationTensor&, const IndexPair&)>& transformation_scheme ):
    m_symmetry_type( CT.m_symmetry_type )
{
    set_direction_pairs();
    std::transform( m_direction_pairs.cbegin(), m_direction_pairs.cend(), std::back_inserter( m_tensor ), [&]( const IndexPair& ipair )
    {
        return transformation_scheme( CT, ipair );
    } );
}


// resize from symmetry type
template<typename Correlation>
void CorrelationTensor<Correlation>::resize( const char symmetry_type )
{
    m_symmetry_type = symmetry_type;
    m_tensor.resize( determine_number_of_correlations(symmetry_type) );
} 

/* fill values from container assuming this to be initialized properly
the container handed over is assumed to have a size of 9 */
template<typename Correlation>
template<typename CorrelationContainer>
void CorrelationTensor<Correlation>::initialized_fill( CorrelationContainer&& CC )
{
    if( CC.size() != 9 )
    {
        error::INVALID_SIZE( CC.size(), __PRETTY_FUNCTION__ );
    }
    switch( m_symmetry_type )
    {
        case 'A':
        {
            this->get_zz() = Correlation{ std::move(CC[8]) };
            break;
        }
        case 'B':
        {
            this->get_xx() = Correlation{ std::move(CC[0]) };
            this->get_zz() = Correlation{ std::move(CC[8]) };
            break;
        }
        case 'C':
        {
            this->get_xx() = Correlation{ std::move(CC[0]) };
            this->get_xy() = Correlation{ std::move(CC[1]) };
            this->get_yx() = Correlation{ std::move(CC[3]) };
            this->get_zz() = Correlation{ std::move(CC[8]) };
            break;
        }
        case 'D':
        {
            this->get_xx() = Correlation{ std::move(CC[0]) };
            this->get_xy() = Correlation{ std::move(CC[1]) };
            this->get_xz() = Correlation{ std::move(CC[2]) };
            this->get_yx() = Correlation{ std::move(CC[3]) };
            this->get_yy() = Correlation{ std::move(CC[4]) };
            this->get_yz() = Correlation{ std::move(CC[5]) };
            this->get_zx() = Correlation{ std::move(CC[6]) };
            this->get_zy() = Correlation{ std::move(CC[7]) };
            this->get_zz() = Correlation{ std::move(CC[8]) };
            break;
        }
        default:
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    }
}

// operate with the correlation tensor and a 3x3-dimensional container
// ThreeTimesThreeElements must be a container with (,) operator and 3x3 elements (users responsibility)
// Operator must be a function eating a tensor element and a container element (users responsibility)
template<typename Correlation>
template<typename ThreeTimesThreeElements, typename ComponentOperator>
void CorrelationTensor<Correlation>::operate( const ThreeTimesThreeElements& container, const ComponentOperator& operation )
{
    switch( m_symmetry_type )
    {
        case 'A':
        {
            operation( this->get_xx(), container(0,0) ); // operation for xx
            break;
        }
        case 'B':
        {
            operation( this->get_xx(), container(0,0) ); // operation for xx
            operation( this->get_zz(), container(2,2) ); // operation for zz
            break;
        }
        case 'C':
        {
            operation( this->get_xx(), container(0,0) ); // operation for xx
            operation( this->get_xy(), container(0,1) ); // operation for xy
            operation( this->get_yx(), container(1,0) ); // operation for yx
            operation( this->get_zz(), container(2,2) ); // operation for zz
            break;
        }
        case 'D':
        {
            operation( this->get_xx(), container(0,0) ); // operation for xx
            operation( this->get_xy(), container(0,1) ); // operation for xy
            operation( this->get_xz(), container(0,2) ); // operation for xz
            operation( this->get_yx(), container(1,0) ); // operation for yx
            operation( this->get_yy(), container(1,1) ); // operation for yy
            operation( this->get_yz(), container(1,2) ); // operation for yz
            operation( this->get_zx(), container(2,0) ); // operation for zx
            operation( this->get_zy(), container(2,1) ); // operation for zy
            operation( this->get_zz(), container(2,2) ); // operation for zz
            break;
        }
        default:
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    }
}

// print : assumes that outstream operator is defined for class Correlation
template<typename Correlation>
void CorrelationTensor<Correlation>::print( const size_t my_rank ) const
{
    switch( m_symmetry_type )
    {
        case 'A':
        {
            std::cout << "g_XX = g_YY = g_ZZ:\n";
            this->get_zz().print(my_rank);
            break;
        }
        case 'B':
        {
            std::cout << "g_XX = g_YY:\n";
            this->get_xx().print(my_rank);
            std::cout << "\ng_ZZ:\n";
            this->get_zz().print(my_rank);
            break;
        }
        case 'C':
        {
            std::cout << "g_XX = g_YY:\n";
            this->get_xx().print(my_rank);
            std::cout << "\ng_XY:\n";
            this->get_xy().print(my_rank);
            std::cout << "\ng_YX:\n";
            this->get_yx().print(my_rank);
            std::cout << "\ng_ZZ:\n";
            this->get_zz().print(my_rank);
            break;
        }
        case 'D':
        {
            std::cout << "g_XX:\n";
            this->get_xx().print(my_rank);
            std::cout << "\ng_XY:\n";
            this->get_xy().print(my_rank);
            std::cout << "\ng_XZ:\n";
            this->get_xz().print(my_rank);
            std::cout << "\ng_YX:\n";
            this->get_yx().print(my_rank);
            std::cout << "\ng_YY:\n";
            this->get_yy().print(my_rank);
            std::cout << "\ng_YZ:\n";
            this->get_yz().print(my_rank);
            std::cout << "\ng_ZX:\n";
            this->get_zx().print(my_rank);
            std::cout << "\ng_ZY:\n";
            this->get_zy().print(my_rank);
            std::cout << "\ng_ZZ:\n";
            this->get_zz().print(my_rank);
            break;
        }
        default:
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
    }
}


// iterate over m_tensor and m_direction_pairs simultaneously using a binary function
template<typename Correlation>
template<typename BinaryFunction>
void CorrelationTensor<Correlation>::iterate( BinaryFunction f )
{
    stda::for_2each( m_tensor.begin(), m_tensor.end(), m_direction_pairs.cbegin(), f );
}

// const iterate over m_tensor and m_direction_pairs simultaneously using a binary function
template<typename Correlation>
template<typename BinaryFunction>
void CorrelationTensor<Correlation>::const_iterate( BinaryFunction f ) const
{
    stda::for_2each( m_tensor.cbegin(), m_tensor.cend(), m_direction_pairs.cbegin(), f );
}

// iterate over m_tensor, m_tensor of another CorrelationTensor and m_direction_pairs simultaneously using a trinary function
template<typename Correlation>
template<typename TrinaryFunction>
void CorrelationTensor<Correlation>::iterate2( CorrelationTensor<Correlation>& other, TrinaryFunction f )
{
    stda::for_3each( m_tensor.begin(), m_tensor.end(), other.begin(), m_direction_pairs.cbegin(), f );
}

// const iterate over m_tensor, m_tensor of another CorrelationTensor and m_direction_pairs simultaneously using a trinary function
template<typename Correlation>
template<typename TrinaryFunction>
void CorrelationTensor<Correlation>::const_iterate2( const CorrelationTensor<Correlation>& other, TrinaryFunction f ) const
{
    stda::for_3each( m_tensor.cbegin(), m_tensor.cend(), other.cbegin(), m_direction_pairs.cbegin(), f );
}

// Operators
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::operator()( const uint a, const uint b )
{   
    if( a > 2 )
    {
        error::INVALID_INDEX( a, __PRETTY_FUNCTION__ );
    }
    else if( b > 2 )
    {
        error::INVALID_INDEX( b, __PRETTY_FUNCTION__ );
    }
    else if( a == 0 )
    {
        if( b == 0 ){       return get_xx(); }
        else if( b == 1 ){  return get_xy(); }
        else{               return get_xz(); }
    }
    else if( a == 1 )
    {
        if( b == 0 ){       return get_yx(); }
        else if( b == 1 ){  return get_yy(); }
        else{               return get_yz(); }
    }
    else
    {
        if( b == 0 ){       return get_zx(); }
        else if( b == 1 ){  return get_zy(); }
        else{               return get_zz(); }
    }
    return Correlation(); // suppress compiler warnings
}

template<typename Correlation>
Correlation CorrelationTensor<Correlation>::operator()( const uint a, const uint b ) const
{   
    if( a > 2 )
    {
        error::INVALID_INDEX( a, __PRETTY_FUNCTION__ );
    }
    else if( b > 2 )
    {
        error::INVALID_INDEX( b, __PRETTY_FUNCTION__ );
    }
    else if( a == 0 )
    {
        if( b == 0 ){       return get_xx(); }
        else if( b == 1 ){  return get_xy(); }
        else{               return get_xz(); }
    }
    else if( a == 1 )
    {
        if( b == 0 ){       return get_yx(); }
        else if( b == 1 ){  return get_yy(); }
        else{               return get_yz(); }
    }
    else
    {
        if( b == 0 ){       return get_zx(); }
        else if( b == 1 ){  return get_zy(); }
        else{               return get_zz(); }
    }
    return Correlation(); // suppress compiler warnings
}

// add two correlation tensors
template<typename Correlation>
CorrelationTensor<Correlation>& CorrelationTensor<Correlation>::operator+=(const CorrelationTensor<Correlation>& other)
{
    std::transform( m_tensor.begin(), m_tensor.end(), other.cbegin(), m_tensor.begin(), 
        [](Correlation& c1, const Correlation& c2){ return c1 + c2; } );
    return *this;
}

// Get Functions
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_xx()
{
    return m_tensor[0];
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_xy()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "xy", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "xy", __PRETTY_FUNCTION__ );
        case 'C':
        return m_tensor[1];
        case 'D':
        return m_tensor[1];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_xz()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "xz", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "xz", __PRETTY_FUNCTION__ );
        case 'C':
        error::IS_ZERO( "xz", __PRETTY_FUNCTION__ );
        case 'D':
        return m_tensor[2];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_yx()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "yx", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "yx", __PRETTY_FUNCTION__ );
        case 'C':
        return m_tensor[2];
        case 'D':
        return m_tensor[3];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_yy()
{
    switch( m_symmetry_type )
    {
        case 'A':
        return m_tensor[0];
        case 'B':
        return m_tensor[0];
        case 'C':
        return m_tensor[0];
        case 'D':
        return m_tensor[4];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_yz()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "yz", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "yz", __PRETTY_FUNCTION__ );
        case 'C':
        error::IS_ZERO( "yz", __PRETTY_FUNCTION__ );
        case 'D':
        return m_tensor[5];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_zx()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "zx", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "zx", __PRETTY_FUNCTION__ );
        case 'C':
        error::IS_ZERO( "zx", __PRETTY_FUNCTION__ );
        case 'D':
        return m_tensor[6];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_zy()
{
    switch( m_symmetry_type )
    {
        case 'A':
        error::IS_ZERO( "zy", __PRETTY_FUNCTION__ );
        case 'B':
        error::IS_ZERO( "zy", __PRETTY_FUNCTION__ );
        case 'C':
        error::IS_ZERO( "zy", __PRETTY_FUNCTION__ );
        case 'D':
        return m_tensor[7];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation& CorrelationTensor<Correlation>::get_zz()
{
    switch( m_symmetry_type )
    {
        case 'A':
        return m_tensor[0];
        case 'B':
        return m_tensor[1];
        case 'C':
        return m_tensor[3];
        case 'D':
        return m_tensor[8];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}

// Get Const Functions
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_xx() const
{
    return m_tensor[0];
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_xy() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "xy" );
        case 'B':
        return zero_correlation( "xy" );
        case 'C':
        return m_tensor[1];
        case 'D':
        return m_tensor[1];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_xz() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "xz" );
        case 'B':
        return zero_correlation( "xz" );
        case 'C':
        return zero_correlation( "xz" );
        case 'D':
        return m_tensor[2];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_yx() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "yx" );
        case 'B':
        return zero_correlation( "yx" );
        case 'C':
        return m_tensor[2];
        case 'D':
        return m_tensor[3];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_yy() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return m_tensor[0];
        case 'B':
        return m_tensor[0];
        case 'C':
        return m_tensor[0];
        case 'D':
        return m_tensor[4];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_yz() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "yz" );
        case 'B':
        return zero_correlation( "yz" );
        case 'C':
        return zero_correlation( "yz" );
        case 'D':
        return m_tensor[5];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_zx() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "zx" );
        case 'B':
        return zero_correlation( "zx" );
        case 'C':
        return zero_correlation( "zx" );
        case 'D':
        return m_tensor[6];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_zy() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return zero_correlation( "zy" );
        case 'B':
        return zero_correlation( "zy" );
        case 'C':
        return zero_correlation( "zy" );
        case 'D':
        return m_tensor[7];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}
template<typename Correlation>
Correlation CorrelationTensor<Correlation>::get_zz() const
{
    switch( m_symmetry_type )
    {
        case 'A':
        return m_tensor[0];
        case 'B':
        return m_tensor[1];
        case 'C':
        return m_tensor[3];
        case 'D':
        return m_tensor[8];
        default:
        error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        return m_tensor[0];
    }
}

template<typename Correlation>
bool CorrelationTensor<Correlation>::is_diagonal( const size_t linear_index ) const
{
    switch( m_symmetry_type )
    {
        case 'A':
        {
            if( linear_index > 0 ){ error::INVALID_INDEX( linear_index, __PRETTY_FUNCTION__ ); }
            return true;            
        }
        case 'B':
        {
            if( linear_index > 1 ){ error::INVALID_INDEX( linear_index, __PRETTY_FUNCTION__ ); }
            return true;
        }
        case 'C':
        {
            if( linear_index > 4 ){ error::INVALID_INDEX( linear_index, __PRETTY_FUNCTION__ ); }
            return (linear_index==0 || linear_index==3) ? true : false;
        }
        case 'D':
        {
            if( linear_index > 8 ){ error::INVALID_INDEX( linear_index, __PRETTY_FUNCTION__ ); }
            return (linear_index==0 || linear_index==4 || linear_index==8) ? true : false;
        }
        default:
        {
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
            return false;
        }
    }
}


template<typename Correlation>
void CorrelationTensor<Correlation>::set_direction_pairs()
{
    switch( m_symmetry_type )
    {
        case 'A':
        {
            m_direction_pairs = IndexPairList{{0,0}};
            break;
        }
        case 'B':
        {
            m_direction_pairs = IndexPairList{{0,0}, {2,2}};
            break;
        }
        case 'C':
        {
            m_direction_pairs = IndexPairList{{0,0}, {0,1}, {1,0}, {2,2}};
            break;
        }
        case 'D':
        {
            m_direction_pairs = IndexPairList{{0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2}};
            break;
        }
        default:
        {
            error::SYMMETRY_TYPE( m_symmetry_type, __PRETTY_FUNCTION__ );
        }
    }
}


template <typename T, typename = int>
struct HasZeroAccordingTo : std::false_type {};

template <typename T>
struct HasZeroAccordingTo<T, decltype(&T::zero_according_to, 0)> : std::true_type {};

template<typename Correlation>
Correlation CorrelationTensor<Correlation>::zero_correlation( const std::string name ) const
{
    if constexpr ( HasZeroAccordingTo<Correlation>::value )
    {
        return m_tensor[0].zero_according_to();
    } 
    else 
    {
        error::IS_ZERO( name, __PRETTY_FUNCTION__ );
        return Correlation(); // avoid compiler warning
    }
}


};