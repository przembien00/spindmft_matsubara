#pragma once

#include<functional>
#include<string>
#include<iostream>
#include<iomanip>
#include<sstream>

namespace Print_Routines
{

// ====================================================
// =============== CAST VALUE TO STRING ===============
// ====================================================
// rounds a value to a string with custom precision
template<typename T>
std::string round_value_to_string( const T& val, const size_t num_Digits )
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision( num_Digits ) << val;
  return ss.str();
}

// rounds an r value to a string with custom precision
template<typename T>
std::string round_value_to_string( const T&& val, const size_t num_Digits )
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision( num_Digits ) << val;
  return ss.str();
}

// converts a bool to a string
inline const std::string bool_to_string( const bool b )
{
  return b ? "true" : "false";
}

// generate a pretty output string for a given variable according to:
// variable_name           : variable_value
// <---num_PreColonSpace-->
template<typename T>
std::string quantity_to_output_line( size_t num_PreColonSpace, const std::string variable_name, const T& variable_value )
{
    int spacing = num_PreColonSpace - variable_name.size();
    if( spacing <= 0 )
    {
      spacing = 1;
    }
    std::stringstream ss;
    ss << variable_name << std::string( spacing, ' ' ) << ": " << variable_value << "\n";
    return ss.str();
}

// see above
template<typename T>
std::string quantity_to_output_line( size_t num_PreColonSpace, const std::string variable_name, const T&& variable_value )
{
    int spacing = num_PreColonSpace - variable_name.size();
    if( spacing <= 0 )
    {
      spacing = 1;
    }
    std::stringstream ss;
    ss << variable_name << std::string( spacing, ' ' ) << ": " << variable_value << "\n";
    return ss.str();
}

// function : generates a pretty output string as python comment for a given variable
// # variable_name           : variable_value
//   <---num_PreColonSpace-->
template<typename T>
inline std::string quantity_to_output_line_ht( size_t num_PreColonSpace, const std::string variable_name, const T& variable_value )
{
    int spacing = num_PreColonSpace - variable_name.size();
    if( spacing <= 0 )
    {
      spacing = 1;
    }
    std::stringstream ss;
    ss << "# " << variable_name << std::string( spacing, ' ' ) << ": " << variable_value << "\n";
    return ss.str();
}

// see above
template<typename T>
inline std::string quantity_to_output_line_ht( size_t num_PreColonSpace, const std::string variable_name, const T&& variable_value )
{
    int spacing = num_PreColonSpace - variable_name.size();
    if( spacing <= 0 )
    {
      spacing = 1;
    }
    std::stringstream ss;
    ss << "# " << variable_name << std::string( spacing, ' ' ) << ": " << variable_value << "\n";
    return ss.str();
}


// ==================================================
// ================= TRIM STRINGS ===================
// ==================================================
// cut a string with a length of more than max_length at max_length-3 adding ... 
inline void cut_if_too_large( std::string& s, const size_t max_length )
{
    if( s.size() > max_length )
    {
        s = s.substr( 0, max_length-3 );
        s += "...";
    }
}

// cut an r value string with a length of more than max_length at max_length-3 adding ... and return the cutted string
inline std::string cut_if_too_large( std::string&& s, const size_t max_length )
{
    if( s.size() > max_length )
    {
        s = s.substr( 0, max_length-3 );
        s += "...";
    }
    return s;
}

// remove zero's from the back of a string
inline void remove_zeros( std::string& s )
{
    if( s.find('.') == std::string::npos )
    {
        return; // number in s is not decimal!
    }
    while( s.back() == '0' && s.size() > 1 )
    {
        s.pop_back();
    }
    if( s.back() == '.' ) // remove the dot
    {
        s.pop_back();
    }
}

// remove zero's from the back of a string and return 
inline std::string remove_zeros( std::string&& s )
{    
    if( s.find('.') == std::string::npos )
    {
        return s; // number in s is not decimal!
    }
    while( s.back() == '0' && s.size() > 1 )
    {
        s.pop_back();
    }
    if( s.back() == '.' ) // remove the dot
    {
        s.pop_back();
    }
    return s;
}


// =================================================
// =============== STRING CHAINS ===================
// =================================================
// split a given string at a given delimiter and write it into a vector
inline std::vector<std::string> split_string_at_delimiter( const std::string& str, const char delimiter )
{
    if( str=="" ){ return std::vector<std::string>{}; }
    std::stringstream ss( str );
    std::vector<std::string> separated_str{};
    while( ss.good() ) // separate the coupling string at each delimiter
    {
      std::string substr{};
      std::getline( ss, substr, delimiter );
      separated_str.push_back( substr );
    }
    return separated_str;
}

// concatenate the strings of a given vector
inline std::string concatenate_string_with_delimiter( const std::vector<std::string>& strs, const char delimiter )
{
    if( strs.size() == 0 ){ return std::string{}; }
    std::string conc_str{};
    for( const auto& str : strs )
    {
        conc_str += str + std::string{delimiter};
    }
    conc_str.pop_back();
    return conc_str;
}

// ==============================================
// ============ PRINT ITERABLE OBJECT ===========
// ==============================================
// prints elements of a vector assuming that the outstream operator is defined for the elements
template<typename T>
void print( const std::vector<T>& v )
{
    for( auto it = v.cbegin(); it != v.cend(); ++it )
    {
        std::cout << *it << "\n";
    }
}

// prints elements of a container assuming that the outstream operator is defined for the elements
template<typename T>
std::ostream& operator<<( std::ostream& os, const std::vector<T>& v )
{
    for( auto it = v.cbegin(); it != v.cend(); ++it )
    {
        os << *it << "\n";
    }
    return os;
}


// ==============================================
// ========= PRINT BLAZE MATRIX OBJECTS =========
// ==============================================
template<typename M>
void print_matrix( const M& matrix, const std::string name = "" )
{
    if( name != "" ){ std::cout << "Matrix " << name << ":\n"; }
    for( size_t i = 0; i < matrix.rows(); ++i )
    {
        for( auto x = matrix.cbegin(i); x != matrix.cend(i); ++x )
        {
            std::cout << std::setw(8) << *x << " ";
        }
        std::cout << "\n";
    }
}

template<typename M>
void print_vector( const M& vector, const std::string name = "" )
{
    if( name != "" ){ std::cout << "Vector " << name << ":\n"; }
    for( auto x : vector )
    {
        std::cout << std::setw(8) << x << " ";
    }
    std::cout << "\n";
}

template<typename M>
void print_matrix_of_matrix( const M& matrix_of_matrix, const std::string name = "" )
{
    if( name != "" ){ std::cout << "MatrixOfMatrix " << name << ":\n"; }
    for( size_t i = 0; i < matrix_of_matrix.rows(); ++i )
    {
        size_t j{};
        for( auto x = matrix_of_matrix.cbegin(i); x != matrix_of_matrix.cend(i); ++x )
        {
            std::cout << "MatrixElement (" << i << "," << j++ << "):\n";
            print_matrix(*x);
        }
        std::cout << "\n";
    }
}

template<typename M>
void print_matrix_of_vector( const M& matrix_of_vector, const std::string name = "" )
{
    if( name != "" ){ std::cout << "MatrixOfVector " << name << ":\n"; }
    for( size_t i = 0; i < matrix_of_vector.rows(); ++i )
    {
        size_t j{};
        for( auto x = matrix_of_vector.cbegin(i); x != matrix_of_vector.cend(i); ++x )
        {
            std::cout << "MatrixElement (" << i << "," << j++ << "):\n";
            print_vector(*x);
        }
        std::cout << "\n";
    }
}


// ========================================
// =============== PRINT IF ===============
// ========================================
// print obj if my_rank is 0
template<typename T>
void print_R0( const size_t my_rank, T& obj )
{
    if( my_rank == 0 )
    {
        std::cout << obj;
    }
}

// print r value obj if my_rank is 0
template<typename T>
void print_R0( const size_t my_rank, T&& obj )
{
    if( my_rank == 0 )
    {
        std::cout << obj;
    }
}

// call print function of class object if my_rank is 0
template<typename T>
void call_print_R0( const size_t my_rank, T& obj )
{
    if( my_rank == 0 )
    {
        obj.print();
    }
}

};