#pragma once

#include<tuple>
#include<utility>
#include<functional>
#include<numeric>

namespace Standard_Algorithms
{

using IndexList = std::vector<size_t>;

// =========================================================
// =============== ITERATING OVER CONTAINERS ===============
// =========================================================
// iterating over 2 containers simultaneously
// only the first range terminates the iteration
template<class It1, class It2, class BinaryFunction>
constexpr BinaryFunction for_2each(It1 first_1, It1 last_1, It2 first_2, BinaryFunction f)
{
    while( first_1 != last_1 )
    {
        f( *first_1, *first_2 );
        ++first_1;
        ++first_2;
    }
    return f;
}

// iterating over 3 containers simultaneously
// only the first range terminates the iteration
template<class It1, class It2, class It3, class TrinaryFunction >
constexpr TrinaryFunction for_3each( It1 first_1, It1 last_1, It2 first_2, It3 first_3, TrinaryFunction f )
{
    while( first_1 != last_1 )
    {
        f( *first_1, *first_2, *first_3 );
        ++first_1;
        ++first_2;
        ++first_3;
    }
    return f;
}

// iterating over n containers simultaneously, should be used with more than one container
template <typename NnaryFunction, class It1, typename... RestIters>
constexpr NnaryFunction for_n_each( NnaryFunction&& f, It1 first, It1 last, RestIters&&... rest ) 
{
  while( first != last ) 
  {
    std::invoke( std::forward<NnaryFunction>( f ), *first, *std::forward<RestIters>( rest )...);
    std::make_tuple( ++rest... );
    ++first;
  }
  return f;
}


// =========================================================
// ==================== CALL FUNCTION IF ===================
// =========================================================
// call function if my_rank is 0
template<typename T>
void call_R0( const size_t my_rank, T& func )
{
    if( my_rank == 0 )
    {
        func();
    }
}

// call rvalue function if my_rank is 0
template<typename T>
void call_R0( const size_t my_rank, T&& func )
{
    if( my_rank == 0 )
    {
        func();
    }
}


// =========================================================
// ==================== SORTING ROUTINES ===================
// =========================================================
/* returns a list of indices corresponding to the sorted order of the inserted vector v
v is not changed in the process */
template<typename T, typename Compare = std::less<T>>
IndexList argsort( const std::vector<T>& v, Compare comp = Compare() )
{
    IndexList indices(v.size()); // vector of indices
    std::iota(indices.begin(), indices.end(), 0); 
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) 
    {
        return comp(v[i],v[j]);
    } );
    return indices;
}

template<typename T, typename Compare = std::less<T>>
IndexList argsort_flip( const std::vector<T>& v, Compare comp = Compare() ) 
{
    IndexList indices(v.size()); // vector of indices
    std::iota(indices.begin(), indices.end(), 0); 
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) 
    {
        return comp(v[j],v[i]); // flipped
    } );
    return indices;
}

/* extended argsort : indices are saved in a list of groups
each group contains indices which correspond to the same value according to the condition "equal" */
template<typename T1, typename T2>
std::vector<IndexList> argsort_extended( const std::vector<T1>& v, const T2 equal ) 
{
    auto s = argsort( v );
    std::vector<IndexList> s_groups{};
    s_groups.emplace_back( IndexList{s[0]} );

    auto it1 = s.cbegin(), it2 = s.cbegin();
    while(++it2 != s.cend())
    {
        if(!equal(v[*it1],v[*it2]))
        {
            s_groups.emplace_back( IndexList{} );
        }
        (--s_groups.end())->emplace_back(*it2);
        ++it1;
    }

    return s_groups;
}

template<typename T1, typename T2>
std::vector<IndexList> argsort_flip_extended( const std::vector<T1>& v, const T2 equal ) 
{
    auto s = argsort_flip( v );
    std::vector<IndexList> s_groups{};
    s_groups.emplace_back( IndexList{s[0]} );

    auto it1 = s.cbegin(), it2 = s.cbegin();
    while(++it2 != s.cend())
    {
        if(!equal(v[*it1],v[*it2]))
        {
            s_groups.emplace_back( IndexList{} );
        }
        (--s_groups.end())->emplace_back(*it2);
        ++it1;
    }

    return s_groups;
}

// returns the (first) index of the minimum element
template<typename T> // T must be iterable
size_t argmin( const T& v ) 
{
    auto min = v[0];
    size_t min_index = 0, current_index = 0;
    for( auto vit = v.cbegin()+1; vit != v.cend(); ++vit )
    {
        if( *vit < min )
        {
            min = *vit;
            min_index = current_index;
        }
        ++current_index;
    }
    return min_index;
}


// =========================================================
// ===================== COMBINATORICS =====================
// =========================================================
/* recursively find all (N over k) possibilities to draw k out of N elements
the elements are always in the initial order */
inline std::vector<std::vector<size_t>> possibilities( const std::vector<size_t>& elements, const size_t k ) 
{
    size_t N = elements.size();
    if(k>N)
    {
        return std::vector<std::vector<size_t>>{};
    }
    else if(k==N || N==0)
    {
        return std::vector<std::vector<size_t>>{elements};
    }
    else
    {
        // separate the first element
        size_t first = elements[0];
        std::vector<size_t> reduced_elements(elements.begin()+1,elements.end());

        // include the first element
        auto with_first = possibilities(reduced_elements,k-1);
        for( auto& wf : with_first )
        {
            wf.insert(wf.begin(),first); 
        }

        // exclude the first element
        auto without_first = possibilities(reduced_elements,k);

        // merge everything:
        std::vector<std::vector<size_t>> all{};
        std::merge(with_first.cbegin(), with_first.cend(), without_first.cbegin(), without_first.cend(), std::back_inserter(all));
        return all;
    }
}


};