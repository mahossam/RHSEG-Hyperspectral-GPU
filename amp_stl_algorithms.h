/*----------------------------------------------------------------------------
* Copyright (c) Microsoft Corp.
*
* Licensed under the Apache License, Version 2.0 (the "License"); you may not 
* use this file except in compliance with the License.  You may obtain a copy 
* of the License at http://www.apache.org/licenses/LICENSE-2.0  
* 
* THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED 
* WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE, 
* MERCHANTABLITY OR NON-INFRINGEMENT. 
*
* See the Apache Version 2.0 License for specific language governing 
* permissions and limitations under the License.
*---------------------------------------------------------------------------
* 
* C++ AMP standard algorithm library.
*
* This file contains the C++ AMP standard algorithms
*---------------------------------------------------------------------------*/
#pragma once

#include <xx_amp_algorithms_impl.h>
#include <amp_iterators.h>

namespace amp_stl_algorithms
{
    //----------------------------------------------------------------------------
    // all_of, any_of, none_of
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator,  typename UnaryPredicate>
    bool all_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p );

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    bool any_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p );

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    bool none_of( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p ); 

    // non-standard versions which store the result on the accelerator
    template<typename ConstRandomAccessIterator,  typename UnaryPredicate, typename OutputIterator>
    void all_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p, OutputIterator dest );

    template<typename ConstRandomAccessIterator,  typename UnaryPredicate, typename OutputIterator>
    void any_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p, OutputIterator dest );

    template<typename ConstRandomAccessIterator,  typename UnaryPredicate, typename OutputIterator>
    void none_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p, OutputIterator dest );

    //----------------------------------------------------------------------------
    // for_each
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename UnaryFunction>
    UnaryFunction for_each( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryFunction f );

    // non-standard: no return
    template<typename ConstRandomAccessIterator, typename UnaryFunction>
    void for_each_no_return( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryFunction f );

    //----------------------------------------------------------------------------
    // count, count_if
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename T>
    typename std::iterator_traits<ConstRandomAccessIterator>::difference_type
        count( ConstRandomAccessIterator first, ConstRandomAccessIterator last, const T &value ); 

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    typename std::iterator_traits<ConstRandomAccessIterator>::difference_type
        count_if( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p ); 

    //----------------------------------------------------------------------------
    // mismatch
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    std::pair<ConstRandomAccessIterator1,ConstRandomAccessIterator2>
        mismatch( ConstRandomAccessIterator1 first1, ConstRandomAccessIterator1 last1, ConstRandomAccessIterator2 first2 );

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename BinaryPredicate>
    std::pair<ConstRandomAccessIterator1,ConstRandomAccessIterator2>
        mismatch( ConstRandomAccessIterator1 first1,
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2,
        BinaryPredicate p );

    //----------------------------------------------------------------------------
    // equal
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    bool equal( ConstRandomAccessIterator1 first1, ConstRandomAccessIterator1 last1, ConstRandomAccessIterator2 first2 );

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename BinaryPredicate>
    bool equal( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1, 
        ConstRandomAccessIterator2 first2, 
        BinaryPredicate p ); 

    //----------------------------------------------------------------------------
    // find, find_if, find_if_not, find_end, find_first_of, adjacent_find
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename T>
    ConstRandomAccessIterator find( ConstRandomAccessIterator first, ConstRandomAccessIterator last, const T& value );

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    ConstRandomAccessIterator find_if(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p );

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    ConstRandomAccessIterator find_if_not( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate q );

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    ConstRandomAccessIterator1 find_end ( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename Predicate>
    ConstRandomAccessIterator1 find_end ( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        Predicate pred);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    ConstRandomAccessIterator1 find_first_of ( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename Predicate>
    ConstRandomAccessIterator1 find_first_of ( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        Predicate pred);

    template<typename ConstRandomAccessIterator>
    ConstRandomAccessIterator adjacent_find (ConstRandomAccessIterator first, ConstRandomAccessIterator last);

    template<typename ConstRandomAccessIterator, typename Predicate>
    ConstRandomAccessIterator adjacent_find (ConstRandomAccessIterator first,  ConstRandomAccessIterator last, Predicate pred);

    //----------------------------------------------------------------------------
    // search, search_n
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    ConstRandomAccessIterator1 search( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1, 
        ConstRandomAccessIterator2 first2,
        ConstRandomAccessIterator2 last2);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename Predicate>
    ConstRandomAccessIterator1 search( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        Predicate pred);

    template<typename ConstRandomAccessIterator, typename Size, typename Type>
    ConstRandomAccessIterator search_n ( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        Size count, 
        const Type& val);

    template<typename ConstRandomAccessIterator, typename Size, typename Type, typename Predicate>
    ConstRandomAccessIterator search_n ( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        Size count, 
        const Type& val, 
        Predicate pred);

    //----------------------------------------------------------------------------
    // copy, copy_if
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename RandomAccessIterator>
    RandomAccessIterator copy( ConstRandomAccessIterator first,  ConstRandomAccessIterator last, RandomAccessIterator dest_beg );

    template<typename ConstRandomAccessIterator, typename RandomAccessIterator, typename BinaryPredicate>
    RandomAccessIterator copy_if( ConstRandomAccessIterator first,  
        ConstRandomAccessIterator last,
        RandomAccessIterator dest,
        BinaryPredicate pred);

    template<typename ConstRandomAccessIterator, typename Size, typename RandomAccessIterator>
    RandomAccessIterator copy_n(ConstRandomAccessIterator first, Size count, RandomAccessIterator result);

    template<typename ConstRandomAccessIterator, typename RandomAccessIterator>
    RandomAccessIterator copy_backward( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_last ); 

    template<typename ConstRandomAccessIterator, typename RandomAccessIterator>
    RandomAccessIterator move( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator, typename RandomAccessIterator>
    RandomAccessIterator move_backward( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_last ); 

    template<typename RandomAccessIterator, typename T>
    void fill( RandomAccessIterator first, RandomAccessIterator last, const T& value );

    template<typename RandomAccessIterator, typename Size, typename T>
    void fill_n( RandomAccessIterator first, Size count, const T& value );

    /*
    template<typename RandomAccessIterator, typename Size, typename T>
    RandomAccessIterator fill_n( RandomAccessIterator first, Size count, const T& value );*/

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename UnaryFunction>
    RandomAccessIterator transform( ConstRandomAccessIterator first1, 
        ConstRandomAccessIterator last1, 
        RandomAccessIterator result,
        UnaryFunction func);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename BinaryFunction>
    RandomAccessIterator transform( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        RandomAccessIterator result,
        BinaryFunction func);

    template<typename RandomAccessIterator, typename Generator>
    void generate( RandomAccessIterator first, RandomAccessIterator last, Generator g );

    template<typename RandomAccessIterator, typename Size, typename Generator>
    void generate_n( RandomAccessIterator first, Size count, Generator g );

    /*
    template<typename RandomAccessIterator, typename Size, typename Generator>
    RandomAccessIterator generate_n( RandomAccessIterator first, Size count, Generator g );*/

    template<typename RandomAccessIterator, typename T>
    RandomAccessIterator remove( RandomAccessIterator first, RandomAccessIterator last, const T& value );

    template<typename RandomAccessIterator, typename UnaryPredicate>
    RandomAccessIterator remove_if( RandomAccessIterator first, RandomAccessIterator last, UnaryPredicate p );

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename T>
    RandomAccessIterator remove_copy( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first,
        const T& value );

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename UnaryPredicate>
    RandomAccessIterator remove_copy_if( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first,
        UnaryPredicate p );

    template<typename RandomAccessIterator, typename T>
    void replace( RandomAccessIterator first, 
        RandomAccessIterator last,
        const T& old_value, 
        const T& new_value );

    template<typename RandomAccessIterator, typename UnaryPredicate, typename T>
    void replace_if( RandomAccessIterator first, 
        RandomAccessIterator last,
        UnaryPredicate p, 
        const T& new_value ); 



    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename T>
    RandomAccessIterator replace_copy( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first,
        const T& old_value,
        const T& new_value );

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename UnaryPredicate, typename T>
    RandomAccessIterator replace_copy_if( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first,
        UnaryPredicate p,
        const T& new_value ); 

    template<typename T>
    void swap( T& a, T& b ) restrict(cpu,amp);

    template<typename T, size_t N>
    void swap( T (&a)[N], T (&b)[N]) restrict(cpu,amp);

    template<typename RandomAccessIterator1, typename RandomAccessIterator2>
    RandomAccessIterator2 swap_ranges( RandomAccessIterator1 first1,
        RandomAccessIterator1 last1, 
        RandomAccessIterator2 first2 );

    template<typename RandomAccessIterator1, typename RandomAccessIterator2>
    void iter_swap( RandomAccessIterator1 a, RandomAccessIterator2 b ) restrict(cpu,amp);

    template<typename RandomAccessIterator>
    void reverse( RandomAccessIterator first, RandomAccessIterator last );

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator>
    RandomAccessIterator reverse_copy( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last, 
        RandomAccessIterator result);

    template<typename RandomAccessIterator>
    void rotate( RandomAccessIterator first, RandomAccessIterator middle, RandomAccessIterator last);

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator>
    RandomAccessIterator rotate_copy( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator n_first,
        ConstRandomAccessIterator last, 
        RandomAccessIterator d_first );

    template<typename RandomAccessIterator>
    void random_shuffle( RandomAccessIterator first, RandomAccessIterator last );

    template<typename RandomAccessIterator, typename RandomNumberGenerator>
    void random_shuffle( RandomAccessIterator first, 
        RandomAccessIterator last,
        RandomNumberGenerator& r );

    template<typename RandomAccessIterator, typename RandomNumberGenerator>
    void random_shuffle( RandomAccessIterator first, 
        RandomAccessIterator last, 
        RandomNumberGenerator&& r ); 

    template<typename RandomAccessIterator, typename UniformRandomNumberGenerator>
    void shuffle( RandomAccessIterator first, 
        RandomAccessIterator last, 
        UniformRandomNumberGenerator&& g ); 

    template<typename ConstRandomAccessIterator>
    ConstRandomAccessIterator unique( ConstRandomAccessIterator first, ConstRandomAccessIterator last);

    template<typename ConstRandomAccessIterator, typename BinaryPredicate>
    ConstRandomAccessIterator unique( ConstRandomAccessIterator first, ConstRandomAccessIterator last, BinaryPredicate comp);

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator>
    ConstRandomAccessIterator unique_copy( ConstRandomAccessIterator first, ConstRandomAccessIterator last, RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename BinaryPredicate>
    ConstRandomAccessIterator unique_copy( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last, 
        RandomAccessIterator d_first, 
        BinaryPredicate p); 

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    bool is_partitioned( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p );

    template<typename RandomAccessIterator, typename UnaryPredicate>
    RandomAccessIterator partition( RandomAccessIterator first, RandomAccessIterator last, UnaryPredicate comp);

    template<typename RandomAccessIterator, typename UnaryPredicate>
    RandomAccessIterator stable_partition( RandomAccessIterator first, RandomAccessIterator last, UnaryPredicate p );

    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    ConstRandomAccessIterator partition_point( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p);



    template<typename ConstRandomAccessIterator>
    bool is_sorted( ConstRandomAccessIterator first, ConstRandomAccessIterator last );

    template<typename ConstRandomAccessIterator, typename Compare>
    bool is_sorted( ConstRandomAccessIterator first, ConstRandomAccessIterator last, Compare comp ); 

    template<typename ConstRandomAccessIterator>
    ConstRandomAccessIterator is_sorted_until( ConstRandomAccessIterator first, ConstRandomAccessIterator last );

    template<typename ConstRandomAccessIterator, typename Compare>
    ConstRandomAccessIterator is_sorted_until( ConstRandomAccessIterator first, ConstRandomAccessIterator last, Compare comp ); 

    template<typename RandomAccessIterator>
    void sort( RandomAccessIterator first, RandomAccessIterator last );

    template<typename RandomAccessIterator, typename Compare>
    void sort( RandomAccessIterator first, RandomAccessIterator last, Compare comp ); 

    template<typename RandomAccessIterator>
    void partial_sort( RandomAccessIterator first, 
        RandomAccessIterator middle, 
        RandomAccessIterator last );

    template<typename RandomAccessIterator, typename Compare>
    void partial_sort( RandomAccessIterator first, 
        RandomAccessIterator middle,
        RandomAccessIterator last, Compare comp );

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator>
    RandomAccessIterator partial_sort_copy( ConstRandomAccessIterator first,
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first, 
        RandomAccessIterator d_last ); 

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename Compare>
    RandomAccessIterator partial_sort_copy( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        RandomAccessIterator d_first, 
        RandomAccessIterator d_last,
        Compare comp ); 

    template<typename RandomAccessIterator>
    void stable_sort( RandomAccessIterator first, RandomAccessIterator last );

    template<typename RandomAccessIterator, typename Compare>
    void stable_sort( RandomAccessIterator first, RandomAccessIterator last, Compare comp ); 

    template<typename RandomAccessIterator>
    void nth_element( RandomAccessIterator first, 
        RandomAccessIterator nth, 
        RandomAccessIterator last ); 

    template<typename RandomAccessIterator, typename Compare>
    void nth_element( RandomAccessIterator first, 
        RandomAccessIterator nth,
        RandomAccessIterator last, Compare comp ); 

    template<typename ConstRandomAccessIterator, typename T>
    ConstRandomAccessIterator lower_bound( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value ); 

    template<typename ConstRandomAccessIterator, typename T, typename Compare>
    ConstRandomAccessIterator lower_bound( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value, Compare comp ); 

    template<typename ConstRandomAccessIterator, typename T>
    ConstRandomAccessIterator upper_bound( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value ); 

    template<typename ConstRandomAccessIterator, typename T, typename Compare>
    ConstRandomAccessIterator upper_bound( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value, Compare comp ); 



    template<typename ConstRandomAccessIterator, typename T>
    bool binary_search( ConstRandomAccessIterator first, ConstRandomAccessIterator last, const T& value ); 

    template<typename ConstRandomAccessIterator, typename T, typename Compare>
    bool binary_search( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value, 
        Compare comp );

    template<typename ConstRandomAccessIterator, typename T>
    std::pair<ConstRandomAccessIterator,ConstRandomAccessIterator> 
        equal_range( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value ); 

    template<typename ConstRandomAccessIterator, typename T, typename Compare>
    std::pair<ConstRandomAccessIterator,ConstRandomAccessIterator> 
        equal_range( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        const T& value, 
        Compare comp ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator>
    RandomAccessIterator merge( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 fast1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        RandomAccessIterator result);

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename BinaryPredicate>
    RandomAccessIterator merge( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        RandomAccessIterator result,
        BinaryPredicate comp);

    template<typename RandomAccessIterator>
    void inplace_merge( RandomAccessIterator first,
        RandomAccessIterator middle,
        RandomAccessIterator last ); 

    template<typename RandomAccessIterator, typename Compare>
    void inplace_merge( RandomAccessIterator first,
        RandomAccessIterator middle,
        RandomAccessIterator last,
        Compare comp ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    bool includes( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2 ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename Compare>
    bool includes( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2, 
        Compare comp );

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator>
    RandomAccessIterator set_difference( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename Compare>
    RandomAccessIterator set_difference( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first, 
        Compare comp ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator>
    RandomAccessIterator set_intersection( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename Compare>
    RandomAccessIterator set_intersection( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first, 
        Compare comp ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator>
    RandomAccessIterator set_symmetric_difference( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename Compare>
    RandomAccessIterator set_symmetric_difference( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first, 
        Compare comp); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator>
    RandomAccessIterator set_union( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename Compare>
    RandomAccessIterator set_union( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        RandomAccessIterator d_first, 
        Compare comp ); 


    // TODO: consider supporting the heap functions (is_heap et al)

#ifdef max
#error amp_stl_algorithms encoutered a definition of the macro max
#endif

    template<typename T> 
    const T& max( const T& a, const T& b ) restrict (cpu,amp);

    template<typename T, typename Compare>
    const T& max( const T& a, const T& b, Compare comp ) restrict (cpu,amp);

#ifdef max
#error amp_stl_algorithms encoutered a definition of the macro max
#endif

    template<typename T> 
    const T& min( const T& a, const T& b ) restrict (cpu,amp);

    template<typename T, typename Compare>
    const T& min( const T& a, const T& b, Compare comp ) restrict (cpu,amp);

    template<typename T> 
    std::pair<const T&,const T&> minmax( const T& a, const T& b ) restrict (cpu,amp);

    template<typename T, typename Compare>
    std::pair<const T&,const T&> minmax( const T& a, const T& b, Compare comp ) restrict (cpu,amp);

    // TODO: enable initializer list in amp restricted code
    //
    // template<typename T>
    // T max( std::initializer_list<T> ilist) restrict(cpu,amp);
    //
    // template<typename T, typename Compare>
    // T max( std::initializer_list<T> ilist, Compare comp ); 
    //
    // template<typename T>
    // T min( std::initializer_list<T> ilist) restrict(cpu,amp);
    //
    // template<typename T, typename Compare>
    // T min( std::initializer_list<T> ilist, Compare comp ); 
    //
    // template<typename T>
    // std::pair<T,T> minmax( std::initializer_list ilist);
    //
    // template<typename T, typename Compare>
    // std::pair<T,T> minmax( std::initializer_list ilist, Compare comp ); 

    template<typename ConstRandomAccessIterator> 
    ConstRandomAccessIterator max_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last );

    template<typename ConstRandomAccessIterator, typename Compare>
    ConstRandomAccessIterator max_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last, Compare comp ); 

    template<typename ConstRandomAccessIterator> 
    ConstRandomAccessIterator min_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last );

    template<typename ConstRandomAccessIterator, typename Compare>
    ConstRandomAccessIterator min_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last, Compare comp ); 

    template<typename ConstRandomAccessIterator> 
    std::pair<ConstRandomAccessIterator,ConstRandomAccessIterator> 
        minmax_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last ); 

    template<typename ConstRandomAccessIterator, typename Compare>
    std::pair<ConstRandomAccessIterator,ConstRandomAccessIterator> 
        minmax_element( ConstRandomAccessIterator first, ConstRandomAccessIterator last, Compare comp ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2>
    bool lexicographical_compare( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2 ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename Compare>
    bool lexicographical_compare( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        ConstRandomAccessIterator2 last2,
        Compare comp );

    // TODO: add decls for is_permuation, next_permutation, prev_permutation

    template<typename RandomAccessIterator, typename T>
    void iota( RandomAccessIterator first, RandomAccessIterator last, T value );

    // non-standard
    template<typename ConstRandomAccessIterator, typename T>
    T reduce( ConstRandomAccessIterator first, ConstRandomAccessIterator last, T init );

    // non-standard
    template<typename ConstRandomAccessIterator, typename T, typename BinaryOperation>
    T reduce( ConstRandomAccessIterator first, ConstRandomAccessIterator last, T init, BinaryOperation op ); 

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2, typename T>
    T inner_product( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2, 
        T value ); 

    template<
        typename ConstRandomAccessIterator1,
        typename ConstRandomAccessIterator2,
        typename T,
        typename BinaryOperation1,
        typename BinaryOperation2>
        T inner_product( ConstRandomAccessIterator1 first1, 
        ConstRandomAccessIterator1 last1,
        ConstRandomAccessIterator2 first2,
        T value,
        BinaryOperation1 op1,
        BinaryOperation2 op2 ); 

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator>
    RandomAccessIterator adjacent_difference( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last, 
        RandomAccessIterator d_first ); 

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename BinaryOperation>
    RandomAccessIterator adjacent_difference( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last, 
        RandomAccessIterator d_first,
        BinaryOperation op );

    template <typename ConstRandomAccessIterator,typename RandomAccessIterator>
    RandomAccessIterator partial_sum( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        RandomAccessIterator result );

    template <typename ConstRandomAccessIterator,typename RandomAccessIterator, typename BinaryOperation>
    RandomAccessIterator partial_sum( ConstRandomAccessIterator first, 
        ConstRandomAccessIterator last,
        RandomAccessIterator result, 
        BinaryOperation binary_op );

}// namespace amp_stl_algorithms

#include <xx_amp_stl_algorithms_impl_inl.h>
