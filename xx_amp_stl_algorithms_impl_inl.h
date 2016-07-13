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
* This file contains the implementation for C++ AMP standard algorithms
*---------------------------------------------------------------------------*/
#pragma once

#include <amp_stl_algorithms.h>
#include <functional>
#include <amp_algorithms.h>

namespace amp_stl_algorithms
{

    namespace _details
    {
        template<class ConstRandomAccessIterator>
        concurrency::array_view<typename std::iterator_traits<ConstRandomAccessIterator>::value_type> 
            create_section(ConstRandomAccessIterator iter, typename std::iterator_traits<ConstRandomAccessIterator>::difference_type distance) 
        {
            typedef std::iterator_traits<ConstRandomAccessIterator>::value_type value_type;
            typedef std::iterator_traits<ConstRandomAccessIterator>::difference_type difference_type;
            auto base_view = _details::array_view_iterator_helper<value_type>::get_base_array_view(iter);
            difference_type start = std::distance(begin(base_view), iter);
            return base_view.section(concurrency::index<1>(start), concurrency::extent<1>(distance));
        }		
    }

    //----------------------------------------------------------------------------
    // for_each
    //----------------------------------------------------------------------------
    template< typename ConstRandomAccessIterator, typename UnaryFunction >
    void for_each_no_return( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryFunction f )
    {	
        typedef std::iterator_traits<ConstRandomAccessIterator>::difference_type difference_type;

        difference_type distance = std::distance(first, last);
        if (distance <= 0) 
        {
            return;
        }

        auto section_view = _details::create_section(first, distance);
        concurrency::parallel_for_each(concurrency::extent<1>(distance), [f,section_view] (concurrency::index<1> idx) restrict(amp) {
            f(section_view[idx]);
        });
    }

    // UnaryFunction CANNOT contain any array, array_view or textures. Needs to be blittable.
    template< typename ConstRandomAccessIterator, typename UnaryFunction >
    UnaryFunction for_each( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryFunction f )
    {
        concurrency::array_view<UnaryFunction> functor_av(concurrency::extent<1>(1), &f);
        for_each_no_return( 
            first, last, 
            [functor_av] (const decltype(*first)& val) restrict (amp)
        {
            functor_av(0)(val);
        }
        );
        functor_av.synchronize();
        return f;
    }

    //----------------------------------------------------------------------------
    // generate, generate_n
    //
    // The "Generator" functor needs to be callable as "g()" and must return a type
    // that is assignable to RandomAccessIterator::value_type.  The functor needs
    // to be blittable and cannot contain any array, array_view, or textures.
    //----------------------------------------------------------------------------

    template<typename RandomAccessIterator, typename Size, typename Generator>
    void generate_n(RandomAccessIterator begin, Size count, Generator g)
    {
        if (count <= 0) 
            return;

        auto section_view = _details::create_section(begin, count);

        concurrency::parallel_for_each(section_view.extent, [g,section_view] (concurrency::index<1> idx) restrict(amp) {
            section_view[idx] = g();
        });
    }

    template <typename RandomAccessIterator, typename Generator>
    void generate(RandomAccessIterator begin, RandomAccessIterator end, Generator g)
    {
        typedef std::iterator_traits<RandomAccessIterator>::difference_type difference_type;

        difference_type distance = std::distance(begin, end);

        amp_stl_algorithms::generate_n(begin, distance, g);
    }

    //----------------------------------------------------------------------------
    // fill, fill_n
    //----------------------------------------------------------------------------

    template<typename RandomAccessIterator, typename T>
    void fill( RandomAccessIterator begin, RandomAccessIterator end, const T& value )
    {
        amp_stl_algorithms::generate(begin, end, [value] () restrict(amp) { return value; });
    }

    template<typename RandomAccessIterator, typename Size, typename T>
    void fill_n( RandomAccessIterator begin, Size count, const T& value )
    {
        amp_stl_algorithms::generate_n(begin, count, [value] () restrict(amp) { return value; });
    }

    //----------------------------------------------------------------------------
    // transform (Unary)
    //
    // The "UnaryFunction" functor needs to be callable as "func(ConstRandomAccessIterator::value_type)".
    // The functor needs to be blittable and cannot contain any array, array_view, or textures.
    //----------------------------------------------------------------------------

    template<typename ConstRandomAccessIterator,typename RandomAccessIterator, typename UnaryFunction>
    RandomAccessIterator transform( ConstRandomAccessIterator begin1, 
        ConstRandomAccessIterator end1, 
        RandomAccessIterator result_begin,
        UnaryFunction func)
    {
        typedef std::iterator_traits<ConstRandomAccessIterator>::difference_type difference_type;

        difference_type distance = std::distance(begin1, end1);
        if (distance <= 0) 
            return result_begin;

        auto input_view = _details::create_section(begin1, distance);
        auto output_view = _details::create_section(result_begin, distance);

        concurrency::parallel_for_each(output_view.extent, [func,input_view,output_view] (concurrency::index<1> idx) restrict(amp) {
            output_view[idx] = func(input_view[idx]);
        });

        return result_begin;
    }


    //----------------------------------------------------------------------------
    // transform (Binary)
    //
    // The "BinaryFunction" functor needs to be callable as "func(ConstRandomAccessIterator1::value_type, ConstRandomAccessIterator2::value_type)".
    // The functor needs to be blittable and cannot contain any array, array_view, or textures.
    //----------------------------------------------------------------------------

    template<typename ConstRandomAccessIterator1, typename ConstRandomAccessIterator2,typename RandomAccessIterator, typename BinaryFunction>
    RandomAccessIterator transform( ConstRandomAccessIterator1 begin1, 
        ConstRandomAccessIterator1 end1,
        ConstRandomAccessIterator2 begin2, 
        RandomAccessIterator result_begin,
        BinaryFunction func)
    {
        typedef std::iterator_traits<ConstRandomAccessIterator1>::difference_type difference_type;

        difference_type distance = std::distance(begin1, end1);
        if (distance <= 0) 
            return result_begin;

        auto input1_view = _details::create_section(begin1, distance);
        auto input2_view = _details::create_section(begin2, distance);
        auto output_view = _details::create_section(result_begin, distance);

        concurrency::parallel_for_each(output_view.extent, [func,input1_view,input2_view,output_view] (concurrency::index<1> idx) restrict(amp) {
            output_view[idx] = func(input1_view[idx], input2_view[idx]);
        });

        return result_begin;
    }


    //----------------------------------------------------------------------------
    // all_of, any_of, none_of
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator,  typename UnaryPredicate >
    bool all_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p ) 
    {
        return !amp_stl_algorithms::any_of(
            first, last, 
            [p] (const decltype(*first)& val) restrict(amp) { return !p(val); }
        );
    }


    // Non-standard, OutputIterator must yield an int referenece, where the result will be
    // stored. This allows the function to eschew synchronization
    template<typename ConstRandomAccessIterator, typename UnaryPredicate, typename OutputIterator >
    void any_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p,  OutputIterator dest)
    {
        auto section_view = _details::create_section(dest, 1);
        amp_stl_algorithms::for_each_no_return(
            first, last, 
            [section_view, p] (const decltype(*first)& val) restrict(amp) 
        {
            int *accumulator = &section_view(0);
            if (*accumulator == 0)
            {
                if (p(val))
                {
                    concurrency::atomic_exchange(accumulator, 1);
                }
            }
        }
        );
    }

    // Standard, builds of top of the non-standard async version above, and adds a sync to
    // materialize the result.
    template<typename ConstRandomAccessIterator, typename UnaryPredicate >
    bool any_of(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p )
    {
        int found_any_storage = 0;
        concurrency::array_view<int> found_any_av(concurrency::extent<1>(1), &found_any_storage);
        amp_stl_algorithms::any_of(first, last, p, amp_stl_algorithms::begin(found_any_av));
        found_any_av.synchronize();
        return found_any_storage == 1;
    }

    template<typename ConstRandomAccessIterator, typename UnaryPredicate >
    bool none_of( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p )
    {
        return !amp_stl_algorithms::any_of(first, last, p);
    }


    //----------------------------------------------------------------------------
    // count
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename T >
    typename std::iterator_traits<ConstRandomAccessIterator>::difference_type
        count( ConstRandomAccessIterator first, ConstRandomAccessIterator last, const T &value )
    {
        return amp_stl_algorithms::count_if(
            first, last, 
            [value] (const decltype(*first)& cur_val) restrict(amp) { return cur_val==value; }
        );
    }

    //----------------------------------------------------------------------------
    // count_if
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename UnaryPredicate >
    typename std::iterator_traits<ConstRandomAccessIterator>::difference_type
        count_if( ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p )
    {
        typedef typename std::iterator_traits<ConstRandomAccessIterator>::difference_type diff_type;

        static const int num_threads = 10 * 1024;
        diff_type count = 0;
        concurrency::array_view<diff_type> count_av(1, &count);

        diff_type element_count = std::distance(first, last);
        auto section_view = _details::create_section(first, element_count);

        concurrency::parallel_for_each(
            concurrency::extent<1>(num_threads),
            [section_view,element_count,p,count_av] (concurrency::index<1> idx) restrict (amp) 
        {
            int tid = idx[0];
            diff_type local_count = 0;
            for (diff_type i=tid; i<element_count; i += num_threads) 
            {
                if (p(section_view(i)))
                {
                    local_count++;
                }
            }
            if (local_count > 0)
            {
                concurrency::atomic_fetch_add(&count_av(0), local_count);
            }
        }
        );

        count_av.synchronize();
        return count;
    }

    //----------------------------------------------------------------------------
    // find_if
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename UnaryPredicate>
    ConstRandomAccessIterator find_if(ConstRandomAccessIterator first, ConstRandomAccessIterator last, UnaryPredicate p )
    {
        typedef std::iterator_traits<ConstRandomAccessIterator>::difference_type difference_type;

        difference_type distance = std::distance(first, last);
        if (distance <= 0) 
        {
            return last;
        }

        difference_type result_position = distance;
        concurrency::array_view<int> result_av(concurrency::extent<1>(1), &result_position);

        auto section_view = _details::create_section(first, distance);

        concurrency::parallel_for_each(concurrency::extent<1>(distance), [=] (concurrency::index<1> idx) restrict(amp) {
            int i = idx[0];
            if (p(section_view[idx]))
            {
                concurrency::atomic_fetch_min(&result_av(0), i);
            }
        });

        result_av.synchronize();
        return first + result_position;
    }

    //----------------------------------------------------------------------------
    // find
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename T>
    ConstRandomAccessIterator find( ConstRandomAccessIterator first, ConstRandomAccessIterator last, const T& value )
    {
        return amp_stl_algorithms::find_if(first, last, [=] (const decltype(*first)& curr_val) restrict(amp) {
            return curr_val == value;
        });
    }

    //----------------------------------------------------------------------------
    // reduce
    //----------------------------------------------------------------------------
    template<typename ConstRandomAccessIterator, typename T, typename BinaryOperation>
    T reduce( ConstRandomAccessIterator first, ConstRandomAccessIterator last, T init, BinaryOperation op )
    {
        typedef typename std::iterator_traits<ConstRandomAccessIterator>::difference_type diff_type;
        diff_type element_count = std::distance(first, last);
        auto section_view = _details::create_section(first, element_count);

        return op(init, amp_algorithms::reduce(section_view, op));
    }

    template<typename ConstRandomAccessIterator, typename T>
    T reduce( ConstRandomAccessIterator first, ConstRandomAccessIterator last, T init )
    {
        return amp_stl_algorithms::reduce(first, last, init, amp_algorithms::sum<std::remove_const<typename std::iterator_traits<ConstRandomAccessIterator>::value_type>::type>());
    }

}// namespace amp_stl_algorithms
