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
* C++ AMP standard algorithms library.
*
* This file contains the iterator classes for C++ AMP containers
*---------------------------------------------------------------------------*/

#pragma once

#include <amp.h>
#include <iterator>
#include <type_traits>
#include <xx_amp_stl_algorithms_impl.h>

namespace amp_stl_algorithms
{
    using namespace concurrency;
    //----------------------------------------------------------------------------
    // array_view_iterator
    //
    // Provides a local iterator for an array_view. Iterators are only comparable
    // for the same array_view (copies of the array_view are shallow, and are OK,
    // i.e., iterators obtained from different copies are comparable).
    // The behavior is undefined on comparing iterators obtained from different 
    // array_views.
    //
    // Note: only array_view<T, 1> is supported, since only array_view<T, 1> could
    // ensure the linear continuous storage.
    //----------------------------------------------------------------------------

    // TODO: array_view_iterator<value_type> should freely be usable as 
    // array_view_iterator<const value_type>. This currently does not work 
    // properly at several places and needs to be fixed.

    namespace _details
    {
        template <typename value_type>
        struct array_view_iterator_helper;
    }

    template <typename value_type>
    class array_view_iterator : public std::iterator<std::random_access_iterator_tag, typename value_type, int>
    {
        template <typename value_type>
        friend struct _details::array_view_iterator_helper;

        template <typename value_type>
        friend array_view_iterator<value_type> begin(const concurrency::array_view<value_type>& arr) restrict(cpu,amp);

        template <typename value_type>
        friend array_view_iterator<value_type> end(const concurrency::array_view<value_type>& arr) restrict(cpu,amp);

    public:
        ~array_view_iterator() restrict(cpu,amp)
        {
        }

        array_view_iterator(const array_view_iterator& other) restrict(cpu,amp)
            : m_base_view(other.m_base_view), 
            m_position(other.m_position)
        {
        }

        // TODO: Should be able to copy construct a const array_view iterator from a non-const array_view iterator

        array_view_iterator() restrict(cpu,amp)
            : m_base_view(_details::empty_array_view_factory<value_type>::create())
        {
        }

        array_view_iterator& operator=(const array_view_iterator& iter) restrict(cpu,amp)
        {
            if (this != &iter)
            {
                m_base_view = iter.m_base_view;
                m_position = iter.m_position;
            }

            return *this;
        }

        // TODO: Should be able to assign a non-const array_view iterator to a const array_view iterator

        // Prefix
        array_view_iterator& operator ++() restrict(cpu,amp)
        {
            m_position++;
            return *this;
        }

        // Postfix
        array_view_iterator operator ++(int) restrict(cpu,amp)
        {
            auto temp = *this;
            m_position++;
            return temp;
        }

        // Prefix
        array_view_iterator& operator --() restrict(cpu,amp)
        {
            m_position--;
            return *this;
        }

        // Postfix
        array_view_iterator operator --(int) restrict(cpu,amp)
        {
            auto temp = *this;
            m_position--;
            return temp;
        }

        bool operator==(const array_view_iterator& rhs) const restrict(cpu,amp)
        {
            return (m_position == rhs.m_position);
        }

        bool operator!=(const array_view_iterator& other) const restrict(cpu,amp)
        {
            return !(*this == other);
        }


        reference operator*() const restrict(cpu,amp)
        {
            return *deref();
        }

        reference operator->() const restrict(cpu,amp)
        {
            return *deref();
        }

        array_view_iterator operator+(difference_type delta) const restrict(cpu,amp)
        {
            auto temp = *this;
            temp.m_position += delta;
            return temp;
        }

        array_view_iterator operator-(difference_type delta) const restrict(cpu,amp)
        {
            return this->operator+(-delta);
        }

        difference_type operator-(const array_view_iterator& other) const restrict(cpu,amp)
        {
            return (m_position - other.m_position);
        }

        friend array_view_iterator operator+(difference_type delta, const array_view_iterator&) restrict(cpu,amp);
        friend array_view_iterator operator-(difference_type delta, const array_view_iterator&) restrict(cpu,amp);

        bool operator < (const array_view_iterator& rhs) const restrict(cpu,amp)
        {
            return (m_position < rhs.m_position);
        }

        bool operator > (const array_view_iterator& rhs) const restrict(cpu,amp)
        {
            return rhs < *this;
        }

        bool operator <= (const array_view_iterator& rhs) const restrict(cpu,amp)
        {
            return (m_position <= rhs.m_position);
        }

        bool operator >= (const array_view_iterator& rhs) const restrict(cpu,amp)
        {
            return rhs <= *this;
        }

        array_view_iterator& operator+=(difference_type delta) restrict(cpu,amp)
        {
            m_position += delta;
            return *this;
        }

        array_view_iterator& operator-=(difference_type delta) restrict(cpu,amp)
        {
            return this->operator+=(-delta);
        }

        value_type& operator[](difference_type idx) const restrict(cpu,amp)
        {
            return *(*this + idx);
        }

    private:
        // TODO: Should be able to construct a const array_view_iterator from a non-const array_view
        array_view_iterator(const array_view<value_type>& arr_view, difference_type position) restrict(cpu,amp)
            : m_base_view(arr_view), 
            m_position(position)
        {
        }

        value_type * deref() const restrict(cpu,amp)
        {
            return &m_base_view[m_position];
        }

        array_view<value_type> m_base_view;
        difference_type m_position;
    };

    namespace _details
    {
        template<typename value_type>
        struct array_view_iterator_helper
        {
            static const array_view<value_type> & get_base_array_view(const array_view_iterator<value_type> & iter)
            {
                return iter.m_base_view;
            }
        };
    }

    // Friends of array_view_iterator
    template <typename value_type>
    array_view_iterator<value_type> operator+(typename array_view_iterator<value_type>::difference_type delta, const array_view_iterator<value_type>& iter) restrict(cpu,amp)
    {
        return iter+delta;
    }

    template <typename value_type>
    array_view_iterator<value_type> operator-(typename array_view_iterator<value_type>::difference_type delta, const array_view_iterator<value_type>& iter) restrict(cpu,amp)
    {
        return iter-delta;
    }

    // end of array_view_iterator
    //----------------------------------------------------------------------------

    ///////////////////////////////////////////////////////////////////////////////
    // iterator_traits
    //
    // Provides typedefs that programmers could use to define iterators. e.g.,
    //
    // iterator_traits<array_view<float,2>>::const_iterator_type myAvIter;
    //
    // In future release, this should be offered as member typedefs of the 
    // respective classes.
    //
    ///////////////////////////////////////////////////////////////////////////////
    template <typename array_type>
    class iterator_traits
    {
        iterator_traits()
        {
            static_assert(false, "This class must be specialized");
        }
    };

    template <typename value_type>
    class iterator_traits<concurrency::array_view<value_type>>
    {
    public:
        typedef array_view_iterator<value_type> iterator_type;
        typedef array_view_iterator<const value_type> const_iterator_type;
    };

    template <typename value_type>
    class iterator_traits<concurrency::array_view<const value_type>>
    {
    public:
        typedef array_view_iterator<const value_type> const_iterator_type;
    };

    //----------------------------------------------------------------------------
    // begin and end iterators for array views
    //----------------------------------------------------------------------------
    template <typename value_type>
    array_view_iterator<value_type> begin(const concurrency::array_view<value_type>& arr) restrict(cpu,amp)
    {
        return array_view_iterator<value_type>(arr, 0);
    }

    template <typename value_type>
    array_view_iterator<value_type> end(const concurrency::array_view<value_type>& arr) restrict(cpu,amp)
    {
        return array_view_iterator<value_type>(arr, arr.get_extent().size());
    }

} // amp_stl_algorithms