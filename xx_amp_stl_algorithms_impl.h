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
* This file contains the helpers classes in amp_stl_algorithms::_details namespace
*---------------------------------------------------------------------------*/

#pragma once

#include <amp.h>
#include <assert.h>

namespace amp_stl_algorithms
{
    namespace _details
    {
        using namespace concurrency;

        void amp_assert(bool cond) restrict(cpu)
        {
            assert(cond);
        }

        void amp_assert(bool cond) restrict(amp)
        {
            // TODO
        }

        template <typename array_type>
        void assert_arrays_are_same_toplevel_resource(const array_type& a1, const array_type& a2) restrict(cpu,amp)
        {
            // TODO
        }

        //////////////////////////////////////////////////////////////////////////
        // Empty array factories
        //
        // Random iterators require a default constructors. Because array_view
        // iterators must hold an array view or equivalent, there must be a way
        // to create a "default" array_view. It is possible to simply not provide
        // a default constructor, but then some algorithms may not work.
        //
        // Instead the following classes and functions make a good effort to create
        // a "real" array_view pointing to at least superficially valid data -- as 
        // long as it is not dereferenced.
        //
        // Note: only array_view<T, 1> is supported.
        /////////////////////////////////////////////////////////////////////////////

        // Specialization for array_views 
        template <typename value_type>
        class empty_array_view_factory
        {
        public:
            // On the CPU, use local static variables as stable storage
            static concurrency::array_view<value_type> create() restrict(cpu)
            {
                static value_type stable_storage;
                return concurrency::array_view<value_type>(1, &stable_storage);
            }

            // On the accelerator, use the stable storage helper.
            static concurrency::array_view<value_type> create() restrict(amp)
            {
                return concurrency::array_view<value_type>(0, nullptr);
            }
        };

        template <class value_type, int rank>
        concurrency::array_view<value_type> make_array_view(concurrency::array<value_type, rank>& arr) restrict(cpu,amp)
        {
            return arr.view_as(concurrency::extent<1>(arr.get_extent().size()));
        }

        template <class value_type, int rank>
        concurrency::array_view<const value_type> make_array_view(const concurrency::array<value_type, rank>& arr) restrict(cpu,amp)
        {
            return arr.view_as(concurrency::extent<1>(arr.get_extent().size()));
        }

    } // namespace _details

} // namespace amp_stl_algorithms
