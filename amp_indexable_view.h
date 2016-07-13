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
* C++ AMP algorithms library.
*
* This file contains the C++ AMP indexable_view concept definition and 
* definition of other indexable_view types
*---------------------------------------------------------------------------*/

#pragma once

namespace amp_algorithms
{
    using namespace concurrency;

    // The type indexable_view_traits defines the traits of any type that
    // conforms to the indexable_view concept. An indexable view should:
    // a. Have a static member 'rank' of type int, which indicates the rank
    //    of the view
    // b. A typedef 'value_type' to indicate the element type of the view.
    // c. Must have an extent property denoting the dense zero based extents
    //    of the view which returns an object convertible to 'const extent<rank>&'
    // d. Must have an 'operator[]' overload with the following properties:
    //      i) Have the 'restrict(cpu, amp)' restriction qualifier
    //     ii) Have a single parameter which 'const index<rank>&' is convertible to
    //    iii) The return type of the 'operator()(const index<rank> &idx)'
    //         must be convertible to 'const value_type&'
    // e. Objects of the type view must be capturable by value in a p_f_e kernel

    template <typename view>
    struct indexable_view_traits
    {
        typedef typename view::value_type value_type;
        static const int rank = view::rank;
        static const bool is_writable = std::is_convertible<typename decltype(std::declval<view>()[std::declval<concurrency::index<rank>>()]), value_type&>::value;

        // TODO: Other traits such as whether the storage is dense in LSD or in which dimension is the storage contiguous

        static_assert(std::is_convertible<typename decltype(std::declval<view>().extent), const concurrency::extent<rank>&>::value, "Not a valid indexable_view. Should have a member 'extent'");

        // TODO: This needs to ensure that the operator() has both the cpu and amp restriction qualifiers. The current check 
        // checks that it has the required parameter and return types and just the "cpu" restriction qualifier.
        static_assert(std::is_convertible<typename decltype(std::declval<view>()[std::declval<concurrency::index<rank>>()]), const value_type&>::value, "Not a valid indexable_view.");

        // TODO: More checks to statically verify the view template parameter and issue appropriate errors.
    };


    // The class functor_view is an example of an indexable_view
    // It can be constructed from an extent<rank> and a functor object
    template <typename Functor, int Rank = 1>
    class functor_view
    {
    public:
        static const int rank = Rank;
        typedef decltype((*((Functor*)NULL))(*((concurrency::index<rank>*)NULL))) functor_return_type;
        typedef typename std::remove_const<typename std::remove_reference<functor_return_type>::type>::type value_type;

        functor_view(const concurrency::extent<rank> &ext, const Functor &functor)
            : _M_extent(ext), _M_functor(functor)
        {
        }

        functor_return_type operator[](const concurrency::index<rank> &idx) const restrict(cpu, amp)
        {
            return _M_functor(idx);
        }

        __declspec(property(get=get_extent)) Concurrency::extent<rank> extent;
        Concurrency::extent<rank> get_extent() const restrict(cpu,amp)
        {
            return _M_extent;
        }

    private:
        concurrency::extent<rank> _M_extent;
        Functor _M_functor;
    };

    // This is a template method for constructing a functor_view object which can be used with the "auto" keyword
    // obviating the need to explicitly specify the functor type which is not very straightforward to obtain for lambdas
    template <typename Functor, int Rank>
    functor_view<Functor, Rank> make_indexable_view(const concurrency::extent<Rank> &ext, const Functor &functor)
    {
        return functor_view<Functor, Rank>(ext, functor);
    }

} // namespace amp_algorithms
