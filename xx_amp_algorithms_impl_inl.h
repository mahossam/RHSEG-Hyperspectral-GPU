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
* This file contains the helpers templates on which amp_algorithms.h depends.
*---------------------------------------------------------------------------*/

#pragma once

namespace amp_algorithms
{

    namespace _details
    {

        // Scan helper that converts binary functions from C++ AMP library to DirectX scan operation codes
        template <typename BinaryFunction>
        struct _dx_scan_op_helper
        {
            static const bool is_op_supported = false;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::sum<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_ADD;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::max<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_MAX;
        };

        // Further specialize amp_algorithms::max for uint and mark as not supported.
        template <>
        struct _dx_scan_op_helper<amp_algorithms::max<unsigned int>>
        {
            // max is not supported for uint, as our implementation is based on int,
            // this will return incorrect results for values that are greater than numeric_limits<int>::max().
            static const bool is_op_supported = false;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::min<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_MIN;
        };

        // Further specialize amp_algorithms::min for uint and mark as not supported.
        template <>
        struct _dx_scan_op_helper<amp_algorithms::min<unsigned int>>
        {
            // min is not supported for uint, as our implementation is based on int,
            // this will return incorrect results for values that are greater than numeric_limits<int>::max().
            static const bool is_op_supported = false;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::mul<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_MUL;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::bit_and<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_AND;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::bit_or<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_OR;
        };

        template <typename T>
        struct _dx_scan_op_helper<amp_algorithms::bit_xor<T>>
        {
            static const bool is_op_supported = true;
            static const D3DX11_SCAN_OPCODE dx_op_type = D3DX11_SCAN_OPCODE_XOR;
        };

    } // namespace amp_algorithms::_details

} // namespace amp_algorithms