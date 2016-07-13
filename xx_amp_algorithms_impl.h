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
* This file contains the helpers classes in amp_algorithms::_details namespace
*---------------------------------------------------------------------------*/

#pragma once

#include <amp.h>
#include <assert.h>
#include <amp_indexable_view.h>
#include <sstream>
#include <d3d11.h>
#include <d3dcsx.h>
#include <wrl\client.h>

#pragma comment(lib, "d3dcsx")

namespace amp_algorithms
{

    namespace _details
    {
        inline concurrency::accelerator_view auto_select_target()
        {
            static concurrency::accelerator_view auto_select_accelerator_view = concurrency::accelerator(concurrency::accelerator::cpu_accelerator).create_view();
            return auto_select_accelerator_view;
        }

        template <int _Rank, typename _Kernel_type>
        void parallel_for_each(const concurrency::accelerator_view &_Accl_view, const concurrency::extent<_Rank>& _Compute_domain, const _Kernel_type &_Kernel)
        {
            _Host_Scheduling_info _SchedulingInfo = { NULL };
            if (_Accl_view != _details::auto_select_target()) 
            {
                _SchedulingInfo._M_accelerator_view = details::_Get_accelerator_view_impl_ptr(_Accl_view);
            }

            details::_Parallel_for_each(&_SchedulingInfo, _Compute_domain, _Kernel);
        }

        template <int _Dim0, int _Dim1, int _Dim2, typename _Kernel_type>
        void parallel_for_each(const concurrency::accelerator_view &_Accl_view, const concurrency::tiled_extent<_Dim0, _Dim1, _Dim2>& _Compute_domain, const _Kernel_type& _Kernel)
        {
            _Host_Scheduling_info _SchedulingInfo = { NULL };
            if (_Accl_view != _details::auto_select_target()) 
            {
                _SchedulingInfo._M_accelerator_view = details::_Get_accelerator_view_impl_ptr(_Accl_view);
            }

            details::_Parallel_for_each(&_SchedulingInfo, _Compute_domain, _Kernel);
        }

        template <int _Dim0, int _Dim1, typename _Kernel_type>
        void parallel_for_each(const concurrency::accelerator_view &_Accl_view, const concurrency::tiled_extent<_Dim0, _Dim1>& _Compute_domain, const _Kernel_type& _Kernel)
        {
            _Host_Scheduling_info _SchedulingInfo = { NULL };
            if (_Accl_view != _details::auto_select_target()) 
            {
                _SchedulingInfo._M_accelerator_view = details::_Get_accelerator_view_impl_ptr(_Accl_view);
            }

            details::_Parallel_for_each(&_SchedulingInfo, _Compute_domain, _Kernel);
        }

        template <int _Dim0, typename _Kernel_type>
        void parallel_for_each(const concurrency::accelerator_view &_Accl_view, const concurrency::tiled_extent<_Dim0>& _Compute_domain, const _Kernel_type& _Kernel)
        {
            _Host_Scheduling_info _SchedulingInfo = { NULL };
            if (_Accl_view != _details::auto_select_target()) 
            {
                _SchedulingInfo._M_accelerator_view = details::_Get_accelerator_view_impl_ptr(_Accl_view);
            }

            details::_Parallel_for_each(&_SchedulingInfo, _Compute_domain, _Kernel);
        }

        // Reduction implementation

        // This function performs an in-place reduction through co-operating threads within a tile.
        // The input data is in the parameter "mem" and is reduced in-place modifying its existing contents
        // The output (reduced result) is contained in "mem[0]" at the end of this function
        // The parameter "partial_data_length" is used to indicate if the size of data in "mem" to be
        // reduced is same as the tile size and if not what is the length of valid data in "mem".
        template <typename T, unsigned int tile_size, typename functor>
        void tile_local_reduction(T* const mem, concurrency::tiled_index<tile_size> tid, const functor& op, int partial_data_length) restrict(amp)
        {
            // local index
            int local = tid.local[0];

            if (partial_data_length < tile_size) 
            {
                // unrolled for performance
                if (partial_data_length >  512) { if (local < (partial_data_length - 512)) { mem[0] = op(mem[0], mem[512]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >  256) { if (local < (partial_data_length - 256)) { mem[0] = op(mem[0], mem[256]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >  128) { if (local < (partial_data_length - 128)) { mem[0] = op(mem[0], mem[128]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >   64) { if (local < (partial_data_length -  64)) { mem[0] = op(mem[0], mem[ 64]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >   32) { if (local < (partial_data_length -  32)) { mem[0] = op(mem[0], mem[ 32]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >   16) { if (local < (partial_data_length -  16)) { mem[0] = op(mem[0], mem[ 16]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >    8) { if (local < (partial_data_length -   8)) { mem[0] = op(mem[0], mem[  8]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >    4) { if (local < (partial_data_length -   4)) { mem[0] = op(mem[0], mem[  4]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >    2) { if (local < (partial_data_length -   2)) { mem[0] = op(mem[0], mem[  2]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (partial_data_length >    1) { if (local < (partial_data_length -   1)) { mem[0] = op(mem[0], mem[  1]); } tid.barrier.wait_with_tile_static_memory_fence(); }
            }
            else
            {
                // unrolled for performance
                if (tile_size >= 1024) { if (local < 512) { mem[0] = op(mem[0], mem[512]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=  512) { if (local < 256) { mem[0] = op(mem[0], mem[256]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=  256) { if (local < 128) { mem[0] = op(mem[0], mem[128]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=  128) { if (local <  64) { mem[0] = op(mem[0], mem[ 64]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=   64) { if (local <  32) { mem[0] = op(mem[0], mem[ 32]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=   32) { if (local <  16) { mem[0] = op(mem[0], mem[ 16]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=   16) { if (local <   8) { mem[0] = op(mem[0], mem[  8]); } tid.barrier.wait_with_tile_static_memory_fence(); }   
                if (tile_size >=    8) { if (local <   4) { mem[0] = op(mem[0], mem[  4]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=    4) { if (local <   2) { mem[0] = op(mem[0], mem[  2]); } tid.barrier.wait_with_tile_static_memory_fence(); }
                if (tile_size >=    2) { if (local <   1) { mem[0] = op(mem[0], mem[  1]); } tid.barrier.wait_with_tile_static_memory_fence(); }
            }
        }

        // Generic reduction of a 1D indexable view with a reduction binary functor
        template<unsigned int tile_size, 
            unsigned int max_tiles,
            typename InputIndexableView,
            typename BinaryFunction>
            typename std::result_of<BinaryFunction(const typename indexable_view_traits<InputIndexableView>::value_type&, const typename indexable_view_traits<InputIndexableView>::value_type&)>::type
            reduce(const concurrency::accelerator_view &accl_view, const InputIndexableView &input_view, const BinaryFunction &binary_op)
        {
            // The input view must be of rank 1
            static_assert(indexable_view_traits<InputIndexableView>::rank == 1, "The input indexable view must be of rank 1");
            typedef typename std::result_of<BinaryFunction(const typename indexable_view_traits<InputIndexableView>::value_type&, const typename indexable_view_traits<InputIndexableView>::value_type&)>::type result_type;

            // runtime sizes
            int n = input_view.extent.size();
            unsigned int tile_count = (n+tile_size-1) / tile_size;
            tile_count = std::min(tile_count, max_tiles);   

            // simultaneous live threads
            const unsigned int thread_count = tile_count * tile_size;

            // global buffer (return type)
            concurrency::array_view<result_type> global_buffer_view(concurrency::array<result_type>(tile_count, concurrency::accelerator(concurrency::accelerator::cpu_accelerator).default_view, accl_view));

            // configuration
            concurrency::extent<1> extent(thread_count);

            _details::parallel_for_each (
                accl_view,
                extent.tile<tile_size>(),
                [=] (concurrency::tiled_index<tile_size> tid) restrict(amp)
            {
                // shared tile buffer
                tile_static result_type local_buffer[tile_size];

                int idx = tid.global[0];

                // this threads's shared memory pointer
                result_type& smem = local_buffer[ tid.local[0] ];

                // this variable is used to test if we are on the edge of data within tile
                int partial_data_length = n - tid.tile[0] * tile_size;

                // initialize local buffer
                smem = input_view[concurrency::index<1>(idx)];
                // next chunk
                idx += thread_count;

                // fold data into local buffer
                while (idx < n)
                {
                    // reduction of smem and X[idx] with results stored in smem
                    smem = binary_op(smem, input_view[concurrency::index<1>(idx)]);

                    // next chunk
                    idx += thread_count;
                }

                // synchronize
                tid.barrier.wait_with_tile_static_memory_fence();

                // reduce all values in this tile
                _details::tile_local_reduction(&smem, tid, binary_op, partial_data_length);

                // only 1 thread per tile does the inter tile communication
                if (tid.local[0] == 0)
                {
                    // write to global buffer in this tiles
                    global_buffer_view[ tid.tile[0] ] = smem;
                }
            });

            // 2nd pass reduction
            result_type *pGlobalBufferViewData = global_buffer_view.data();
            result_type retVal = pGlobalBufferViewData[0];
            for (unsigned int i = 1; i < tile_count; ++i) {
                retVal = binary_op(retVal, pGlobalBufferViewData[i]);
            }

            return retVal;
        }

        inline void _check_hresult(HRESULT _hr, std::string _exception_msg = "")
        {
            if (FAILED(_hr))
            { 
                std::stringstream _out;
                _out << _exception_msg << " 0x" << std::hex << _hr << std::endl; 
                throw std::runtime_error(_out.str());
            }
        }

        inline Microsoft::WRL::ComPtr<ID3D11Device> _get_d3d11_device_ptr(const concurrency::accelerator_view &av)
        {
            IUnknown *u = direct3d::get_device(av);
            Microsoft::WRL::ComPtr<ID3D11Device> dev_ptr;
            auto hr = u->QueryInterface(__uuidof(ID3D11Device), reinterpret_cast<void**>(dev_ptr.GetAddressOf()));
            u->Release();
            _check_hresult(hr);
            return dev_ptr;
        }

        template<typename T, unsigned int Rank>
        inline Microsoft::WRL::ComPtr<ID3D11Buffer> _get_d3d11_buffer_ptr(const array<T, Rank> &a)
        {
            IUnknown *u = direct3d::get_buffer(a);
            Microsoft::WRL::ComPtr<ID3D11Buffer> buf_ptr;
            auto hr = u->QueryInterface(__uuidof(ID3D11Buffer), reinterpret_cast<void**>(buf_ptr.GetAddressOf()));
            u->Release();
            _check_hresult(hr);
            return buf_ptr;
        }

        inline Microsoft::WRL::ComPtr<ID3D11UnorderedAccessView> _create_d3d11_uav(Microsoft::WRL::ComPtr<ID3D11Device> &device, Microsoft::WRL::ComPtr<ID3D11Buffer> &pSrcBuff, DXGI_FORMAT view_format)
        {
            D3D11_UNORDERED_ACCESS_VIEW_DESC desc;
            ZeroMemory(&desc, sizeof(desc));
            desc.Format = view_format;
            desc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;

            D3D11_BUFFER_DESC descBuff;
            pSrcBuff->GetDesc(&descBuff);
            desc.Buffer.FirstElement = 0;
            desc.Buffer.NumElements  = descBuff.ByteWidth / sizeof(int);

            Microsoft::WRL::ComPtr<ID3D11UnorderedAccessView> pUAView;
            _check_hresult(device->CreateUnorderedAccessView(reinterpret_cast<ID3D11Resource*>(pSrcBuff.Get()), &desc, pUAView.GetAddressOf()), "Failed to create view");

            return pUAView;
        }

        template <typename _Type>
        struct _dx_scan_type_helper
        {
            static const bool is_type_supported = false;
        };

        template <>
        struct _dx_scan_type_helper<int>
        {
            static const bool is_type_supported = true;
            static const D3DX11_SCAN_DATA_TYPE dx_scan_type = D3DX11_SCAN_DATA_TYPE_INT;
            static const DXGI_FORMAT dx_view_type = DXGI_FORMAT_R32_SINT;
        };

        template <>
        struct _dx_scan_type_helper<unsigned int>
        {
            // Note: Despite what the MSDN says D3DCSX does not support uint, 
            // we can partially support it by treating it as int.
            static const bool is_type_supported = true;
            static const D3DX11_SCAN_DATA_TYPE dx_scan_type = D3DX11_SCAN_DATA_TYPE_INT;
            static const DXGI_FORMAT dx_view_type = DXGI_FORMAT_R32_SINT;
        };

        template <>
        struct _dx_scan_type_helper<float>
        {
            static const bool is_type_supported = true;
            static const D3DX11_SCAN_DATA_TYPE dx_scan_type = D3DX11_SCAN_DATA_TYPE_FLOAT;
            static const DXGI_FORMAT dx_view_type = DXGI_FORMAT_R32_FLOAT;
        };

        struct _dx_state_cleaner 
        {
            _dx_state_cleaner(Microsoft::WRL::ComPtr<ID3D11DeviceContext> &context) : m_immediate_context(context)
            {
                memset(m_uavs, 0, D3D11_PS_CS_UAV_REGISTER_COUNT * sizeof(ID3D11UnorderedAccessView*));
                m_immediate_context->CSGetUnorderedAccessViews(0, D3D11_PS_CS_UAV_REGISTER_COUNT, m_uavs);
            }

            ~_dx_state_cleaner()
            {
                m_immediate_context->CSSetUnorderedAccessViews(0, D3D11_PS_CS_UAV_REGISTER_COUNT, m_uavs, nullptr);
                for(unsigned int i = 0; i < D3D11_PS_CS_UAV_REGISTER_COUNT; ++i)
                {
                    if (m_uavs[i] != nullptr)
                    {
                        m_uavs[i]->Release();
                    }
                }
            }

        private:
            ID3D11UnorderedAccessView *m_uavs[D3D11_PS_CS_UAV_REGISTER_COUNT];
            Microsoft::WRL::ComPtr<ID3D11DeviceContext> m_immediate_context;
        };

        // Declaration of Scan helper that converts binary functions from C++ AMP library to DirectX scan operation codes
        template <typename BinaryFunction>
        struct _dx_scan_op_helper;

    } // namespace amp_algorithms::_details

} // namespace amp_algorithms
