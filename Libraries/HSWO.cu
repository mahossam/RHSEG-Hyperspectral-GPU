#include <map>
#include <string>

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/extrema.h>
#include <thrust/remove.h>

#include <vector_types.h>
#include <cuda.h>
#if CUDA_42
#include <cutil_inline.h>
#endif
#include <cfloat>
#include <time.h>

#include "HSWO.h"
#include "HSWODevice.h"

#include <sstream>
#include "Libraries/Logger.h"

#ifdef LOG_ENABLED
	#define stream_to_log(theStream) {stream_to_log_inner(theStream)}
	#define stream_to_log_sameLine(theStream) {stream_to_log_inner_sameLine(theStream)}
#else
	#define stream_to_log(theStream)
	#define stream_to_log_sameLine(theStream)
#endif

#define stream_to_log_inner(theStream) { using namespace HyperSpectralToolbox; \
					auto sp_ofstream = Logger::openStream(); \
					(*sp_ofstream) << theStream << "\n";	\
					sp_ofstream->close(); }					
			
#define stream_to_log_inner_sameLine(theStream) { using namespace HyperSpectralToolbox; \
					auto sp_ofstream = Logger::openStream(); \
					(*sp_ofstream) << theStream;	\
					sp_ofstream->close(); }	

__global__ void atmoic_exp(int* mutex)
{
	int indexi = blockIdx.x * blockDim.x + threadIdx.x;
	int indexj = blockIdx.y * blockDim.y + threadIdx.y;

	for (int i=0; i < 32; i++) {
		for (int j=0; j < 32; j++) {
			if ((indexi % 32 == i) && (indexj % 32 == j)) {

				while(atomicCAS(mutex,0, 1) == 1)
				{
					printf("Thread %d, %d - cas = %d\n", indexi, indexj, 1);
					//printf("Thread %d - cas = %d\n", indexi, 1);
				}

				printf("Thread %d, %d - inside section \n", indexi, indexj);
				//printf("Thread %d - inside section\n", indexi);

				atomicExch(mutex, 0 );	

			}		
		}
	}
	
}
namespace HyperSpectralToolbox
{
	void run_expr()
	{	
		//mutex initialization (for solution2)
		//int N = 64;
		printf( "expr atomic =============\n" );
		int* dev_mutex2;
		cudaMalloc((void**)&dev_mutex2, sizeof(int));		
		cudaMemset(dev_mutex2, 0, sizeof(int));

		// Kernel invocation 
		dim3 numBlocks2(2,2); //
		//dim3 numBlocks2(3,1); //

		dim3 threadsPerBlock2(8, 4);  //
		//dim3 threadsPerBlock2(32, 1);  //
		atmoic_exp<<<numBlocks2 , threadsPerBlock2>>>(dev_mutex2);

		cudaThreadSynchronize();
		printf("kernel finished\n");
		cudaFree( dev_mutex2 ); //for solution 2
	}
}

//#define MY_FLT_MAX         3.402823466e+38F        /* max value */
//#define MY_FLT_MAX         (3.402823466*1000000.0f)        /* max value */
#define MY_FLT_MAX         (FLT_MAX)        /* max value */

__global__ void kernel_compute_dissim(int* keysFilter ,int* adj, int nRegions, const int nBands
	, float* regionsSums, int* regionsPixelsCount
	, float* regionsAdjMinimums
	, int* regionsAdjMinimumsLabels
	, int* needsAdjRecomputation) 
{ 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	//printf("[Thread %d]: blockIdx.x = %d, blockDim.x = %d, threadIdx.x = %d: Hi ! \n", index, blockIdx.x , blockDim.x , threadIdx.x);

	if(index < nRegions)
	{
		int regionID = keysFilter[index]-1;
		if (regionID == -2) return;

		float SHregionsAdjMinimums;
		SHregionsAdjMinimums = MY_FLT_MAX;
		int SHregionsAdjMinimumsLabels;
		SHregionsAdjMinimumsLabels = -1;

		//printf(" - [Thread %d] region %d\n", index, regionID+1 );

		// ****** Optimization ****************
		if(needsAdjRecomputation[regionID] == 1)
		//if(true)
		{
			///////////////////////////////////////
			int offset_regionID = regionID*MAX_REGION_ADJ;	
			float r1nPixels = regionsPixelsCount[regionID];

			//printf(" - [Thread %d] region %d pixel count = %f \n", index, regionID+1 , r1nPixels);

			for(int a=0; a<MAX_REGION_ADJ; a++)
			{				
				int regionIDadj = adj[offset_regionID+a]-1;				

				if (regionIDadj <= -1) 
				{
					break;
				}
				//if(regionID > regionIDadj) continue;  //optimization			
				float sum_eculidean = 0.0f;
				float r2nPixels = regionsPixelsCount[regionIDadj];
				
				//printf(" - [Thread %d] ADJ region %d pixel count = %f \n", index, regionIDadj+1, r2nPixels);

				float temp;

				int band = 0;

				for (; band < nBands; band++)
				{
					temp = regionsSums[regionID*nBands+band]/(r1nPixels) - regionsSums[regionIDadj*nBands+band]/(r2nPixels);
					//					printf("	--- [Thread %d] region %d sum at band %d = %f \n", index, regionID, band, regionsSums[regionID*nBands+band]);
					//					printf("	--- [Thread %d] region %d sum at band %d = %f \n", index, regionIDadj, band, regionsSums[regionIDadj*nBands+band]);
					sum_eculidean += temp * temp;
				}

				//output[offset_regionID+a] = sqrtf(sum_eculidean) / float(nBands);
				//output[offset_regionID+a] = sqrtf(sum_eculidean*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));
				
				temp = sqrtf(sum_eculidean*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));
				//printf("temp = %f\n", temp);

				if(temp < SHregionsAdjMinimums)  
				{
					SHregionsAdjMinimums = temp;
					SHregionsAdjMinimumsLabels = regionIDadj+1;
					//printf(" - [Thread %d] region %d min distance to neighbor %d = %f \n", index, regionID+1, regionsAdjMinimumsLabels[regionID], regionsAdjMinimums[regionID]);
				}				

				//if(temp < regionsAdjMinimums[regionID])  
				//{
					//regionsAdjMinimums[regionID] = temp;
					//regionsAdjMinimumsLabels[regionID] = regionIDadj+1;
					//printf(" - [Thread %d] region %d min distance to neighbor %d = %f \n", index, regionID+1, regionsAdjMinimumsLabels[regionID], regionsAdjMinimums[regionID]);
				//}				

			}
			//}
			regionsAdjMinimums[regionID] = SHregionsAdjMinimums;
			regionsAdjMinimumsLabels[regionID] = SHregionsAdjMinimumsLabels;
			needsAdjRecomputation[regionID] = 0;
		}
		else
		{
			//no recomputation needed: do nothing
		}
	}
} 

#if ONE_DIMENSIONAL_SPEC_KERNEL
__global__ void kernel_compute_regions_dissim(const int* keysFilter ,const int* adj, int nRegions, const int nBands
	, const float* regionsSums, const int* regionsPixelsCount
	, float* regionsMinimums
	, int* regionsMinimumsLabels
	, int* needsRecomputation) 
{ 
	int index1 = blockIdx.x * blockDim.x + threadIdx.x;	
	int canCompute = 1;

	//if(index1 < nRegions)
	//{
		int regionID1 = keysFilter[index1]-1;
		if (regionID1 == -2) 
		{
			canCompute = 0;
		}

		//__shared__ float r1nPixels;		
		float r1nPixels = regionsPixelsCount[regionID1];
		float SHregionsMinimums = MY_FLT_MAX;
		int SHregionsMinimumsLabels = -1;		

		// ****** Optimization ****************
		//if(true)			
		//if(needsRecomputation[regionID1] == 1)			
		if(canCompute)
		{
			for(int i=0; i < nRegions; ++i)
			{
				int regionID2 = keysFilter[i]-1;		
				int canComputePair = 1;
				if (regionID2 == -2)
				{
					canComputePair = 0;
				}

				if (regionID1 == regionID2)
				{
					canComputePair = 0;
				}

				//__syncthreads();  

				if(canComputePair)
				{
					int offset_regionID = regionID1*MAX_REGION_ADJ;	
					int foundAdjacent = 0;
					//int runFlag = 1;
					//__syncthreads(); 
					for(int a=0; a<MAX_REGION_ADJ; ++a)
					{
						//if(runFlag)
						//{
							int regionIDadj = adj[offset_regionID+a]-1;				

							if (regionIDadj <= -1) 
							{
								break;
								//runFlag = 0;
							}
							if(regionIDadj == regionID2) {
								foundAdjacent = 1;	//after adding this line, the kernel slowed down to 489 seconds!
								// seems pretty nasty thread divergence happned !!
								// Fixed: by adding the __syncthreads(); line after the loop
								break;	
								//runFlag = 0;
							}
						//}
					}

					__syncthreads();  //The fix to thread divergence problem

					if(foundAdjacent == 0)
					{
						float r2nPixels = regionsPixelsCount[regionID2];

						//if(regionID > regionIDadj) continue;  //optimization			
						float sum_dissim = 0.0f;
						float temp;
						for (int band = 0; band < nBands; band++)
						{
							temp = regionsSums[regionID1*nBands+band]/(r1nPixels) - regionsSums[regionID2*nBands+band]/(r2nPixels);
							sum_dissim += temp * temp;
						}

						temp = sqrtf(sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));

						if(temp < SHregionsMinimums) 
						{
							SHregionsMinimums = temp;
							SHregionsMinimumsLabels = regionID2+1;
						}
					}
					
				}				
				//if(temp < regionsMinimums[regionID1]) 
				//{
				//	regionsMinimums[regionID1] = temp;
				//	regionsMinimumsLabels[regionID1] = regionID2+1;
				//}				
			}	
			//__syncthreads();
			regionsMinimums[regionID1] = SHregionsMinimums;
			regionsMinimumsLabels[regionID1] = SHregionsMinimumsLabels;
			needsRecomputation[regionID1] = 0;
		}		
	//}
} 
#endif

__global__ void printTestArr(int* testArr)
{		
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	//for(int t=0; t< 1024;++t)
	printf("[%d] = %d\n", index,testArr[index]);
}

#define BLOCK_BANDS 16
#define BLOCK_WIDTH 16

__global__ void kernel_compute_regions_dissim_2dim(/*int *testArr,*/ const int* keysFilter ,const int* adj, int nRegions, const int nBands
	, const float* regionsSums, const int* regionsPixelsCount
	/*
	, float* regionsMinimums
	, int* regionsMinimumsLabels
	*/
	, int* needsRecomputation	
	, float* _temp_all_dissims) 
{ 
	int region1Index = blockIdx.y * blockDim.y + threadIdx.y;
	int region2Index = blockIdx.x * blockDim.x + threadIdx.x;
	

	//__shared__ int shared_regions1_pixels_count[BLOCK_WIDTH];
	//__shared__ int shared_regions2_pixels_count[BLOCK_WIDTH];

	__shared__ float shared_regions1_sums[BLOCK_WIDTH][BLOCK_BANDS];
	__shared__ float shared_regions2_sums[BLOCK_WIDTH][BLOCK_BANDS];
	if(threadIdx.x == 0) 
	{
		//shared_regions1_pixels_count[threadIdx.y] = regionsPixelsCount[region1Index];
		for(int b=0; b<BLOCK_BANDS; ++b)
		{
			shared_regions1_sums[threadIdx.y][b] = regionsSums[region1Index*nBands+b];
		}
	}

	if(threadIdx.y == 0) 
	{
		//shared_regions2_pixels_count[threadIdx.x] = regionsPixelsCount[region2Index];
		for(int b=0; b<BLOCK_BANDS; ++b)
		{
			shared_regions2_sums[threadIdx.x][b] = regionsSums[region2Index*nBands+b];
		}
	}
	
	__syncthreads();


	float sum_dissim = 0.0f;
	float temp = MY_FLT_MAX;
	int can_compute = 1;

	int region1ID = keysFilter[region1Index]-1;				
	if (region1ID == -2) 
	{
		can_compute = 0;
	}

#if SPECTRAL_DYNAMIC_PROGRAMMING
	if( region1ID != -2)
	{
		if( needsRecomputation[region1ID] == 0 )
		{
			can_compute = 0;
		}
	}
#endif

	int region2ID = keysFilter[region2Index]-1;				
	if (region2ID == -2) 
	{
		can_compute = 0;
	}

	if (region1ID == region2ID) 
	{
		can_compute = 0;
	}	

	//__syncthreads();

	if(can_compute)
	{
		float r1nPixels = regionsPixelsCount[region1ID];
		//float r1nPixels = shared_regions1_pixels_count[threadIdx.y];

		// ****** Optimization ****************			
		//if(needsRecomputation[regionID1] == 1)			
		int offset_regionID = region1ID*MAX_REGION_ADJ;	
		int foundAdjacent = 0;
		for(int a=0; a<MAX_REGION_ADJ; ++a)
		{				
			int regionIDadj = adj[offset_regionID+a]-1;				

			if (regionIDadj <= -1) 
			{
				break;
			}
			if(regionIDadj == region2ID) {						
				foundAdjacent = 1;
				break;
			}
		}

		__syncthreads();

		if(foundAdjacent == 0)
		{
			float r2nPixels = regionsPixelsCount[region2ID];
			//float r2nPixels = shared_regions2_pixels_count[threadIdx.x];

			//if(region1ID > region2ID) return;  //optimization			
			int band = 0;
			for (band = 0; band < BLOCK_BANDS; ++band)
			{
				temp = (shared_regions1_sums[threadIdx.y][band])/(r1nPixels) - (shared_regions2_sums[threadIdx.x][band])/(r2nPixels);
				sum_dissim += temp * temp;
			}

			for (band = band; band < nBands; ++band)
			{
				temp = regionsSums[region1ID*nBands+band]/(r1nPixels) - regionsSums[region2ID*nBands+band]/(r2nPixels);
				sum_dissim += temp * temp;
			}


			temp = sqrtf(sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));

			//if(temp < LocalregionsMinimums) 
			{
				//LocalregionsMinimums = temp;
				//LocalregionsMinimumsLabels = region2ID+1;
				//local_minimums_store[double_idx.local[0]][double_idx.local[1]] = temp;
				//local_minimums_store_labels[double_idx.local[0]][double_idx.local[1]] = region2ID+1;

			}			
		}

		_temp_all_dissims[region1ID*nRegions + region2ID] = temp;	//if region1, region2 pair is to be computed but adjacent to each other, temp will equal MY_FLT_MAX
																	//and the old value in _temp_all_dissims will be rest. (this is a must for correctness)
																	//This solves the changing adjacents problem with _temp_all_dissims array
		//if(temp < regionsMinimums[regionID1]) 
		//{
		//	regionsMinimums[regionID1] = temp;
		//	regionsMinimumsLabels[regionID1] = regionID2+1;
		//}				
		//}			


		//__syncthreads(); //adding this syncthreads increased the time by 4 seconds, so I removed it
		//race (global minimums update)

		////if(region1ID >= 0)
		////{
		//	//mutex		(for CUDA, must serialize the wraps (using these loops and ifs) to avoid deadlock caused by wrap scheduler when he detects path divergence)
		//	//int indexi = blockIdx.x * blockDim.x + threadIdx.x;
		//	//int indexj = blockIdx.y * blockDim.y + threadIdx.y;
		//	
		//	for (int i=0; i < 32; ++i) {
		//		for (int j=0; j < 32; ++j) {
		//			//if ((indexi % 32 == i) && (indexj % 32 == j)) {
		//			if ((region1Index % 32 == i) && (region2Index % 32 == j)) {

		//				//while(atomicCAS(mutex,0, 1) == 1)
		//				//{
		//				//};

		//				///printf("Thread %d, %d - inside section \n", indexi, indexj);
		//				//critical section
		//				//__threadfence();
		//				//testing if wrap threads are serialized or not
		//				//testArr[region1ID] = testArr[region1ID] + 1;

		//				float glbl = regionsMinimums[region1ID] ;

		//				//if(local_minimums_final[double_idx.local[0]] < glbl)
		//				if(temp < glbl)
		//				{
		//					regionsMinimums[region1ID] = temp;
		//					regionsMinimumsLabels[region1ID] = region2ID+1;
		//					//needsRecomputation[regionID1] = 0;
		//				}

		//				//end critical section		
		//				//atomicExch(mutex, 0 );
		//			}
		//		}
		//	}
		////}
	}
	

} 

__global__ void update_all_spectral_minimums(const int maxRegions, const float* _temp_all_dissims, float* regionsMinimums, int* regionsMinimumsLabels)
{
	int regionIndex = blockIdx.x * blockDim.x + threadIdx.x;
	
	for(int i=0; i<maxRegions; ++i)
	{
		float val = _temp_all_dissims[regionIndex*maxRegions+i];
		if(val < regionsMinimums[regionIndex])  
		{
			regionsMinimums[regionIndex] = val;
			regionsMinimumsLabels[regionIndex] = i+1;			
		}
	}

}

__global__ void recompute_others_best_region_dissimilarity_to(int theNewMergedRegionLabel, const int nMaxRegions, int* keysFilter
	, int* needsRecomputation, int* regionsMinimumsLabels, float* regionsMinimums, const int nBands
	, float* regionsSums, int* regionsPixelsCount) 
{ 
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int canCompute = 1;

	if(index < nMaxRegions)
	{
		int regionID = keysFilter[index]-1;
		int theNewMergedRegionID = theNewMergedRegionLabel-1;
		if(regionID == -2) 
		{
			canCompute = 0;
		}
		if(regionID == (theNewMergedRegionID))
		{
			canCompute = 0;
		}

		__syncthreads();
		if(canCompute == 1)
		{
			if(needsRecomputation[regionID] == 0) //compute only for the "no need for computaion" regions
			{
				float temp; float sum_distance=0.0;

				float r_nPixels = regionsPixelsCount[regionID];
				float theNewMergedRegionNPixels = regionsPixelsCount[theNewMergedRegionID];

				for (int band = 0; band < nBands; band++)
				{
					temp = regionsSums[regionID*nBands+band]/(r_nPixels) - regionsSums[theNewMergedRegionID*nBands+band]/(theNewMergedRegionNPixels);
					sum_distance += temp * temp;
				}

				temp = sqrtf(sum_distance*((r_nPixels*theNewMergedRegionNPixels)/(r_nPixels+theNewMergedRegionNPixels)));

				if(temp < regionsMinimums[regionID])  
				{
					regionsMinimums[regionID] = temp;
					regionsMinimumsLabels[regionID] = theNewMergedRegionLabel;			
				}
			}
		}
	}
}


//__global__ void kernel_compute_regions_dissim(int* keysFilter ,int* adj, int nRegions, const int nBands
//	, float* regionsSums, float* output, int* regionsPixelsCount
//	, float* regionsMinimums
//	, int* regionsMinimumsLabels
//	, int* needsRecomputation
//	, float lastMinDissim
//	, float spclustWeight) 
//{ 
//	int index1 = blockIdx.x * blockDim.x + threadIdx.x;
//	int index2 = blockIdx.y * blockDim.y + threadIdx.y;
//	
//	if(index1 < nRegions && index2 < nRegions)
//	{
//		int regionID1 = keysFilter[index1]-1;
//		int regionID2 = keysFilter[index2]-1;
//		if (regionID1 == -2) return;
//		if (regionID2 == -2) return;
//		if (regionID1 == regionID2) return;
//
//		int offset_regionID = regionID1*MAX_REGION_ADJ;	
//		for(int a=0; a<MAX_REGION_ADJ; a++)
//		{				
//				int regionIDadj = adj[offset_regionID+a]-1;				
//
//				if (regionIDadj <= -1) 
//				{
//					break;
//				}
//				if(regionIDadj == regionID2) return;
//		}
//		//printf(" - [Thread %d] region %d\n", index, regionID+1 );
//
//		// ****** Optimization ****************
//		//if(needsAdjRecomputation[regionID] == 1)
//		if(true)
//		{
//			///////////////////////////////////////
//			float r1nPixels = regionsPixelsCount[regionID1];
//
//			//printf(" - [Thread %d] region %d pixel count = %f \n", index, regionID+1 , r1nPixels);
//
//			//if(regionID > regionIDadj) continue;  //optimization			
//			float sum_dissim = 0.0f;
//			float r2nPixels = regionsPixelsCount[regionID2];
//
//			//printf(" - [Thread %d] ADJ region %d pixel count = %f \n", index, regionIDadj+1, r2nPixels);
//
//			float temp; int band = 0;
//
//			for (; band < nBands; band++)
//			{
//				temp = regionsSums[regionID1*nBands+band]/(r1nPixels) - regionsSums[regionID2*nBands+band]/(r2nPixels);
//				sum_dissim += temp * temp;
//			}
//
//			temp = sqrtf(sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));
//			//printf("temp = %f\n", temp);
//
//			if(temp < regionsMinimums[regionID1] && temp < (lastMinDissim*spclustWeight) ) 
//			{
//				regionsMinimums[regionID1] = temp;
//				regionsMinimumsLabels[regionID1] = regionID2+1;
//				//printf(" - [Thread %d] region %d min distance to neighbor %d = %f \n", index, regionID+1, regionsAdjMinimumsLabels[regionID], regionsAdjMinimums[regionID]);
//			}				
//
//
//			needsRecomputation[regionID1] = 0;
//		}
//		else
//		{
//			//no recomputation needed: do nothing
//		}
//	}
//} 

__global__ void reset_others_best_region_computation_flags_from(int region1Label, int region2Label, const int nMaxRegions, int* keysFilter
	, int* needsRecomputation, int* regionsMinimumsLabels, float* regionsMinimums)  
{ 
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index < nMaxRegions)
	{
		int regionID = keysFilter[index]-1;
		if(regionID == -2) return;
		if(regionsMinimumsLabels[regionID] == region1Label || regionsMinimumsLabels[regionID] == region2Label) 
		{	
			needsRecomputation[regionID] = 1;
			regionsMinimumsLabels[regionID] = -1;
			regionsMinimums[regionID] = MY_FLT_MAX;			
		}
	}
}

__global__ void kernel_fill_empty_dissims(const int nMaxRegions, float* output) 
{ 
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index < nMaxRegions)
	{
		int regionID = index;
		int offset_regionID = regionID*MAX_REGION_ADJ;	

		for(int a=0; a<MAX_REGION_ADJ; a++)
		{
			output[offset_regionID+a] = MY_FLT_MAX;
		}		
	}
}

template <typename T>
__global__ void kernelInitializeArray(T* __restrict__ a, const T value, const size_t n) {
      int tid = threadIdx.x + blockDim.x * blockIdx.x;
      if (tid < n) {
           a[tid] = value;
       }
}

namespace HyperSpectralToolbox
{
	void DeviceInitializeData(HSWO* object)
	{	
		//if (cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual) != cudaSuccess)
		{
			//printf("failed config\n");
		}

		// Allocate whole data cube in device memory
		int N = object->m_image_width*object->m_image_height*object->m_nBands;
		object->m_deviceData.m_dev_nRegions = object->m_image_width*object->m_image_height;
		
		//size_t size = N*sizeof(float);
		//CUDA_SAFE_CALL( cudaMalloc((void**)&object->m_deviceData.m_dev_pDataCube, size) );
		
		// Copy data cube from host memory to device memory
		//CUDA_SAFE_CALL( cudaMemcpy(object->m_deviceData.m_dev_pDataCube, object->m_pDataCube.get(), size, cudaMemcpyHostToDevice) );
		
		//regions keys
		size_t regionKeysSize = object->m_image_height*object->m_image_width*sizeof(int);
		cudaMalloc((void**)&(object->m_deviceData.m_dev_RegionKeys), regionKeysSize);		
		
		
		//adjacents
		size_t regionAdjSize = object->m_image_height*object->m_image_width*MAX_REGION_ADJ*sizeof(int);
		cudaMalloc((void**)&(object->m_deviceData.m_dev_pRegionAdjancencies), regionAdjSize);		
		
		
		//means
		int nMeans = object->m_image_height*object->m_image_width *object->m_nBands;
		
		size_t regionMeans = nMeans * sizeof(float);
		cudaMalloc((void**)&(object->m_deviceData.m_dev_pRegionMeans), regionMeans);
		
		
		//TODO can be removed, I now don't use the dissims ....
		//allocate device memory for regions dissimilarities 
		//size_t regionDissimsAdjSize = object->m_image_height*object->m_image_width*MAX_REGION_ADJ*sizeof(float);
		//CUDA_SAFE_CALL( cudaMalloc((void**)&(object->m_deviceData.m_dev_pDissims), regionDissimsAdjSize) );

		//regions pixel count
		{
			size_t regionPixelCountsSize = object->m_image_height*object->m_image_width*sizeof(int);
			cudaMalloc((void**)&(object->m_deviceData.m_dev_RegionPixelsCount), regionPixelCountsSize);
		}

		//m_dev_regionsAdjMinimums
		size_t regionsAdjMinimums = object->m_image_height*object->m_image_width;
		cudaMalloc((void**)&(object->m_deviceData.m_dev_regionsAdjMinimums), regionsAdjMinimums*sizeof(float));
		//m_dev_regionsAdjMinimumsLabels
		cudaMalloc((void**)&(object->m_deviceData.m_dev_regionsAdjMinimumsLabels), regionsAdjMinimums*sizeof(int)) ;
		//m_dev_needsAdjRecomputation
		cudaMalloc((void**)&(object->m_deviceData.m_dev_needsAdjRecomputation), regionsAdjMinimums*sizeof(int)) ;

		//m_dev_regionsMinimums		
		cudaMalloc((void**)&(object->m_deviceData.m_dev_regionsMinimums), regionsAdjMinimums*sizeof(float)) ;
		//m_dev_regionsMinimumsLabels
		cudaMalloc((void**)&(object->m_deviceData.m_dev_regionsMinimumsLabels), regionsAdjMinimums*sizeof(int)) ;
		//m_dev_needsRecomputation
		cudaMalloc((void**)&(object->m_deviceData.m_dev_needsRecomputation), regionsAdjMinimums*sizeof(int)) ;

		//new, for the 2dim spectral kernel
		size_t temp_all_dissims_size = object->m_deviceData._temp_all_dissims_size * sizeof(float);
		CUDA_SAFE_CALL (cudaMalloc((void**)&(object->m_deviceData._temp_all_dissims), temp_all_dissims_size); );
		
#if SPECTRAL_DYNAMIC_PROGRAMMING
		//int	threadsPerBlock = object->m_deviceData.m_threadsPerBlock;
		//int number_of_blocks =  object->m_deviceData._temp_all_dissims_size / (threadsPerBlock*threadsPerBlock);
		//dim3 gridDim(number_of_blocks, 1);
		//dim3 blockDim(threadsPerBlock*threadsPerBlock, 1);
		//kernelInitializeArray<float> <<<gridDim, blockDim>>>(object->m_deviceData._temp_all_dissims, MY_FLT_MAX, object->m_deviceData._temp_all_dissims_size);
		//cudaDeviceSynchronize();
#endif
		printf("temp_all_dissims_size = %d\n",temp_all_dissims_size);


	}
	
	void DeviceExit(HSWO* object)
	{
		// Free used device memory
		//if(object->m_deviceData.m_dev_pDataCube != NULL) 
		//{
			//CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_pDataCube) );;			
		//}
		//object->m_deviceData.m_dev_pDataCube = NULL;
		
		//========
		if(object->m_deviceData.m_dev_pRegionMeans != NULL) 
		{
			CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_pRegionMeans) );;			
		}
		object->m_deviceData.m_dev_pRegionMeans = NULL;
		
		//========
		//if(object->m_deviceData.m_dev_pDissims != NULL) 
		//{
		//		CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_pDissims) );;			
		//}
		//object->m_deviceData.m_dev_pDissims = NULL;
		
		//========
		if(object->m_deviceData.m_dev_pRegionAdjancencies != NULL) 
		{
			CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_pRegionAdjancencies) );;			
		}
		object->m_deviceData.m_dev_pRegionAdjancencies = NULL;
		
		//========
		if(object->m_deviceData.m_dev_RegionKeys != NULL) 
		{
			CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_RegionKeys) );;			
		}
		object->m_deviceData.m_dev_RegionKeys = NULL;

		//========
		if(object->m_deviceData.m_dev_RegionPixelsCount != NULL) 
		{
			CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_RegionPixelsCount) );;			
		}
		object->m_deviceData.m_dev_RegionPixelsCount = NULL;

		//========		
		{
			if(object->m_deviceData.m_dev_regionsAdjMinimums != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_regionsAdjMinimums) );;			
			}
			object->m_deviceData.m_dev_regionsAdjMinimums = NULL;

			if(object->m_deviceData.m_dev_regionsAdjMinimumsLabels != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_regionsAdjMinimumsLabels) );;			
			}
			object->m_deviceData.m_dev_regionsAdjMinimumsLabels = NULL;

			if(object->m_deviceData.m_dev_needsAdjRecomputation != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_needsAdjRecomputation) );;			
			}
			object->m_deviceData.m_dev_needsAdjRecomputation = NULL;
		}

		//========		
		{
			if(object->m_deviceData.m_dev_regionsMinimums != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_regionsMinimums) );;			
			}
			object->m_deviceData.m_dev_regionsMinimums = NULL;

			if(object->m_deviceData.m_dev_regionsMinimumsLabels != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_regionsMinimumsLabels) );;			
			}
			object->m_deviceData.m_dev_regionsMinimumsLabels = NULL;

			if(object->m_deviceData.m_dev_needsRecomputation != NULL) 
			{
				CUDA_SAFE_CALL( cudaFree(object->m_deviceData.m_dev_needsRecomputation) );;			
			}
			object->m_deviceData.m_dev_needsRecomputation = NULL;
		}

		//=====
		if(object->m_deviceData._temp_all_dissims != NULL) 
		{
			CUDA_SAFE_CALL( cudaFree(object->m_deviceData._temp_all_dissims) );;			
		}
		object->m_deviceData._temp_all_dissims = NULL;
	}
	
	void DeviceInitStep(HSWO* object)
	{
		//copy the following data to device
		//1- region adjacencies
		//2- region means
		//3- keys (regionID filter)
		
		//size_t regionKeysSize = object->m_image_height*object->m_image_width*sizeof(int);
		//cudaMemcpy(object->m_deviceData.m_dev_RegionKeys, object->m_deviceData.m_pRegionKeys, regionKeysSize , cudaMemcpyHostToDevice) ;
		
		//cudaError_t err = cudaGetLastError(); 
		//std::string serr = std::string(cudaGetErrorString(err)); 		
	}

	void DeviceUploadInitialMeans(HSWO* object)
	{
		//========================================================================
		//upload the means data				
		/*
		//int nMeans = object->m_image_height*object->m_image_width *object->m_nBands;
		//size_t regionMeans = nMeans * sizeof(float);
		cudaMemcpy(object->m_deviceData.m_dev_pRegionMeans,
			object->m_deviceData.m_pRegionsMeans,
			regionMeans, cudaMemcpyHostToDevice) ;
			*/
	
		std::hash_map<int, HSWO::Region*>::iterator it = object->m_regions.begin();
		for (;it != object->m_regions.end(); it++)
		{			
			int currentRegionLabel = it->first;
 			int device_pointer_offset = (currentRegionLabel-1)*object->m_nBands; 			
			cudaMemcpy(object->m_deviceData.m_dev_pRegionMeans+device_pointer_offset,
				it->second->sumOfPixels,
				object->m_nBands*sizeof(float), cudaMemcpyHostToDevice) ;


			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << std::endl;}
		}
	}
	
	void DeviceUploadInitialAdjacents(HSWO* object)
	{
		//upload the initial adjaccents data		
		size_t regionAdjSize = object->m_image_height*object->m_image_width*MAX_REGION_ADJ*sizeof(int);
		cudaMemcpy(object->m_deviceData.m_dev_pRegionAdjancencies, object->m_deviceData.m_pRegionAdjancencies
			, regionAdjSize, cudaMemcpyHostToDevice) ;
		
		cudaError_t err = cudaGetLastError(); 
		if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << std::endl;}
	}	

	void DeviceUpdateRegionMean(HSWO* object, int lastChangedRegionLabel, HSWO::Region* ptrLastChangedRegion)
	{
		int nMeans = object->m_nBands;
		size_t regionMeans = nMeans * sizeof(float);
		int device_pointer_offset = (lastChangedRegionLabel-1)*object->m_nBands; 
		
		/*
		cudaMemcpy(object->m_deviceData.m_dev_pRegionMeans+device_pointer_offset,
			&(object->m_deviceData.m_pRegionsMeans[(lastChangedRegionLabel-1)*object->m_nBands]),
			regionMeans, cudaMemcpyHostToDevice) ;
			*/

		cudaMemcpyAsync(object->m_deviceData.m_dev_pRegionMeans+device_pointer_offset,
			ptrLastChangedRegion->sumOfPixels,
			regionMeans, cudaMemcpyHostToDevice) ;

		cudaError_t err = cudaGetLastError(); 
		if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << std::endl;}
	}

	void DeviceExitStep(HSWO* object)
	{
	}
	
	typedef thrust::tuple<float,int> Tuple; 
	struct min_index 
	{ 
		__host__ __device__ 
			Tuple operator()(Tuple a, Tuple b) 
		{ 
			if (thrust::get<0>(a) < thrust::get<0>(b)) 
				return a; 
			else 
				return b; 
		} 
	}; 

	void DeviceCalcAllDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim)
	{
		int N = object->_max_nRegions;
		//int N = object->m_deviceData.m_regions.size();
		int maxRegions = object->_max_nRegions;
		
		//emptying all dissims values from previous computation
		//for
		//thrust::device_ptr<float> pdev1(object->m_deviceData.m_dev_pDissims);
		//thrust::fill_n(pdev1, maxRegions*MAX_REGION_ADJ, MY_FLT_MAX);
		
		/*dim3 numBlocks1(256,1);
		int threadsForEveryBlock1 = 256;
		dim3 threadsPerBlock1(threadsForEveryBlock1, 1); 
		kernel_fill_empty_dissims<<<numBlocks1, threadsPerBlock1>>>(maxRegions,object->m_deviceData.m_dev_pDissims);
		*/

		// Kernel invocation 
		//int threadsForEveryBlock2 = threadsPerBlock;
		int threadsForEveryBlock2 = HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH;	
		//printf("threads per block = %d \n", threadsForEveryBlock2);

		int nBlocks = N/threadsForEveryBlock2;
		dim3 numBlocks2(nBlocks,1);
		//printf("n blocks = %d \n", numBlocks2.x);
		
		dim3 threadsPerBlock2(threadsForEveryBlock2, 1); 
		//kernel_compute_dissim<<<numBlocks2, threadsPerBlock2, object->m_nBands*threadsForEveryBlock2* sizeof(float) >>>(object->m_deviceData.m_dev_RegionKeys
		kernel_compute_dissim<<<numBlocks2, threadsPerBlock2>>>(object->m_deviceData.m_dev_RegionKeys
															  ,object->m_deviceData.m_dev_pRegionAdjancencies
															  ,N		
															  ,object->m_nBands
															  ,object->m_deviceData.m_dev_pRegionMeans															  
															  ,object->m_deviceData.m_dev_RegionPixelsCount
															  ,object->m_deviceData.m_dev_regionsAdjMinimums
															  ,object->m_deviceData.m_dev_regionsAdjMinimumsLabels
															  ,object->m_deviceData.m_dev_needsAdjRecomputation); 
		
		
#if CUDA_42
		cudaThreadSynchronize();
#else
		cudaDeviceSynchronize();
#endif

		//printf("done kernel\n");
		
		// reduction step: get the minimum dissim measure and its region pair IDs
		// so that we can merge them on the host		
		//thrust::device_ptr<float> pdev(object->m_deviceData.m_dev_pDissims);
		//int maxDissims = maxRegions*MAX_REGION_ADJ;
		//thrust::device_ptr<float> result = thrust::min_element(pdev,pdev+(maxDissims) );
		//int index = (result - pdev);	

		thrust::device_ptr<float> pdev(object->m_deviceData.m_dev_regionsAdjMinimums);
		int maxDissims = maxRegions;
		
		thrust::device_ptr<float> result = thrust::min_element(pdev,pdev+(maxDissims) );
		minDissim = *result;		
		
		label1 = (result - pdev) + 1;		
		cudaMemcpy(&label2, object->m_deviceData.m_dev_regionsAdjMinimumsLabels+(label1-1), sizeof(int), cudaMemcpyDeviceToHost);
		

		//thrust::counting_iterator<int> Y(0);  
		//Tuple init(pdev[0],Y[0]); 
		//Tuple result = thrust::reduce 
		//	(thrust::make_zip_iterator(thrust::make_tuple(pdev, Y)), 
		//	thrust::make_zip_iterator(thrust::make_tuple(pdev+(maxDissims),   Y +(maxDissims))), 
		//	init, 
		//	min_index());  
		//float value; int index;  thrust::tie(value,index) = result; 
 
		////int index = result - object->m_deviceData.m_dev_pDissims;
				
		object->_last_merge_dissim_value = minDissim;
		object->_last_block_count = nBlocks;
		
		//return index;
		
	}

	
	void DeviceCalcRegionDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim)
	{
		int N = object->_max_nRegions;
		int maxRegions = object->_max_nRegions;

		
		// Kernel invocation parameters
		int threadsForEveryBlock2 = threadsPerBlock;	
		int nBlocks = N/threadsForEveryBlock2;
#if ONE_DIMENSIONAL_SPEC_KERNEL
		dim3 numBlocks2(nBlocks,1);
		dim3 threadsPerBlock2(threadsForEveryBlock2, 1); 
#else		
		//int *testArr = NULL;

		//cudaMalloc((void**)&testArr, 1024*sizeof(int));		
		//cudaMemset(testArr, 0, 1024*sizeof(int));
		
		
		//mutex initialization (for solution2)
		//int* dev_mutex2;
		//cudaMalloc((void**)&dev_mutex2, sizeof(int));		
		//cudaMemset(dev_mutex2, 0, sizeof(int));

		dim3 numBlocks2(nBlocks,nBlocks); //for solution 2
		dim3 threadsPerBlock2(threadsForEveryBlock2, threadsForEveryBlock2);  //for solution 2
#endif		
		
		//clock_t startClock2 = clock();

#if ONE_DIMENSIONAL_SPEC_KERNEL
#else
		//emptying all dissims values from previous computation		
		int number_of_blocks =  object->m_deviceData._temp_all_dissims_size / (threadsPerBlock*threadsPerBlock);
		dim3 gridDim(number_of_blocks, 1);
		dim3 blockDim(threadsPerBlock*threadsPerBlock, 1);
		kernelInitializeArray<float> <<<gridDim, blockDim>>>(object->m_deviceData._temp_all_dissims, MY_FLT_MAX, object->m_deviceData._temp_all_dissims_size);
		cudaDeviceSynchronize();

#endif
		//clock_t endClock2 = clock();
		//HSWO::time_count += (endClock2 - startClock2);
		
		//clock_t startClock2 = clock();

#if ONE_DIMENSIONAL_SPEC_KERNEL
		kernel_compute_regions_dissim<<<numBlocks2, threadsPerBlock2>>>(object->m_deviceData.m_dev_RegionKeys
																,object->m_deviceData.m_dev_pRegionAdjancencies
																,N		
																,object->m_nBands
																,object->m_deviceData.m_dev_pRegionMeans
																,object->m_deviceData.m_dev_RegionPixelsCount														  
																,object->m_deviceData.m_dev_regionsMinimums
																,object->m_deviceData.m_dev_regionsMinimumsLabels
																,object->m_deviceData.m_dev_needsRecomputation															  
																);
#else
		kernel_compute_regions_dissim_2dim<<<numBlocks2, threadsPerBlock2>>>(/*testArr, */object->m_deviceData.m_dev_RegionKeys //for solution 2		
															  ,object->m_deviceData.m_dev_pRegionAdjancencies
															  ,N		
															  ,object->m_nBands
															  ,object->m_deviceData.m_dev_pRegionMeans
															  ,object->m_deviceData.m_dev_RegionPixelsCount
															  /*
															  ,object->m_deviceData.m_dev_regionsMinimums
															  ,object->m_deviceData.m_dev_regionsMinimumsLabels															  
															  */
															  ,object->m_deviceData.m_dev_needsRecomputation															  
															  ,object->m_deviceData._temp_all_dissims); 
#endif
		
#if CUDA_42
		cudaThreadSynchronize();
#else
		cudaDeviceSynchronize();
#endif
		//printf("\nreached here 2\n");

		//clock_t endClock2 = clock();
		//HSWO::time_count += (endClock2 - startClock2);

		//printf("Test=======================\n");
		//printTestArr<<< 1024/128, 128>>>(testArr);

		//// First reduction
		//int number_of_blocks_reduction_1 =  object->_max_nRegions / (HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH);
		//dim3 gridDim_reduction_1(number_of_blocks_reduction_1, 1);
		//dim3 blockDim_reduction_1(HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH, 1);
		//update_all_spectral_minimums<<< gridDim_reduction_1, blockDim_reduction_1>>>(object->_max_nRegions
		//													,object->m_deviceData._temp_all_dissims
		//													,object->m_deviceData.m_dev_regionsMinimums
		//													,object->m_deviceData.m_dev_regionsMinimumsLabels);
		//cudaDeviceSynchronize();

		////Final reduction
		//thrust::device_ptr<float> pdev(object->m_deviceData.m_dev_regionsMinimums);
		//int maxDissims = maxRegions;

		//thrust::device_ptr<float> result = thrust::min_element(pdev,pdev+(maxDissims) );
		//label1 = (result - pdev) + 1;		
		//cudaMemcpy(&label2, object->m_deviceData.m_dev_regionsMinimumsLabels+(label1-1), sizeof(int), cudaMemcpyDeviceToHost);
		//minDissim = *result;

#if ONE_DIMENSIONAL_SPEC_KERNEL
		//Final reduction
		thrust::device_ptr<float> pdev(object->m_deviceData.m_dev_regionsMinimums);
		int maxDissims = maxRegions;

		thrust::device_ptr<float> result = thrust::min_element(pdev,pdev+(maxDissims) );
		label1 = (result - pdev) + 1;		
		cudaMemcpy(&label2, object->m_deviceData.m_dev_regionsMinimumsLabels+(label1-1), sizeof(int), cudaMemcpyDeviceToHost);
		minDissim = *result;

#else

		thrust::device_ptr<float> pdev(object->m_deviceData._temp_all_dissims);
		int maxDissims = object->m_deviceData._temp_all_dissims_size;

		thrust::device_ptr<float> result = thrust::min_element(pdev,pdev+(maxDissims) );
		int offset = (result - pdev);		
		label1 = offset / maxRegions +1;
		label2 = offset % maxRegions +1 ;

		minDissim = *result;

#endif		
		//printf("last mindissim = %f, new dissim = %f, %d, %d\n", object->_last_merge_dissim_value, minDissim , label1, label2);
		//printf("\nreached here 2.2\n");

		if(minDissim < (object->_last_merge_dissim_value*object->_spclustWeight) )
		{
			//printf("found spclust merge\n");			
		}
		else
		{
			minDissim = FLT_MAX;
			label1 = -1;
			label2 = -1;
		}
		object->_last_block_count = nBlocks;

#if ONE_DIMENSIONAL_SPEC_KERNEL
#else
		//printf("\nreached here 2.3\n");
		//cudaFree( dev_mutex2 ); //for solution 2
		//cudaFree( testArr);
#endif
	}

	void DeviceResetOthersBestRegionComputationFlagsFrom(HSWO* object, int region1Label, int region2Label ) 
	{		
		int maxRegions = object->_max_nRegions;
		
		// Kernel invocation 
		//int threadsForEveryBlock2 = object->m_deviceData.m_threadsPerBlock;	
		int threadsForEveryBlock2 = HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH;	
		int nBlocks = maxRegions/threadsForEveryBlock2;
		dim3 numBlocks2(nBlocks,1);
		
		dim3 threadsPerBlock2(threadsForEveryBlock2, 1); 
		reset_others_best_region_computation_flags_from<<<numBlocks2, threadsPerBlock2>>>(region1Label, region2Label
																, maxRegions
																, object->m_deviceData.m_dev_RegionKeys
															  ,object->m_deviceData.m_dev_needsRecomputation													
															  ,object->m_deviceData.m_dev_regionsMinimumsLabels
															  ,object->m_deviceData.m_dev_regionsMinimums); 
		
		cudaThreadSynchronize();
	}

	void DeviceRecomputeOthersBestRegionsDissimilarityTo(HSWO* object, int regionLabel)
	{		
		int maxRegions = object->_max_nRegions;
		
		// Kernel invocation 
		//int threadsForEveryBlock2 = object->m_deviceData.m_threadsPerBlock;	
		int threadsForEveryBlock2 = HARD_CODED_1_DIM_CUDA_KERNELS_BLOCK_WIDTH;	
		int nBlocks = maxRegions/threadsForEveryBlock2;
		dim3 numBlocks2(nBlocks,1);
		
		dim3 threadsPerBlock2(threadsForEveryBlock2, 1); 
		recompute_others_best_region_dissimilarity_to<<<numBlocks2, threadsPerBlock2>>>(regionLabel
																, maxRegions
																, object->m_deviceData.m_dev_RegionKeys
															  ,object->m_deviceData.m_dev_needsRecomputation													
															  ,object->m_deviceData.m_dev_regionsMinimumsLabels
															  ,object->m_deviceData.m_dev_regionsMinimums
															  ,object->m_nBands
															  ,object->m_deviceData.m_dev_pRegionMeans
															  ,object->m_deviceData.m_dev_RegionPixelsCount);  
		
		cudaThreadSynchronize();
	}

	void DeviceCalcAllRegionsMeanVectors(HSWO* object)
	{
		/*
		//allocate region means matrix of nRegions rows x nBands columns
		// rows : regions
		// columns : nBands columns with each columns value 
		// represent the mean value of the region at this 
		// band
		 
		int N = object->m_dev_nRegions * object->m_nBands;
		size_t size = N*sizeof(float);
		CUDA_SAFE_CALL( cudaMalloc((void**)&object->m_dev_pRegionMeans, size) );
		
		//calculate
		
		// Copy data cube from host memory to device memory
		CUDA_SAFE_CALL( cudaMemcpy(object->m_dev_pDataCube, object->m_pDataCube, size, cudaMemcpyHostToDevice) );
*/		
	}
	
	void DeviceThrustRemove(int* pArray, int size, int valueToRemove)
	{
		thrust::device_ptr<int> pdevArray(pArray);
		thrust::remove(pdevArray, pdevArray+size, valueToRemove);
	}

	//template <typename T> 
	void DeviceThrustFill(int* pdeviceArray, int size, int newValue)
	{ 
		thrust::device_ptr<int> pdevArray(pdeviceArray);
		thrust::fill_n(pdevArray, size, newValue);
	}
	void DeviceThrustFillFloat(float* pdeviceArray, int size, float newValue)
	{ 
		thrust::device_ptr<float> pdevArray(pdeviceArray);
		thrust::fill_n(pdevArray, size, newValue);
	}

	void DeviceThrustResetAdjacencyOfMergedRegion(int* m_deviceDatam_dev_pRegionAdjancencies, int lastMergedRegion_OffsetInMatrix)
	{
		//emptying all adjacencies on the device of the merged region
		thrust::device_ptr<int> pdevArray(m_deviceDatam_dev_pRegionAdjancencies+lastMergedRegion_OffsetInMatrix);
		thrust::fill_n(pdevArray, MAX_REGION_ADJ, -1);
	}

	void DeviceThrustReplace(int* m_deviceData_DOT_m_dev_pRegionAdjancencies_PLUS_currentAdjacentRegion_OffsetInMatrix
		, int nAdjacentsOfAdjacent, int adjRegionLabel, int regionLabel )
	{		
		thrust::device_ptr<int> pdevArray(m_deviceData_DOT_m_dev_pRegionAdjancencies_PLUS_currentAdjacentRegion_OffsetInMatrix);
		thrust::replace(pdevArray, pdevArray+nAdjacentsOfAdjacent, adjRegionLabel, regionLabel);		
	}

	float computeMax(HSWO* object)
	{
		if(object->m_deviceData.m_dev_pDataCube == NULL) exit(0);
		
		int N = object->m_image_width*object->m_image_height*object->m_nBands;
		
		thrust::device_ptr<float> dev_ptr(object->m_deviceData.m_dev_pDataCube);	
		
		thrust::device_ptr<float> result = thrust::max_element(dev_ptr,dev_ptr+N);
		
		return *result;
	}
	
	float computeMin(HSWO* object)
	{
		if(object->m_deviceData.m_dev_pDataCube == NULL) exit(0);
		
		int N = object->m_image_width*object->m_image_height*object->m_nBands;
		
		// Use already copied data to device
		thrust::device_ptr<float> dev_ptr(object->m_deviceData.m_dev_pDataCube);	
		
		thrust::device_ptr<float> result = thrust::min_element(dev_ptr,dev_ptr+N);
		
		return *result;
	}
	
	void runThrustExpr(thrust::host_vector<int>& h_vec, int N)
	{
		// transfer data to the device
		thrust::device_vector<int> d_vec = h_vec;
		
		// sort data on the device
		
		thrust::sort(d_vec.begin(), d_vec.end());
		
		// transfer data back to host
		thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
	}
}