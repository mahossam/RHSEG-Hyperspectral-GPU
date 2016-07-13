#if USE_AMP
#define NOMINMAX

#include "HSWO.h"
#include <cmath>
#include <climits>
#include <algorithm>
#include <iostream>
#include "./Libraries/definitions.h"
#include "./Libraries/Logger.h"
#include "./Libraries/DebuggingUtilities.h"

#include <amp_math.h>
#include <amp.h> 
#include <amp_stl_algorithms.h>


#include "HSWODevice.h"

#include "myglobals.h"

using namespace concurrency; 
using namespace HyperSpectralToolbox;
//using namespace amp_stl_algorithms;
using namespace std;

#define PI 3.1415926535897932384626433832795
#define MY_FLT_MAX         3.402823466e+38F        /* max value */

void HyperSpectralToolbox::DeviceInitializeData(HSWO* object){}
void HyperSpectralToolbox::DeviceExit(HSWO* object){}

void HyperSpectralToolbox::DeviceInitStep(HSWO* object){}
void HyperSpectralToolbox::DeviceExitStep(HSWO* object){}

void HyperSpectralToolbox::DeviceUpdateRegionMean(HSWO* object, int lastChangedRegionID, HSWO::Region* ptrLastChangedRegionID){}
void HyperSpectralToolbox::DeviceUploadInitialMeans(HSWO* object){}
void HyperSpectralToolbox::DeviceUploadInitialAdjacents(HSWO* object){}


void HyperSpectralToolbox::DeviceCalcAllDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim){}
void HyperSpectralToolbox::DeviceCalcRegionDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim){}
void HyperSpectralToolbox::DeviceResetOthersBestRegionComputationFlagsFrom(HSWO* object, int region1Label, int region2Label ) {}
void HyperSpectralToolbox::DeviceRecomputeOthersBestRegionsDissimilarityTo(HSWO* object, int regionLabel){}

void HyperSpectralToolbox::DeviceCalcAllRegionsMeanVectors(HSWO* object){}

//template <typename T> 
void HyperSpectralToolbox::DeviceThrustFill(int* pdeviceArray, int size, int newValue){}
void HyperSpectralToolbox::DeviceThrustFillFloat(float* pdeviceArray, int size, float newValue){}
void HyperSpectralToolbox::DeviceThrustRemove(int* pdeviceArray, int size, int valueToRemove){}

void HyperSpectralToolbox::DeviceThrustResetAdjacencyOfMergedRegion(int* m_deviceDatam_dev_pRegionAdjancencies, int lastMergedRegion_OffsetInMatrix){}
void HyperSpectralToolbox::DeviceThrustReplace(int* m_deviceData_DOT_m_dev_pRegionAdjancencies_PLUS_currentAdjacentRegion_OffsetInMatrix, int nAdjacentsOfAdjacent, int adjRegionLabel, int regionLabel ){}


//----------------------------------------------------------------------------
// This is an implementation of the reduction algorithm using a simple
// parallel_for_each. Multiple kernel launches are required to synchronize
// memory access among threads in separate tiles.
//----------------------------------------------------------------------------
float HyperSpectralToolbox::reduction_simple_1_min_INPLACE(Concurrency::array<float, 1>& source, int& minIndex)
{
	// Using array, as we mostly need just temporary memory to store
	// the algorithm state between iterations and in the end we have to copy
	// back only the first element.

	unsigned element_count = source.extent.size();
	auto av = source.view_as(Concurrency::extent<1>(element_count));

	if (element_count == 1)
	{
		return av[0];
	}

	array<int, 1> indexplaceholder(element_count);
	auto av_indexplaceholder = indexplaceholder.view_as(Concurrency::extent<1>(element_count));
	parallel_for_each(Concurrency::extent<1>(element_count), [&indexplaceholder] (index<1> idx) restrict(amp)
	{
		indexplaceholder[idx] = idx[0];			
	});

	// Takes care of odd input elements we could completely avoid tail sum
	// if we would require source to have even number of elements.
	// Each thread reduces two elements.
	int n = element_count;
	for (unsigned s = n / 2; s > 0; s /= 2)
	{
		parallel_for_each(Concurrency::extent<1>(s), [s, &source, &indexplaceholder] (index<1> idx) restrict(amp)
		{
			if(source[idx] > source[idx+s])
			{
				source[idx] = source[idx+s];
				indexplaceholder[idx] = indexplaceholder[idx+s];
			}

		});

		// Reduce the tail in cases where the number of elements is odd.	            
		//barrier
		if ((n & 0x1) && (n != 1))
		{
			float maxy = av[0];
			int maxyIndex = 0;
			for (int i = 1; i< s-1; i++)
			{
				if(av[i] > maxy)
				{
					maxy = av[i];
					maxyIndex = i;
				}
			}

			if(av[n - 1] < maxy)
			{
				av[maxyIndex] = av[n-1];
				av_indexplaceholder[maxyIndex] = av_indexplaceholder[n-1];
			}                
		}		
		n = s;
	}

	minIndex = av_indexplaceholder[0];
	return av[0];
}

void HSWO::Device_doStep()
{
#define tile_width 16
#define MAX_BANDS 220
	//#define FLAG_VALUE (MY_FLT_MAX-9999.0f)

	int regionLabel = -1, adjRegionLabel = -1; float minDissim=FLT_MAX;
	int maxRegions = this->m_image_height*this->m_image_width;

	int nBands = this->m_nBands;
	//int nRegions = this->m_regions.size();

	Concurrency::array<int, 1>& KeysFilter = *m_deviceData.amp_keysFilter;
	Concurrency::array<int, 1>& adj = *m_deviceData.amp_adj;
	Concurrency::array<float, 1>& regionsSums= *m_deviceData.amp_regionsSums;
	Concurrency::array<int, 1>& regionsPixelsCount= *m_deviceData.amp_regionsPixelsCount;

	{
		Concurrency::array<float, 1>& regionsAdjMinimums= *m_deviceData.amp_regionsAdjMinimums;
		Concurrency::array<int, 1>& regionsAdjMinimumsLabels= *m_deviceData.amp_regionsAdjMinimumsLabels;	
		Concurrency::array<int, 1>& needsAdjRecomputation = *m_deviceData.amp_needsAdjRecomputation;	

		parallel_for_each(KeysFilter.extent, [nBands, &KeysFilter, &adj, &regionsSums, &regionsPixelsCount, &regionsAdjMinimums, &regionsAdjMinimumsLabels, &needsAdjRecomputation] (index<1> idx) restrict(amp)
		{
			int indexval = idx[0];
			int regionID = KeysFilter[indexval]-1;
			if (regionID == -2) return;

			amp_float SHregionsAdjMinimums;
			SHregionsAdjMinimums = MY_FLT_MAX;
			int SHregionsAdjMinimumsLabels;
			SHregionsAdjMinimumsLabels = -1;

			//printf(" - [Thread %d] region %d\n", index, regionID+1 );

			// ****** Optimization ****************
			//if(needsAdjRecomputation[regionID] == 1) //TODO: VIP dynamic programming is not available 
													// as long as the reduction step is using inplace reduction
													//check the device_merge_region statements that resets the WHOLE minumums array
													//and fix it before trying the dynamic solution
			if(true)
			{
				///////////////////////////////////////
				int offset_regionID = regionID*MAX_REGION_ADJ;	
				amp_float r1nPixels = regionsPixelsCount[regionID];

				//printf(" - [Thread %d] region %d pixel count = %f \n", index, regionID+1 , r1nPixels);

				for(int a=0; a<MAX_REGION_ADJ; a++)
				{				
					int regionIDadj = adj[offset_regionID+a]-1;				

					if (regionIDadj <= -1) 
					{
						break;
					}
					//if(regionID > regionIDadj) continue;  //optimization			
					amp_float sum_eculidean = 0.0f;
					amp_float r2nPixels = regionsPixelsCount[regionIDadj];

					//printf(" - [Thread %d] ADJ region %d pixel count = %f \n", index, regionIDadj+1, r2nPixels);

					amp_float temp;

					int band = 0;

					for (; band < nBands; ++band)
					{
						temp = regionsSums[regionID*nBands+band]/(r1nPixels) - regionsSums[regionIDadj*nBands+band]/(r2nPixels);
						//					printf("	--- [Thread %d] region %d sum at band %d = %f \n", index, regionID, band, regionsSums[regionID*nBands+band]);
						//					printf("	--- [Thread %d] region %d sum at band %d = %f \n", index, regionIDadj, band, regionsSums[regionIDadj*nBands+band]);
						sum_eculidean += temp * temp;
					}

					//output[offset_regionID+a] = sqrtf(sum_eculidean) / float(nBands);
					//output[offset_regionID+a] = sqrtf(sum_eculidean*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));
					double d = sum_eculidean*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels));
					temp =  amp_math::sqrt(d);
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
				//needsAdjRecomputation[regionID] = 0;
			}
			else
			{
				//no recomputation needed: do nothing
			}
		}
		);
		

		/*
		array<int, 1> mutex(1);
		auto av_mutex = mutex.view_as(mutex.extent);
		av_mutex[0] = 0;
		av_mutex.synchronize();

		parallel_for_each(Concurrency::extent<2>(KeysFilter.extent.size(),KeysFilter.extent.size()).tile<tile_width,tile_width>()
			, [&mutex, nBands, &KeysFilter, &adj, &regionsSums, &regionsPixelsCount, &regionsAdjMinimums, &regionsAdjMinimumsLabels] (tiled_index<tile_width, tile_width> double_idx) restrict(amp)
		{
			int region1Index = double_idx.global[0];
			int region2Index = double_idx.global[1];

			amp_float sum_dissim = 0.0f;
			amp_float temp = MY_FLT_MAX;
			int can_compute = 1;

			int region1ID = KeysFilter[region1Index]-1;				
			if (region1ID == -2) 
			{
				//return;
				can_compute = 0;
			}

			int region2ID = KeysFilter[region2Index]-1;				
			if (region2ID == -2) 
			{
				//return;
				can_compute = 0;
			}

			if (region1ID == region2ID) 
			{
				//return;
				can_compute = 0;
			}


			double_idx.barrier.wait();

			if(can_compute)
			{
				amp_float r1nPixels = regionsPixelsCount[region1ID];

				// ****** Optimization ****************			
				//if(needsRecomputation[regionID1] == 1)			
				//if(true)		
				//{
				//for(int i=0; i < maxRegions; i++)
				//{
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

				if(foundAdjacent == 1)
				{
					amp_float r2nPixels = regionsPixelsCount[region2ID];
					//float r2nPixels = local_region2_nPixels[double_idx.local[1]];

					//if(region1ID > region2ID) return;  //optimization			
					for (int band = 0; band < nBands; ++band)
					{
						temp = regionsSums[region1ID*nBands+band]/(r1nPixels) - regionsSums[region2ID*nBands+band]/(r2nPixels);
						sum_dissim += temp * temp;
					}

					double d = sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels));
					temp =  amp_math::sqrt(d);
				}
			}

			double_idx.barrier.wait();

			//race (global minimums update)
			//if(double_idx.local[1] == 0)
			if(region1ID >= 0)				
			{
				//mutex			
				int ev1 = 0;
				while(Concurrency::atomic_compare_exchange(&mutex[0],&ev1, 1) != true)
				{
				};
				//critical section

				amp_float glbl = regionsAdjMinimums[region1ID] ;

				//if(local_minimums_final[double_idx.local[0]] < glbl)
				if(temp < glbl)
				{
					regionsAdjMinimums[region1ID] = temp;
					regionsAdjMinimumsLabels[region1ID] = region2ID+1;
				}

				//end critical section
				int ev2 = 1;
				Concurrency::atomic_compare_exchange(&mutex[0], &ev2, 0);

				//regionsMinimums[region1ID] = LocalregionsMinimums;
				//regionsMinimumsLabels[region1ID] = LocalregionsMinimumsLabels;
				//needsRecomputation[regionID1] = 0;
			}		
		}
		);
		*/
			
		//Reduction step (finding minimum)
		minDissim = reduction_simple_1_min_INPLACE(*(m_deviceData.amp_regionsAdjMinimums), regionLabel);
		++regionLabel;
		auto av_minimumsLabels = m_deviceData.amp_regionsAdjMinimumsLabels->view_as(m_deviceData.amp_regionsAdjMinimumsLabels->extent);
		adjRegionLabel = av_minimumsLabels[regionLabel-1];		

		_last_merge_dissim_value = minDissim;
		m_deviceData.currentStepRegionID_to_merge_1 = regionLabel;
		m_deviceData.currentStepRegionID_to_merge_2 = adjRegionLabel;
		_last_block_count = -17;

		//if (adjRegionLabel == -1)
		//{
		//	std::cout<<"yanhaaaaaaaaaaaaaaar !\n";
		//}

		Region* region1 = this->m_deviceData.allRegions[regionLabel-1];		
		Region* region2 = this->m_deviceData.allRegions[adjRegionLabel-1];

		//cout << "AMP step done\n";

		//now we can merge the minimum dissimilarity pair
		Device_merge_regions(regionLabel, region1, adjRegionLabel, region2);
	}

	///========================================================================
	/// spectral clustering step
	///========================================================================
	//if(false)
	if(m_regions.size() > m_min_no_of_clusters)	
	{
		regionLabel = -1, adjRegionLabel = -1; minDissim=FLT_MAX;
		Concurrency::array<int, 1>& needsRecomputation= *m_deviceData.amp_needsRecomputation;							
		Concurrency::array<int, 1>& regionsMinimumsLabels= *m_deviceData.amp_regionsMinimumsLabels;	
		Concurrency::array<float, 1>& regionsMinimums= *m_deviceData.amp_regionsMinimums;

		
		/*
		parallel_for_each(KeysFilter.extent, [nBands, maxRegions, &KeysFilter, &adj, &regionsSums, &regionsPixelsCount, &regionsMinimums, &regionsMinimumsLabels] (index<1> idx) restrict(amp)
		{
			int index1 = idx[0];

			int regionID1 = KeysFilter[index1]-1;				
			if (regionID1 == -2) return;

			float r1nPixels;
			r1nPixels = regionsPixelsCount[regionID1];

			float LocalregionsMinimums;
			LocalregionsMinimums = MY_FLT_MAX;
			int LocalregionsMinimumsLabels;
			LocalregionsMinimumsLabels = -1;		

			// ****** Optimization ****************			
			//if(needsRecomputation[regionID1] == 1)			
			if(true)		
			{
				for(int i=0; i < maxRegions; i++)
				{
					int regionID2 = KeysFilter[i]-1;		
					if (regionID2 == -2) continue;

					if (regionID1 == regionID2) continue;

					int offset_regionID = regionID1*MAX_REGION_ADJ;	
					int foundAdjacent = 0;
					for(int a=0; a<MAX_REGION_ADJ; a++)
					{				
						int regionIDadj = adj[offset_regionID+a]-1;				

						if (regionIDadj <= -1) 
						{
							break;
						}
						if(regionIDadj == regionID2) {
							foundAdjacent = 1;
							break;							
						}
					}
					
					// TODO: VIP as learned from CUDA, it is better to put a barrier here to 
					// prevent accidental thread divergence (if  AMP allows it or using tiles or something)
					

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

						temp = amp_math::sqrt(sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels)));

						if(temp < LocalregionsMinimums) 
						{
							LocalregionsMinimums = temp;
							LocalregionsMinimumsLabels = regionID2+1;
						}
					}
					//if(temp < regionsMinimums[regionID1]) 
					//{
					//	regionsMinimums[regionID1] = temp;
					//	regionsMinimumsLabels[regionID1] = regionID2+1;
					//}				
				}		
				regionsMinimums[regionID1] = LocalregionsMinimums;
				regionsMinimumsLabels[regionID1] = LocalregionsMinimumsLabels;
				//needsRecomputation[regionID1] = 0;
			}		
		}
		); 
		*/

		array<int, 1> mutex(1);
		auto av_mutex = mutex.view_as(mutex.extent);
		av_mutex[0] = 0;
		av_mutex.synchronize();


		parallel_for_each(Concurrency::extent<2>(KeysFilter.extent.size(),KeysFilter.extent.size()).tile<tile_width,tile_width>() ,
			[&mutex, nBands, maxRegions, &KeysFilter, &adj, &regionsSums, &regionsPixelsCount, &regionsMinimums, &regionsMinimumsLabels] (tiled_index<tile_width, tile_width> double_idx) restrict(amp)
		{
			int region1Index = double_idx.global[0];
			int region2Index = double_idx.global[1];

			amp_float sum_dissim = 0.0f;
			amp_float temp = MY_FLT_MAX;
			int can_compute = 1;

			int region1ID = KeysFilter[region1Index]-1;				
			if (region1ID == -2) 
			{
				//return;
				can_compute = 0;
			}

			int region2ID = KeysFilter[region2Index]-1;				
			if (region2ID == -2) 
			{
				//return;
				can_compute = 0;
			}

			if (region1ID == region2ID) 
			{
				//return;
				can_compute = 0;
			}

			//local npixels, adj, sums, 
			//tile_static int local_region2_keys[tile_width];

			//tile_static int local_region1_nPixels[tile_width];			
			//tile_static int local_region2_nPixels[tile_width];						
			//tile_static int local_adj[tile_width][MAX_REGION_ADJ];
			//tile_static float local_sums1[tile_width][MAX_BANDS];
			//tile_static float local_sums2[tile_width][MAX_BANDS];

			//tile_static float local_minimums_store[tile_width][tile_width];
			//tile_static int local_minimums_store_labels[tile_width][tile_width];

			//tile_static float local_minimums_final[tile_width];
			//tile_static int local_minimums_final_labels[tile_width];

			//local_minimums_store[double_idx.local[0]][double_idx.local[1]] = MY_FLT_MAX;
			//local_minimums_store_labels[double_idx.local[0]][double_idx.local[1]] = -1;

			//if(double_idx.local[1] == 0)
			//{			
			//	local_minimums_final[double_idx.local[0]] = MY_FLT_MAX;
			//	local_minimums_final_labels[double_idx.local[0]] = -1;
			//}

			//the first thread in the tile(block)			
			//if(double_idx.local[0] == 0 && double_idx.local[1] == 0)
			//{
			//	//load the local arrays from global memory
			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		//local_region2_keys[t] = KeysFilter[region2ID+t];
			//		local_region1_nPixels[t] = regionsPixelsCount[region1ID+t];
			//	}
			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		local_region2_nPixels[t] = regionsPixelsCount[region2ID+t];
			//	}
			//	
			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		for (int a = 0; a < MAX_REGION_ADJ; ++a)
			//		{
			//			int offset_regionID = (region1ID+t)*MAX_REGION_ADJ;	
			//			int adjVal = adj[offset_regionID + a];
			//			local_adj[t][a] = adjVal;
			//			if (adjVal <= -1) 
			//			{
			//				break;
			//			}						
			//		}
			//	}

			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		for (int b = 0; b < nBands; ++b)
			//		{
			//			int offset_regionID = (region1ID+t)*nBands;	
			//			local_sums1[t][b] = regionsSums[offset_regionID + b];
			//		}
			//	}
			//	
			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		for (int b = 0; b < nBands; ++b)
			//		{
			//			int offset_regionID = (region2ID+t)*nBands;	
			//			local_sums2[t][b] = regionsSums[offset_regionID + b];
			//		}
			//	}

			//}

			double_idx.barrier.wait();

			//float LocalregionsMinimums;
			//LocalregionsMinimums = MY_FLT_MAX;
			//int LocalregionsMinimumsLabels;
			//LocalregionsMinimumsLabels = -1;		

			if(can_compute)
			{
				amp_float r1nPixels = regionsPixelsCount[region1ID];

				// ****** Optimization ****************			
				//if(needsRecomputation[regionID1] == 1)			
				//if(true)		
				//{
				//for(int i=0; i < maxRegions; i++)
				//{
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

				// TODO: VIP as learned from CUDA, it is better to put a barrier here to 
				// prevent accidental thread divergence (if AMP allows it)				

				if(foundAdjacent == 0)
				{
					amp_float r2nPixels = regionsPixelsCount[region2ID];
					//float r2nPixels = local_region2_nPixels[double_idx.local[1]];

					//if(region1ID > region2ID) return;  //optimization			
					for (int band = 0; band < nBands; ++band)
					{
						temp = regionsSums[region1ID*nBands+band]/(r1nPixels) - regionsSums[region2ID*nBands+band]/(r2nPixels);
						sum_dissim += temp * temp;
					}

					double d = sum_dissim*((r1nPixels*r2nPixels)/(r1nPixels+r2nPixels));
					temp =  amp_math::sqrt(d);

					//if(temp < LocalregionsMinimums) 
					{
						//LocalregionsMinimums = temp;
						//LocalregionsMinimumsLabels = region2ID+1;
						//local_minimums_store[double_idx.local[0]][double_idx.local[1]] = temp;
						//local_minimums_store_labels[double_idx.local[0]][double_idx.local[1]] = region2ID+1;

					}
				}
				//if(temp < regionsMinimums[regionID1]) 
				//{
				//	regionsMinimums[regionID1] = temp;
				//	regionsMinimumsLabels[regionID1] = regionID2+1;
				//}				
				//}			
			}

			//double_idx.barrier.wait();
			//get the minimum(s) for current tile

			//if(double_idx.local[1] == 0)
			//{		
			//	float miny = MY_FLT_MAX;
			//	int minyLabel = -1;

			//	for (int t = 0; t < tile_width; ++t)
			//	{
			//		if(local_minimums_store[double_idx.local[0]][t] < miny)
			//		{
			//			miny = local_minimums_store[double_idx.local[0]][t];
			//			minyLabel = local_minimums_store_labels[double_idx.local[0]][t];
			//		}
			//	}
			//	local_minimums_final[double_idx.local[0]] = miny;
			//	local_minimums_final_labels[double_idx.local[0]] = minyLabel;
			//}

			double_idx.barrier.wait();

			//race (global minimums update)
			//if(double_idx.local[1] == 0)
			if(region1ID >= 0)				
			{
				//mutex			
				int ev1 = 0;
				while(Concurrency::atomic_compare_exchange(&mutex[0],&ev1, 1) != true)
				{
				};
				//critical section

				amp_float glbl = regionsMinimums[region1ID] ;

				//if(local_minimums_final[double_idx.local[0]] < glbl)
				if(temp < glbl)
				{
					regionsMinimums[region1ID] = temp;
					regionsMinimumsLabels[region1ID] = region2ID+1;
					//TODO needs recomputation flag, should it be updated here or where ? revise with CPU implementation
				}

				//end critical section
				//int ev2 = 1;
				//Concurrency::atomic_compare_exchange(&mutex[0], &ev2, 0);
				Concurrency::atomic_exchange(&mutex[0], 0);

				//regionsMinimums[region1ID] = LocalregionsMinimums;
				//regionsMinimumsLabels[region1ID] = LocalregionsMinimumsLabels;
				//needsRecomputation[regionID1] = 0;
			}		
		}
		); 		
		

		//Reduction step (finding minimum)
		minDissim = reduction_simple_1_min_INPLACE(*(m_deviceData.amp_regionsMinimums), regionLabel);
		++regionLabel;
		auto av_minimumsLabels = m_deviceData.amp_regionsMinimumsLabels->view_as(m_deviceData.amp_regionsMinimumsLabels->extent);
		adjRegionLabel = av_minimumsLabels[regionLabel-1];		

		//std::cout<<
		//	"SPEC last dissim value = " + ts(minDissim)
		//	+ ", r1 = " + ts(regionLabel)
		//	+ ", r2 = " + ts(adjRegionLabel)						
		//	<<endl;

		if(minDissim < (_last_merge_dissim_value*_spclustWeight) )
		{
			//printf("found spclust merge\n");			
			_last_merge_dissim_value = minDissim;
			m_deviceData.currentStepRegionID_to_merge_1 = regionLabel;
			m_deviceData.currentStepRegionID_to_merge_2 = adjRegionLabel;
		}
		else
		{
			minDissim = FLT_MAX;
			regionLabel = -1;
			adjRegionLabel = -1;
		}

		if(regionLabel != -1 && adjRegionLabel != -1)
		{
			Region* region1 = this->m_deviceData.allRegions[regionLabel-1];		
			Region* region2 = this->m_deviceData.allRegions[adjRegionLabel-1];

			Device_merge_regions(regionLabel, region1, adjRegionLabel, region2);
		}
	}
}


void HSWO::Device_Initialize()
{	
	try
	{
		vector<Region>& regionsObjects = *(m_regionsObjects);

		/// Build initial set of regions (each pixel is a region)
		int nPixels = m_image_height * m_image_width;

		this->m_deviceData.allRegions = new Region* [nPixels];		

		for(int y = 0;y < m_image_height;y++){
			for(int x = 0;x < m_image_width;x++)
			{				
				auto createNewRegion = [&](int newRegionLabel)->Region* {
					regionsObjects[newRegionLabel-1].Construct(newRegionLabel, m_nBands);

					//Region* ptrNewRegion = new Region(newRegionLabel, m_nBands);
					Region* ptrNewRegion = &(regionsObjects[newRegionLabel-1]);
					this->m_regions[newRegionLabel] = ptrNewRegion ;

					this->m_deviceData.allRegions[newRegionLabel-1] = ptrNewRegion;

					ptrNewRegion ->pixels_IDs.push_back(newRegionLabel);
					ptrNewRegion->_pixelCount += 1;
					//m_labeled_segmented_image[RegionLabel] = RegionLabel;

					//calculate initial regions means
					for (int band = 0; band < m_nBands; band++)
					{
						(ptrNewRegion->sumOfPixels)[band] = (m_pDataCube.get())[(newRegionLabel-1)+band*m_image_height*m_image_width];
						//this->m_pRegionsSums[(RegionLabel-1)*m_nBands + band] = (m_pDataCube.get())[(RegionLabel-1)+band*m_image_height*m_image_width];
					}
					return ptrNewRegion;
				};

				// new logic, I decided to make the regions and pixels IDs 1 based to easily match 
				// the original RHSEG author logic
				int RegionLabel = x + y * m_image_width +1;

				Region* ptrCurrentRegion = nullptr;
				if(m_regions.find(RegionLabel) == m_regions.end())
				{
					ptrCurrentRegion = createNewRegion(RegionLabel);
				}
				else
				{
					ptrCurrentRegion = m_regions[RegionLabel];
				}


				// add boundary regions (which are the 8 surrounding pixels) for this current region(x,y)
				//topleft
				if(x > 0 && y > 0){
					int pixelIndex = (x - 1) + (y - 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//top
				if(y > 0){
					int pixelIndex = (x) + (y - 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//topright
				if(x < (m_image_width - 1) && y > 0){
					int pixelIndex = (x + 1) + (y - 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//left
				if(x > 0){
					int pixelIndex = (x - 1) + (y) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//right
				if(x < (m_image_width - 1)){
					int pixelIndex = (x + 1) + (y) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//bottomleft
				if(x > 0 && y < (m_image_height - 1)){
					int pixelIndex = (x - 1) + (y + 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//bottom
				if(y < (m_image_height - 1)){
					int pixelIndex = (x) + (y + 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
				//bottomright
				if(x < (m_image_width - 1) && y < (m_image_height - 1)){
					int pixelIndex = (x + 1) + (y + 1) * m_image_width;
					if(m_regions.find(pixelIndex+1) == m_regions.end())
					{
						createNewRegion(pixelIndex+1);
					}
					ptrCurrentRegion->adjacentRegions.insert(m_regions[pixelIndex+1]);
				}
			}

		}

		//DeviceUploadInitialMeans(this);		
		{			
			m_deviceData.amp_regionsSums = new Concurrency::array<float, 1>(getRegionsCount()*m_nBands);
			auto& v = *(m_deviceData.amp_regionsSums);		
			auto av = v.view_as(v.extent);

			for (const auto& r : m_regions)
			{			
				int currentRegionLabel = r.first;
				int index = (currentRegionLabel-1)*m_nBands; 			
				for (int b=0; b< m_nBands; ++b)
				{
					av[index+b] = r.second->sumOfPixels[b];
				}					
			}
			av.synchronize();
		}

		//DeviceUploadInitialAdjacents(this);
		{
			m_deviceData.amp_adj = new Concurrency::array<int, 1>(getRegionsCount()*MAX_REGION_ADJ);
			auto& v = *(m_deviceData.amp_adj);	
			auto av = v.view_as(v.extent);

			for (const auto& itCurrentRegion : m_regions)
			{
				int i=0;
				int currentRegionLabel = itCurrentRegion.first;
				int index = (currentRegionLabel-1)*MAX_REGION_ADJ;

				for (const auto& neighbor : itCurrentRegion.second->adjacentRegions)
				{					
					av[index+i] = neighbor->_regionLabel;
					++i;
				}
				if(i<MAX_REGION_ADJ)
				{
					//fill the remaining indexes = -1
					for(int c=i;c<MAX_REGION_ADJ;++c)
					{
						av[index+c] = -1;
					}
				}		
			}
			av.synchronize();
		}



		// initialize the keys
		{
			m_deviceData.amp_keysFilter = new Concurrency::array<int, 1>(getRegionsCount());
			auto& v = *(m_deviceData.amp_keysFilter);
			auto av = v.view_as(v.extent);

			int nRegions = m_image_height*m_image_width;

			for(int r=0; r<nRegions; ++r)
			{
				av[r] = r+1;
			}			
			av.synchronize();
		}


		// initialize the regions pixels count
		{		
			m_deviceData.amp_regionsPixelsCount = new Concurrency::array<int, 1>(getRegionsCount());
			auto& v = *(m_deviceData.amp_regionsPixelsCount);	
			auto av = v.view_as(v.extent);

			//int nRegions = m_image_height*m_image_width;

			//for(int r=0; r<nRegions; r++)
			//{
			//	v[r] = 1;
			//}

			amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), 1);
			av.synchronize();
		}


		{
			//amp_needsAdjRecomputation
			//amp_regionsAdjMinimums
			//amp_regionsAdjMinimumsLabels
			int nRegions = m_image_height*m_image_width;

			{
				m_deviceData.amp_regionsAdjMinimums = new Concurrency::array<float, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_regionsAdjMinimums);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), FLT_MAX);
				av.synchronize();
			}

			{
				m_deviceData.amp_regionsAdjMinimumsLabels = new Concurrency::array<int, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_regionsAdjMinimumsLabels);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), -1);
				av.synchronize();
			}

			{
				m_deviceData.amp_needsAdjRecomputation = new Concurrency::array<int, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_needsAdjRecomputation);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), 1);
				av.synchronize();
			}

			//TODO do amp HSEG optimization
			//DeviceThrustFill(m_deviceData.m_dev_needsAdjRecomputation, nRegions, 1);
		}


		{
			//m_dev_needsRecomputation
			//m_dev_regionsMinimums
			//m_dev_regionsMinimumsLabels

			int nRegions = m_image_height*m_image_width;

			{
				m_deviceData.amp_regionsMinimums = new Concurrency::array<float, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_regionsMinimums);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), FLT_MAX);
				av.synchronize();
			}

			{
				m_deviceData.amp_regionsMinimumsLabels = new Concurrency::array<int, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_regionsMinimumsLabels);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), -1);
				av.synchronize();
			}

			{
				m_deviceData.amp_needsRecomputation = new Concurrency::array<int, 1>(getRegionsCount());
				auto& v = *(m_deviceData.amp_needsRecomputation);	
				auto av = v.view_as(v.extent);
				amp_stl_algorithms::fill(amp_stl_algorithms::begin(av), amp_stl_algorithms::end(av), 1);
				av.synchronize();
			}

			//TODO do amp HSEG optimization
			//DeviceThrustFill(m_deviceData.m_dev_needsRecomputation, nRegions, 1);
		}		

		//cout<<"Done amp intialization"<<endl;
	}
	catch(std::exception e)
	{
		std::cout<<"Exception in AMP device initialize : "<<e.what() <<endl;
		throw e;
	}
	catch(...)
	{
		throw std::exception("Unknown exception in AMP device initialize");
	}
}


HSWO::Region* HSWO::Device_merge_regions(int region1Label, Region* r1, int region2Label, Region* r2)
{
	try
	{
		//we want to delete any one of the regions, and the deleted regions pixels will get the label of the remaining region
		//, so we will choose to keep region 1, and delete region 2
		// r2:  will be deallocated from memory at the end of the function		
		m_regions.erase(region2Label);

		// 1 - NOW, VERY IMPORTANT, replace region 2 ID at its adjacent regions by region 1 ID		
		for (const auto& currentRegion : r2->adjacentRegions)
		{
			currentRegion->adjacentRegions.erase(r2);

			if(currentRegion->_regionLabel != region1Label)
			{
				currentRegion->adjacentRegions.insert(r1);

				// 2- NOW, VERY IMPORTANT, make every neighbor to the deleted region a neighbor for the region 1 (except the new region itself)
				r1->adjacentRegions.insert(currentRegion);
			}		
		}

		/*
		// 2- NOW, VERY IMPORTANT, make every neighbor to the deleted region a neighbor for the region 1 (except the new region itself)
		for each(int currentRegionLabel in r2->adjacentRegionsLabels)
		{
		if(currentRegionLabel != region1Label)
		{
		m_regions[region1Label]->adjacentRegionsLabels.insert(currentRegionLabel);
		}
		}
		*/

		//DON'T forget to add the deleted region pixels to the new region
		r1->pixels_IDs.merge(r2->pixels_IDs);
		r1->_pixelCount += r2->_pixelCount;

		//update new region band sums
		r1->updateRegionBandsSumsByMergingWith(r2, this->m_nBands);

		//update the mean of the ONLY changed region region1Label
		//r1->calcMeanVector(&(this->m_deviceData.m_pRegionsMeans[(region1Label-1)*m_nBands]), m_pDataCube.get(), m_nBands,
		//m_image_width, m_image_height);

		///==================================================================
		//		printf("GPU merge_regions step A done ()\n");
		stream_to_log("step 2 done");

		int CurrentNRegions= getRegionsCount();
		HSWO::Region* lastMergedRegion = r1;

		//DeviceUpdateRegionMean
		{			
			//DeviceUpdateRegionMean(this, region1Label, lastMergedRegion);
			int nMeans = m_nBands;			
			int index = (region1Label-1)*m_nBands; 

			array_view<float, 1> avSource(nMeans, lastMergedRegion->sumOfPixels);
			auto avDest = m_deviceData.amp_regionsSums->view_as(m_deviceData.amp_regionsSums->extent);

			Concurrency::copy(avSource, amp_stl_algorithms::begin(avDest)+index);			

			//TODO , try sync_async for better performance
			avDest.synchronize();
		}

		stream_to_log("step 2.5 done");


		//--------------
		{
			auto av_regionsPixelsCount = m_deviceData.amp_regionsPixelsCount->view_as(m_deviceData.amp_regionsPixelsCount->extent);
			av_regionsPixelsCount[region1Label-1] = r1->_pixelCount;
			av_regionsPixelsCount.synchronize();
		}

		//		printf("GPU merge_regions step B done (update pixel count)\n");
		stream_to_log("step 3 done");


		// updating the adjacency of the merged region only ===========================================	
		//size_t updateRegionAdjSize = MAX_REGION_ADJ*sizeof(int);
		{
			int lastMergedRegion_OffsetInMatrix = (region1Label-1)*MAX_REGION_ADJ;
			auto& v = *(m_deviceData.amp_adj);	
			auto av = v.view_as(v.extent);

			//emptying all adjacencies on the device of the merged region
			//DeviceThrustResetAdjacencyOfMergedRegion(m_deviceData.m_dev_pRegionAdjancencies, lastMergedRegion_OffsetInMatrix);
			amp_stl_algorithms::fill_n(amp_stl_algorithms::begin(av)+lastMergedRegion_OffsetInMatrix, MAX_REGION_ADJ, -1);

			int aa=0;
			for (const auto& aaRegion : lastMergedRegion->adjacentRegions)
			{
				av[lastMergedRegion_OffsetInMatrix+aa] = aaRegion->_regionLabel;
				++aa;
			}

			av.synchronize();
		}

		//--------------------------------------
		//printf("GPU merge_regions step C done ()\n");
		stream_to_log("step 4 done");

		//loop only changed regions and update their adjacencies data on the device (the adjacent regions of the result merged region)		
		{
			auto& v = *(m_deviceData.amp_adj);	
			auto av = v.view_as(v.extent);

			for (const auto& lastRegionAdj : lastMergedRegion->adjacentRegions)
			{				
				HSWO::Region* currAdjacentRegion = lastRegionAdj;

				int currentAdjacentRegion_OffsetInMatrix = (lastRegionAdj->_regionLabel-1)*MAX_REGION_ADJ;

				//reset then copy current updated adjacency data to device								
				amp_stl_algorithms::fill_n(amp_stl_algorithms::begin(av)+currentAdjacentRegion_OffsetInMatrix, MAX_REGION_ADJ, -1);


				int aa=0;			
				for (const auto& aaRegion : currAdjacentRegion->adjacentRegions)
				{
					av[currentAdjacentRegion_OffsetInMatrix+aa] = aaRegion->_regionLabel;
					++aa;
				}					
			}	

			av.synchronize();
		}

		//		printf("GPU merge_regions step D done ()\n");
		stream_to_log("step 5 done");


		{
			Concurrency::array<int, 1>& KeysFilter = *m_deviceData.amp_keysFilter;
			auto av_KeysFilter = KeysFilter.view_as(KeysFilter.extent);

			Concurrency::array<int, 1>& needsAdjRecomputation= *m_deviceData.amp_needsAdjRecomputation;		
			auto av_needsAdjRecomputation = needsAdjRecomputation.view_as(needsAdjRecomputation.extent);

			Concurrency::array<int, 1>& regionsAdjMinimumsLabels= *m_deviceData.amp_regionsAdjMinimumsLabels;
			auto av_regionsAdjMinimumsLabels = regionsAdjMinimumsLabels.view_as(regionsAdjMinimumsLabels.extent);

			Concurrency::array<float, 1>& regionsAdjMinimums= *m_deviceData.amp_regionsAdjMinimums;
			auto av_regionsAdjMinimums = regionsAdjMinimums.view_as(regionsAdjMinimums.extent);

			///=========///
			Concurrency::array<int, 1>& needsRecomputation= *m_deviceData.amp_needsRecomputation;							
			auto av_needsRecomputation = needsRecomputation.view_as(needsRecomputation.extent);

			Concurrency::array<int, 1>& regionsMinimumsLabels= *m_deviceData.amp_regionsMinimumsLabels;	
			auto av_regionsMinimumsLabels = regionsMinimumsLabels.view_as(regionsMinimumsLabels.extent);

			Concurrency::array<float, 1>& regionsMinimums= *m_deviceData.amp_regionsMinimums;
			auto av_regionsMinimums = regionsMinimums.view_as(regionsMinimums.extent);

			//TODO (100) for optimization, try to just reset the changed region
			//, and let the reduction COPY the vector, not INPLACEing the vector
			amp_stl_algorithms::fill(amp_stl_algorithms::begin(av_regionsAdjMinimums), amp_stl_algorithms::end(av_regionsAdjMinimums),FLT_MAX);
			amp_stl_algorithms::fill(amp_stl_algorithms::begin(av_regionsMinimums), amp_stl_algorithms::end(av_regionsMinimums),FLT_MAX);

			//reset all computational flags for the new region and its adjacents
			{
				r1->resetAllComputationFlags(); //shouldn't be needed in GPU code

				int invalid = 1; int invalidLabel = -1;
				float maxf = FLT_MAX;
				av_needsAdjRecomputation[(region1Label-1)] = invalid;
				av_regionsAdjMinimumsLabels[(region1Label-1)]= invalidLabel;
				//TODO (100) for optimization, try to just reset the changed region
				//, and let the reduction COPY the vector, not INPLACEing the vector
				//av_regionsAdjMinimums[(region1Label-1)]= maxf;

				av_needsRecomputation[(region1Label-1)]= invalid;
				av_regionsMinimumsLabels[(region1Label-1)] = invalidLabel;
				//TODO (100) for optimization, try to just reset the changed region
				//, and let the reduction COPY the vector, not INPLACEing the vector
				//av_regionsMinimums[(region1Label-1)]= maxf;

				for (const auto& currentAdjRegion : r1->adjacentRegions)
				{		
					av_needsAdjRecomputation[currentAdjRegion->_regionLabel-1] = invalid;
					av_regionsAdjMinimumsLabels[currentAdjRegion->_regionLabel-1]= invalidLabel;
					//TODO (100) for optimization, try to just reset the changed region
					//, and let the reduction COPY the vector, not INPLACEing the vector
					//av_regionsAdjMinimums[currentAdjRegion->_regionLabel-1]= maxf;
				}
			}


#if SPECTRAL_DYNAMIC_PROGRAMMING
			//DeviceResetOthersBestRegionComputationFlagsFrom(this, region1Label, region2Label);
			//int nMaxRegions = m_image_height*m_image_width;	
			parallel_for_each(m_deviceData.amp_keysFilter->extent, [=, &KeysFilter
				,&needsRecomputation, &regionsMinimumsLabels, &regionsMinimums](index<1> idx)  restrict(amp)
			{ 
				int index = idx[0];
				int regionID = KeysFilter[index]-1;
				if(regionID == -2) return;
				if(regionsMinimumsLabels[regionID] == region1Label || regionsMinimumsLabels[regionID] == region2Label) 
				{	
					needsRecomputation[regionID] = 1;
					regionsMinimumsLabels[regionID] = -1;
					//TODO (100) for optimization, try to just reset the changed region
					//, and let the reduction COPY the vector, not INPLACEing the vector
					//regionsMinimums[regionID] = FLT_MAX;			
				}
			});
#endif
			//deallocate the deleted region from the memory
			{
				//reset all computational flags for the deleted region
				int invalid = 1; int invalidLabel = -1;
				float maxf = FLT_MAX;
				av_needsAdjRecomputation[(region2Label-1)] = invalid;
				av_regionsAdjMinimumsLabels[(region2Label-1)]= invalidLabel;
				//TODO (100) for optimization, try to just reset the changed region
				//, and let the reduction COPY the vector, not INPLACEing the vector
				//av_regionsAdjMinimums[(region2Label-1)]= maxf;

				av_needsRecomputation[(region2Label-1)]= invalid;
				av_regionsMinimumsLabels[(region2Label-1)] = invalidLabel;
				//TODO (100) for optimization, try to just reset the changed region
				//, and let the reduction COPY the vector, not INPLACEing the vector
				//av_regionsMinimums[(region2Label-1)]= maxf;

				//TODO very important, delete this Region* object if m_regions is heap allocated
				//delete r2;	
				r2->region_deleted = true;

				//update invalid region IDs
				m_deviceData.amp_keysFilter->view_as(m_deviceData.amp_keysFilter->extent)[region2Label-1] = -1;				
			}

			av_KeysFilter.synchronize();

			av_needsAdjRecomputation.synchronize();
			av_regionsAdjMinimumsLabels.synchronize();
			av_regionsAdjMinimums.synchronize();

			av_needsRecomputation.synchronize();
			av_regionsMinimumsLabels.synchronize();
			av_regionsMinimums.synchronize();
		}

		//		printf("GPU merge_regions step E done (invalidation of computations)\n");
		return r1;
	}
	catch(std::runtime_error re)
	{
		std::cout<<"Runtime exception in AMP merge regions : "<<re.what() <<endl;
		throw re;
	}
	catch(std::exception e)
	{
		std::cout<<"Exception in AMP merge regions : "<<e.what() <<endl;
		throw e;
	}
	catch(...)
	{		
		throw std::exception("Unknown exception in AMP merge regions");
	}

	return nullptr;
}


#endif