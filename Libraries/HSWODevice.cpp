#include "HSWO.h"
#include <cmath>
#include <climits>
#include <algorithm>
#include <iostream>
#include "./Libraries/definitions.h"
#include "./Libraries/Logger.h"
#include "./Libraries/DebuggingUtilities.h"

#include "HSWODevice.h"

#include "myglobals.h"

using namespace HyperSpectralToolbox;
using namespace std;

#define PI 3.1415926535897932384626433832795
#define MY_FLT_MAX         3.402823466e+38F        /* max value */

void HyperSpectralToolbox::HSWO::TO_GPU_HSWO(unsigned int threadsPerBlock)
{
#if USE_CUDA
	m_deviceData.m_threadsPerBlock = threadsPerBlock;	
	this->m_deviceData.allRegions = nullptr;
	
	try
	{
		DeviceInitializeData(this);
		
		int nPixels = m_image_height * m_image_width;
		//this->m_deviceData.allRegions.reset(new Region* [nPixels]);		
		this->m_deviceData.allRegions = new Region* [nPixels];		

		/// iterate all current regions :
		
		for each(auto itRegion in m_regions)
		{
			auto currentRegionLabel = itRegion.first;
			auto currentRegion = itRegion.second;
			this->m_deviceData.allRegions[currentRegionLabel-1] = currentRegion;
		}
		

		DeviceUploadInitialMeans(this);

		int N = m_image_height*m_image_width;

		for each(auto itCurrentRegion in m_regions)
		{
			int i=0;
			int regionLabel = itCurrentRegion.first;

			for each(auto neighbor in itCurrentRegion.second->adjacentRegions)
			{
				int index = (regionLabel-1)*MAX_REGION_ADJ + i;
				m_deviceData.m_pRegionAdjancencies[index] = neighbor->_regionLabel;

				i++;
			}
			if(i<MAX_REGION_ADJ)
			{
				//fill the remaining indexes = -1
				for(int c=i;c<MAX_REGION_ADJ;c++)
				{
					m_deviceData.m_pRegionAdjancencies[(itCurrentRegion.first-1)*MAX_REGION_ADJ + c]	= -1;
				}
			}		
		}

		DeviceUploadInitialAdjacents(this);

		// initialize the keys
		int nRegions = m_image_height*m_image_width;

		for(int r=0; r<nRegions; r++)
		{
			m_deviceData.m_pRegionKeys[r] = -1;
		}		
		for each(auto itRegion in m_regions)
		{			
			auto currentRegionLabel = itRegion.first;			
			m_deviceData.m_pRegionKeys[currentRegionLabel-1] = currentRegionLabel;			
		}

		size_t regionKeysSize = nRegions*sizeof(int);
		cudaMemcpy(m_deviceData.m_dev_RegionKeys, m_deviceData.m_pRegionKeys, regionKeysSize , cudaMemcpyHostToDevice) ;
		cudaError_t err = cudaGetLastError(); 
		if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}

		{		
			// initialize the regions pixels count
			int nRegions = m_image_height*m_image_width;

			for(int r=0; r<nRegions; r++)
			{
				m_deviceData.m_pRegionPixelsCount[r] = 1;
			}
			for each(auto itRegion in m_regions)
			{			
				auto currentRegionLabel = itRegion.first;			
				m_deviceData.m_pRegionPixelsCount[currentRegionLabel-1] = itRegion.second->_pixelCount;			
			}

			size_t regionsPixelsCountSize = nRegions*sizeof(int);
			cudaMemcpy(m_deviceData.m_dev_RegionPixelsCount, m_deviceData.m_pRegionPixelsCount, regionsPixelsCountSize , cudaMemcpyHostToDevice) ;

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

		{
			//m_dev_needsAdjRecomputation
			//m_dev_regionsAdjMinimums
			//m_dev_regionsAdjMinimumsLabels

			int nRegions = m_image_height*m_image_width;

			DeviceThrustFillFloat(m_deviceData.m_dev_regionsAdjMinimums, nRegions, FLT_MAX);
			DeviceThrustFill(m_deviceData.m_dev_regionsAdjMinimumsLabels, nRegions, -1);
			DeviceThrustFill(m_deviceData.m_dev_needsAdjRecomputation, nRegions, 1);

			for each(auto itRegion in m_regions)
			{			
				auto currentRegionLabel = itRegion.first;			
				auto currentRegion = itRegion.second;

				int currentRegionNeedsRecomputation = currentRegion->_needsAdjRecomputation == true ? 1: 0;
				if(currentRegionNeedsRecomputation == 0)
				{
					cudaMemcpy(m_deviceData.m_dev_needsAdjRecomputation+(currentRegionLabel-1), &(currentRegionNeedsRecomputation), sizeof(int), cudaMemcpyHostToDevice) ;
					cudaMemcpy(m_deviceData.m_dev_regionsAdjMinimumsLabels+(currentRegionLabel-1), &(currentRegion->_bestAdjLabel), sizeof(int), cudaMemcpyHostToDevice) ;
					cudaMemcpy(m_deviceData.m_dev_regionsAdjMinimums+(currentRegionLabel-1), &(currentRegion->_bestAdjDissim), sizeof(float), cudaMemcpyHostToDevice) ;
				}
			}

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

		{
			//m_dev_needsRecomputation
			//m_dev_regionsMinimums
			//m_dev_regionsMinimumsLabels

			int nRegions = m_image_height*m_image_width;

			DeviceThrustFillFloat(m_deviceData.m_dev_regionsMinimums, nRegions, FLT_MAX);
			DeviceThrustFill(m_deviceData.m_dev_regionsMinimumsLabels, nRegions, -1);
			DeviceThrustFill(m_deviceData.m_dev_needsRecomputation, nRegions, 1);

			for each(auto itRegion in m_regions)
			{			
				auto currentRegionLabel = itRegion.first;			
				auto currentRegion = itRegion.second;
				
				int currentRegionNeedsRecomputation = currentRegion->_needsRegionRecomputation == true ? 1: 0;
				if(currentRegionNeedsRecomputation == 0)
				{
					cudaMemcpy(m_deviceData.m_dev_needsRecomputation+(currentRegionLabel-1), &(currentRegionNeedsRecomputation), sizeof(int), cudaMemcpyHostToDevice) ;
					cudaMemcpy(m_deviceData.m_dev_regionsMinimumsLabels+(currentRegionLabel-1), &(currentRegion->_bestRegionLabel), sizeof(int), cudaMemcpyHostToDevice) ;
					cudaMemcpy(m_deviceData.m_dev_regionsMinimums+(currentRegionLabel-1), &(currentRegion->_bestRegionDissim), sizeof(float), cudaMemcpyHostToDevice) ;
				}
			}

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

	}
	catch(std::exception e)
	{
		std::cout<<"Exception in TO GPU HSWO : "<<e.what() <<endl;
	}
	catch(...)
	{
		throw std::exception("Unknown exception in TO GPU HSWO !");		
		std::cout<<"Unknown exception in TO GPU HSWO  : "<<endl;
	}
#endif
}

#if USE_CUDA
void HSWO::Device_doStep()
{
	try
	{
		/*
		int min_Region_to_RegionAdj_index =  DeviceCalcAllDissims(this, m_threadsPerBlock);		
		int regionLabel = (min_Region_to_RegionAdj_index  / MAX_REGION_ADJ) + 1 ; //for one based labels
		Region* region1 = m_regions[regionLabel];

		int adjRegionIndex = min_Region_to_RegionAdj_index  % MAX_REGION_ADJ;
		int adjRegionLabel = m_deviceData.m_pRegionAdjancencies[min_Region_to_RegionAdj_index];
		
		if(adjRegionIndex > (region1->adjacentRegions.size()-1) )
		throw exception("HSWO: Must not happen adj region does not exist label = -1");

		*/
		{
			int regionLabel = -1, adjRegionLabel = -1; float minDissim=FLT_MAX;
			DeviceCalcAllDissims(this, m_deviceData.m_threadsPerBlock, regionLabel, adjRegionLabel, minDissim);

			Region* region1 = this->m_deviceData.allRegions[regionLabel-1];		
			Region* region2 = this->m_deviceData.allRegions[adjRegionLabel-1];

			//stream_to_log("regionLabel= "<<regionLabel <<", n_adjs="<<m_regions[regionLabel]->adjacentRegions.size()
			//<<"\t"<< "adjRegionIndex= " <<adjRegionIndex <<"\t"<< "adjLabel="<<adjRegionLabel);		

			//printf("region %d and %d are min with dissim = %f\n", regionLabel, adjRegionLabel, minDissim);
			//printf("GPU dostep A done\n");
			stream_to_log("step 1 done");

			//	_last_merge_dissim_value = minDissim;

			//now we can merge the minimum dissimilarity pair
			Device_merge_regions(regionLabel, region1, adjRegionLabel, region2);
		}
		//========================================================================================
		//========================================================================================
		if(m_regions.size() > m_min_no_of_clusters)
		{
			clock_t startClock2 = clock();

			int region1Label = -1, region2Label = -1; float minDissim=FLT_MAX;
			DeviceCalcRegionDissims(this, m_deviceData.m_threadsPerBlock, region1Label, region2Label, minDissim);			

#if SPECTRAL_DYNAMIC_PROGRAMMING
			DeviceThrustFill(m_deviceData.m_dev_needsRecomputation, _max_nRegions, 0); //for dynamic programming optimization to work, this step is necessery like the one inside the loop in CPU code
#endif

			if(region1Label != -1 && region2Label != -1)
			{
				Region* region1 = this->m_deviceData.allRegions[region1Label-1];		
				Region* region2 = this->m_deviceData.allRegions[region2Label-1];

				Device_merge_regions(region1Label, region1, region2Label, region2);
				//cout <<"spectral merge at " <<m_regions.size()<<endl;

#if SPECTRAL_DYNAMIC_PROGRAMMING
				DeviceRecomputeOthersBestRegionsDissimilarityTo(this, region1Label);
#endif
			}


			clock_t endClock2 = clock();
			HSWO::time_count += (endClock2 - startClock2);

			//printf("\nreached here 4\n");
		}

		//printf("GPU dostep B done (merging)\n");
		//======================================
		
		//printf("GPU dostep C done (copy updated key)\n");
		/*
		int regionCounter=0;
		for each(auto regionIt in m_regions)
		{
			m_deviceData.m_pRegionKeys[regionCounter] = regionIt.first;
			regionCounter++;
		}
		
		size_t regionIDsSize = m_regions.size()*sizeof(int);
		cudaMemcpy(m_deviceData.m_dev_RegionKeys, m_deviceData.m_pRegionKeys, regionIDsSize , cudaMemcpyHostToDevice) ;
		*/

		//DeviceThrustRemove(m_deviceData.m_dev_RegionKeys,m_regions.size(), adjRegionLabel) ;
		
		//=======================================
		/*
		auto dumpHswoImage = [&]()->void{
			stream_to_log("dumping the GPU hswo image");
			//for (int b =0; b<1;b++)
			auto device_cluster_map = Device_getClustersLabelsImage();
			{		
				for (int y=0; y<m_image_height; y++)			
				{			
					string str;
					for (int x=0; x<m_image_width; x++)
					{				
						str += ts(device_cluster_map.get()[x  + y * m_image_width]) + "\t";
					}
					stream_to_log(str.c_str());
				}
				stream_to_log("");
			}		
		};
		*/
	}
	catch(std::runtime_error re)
	{
		stream_to_log(string("General runtime error exception (Device_doStep) : ") + string(re.what()));
		std::cout<<"General runtime error exception (Device_doStep) : "<<re.what() <<endl;
		throw re;
	}
	catch(std::exception e)
	{
		stream_to_log(string("General exception (Device_doStep) : ") + e.what());
		std::cout<<"General exception (Device_doStep) : "<<e.what() <<endl;
		throw e;
	}
	catch(...)
	{
		stream_to_log("Unknown exception (Device_doStep)");
		std::cout<<"Unknown exception (Device_doStep)" <<endl;
		throw std::exception("Unknown exception (Device_doStep)");
	}
}
#endif


#if USE_CUDA
void HSWO::Device_Initialize()
{
	//int threadsPerBlock = 512;
	//int blocksPerGrid =	1;
	//VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C);
	try
	{
		vector<Region>& regionsObjects = (*m_regionsObjects);
		DeviceInitializeData(this);
		
		/// Build initial set of regions (each pixel is a region)
		int nPixels = m_image_height * m_image_width;

		//this->m_deviceData.allRegions.reset(new Region* [nPixels]);		
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

					ptrNewRegion->pixels_IDs.push_back(newRegionLabel);
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

		DeviceUploadInitialMeans(this);


		int N = m_image_height*m_image_width;

		for each(auto itCurrentRegion in m_regions)
		{
			int i=0;
			int regionID = itCurrentRegion.first;

			for each(auto neighbor in itCurrentRegion.second->adjacentRegions)
			{
				int index = (itCurrentRegion.first-1)*MAX_REGION_ADJ + i;
				m_deviceData.m_pRegionAdjancencies[index] = neighbor->_regionLabel;

				i++;
			}
			if(i<MAX_REGION_ADJ)
			{
				//fill the remaining indexes = -1
				for(int c=i;c<MAX_REGION_ADJ;c++)
				{
					m_deviceData.m_pRegionAdjancencies[(itCurrentRegion.first-1)*MAX_REGION_ADJ + c]	= -1;
				}
			}		
		}

		DeviceUploadInitialAdjacents(this);

		// initialize the keys
		int nRegions = m_image_height*m_image_width;

		for(int r=0; r<nRegions; r++)
		{
			m_deviceData.m_pRegionKeys[r] = r+1;
		}

		size_t regionKeysSize = nRegions*sizeof(int);
		cudaMemcpy(m_deviceData.m_dev_RegionKeys, m_deviceData.m_pRegionKeys, regionKeysSize , cudaMemcpyHostToDevice) ;
		cudaError_t err = cudaGetLastError(); 
		if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		
		{		
			// initialize the regions pixels count
			int nRegions = m_image_height*m_image_width;

			for(int r=0; r<nRegions; r++)
			{
				m_deviceData.m_pRegionPixelsCount[r] = 1;
			}

			size_t regionsPixelsCountSize = nRegions*sizeof(int);
			cudaMemcpy(m_deviceData.m_dev_RegionPixelsCount, m_deviceData.m_pRegionPixelsCount, regionsPixelsCountSize , cudaMemcpyHostToDevice) ;

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

		{
			//m_dev_needsAdjRecomputation
			//m_dev_regionsAdjMinimums
			//m_dev_regionsAdjMinimumsLabels
			
			int nRegions = m_image_height*m_image_width;

			DeviceThrustFillFloat(m_deviceData.m_dev_regionsAdjMinimums, nRegions, FLT_MAX);
			DeviceThrustFill(m_deviceData.m_dev_regionsAdjMinimumsLabels, nRegions, -1);
			DeviceThrustFill(m_deviceData.m_dev_needsAdjRecomputation, nRegions, 1);

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

		{
			//m_dev_needsRecomputation
			//m_dev_regionsMinimums
			//m_dev_regionsMinimumsLabels

			int nRegions = m_image_height*m_image_width;

			DeviceThrustFillFloat(m_deviceData.m_dev_regionsMinimums, nRegions, FLT_MAX);
			DeviceThrustFill(m_deviceData.m_dev_regionsMinimumsLabels, nRegions, -1);
			DeviceThrustFill(m_deviceData.m_dev_needsRecomputation, nRegions, 1);

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}

	}
	catch(...)
	{
		std::cout<<"Exception in (Device_Initialize)"<<endl;
		throw std::exception("Exception in (Device_Initialize)");
	}

}
#endif

void HSWO::Device_run()
{
	cout << "Started Whole GPU thread # "<< myThreadID <<endl;

	Device_Initialize();
	while(getRegionsCount() > m_min_no_of_clusters)
	{
		Device_doStep();		
	}
	DeviceExit(this);

	cout << "Finished Whole GPU thread # "<< myThreadID ;
	if(nextToTryGPU != nullptr)
	{
		cout <<", trying next thread # "<< nextToTryGPU->myThreadID <<"...\n";
		nextToTryGPU->canWork = false;   //TODO Why ???
	}
	cout<<"\n";
}

// this function executes HSEG on the CPU for regions created by on device while executing GPU RHSEG
// this function I use ONLY for top RHSEG level when running HSEG for the final image,
// since usually the regions at this stage are small, CPU will be faster.
/*
void HSWO::RunOnHost()
{	
	while(getRegionsCount() > m_min_no_of_clusters)
	{
		{	
			int minDissimRegion1Label, minDissimRegion2Label;
			float minDissim = MY_FLT_MAX;
			Region *ptrR1=nullptr, *ptrR2=nullptr;

			// region growing step
			// compute the dissim among all pairs of adjacent regions

			for each(auto itCurrentRegion in m_regions)
			{
				for each(auto neighbor in itCurrentRegion.second->adjacentRegions)
				{
					if((itCurrentRegion.first) < neighbor->_regionLabel)
					{				
						float dissim = measureDissimilarity(itCurrentRegion.second ,neighbor);
						if(dissim < minDissim)
						{
							minDissim = dissim;
							minDissimRegion1Label = itCurrentRegion.first;
							minDissimRegion2Label = neighbor->_regionLabel ;
							ptrR1 = itCurrentRegion.second;		
							ptrR2 = neighbor;
						}				
					}
				}		
			}

			_last_merge_dissim_value = minDissim;
			
			{				
				//now we can merge the minimum dissimilarity pair	
				//merge_regions(ptrR1, minDissimRegion1Label, ptrR2, minDissimRegion2Label);

				///I was forced to unfold the merge_regions code here

				m_regions.erase(minDissimRegion2Label);

				// 1 - NOW, VERY IMPORTANT, replace region 2 ID at its adjacent regions by region 1 ID
				for each(auto currentRegion in ptrR2->adjacentRegions)
				{
					currentRegion->adjacentRegions.erase(ptrR2);		

					if(currentRegion->_regionLabel != minDissimRegion1Label)
					{
						currentRegion->adjacentRegions.insert(ptrR1);

						// 2- NOW, VERY IMPORTANT, make every neighbor to the deleted region a neighbor for the region 1 (except the new region itself)
						ptrR1->adjacentRegions.insert(currentRegion);
					}		
				}

				//DON'T forget to add the deleted region pixels to the new region
				for each(int itPixels in ptrR2->pixels_IDs)
				{		
					ptrR1->pixels_IDs.push_back(itPixels);
					ptrR1->_pixelCount += 1;
				}

				//update new region band sums
				ptrR1->updateRegionBandsSumsByMergingWith(ptrR2, this->m_nBands);

				//deallocate the deleted region from the memory
				delete ptrR2;

				ptrR1->resetAllComputationFlags();
			}
		}

		// spectral clustering step
		// compute the dissim among all pairs of non-adjacent regions
		// them merge all pairs that has dissim value less than or equal to _last_merge_dissim_value
		if(m_regions.size() > m_min_no_of_clusters)	
		{
			int minDissimRegion1Label, minDissimRegion2Label;
			float minDissim = _last_merge_dissim_value*_spclustWeight;
			Region *ptrR1=nullptr, *ptrR2=nullptr;

			for each(auto itCurrentRegion1 in m_regions)        
			{
				int region1Label = itCurrentRegion1.first;

				float bestRegionDissim = MY_FLT_MAX;
				int bestRegionLabel=-1;
				Region *bestRegion=nullptr;	
				auto region1End = itCurrentRegion1.second->adjacentRegions.end();

				//if(itCurrentRegion1.second->_needsRegionRecomputation)
				if(true)
				{				
					for each(auto itCurrentRegion2 in m_regions)
					{
						int region2Label = itCurrentRegion2.first;
						{
							if(itCurrentRegion1.first == itCurrentRegion2.first) continue;
							if (itCurrentRegion1.second->adjacentRegions.find(itCurrentRegion2.second) == region1End)							
							{
								float dissim = measureDissimilarity(itCurrentRegion1.second ,itCurrentRegion2.second);
								if(dissim < minDissim)
								{
									minDissim = dissim;
									minDissimRegion1Label = region1Label;
									minDissimRegion2Label = region2Label;
									ptrR1 = itCurrentRegion1.second;
									ptrR2 = itCurrentRegion2.second;							
								}

								//if(dissim < bestRegionDissim)
								//{
								//	bestRegionDissim = dissim;
								//	bestRegionLabel = itCurrentRegion2.second->_regionLabel;
								//	bestRegion = itCurrentRegion2.second;
								//}
							}
						}
					}

					//itCurrentRegion1.second->_needsRegionRecomputation = false;
					//itCurrentRegion1.second->_bestRegionDissim = bestRegionDissim;
					//itCurrentRegion1.second->_bestRegionLabel = bestRegionLabel;
					//itCurrentRegion1.second->_bestRegion = bestRegion;
				}
				else
				{	
					//if(itCurrentRegion1.second->_bestRegionDissim < minDissim)
					//{
					//	minDissim = itCurrentRegion1.second->_bestRegionDissim;
					//	minDissimRegion1Label = itCurrentRegion1.first;
					//	minDissimRegion2Label = itCurrentRegion1.second->_bestRegionLabel ;
					//	ptrR1 = itCurrentRegion1.second;
					//	ptrR2 = itCurrentRegion1.second->_bestRegion;
					//}
				}						
			}

			//MERGE!
			if(ptrR1 != nullptr)
			{			
				//now we can merge the minimum dissimilarity pair	
				//merge_regions(ptrR1, minDissimRegion1Label, ptrR2, minDissimRegion2Label);

				printf("doing spclust merge after RHSEG\n");
				
				///I was forced to unfold the merge_regions code here
				//Region* spMergedRegion = merge_regions(ptrR1, minDissimRegion1Label, ptrR2, minDissimRegion2Label);

				{
					m_regions.erase(minDissimRegion2Label);

					// 1 - NOW, VERY IMPORTANT, replace region 2 ID at its adjacent regions by region 1 ID
					for each(auto currentRegion in ptrR2->adjacentRegions)
					{
						currentRegion->adjacentRegions.erase(ptrR2);		

						if(currentRegion->_regionLabel != minDissimRegion1Label)
						{
							currentRegion->adjacentRegions.insert(ptrR1);

							// 2- NOW, VERY IMPORTANT, make every neighbor to the deleted region a neighbor for the region 1 (except the new region itself)
							ptrR1->adjacentRegions.insert(currentRegion);
						}		
					}

					//DON'T forget to add the deleted region pixels to the new region
					for each(int itPixels in ptrR2->pixels_IDs)
					{		
						ptrR1->pixels_IDs.push_back(itPixels);
						ptrR1->_pixelCount += 1;
					}

					//update new region band sums
					ptrR1->updateRegionBandsSumsByMergingWith(ptrR2, this->m_nBands);

					//deallocate the deleted region from the memory
					delete ptrR2;

					ptrR1->resetAllComputationFlags();
				}

			}	
		}	
	}	
}
*/
	
#if USE_CUDA
HSWO::Region* HSWO::Device_merge_regions(int region1Label, Region* r1, int region2Label, Region* r2)
{

	try
	{
		//we want to delete any one of the regions, and the deleted regions pixels will get the label of the remaining region
		//, so we will choose to keep region 1, and delete region 2
		// r2:  will be deallocated from memory at the end of the function
		m_regions.erase(region2Label);

		// 1 - NOW, VERY IMPORTANT, replace region 2 ID at its adjacent regions by region 1 ID		
		for each(const auto& currentRegion in r2->adjacentRegions)
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
		
		HSWO::Region* lastMergedRegion = r1;
		DeviceUpdateRegionMean(this, region1Label, lastMergedRegion);

		stream_to_log("step 2.5 done");

		//--------------
		m_deviceData.m_pRegionPixelsCount[region1Label-1] = r1->_pixelCount;

		size_t regionPixelsCountSize = 1*sizeof(int);
		cudaMemcpyAsync(m_deviceData.m_dev_RegionPixelsCount+(region1Label-1), &(m_deviceData.m_pRegionPixelsCount[region1Label-1]), regionPixelsCountSize , cudaMemcpyHostToDevice) ;

//		printf("GPU merge_regions step B done (update pixel count)\n");
		stream_to_log("step 3 done");


		// updating the adjacency of the merged region only ===========================================	
		//size_t updateRegionAdjSize = MAX_REGION_ADJ*sizeof(int);
		{
			int lastMergedRegion_OffsetInMatrix = (region1Label-1)*MAX_REGION_ADJ;

			int aa=0;
			for each(const auto& aaRegion in lastMergedRegion->adjacentRegions)
			{
				m_deviceData.m_pRegionAdjancencies[lastMergedRegion_OffsetInMatrix+aa] = aaRegion->_regionLabel;
				++aa;
			}

			size_t updateRegionAdjSize = aa*sizeof(int);
			
			/*
			//TODO OPT: can be removed this if statement
			if(aa<MAX_REGION_ADJ)
			{
				//fill the remaining indexes = -1
				for(int c=aa;c<MAX_REGION_ADJ;c++)
				{
					m_deviceData.m_pRegionAdjancencies[lastMergedRegion_OffsetInMatrix+c] = -1;					
				}
			}
			*/
			
			
			//emptying all adjacencies on the device of the merged region
			DeviceThrustResetAdjacencyOfMergedRegion(m_deviceData.m_dev_pRegionAdjancencies, lastMergedRegion_OffsetInMatrix);

			//copy last merged region updated adjacency data to device
			cudaMemcpyAsync(m_deviceData.m_dev_pRegionAdjancencies+lastMergedRegion_OffsetInMatrix
				, &(m_deviceData.m_pRegionAdjancencies[lastMergedRegion_OffsetInMatrix]), updateRegionAdjSize, cudaMemcpyHostToDevice) ;			
			
			/*
			//copy last merged region updated adjacency data to device
			cudaMemcpy(m_deviceData.m_dev_pRegionAdjancencies+lastMergedRegion_OffsetInMatrix
				, &(m_deviceData.m_pRegionAdjancencies[lastMergedRegion_OffsetInMatrix]), updateRegionAdjSize, cudaMemcpyHostToDevice) ;	
			*/

			cudaError_t err = cudaGetLastError(); 
			if(err != cudaError::cudaSuccess) {std::cout << std::string(cudaGetErrorString(err)) << endl;}
		}
		//--------------------------------------
//		printf("GPU merge_regions step C done ()\n");
		stream_to_log("step 4 done");

		//loop only changed regions and update their adjacencies data on the device (the adjacent regions of the result merged region)		
		{
			for each(const auto& lastRegionAdj in lastMergedRegion->adjacentRegions)
			{				
				HSWO::Region* currAdjacentRegion = lastRegionAdj;

				int currentAdjacentRegion_OffsetInMatrix = (lastRegionAdj->_regionLabel-1)*MAX_REGION_ADJ;
				
				/*
				int nAdjacentsOfAdjacent = currAdjacentRegion->adjacentRegionsLabels.size();	
				
				// foreach element in host structure, replace each removed region label with the new merged label
				std::replace(m_deviceData.m_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix
					, m_deviceData.m_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix+nAdjacentsOfAdjacent
					, adjRegionLabel, regionLabel);
					*/
				/*
				stream_to_log("after replacing "<<adjRegionLabel<<" with " << regionLabel <<":");
				stream_to_log("adjofadj = ");
				for(int i=0; i<MAX_REGION_ADJ;i++)
				{
					stream_to_log_sameLine(", "<< *(m_deviceData.m_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix+i))
				}

				stream_to_log("\n");
				*/
					
				
				int aa=0;			
				for each(const auto& aaRegion in currAdjacentRegion->adjacentRegions)
				{
					m_deviceData.m_pRegionAdjancencies[currentAdjacentRegion_OffsetInMatrix+aa] = aaRegion->_regionLabel;
					++aa;
				}
					
				size_t updateRegionAdjSize = aa*sizeof(int);
				

				/*
				if(aa<MAX_REGION_ADJ)
				{
					//fill the remaining indexes = -1
					for(int c=aa;c<MAX_REGION_ADJ;c++)
					{
						m_deviceData.m_pRegionAdjancencies[currentAdjacentRegion_OffsetInMatrix+c] = -1;					
					}
				}					
				*/

				
				//reset then copy current updated adjacency data to device								
				DeviceThrustFill(m_deviceData.m_dev_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix, MAX_REGION_ADJ, -1);			

				cudaMemcpyAsync(m_deviceData.m_dev_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix
					, &(m_deviceData.m_pRegionAdjancencies[currentAdjacentRegion_OffsetInMatrix]), updateRegionAdjSize, cudaMemcpyHostToDevice) ;			
				

				//DeviceThrustReplace(m_deviceData.m_dev_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix, nAdjacentsOfAdjacent, adjRegionLabel, regionLabel);				

				/*
				//copy current updated adjacency data to device
				cudaMemcpy(m_deviceData.m_dev_pRegionAdjancencies+currentAdjacentRegion_OffsetInMatrix
					, &(m_deviceData.m_pRegionAdjancencies[currentAdjacentRegion_OffsetInMatrix]), updateRegionAdjSize, cudaMemcpyHostToDevice) ;			
				*/

				//REMOVE ALL CUDEGETLASTERROR FOR PERFORMANCE !!
				//cudaError_t err = cudaGetLastError(); 
				//if(err != cudaError::cudaSuccess) {
				//	std::string serr = std::string(cudaGetErrorString(err));
				//	std::cout << serr << endl;
				//	stream_to_log("Last CUDA error in step 5 " << serr);
				//}
			}	
		}
		
//		printf("GPU merge_regions step D done ()\n");
		stream_to_log("step 5 done");


		//reset all computational flags for the new region and its adjacents
		{
			r1->resetAllComputationFlags(); //shouldn't be needed in GPU code			

			int invalid = 1; int invalidLabel = -1;
			float maxf = FLT_MAX;
			cudaMemcpyAsync(m_deviceData.m_dev_needsAdjRecomputation+(region1Label-1), &(invalid), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimumsLabels+(region1Label-1), &(invalidLabel), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimums+(region1Label-1), &(maxf), sizeof(float), cudaMemcpyHostToDevice) ;			

			invalid = 1; invalidLabel = -1;
			maxf = FLT_MAX;
			cudaMemcpyAsync(m_deviceData.m_dev_needsRecomputation+(region1Label-1), &(invalid), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsMinimumsLabels+(region1Label-1), &(invalidLabel), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsMinimums+(region1Label-1), &(maxf), sizeof(float), cudaMemcpyHostToDevice) ;			
			
			for each(const auto& currentAdjRegion in r1->adjacentRegions)
			{		
				int invalid = 1; int invalidLabel = -1;
				float maxf = FLT_MAX;
				cudaMemcpyAsync(m_deviceData.m_dev_needsAdjRecomputation+(currentAdjRegion->_regionLabel-1), &(invalid), sizeof(int), cudaMemcpyHostToDevice) ;
				cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimumsLabels+(currentAdjRegion->_regionLabel-1), &(invalidLabel), sizeof(int), cudaMemcpyHostToDevice) ;
				cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimums+(currentAdjRegion->_regionLabel-1), &(maxf), sizeof(float), cudaMemcpyHostToDevice) ;			
			}

			//TODO add flag reset logic here for non-adjacent regions
		}

#if SPECTRAL_DYNAMIC_PROGRAMMING
		DeviceResetOthersBestRegionComputationFlagsFrom(this, region1Label, region2Label);
#endif
		//deallocate the deleted region from the memory
		{
			//reset all computational flags for the deleted region
			int invalid = 1; int invalidLabel = -1;
			float maxf = FLT_MAX;
			cudaMemcpyAsync(m_deviceData.m_dev_needsAdjRecomputation+(region2Label-1), &(invalid), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimumsLabels+(region2Label-1), &(invalidLabel), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsAdjMinimums+(region2Label-1), &(maxf), sizeof(float), cudaMemcpyHostToDevice) ;			

			invalid = 1; invalidLabel = -1;
			maxf = FLT_MAX;
			cudaMemcpyAsync(m_deviceData.m_dev_needsRecomputation+(region2Label-1), &(invalid), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsMinimumsLabels+(region2Label-1), &(invalidLabel), sizeof(int), cudaMemcpyHostToDevice) ;
			cudaMemcpyAsync(m_deviceData.m_dev_regionsMinimums+(region2Label-1), &(maxf), sizeof(float), cudaMemcpyHostToDevice) ;			

			//TODO very important, delete this Region* object if m_regions is heap allocated
			//delete r2;	
			r2->region_deleted = true;

			//update invalid region IDs
			m_deviceData.m_pRegionKeys[region2Label-1] = -1;

			size_t regionKeysSize = 1*sizeof(int);
			cudaMemcpyAsync(m_deviceData.m_dev_RegionKeys+(region2Label-1), &(m_deviceData.m_pRegionKeys[region2Label-1]), regionKeysSize , cudaMemcpyHostToDevice) ;

		}

		cudaStreamSynchronize(0);
//		printf("GPU merge_regions step E done (invalidation of computations)\n");

		return r1;
	}
	catch(...)
	{
		std::cout<<"Exception in device merge regions"<<endl;
		throw std::exception("Exception in device merge regions");
	}

	return nullptr;
}

#endif