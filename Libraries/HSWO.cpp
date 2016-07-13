#include "HSWO.h"
#include <cmath>
#include <climits>
#include "./Libraries/definitions.h"
#include "./Libraries/Logger.h"
#include "./Libraries/DebuggingUtilities.h"

#include "HSWODevice.h"

#include "myglobals.h"

using namespace HyperSpectralToolbox;
using namespace std;

#define PI 3.1415926535897932384626433832795f
#define MY_FLT_MAX         FLT_MAX        /* max value */


//float HSWO::all_sums_of_pixels[1024][220];
int HSWO::time_count = 0;
bool HSWO::silentMode = false;

#pragma region Constructors and Distructors

void HSWO::Region::resetAdjComputationFlag()
{
	_bestAdjLabel = -1;
	_bestAdjDissim = MY_FLT_MAX;	//TODO: OOPS!, VIP this should have been MAX_FLT from the beginning !, 
									// discovered and fixed after the mem blocking bug
	_bestAdj = nullptr;
	_needsAdjRecomputation = true;
}

void HSWO::Region::resetRegionComputationFlag()
{
	_bestRegionLabel = -1;
	_bestRegionDissim = MY_FLT_MAX;	//TODO: OOPS!, VIP this should have been MAX_FLT from the beginning !, 
									// discovered and fixed after the mem blocking bug
	_bestRegion = nullptr;
	_needsRegionRecomputation =true;
}

HSWO::HSWO(int nCols, int nRows, int nBands, float spclustWeight, int min_no_of_clusters, shared_ptr<float> spDataCube, unsigned int threadsPerBlock) :
m_nBands(nBands)
{
	//HSWO::time_count = 0;
    //m_nBands = nBands;
	_max_nRegions = nCols * nRows;
	m_deviceData._temp_all_dissims_size = (_max_nRegions) * (_max_nRegions);

    m_pDataCube = spDataCube;
    m_deviceData.m_dev_pDataCube = NULL;

    m_min_no_of_clusters=min_no_of_clusters;

    m_labeled_segmented_image = new int[nCols * nRows];
    m_deviceData.m_labeled_segmented_image = new int[nCols * nRows];
    
    m_tempPixelVector3 = new float[m_nBands];

    m_minDataValue = new float[m_nBands];
    m_meanDataValue = new float[m_nBands];
    m_varDataValue = new float[m_nBands];

    m_image_width = nCols;
    m_image_height = nRows;

    for (int y = 0; y < m_image_height; y++)
    {
        for (int x = 0; x < m_image_width; x++)
        {
            int index = x + y * m_image_width;

            m_labeled_segmented_image[index] = index;
            m_deviceData.m_labeled_segmented_image[index] = index;
        }
    }

    //m_pRegionsSums = new float[nCols* nRows* m_nBands];
    
    //allocate device memory for regions adjacencies and means	
        //allocate device memory for regions adjacencies and means	
    //I commented this file only for testing large files > 400x400, where this allocation is very high and fails with bad_alloc
	// and replaced it with a dummy line
	// this solution works when I use only one core, but doen't work when I use more than one core in RHSEG 
	// (still don't know why)
	m_deviceData.m_pRegionAdjancencies = new int[nCols* nRows* MAX_REGION_ADJ]; 
	//m_deviceData.m_pRegionAdjancencies = new int[4*4*10];  //dummy

    m_deviceData.m_pRegionKeys = new int[nCols* nRows];
    m_deviceData.m_pRegionPixelsCount = new int[nCols* nRows];
    //m_deviceData.m_pRegionsMeans = new float[nCols* nRows* m_nBands];

    m_deviceData.m_threadsPerBlock = threadsPerBlock;	

	this->m_deviceData.allRegions = nullptr;

	_spclustWeight = spclustWeight;
	
	canWork = true;
	myThreadID = -1;
	nextToTryGPU = nullptr;	

	m_regionsObjects = nullptr;
	m_regionsObjects = new vector<Region>(m_image_width*m_image_height);
	cout <<"(HSWO::Region) c++ class size = " << sizeof(HSWO::Region) <<endl;
}

HSWO::~HSWO(void)
{
    if (m_labeled_segmented_image != NULL)
        delete[] m_labeled_segmented_image;
    m_labeled_segmented_image = NULL;

    if (m_tempPixelVector3 != NULL)
        delete[] m_tempPixelVector3;
    m_tempPixelVector3 = NULL;

    if (m_minDataValue != NULL)
        delete[] m_minDataValue;
    m_minDataValue = NULL;

    if (m_meanDataValue != NULL)
        delete[] m_meanDataValue;
    m_meanDataValue = NULL;

    if (m_varDataValue != NULL)
        delete[] m_varDataValue;
    m_varDataValue = NULL;

    /*
    if (m_pRegionsSums  != NULL)
        delete[] m_pRegionsSums ;
    m_pRegionsSums  = NULL;
    */

	auto itCurrentRegion = m_regions.begin();
    for(;itCurrentRegion != m_regions.end(); itCurrentRegion++)
    {
		//TODO very important, delete this Region* object if m_regions is heap allocated
        //delete itCurrentRegion->second;
    }

    /// Device host data deallocation ===========================

    if (m_deviceData.m_labeled_segmented_image != NULL)
        delete[] m_deviceData.m_labeled_segmented_image;
    m_deviceData.m_labeled_segmented_image = NULL;

    //if (m_deviceData.m_pRegionsMeans  != NULL)
        //delete[] m_deviceData.m_pRegionsMeans ;
    //m_deviceData.m_pRegionsMeans  = NULL;
    
    if (m_deviceData.m_pRegionKeys  != NULL)
        delete[] m_deviceData.m_pRegionKeys ;
    m_deviceData.m_pRegionKeys  = NULL;
    
    if (m_deviceData.m_pRegionAdjancencies  != NULL)
        delete[] m_deviceData.m_pRegionAdjancencies ;
    m_deviceData.m_pRegionAdjancencies  = NULL;

    /*
	hash_map<int, Region*>::iterator DeviceItCurrentRegion = m_deviceData.m_regions.begin();
    for(;DeviceItCurrentRegion  != m_deviceData.m_regions.end(); DeviceItCurrentRegion ++)
    {
        delete DeviceItCurrentRegion->second;
    }
	*/

    if (m_deviceData.m_pRegionPixelsCount  != NULL)
        delete[] m_deviceData.m_pRegionPixelsCount ;
    m_deviceData.m_pRegionPixelsCount  = NULL;

	if(this->m_deviceData.allRegions != nullptr)
	{
		delete[] this->m_deviceData.allRegions;
	}

	if(m_regionsObjects != nullptr)
	{
		delete m_regionsObjects;
	}
}

#pragma endregion Constructors and Distructors

inline float HSWO::measureDissimilarity(Region* region1, Region* region2)
//inline float HSWO::measureDissimilarity(Region* region1, int region1Label, Region* region2, int region2Label)
{
//	float* mean_vectorR1 = m_tempPixelVector1;
//	memset(mean_vectorR1, 0.0, m_nBands*sizeof(float));

//	float* mean_vectorR2 = m_tempPixelVector2;
//	memset(mean_vectorR2, 0.0, m_nBands*sizeof(float));

//	r1.calcMeanVector(mean_vectorR1, m_pDataCube, m_nBands, m_image_width, m_image_height);
//	r2.calcMeanVector(mean_vectorR2, m_pDataCube, m_nBands, m_image_width, m_image_height);
	//clock_t startClock2 = clock();

    float temp_sum1 = 0;
	float temp_val = 0;
    
	float r1nPixles = region1->_pixelCount;
    float r2nPixles = region2->_pixelCount;

    //float* r1Sums = region1->sumOfPixels.get();
	float* r1Sums = region1->sumOfPixels;
    //float* r2Sums = region2->sumOfPixels.get();
	float* r2Sums = region2->sumOfPixels;

    //float maxDiff = 0.0;
	
	//int BlockWidth = 20;
	//int nBlocks = HARD_CODED_N_BANDS/BlockWidth;
	

//	for (int B_ind = 0; B_ind < nBlocks; ++B_ind)
	{
		//int start_offset = B_ind*BlockWidth;
		//int end_offset = start_offset + BlockWidth;
		//for (int band = start_offset; band < (end_offset); ++band)
		for (int band = 0; band < (m_nBands); ++band)
		{
			//m_tempPixelVector3[band] = (r1Sums[band])/(r1nPixles) - (r2Sums[band])/(r2nPixles);		
			temp_val = (r1Sums[band])/(r1nPixles) - (r2Sums[band])/(r2nPixles);		
			//temp_sum1 += (m_tempPixelVector3[band]) * (m_tempPixelVector3[band]);
			temp_sum1 += (temp_val) * (temp_val);

			//float diff = fabs(r1Sums[band]/r1nPixles - r2Sums[band]/r2nPixles);
			//if(diff > maxDiff) maxDiff = diff;
			//cout << band << endl;
		}
	}

	//cout << "fin B" << endl;

	float d = temp_sum1*( (r1nPixles * r2nPixles) 	/  (r1nPixles + r2nPixles) );

    float dissim = sqrt( d );

	//clock_t endClock2 = clock();

	//HSWO::time_count += (endClock2 - startClock2);

    //float dissim = sqrt(temp_sum1) / float(m_nBands);
    //float dissim = maxDiff;
    return dissim;
}

inline void HSWO::measureMeanVector_for_merging_candidates(Region & r1, Region & r2,
                                                           float*& mean_vector_for_merging_candids)
{
    /*
    list<Pixel>::iterator it1 = r1.pixels.begin();
    list<Pixel>::iterator it2 = r2.pixels.begin();

    memset(m_tempPixelVector3, 0.0, m_nBands);

    int nPixels1 = r1.pixels.size();
    int nPixels2 = r2.pixels.size();

    {
    while(it1 != r1.pixels.end())
    {
    int x = it->x;
    int y = it->y;

    {
    for(int i=0; i< m_nBands; i++) meanVector3[i] += pDataCube[x + y*image_width];
    }

    it1++;
    }

    //proceed to the 2nd region
    while(it2 != r2.pixels.end())
    {
    int x = it->x;
    int y = it->y;

    {
    for(int i=0; i< m_nBands; i++) meanVector3[i] += pDataCube[x + y*image_width];
    }

    it2++;
    }

    for(int i=0; i<m_nBands; i++) meanVector3[i] = meanVector3[i]/float(nPixels1+nPixels2);
    }

    mean_vector_for_merging_candids = mean_vector3;
    */
}

inline void HSWO::Region::updateRegionBandsSumsByMergingWith(Region* region2, int nBands)
//inline void HSWO::Region::updateRegionBandsSumsByMergingWith(float* ptrSums, int region1Label, int region2Label)
{	
	//auto ptr1 = sumOfPixels.get();
	auto ptr1 = sumOfPixels;
	//auto ptr2 = region2->sumOfPixels.get();
	auto ptr2 = region2->sumOfPixels;

    for (int band = 0; band < nBands; ++band)
    {
        ptr1[band] += ptr2[band];
        //ptrSums[(region1Label-1)*_nBands+band] += ptrSums[(region2Label-1)*_nBands+band];
    }
    
    /*
    int nRegionPixels = pixels_IDs.size();
    for each (int it in pixels_IDs)
    {
        //int x = it->x;
        //int y = it->y;
        int index = it;
        {
            for (int band = 0; band < nBands; band++)
            {
                meanVector[band] += pDataCube[(index-1)+band*nImagePixels];
            }
        }
    }
    for (int band = 0; band < _nBands; band++)
    {
        meanVector[band] = meanVector[band] / float(nRegionPixels);
    }
    */
}

void HSWO::computeDataStatistics()
{
    float value;
    int nPixels = m_image_height * m_image_width;
    for(int band=0; band<m_nBands; band++) m_minDataValue[band]= MY_FLT_MAX;
    memset(m_meanDataValue, 0, m_nBands*sizeof(float));
    memset(m_varDataValue, 0, m_nBands*sizeof(float));

    float* sums = new float[m_nBands];
    float* sums_sq = new float[m_nBands];

    memset(sums, 0, m_nBands*sizeof(float));
    memset(sums_sq, 0, m_nBands*sizeof(float));

    for (int y = 0; y < m_image_height; y++)
    {
        for (int x = 0; x < m_image_width; x++)
        {
            int index = x + y * m_image_width;

            for (int  band = 0; band < m_nBands; band++)
            {
                value = (float) m_pDataCube.get()[(index-1) + band * nPixels];
                if (m_minDataValue[band] > value)
                    m_minDataValue[band] = value;
                sums[band] += value;
                sums_sq[band] += value * value;
            }
        }
    }

    for (int band = 0; band < m_nBands; band++)
    {
        m_meanDataValue[band] = sums[band] / ((float) nPixels);
        m_varDataValue[band] = (sums[band] * sums[band]) / ((float) nPixels);
        m_varDataValue[band] = sums_sq[band] - m_varDataValue[band];
        m_varDataValue[band] = m_varDataValue[band] / (((float) nPixels) - 1.0);
        //		if (((params.debug > 0) && (params.nbands < 15)) || (params.debug > 2))
        //		{
        //			params.log_fs << "For band " << band << ", minimum = " << min_stat[band];
        //			params.log_fs << ", mean = " << mean_stat[band];
        //			params.log_fs << ", and variance = " << var_stat[band] << endl;
        //		}
    }

    delete[] sums;
    delete[] sums_sq;
}

inline bool HSWO::getPixelBoundryRegions(const int x, const int y, const int pixelLabel, int& topleftRegionLabel, int& topRegionLabel, int& toprightRegionLabel
                                         , int& leftRegionLabel, int& rightRegionLabel, int& bottomleftRegionLabel, int& bottomRegionLabel, int& bottomrightRegionLabel)
{
    /*
    int index;

    bool isBoundryPixel = false;	
    //topleft
    if(x > 0 && y > 0 )
    {
        index = (x-1) + (y-1)*m_image_width;
        topleftRegionLabel = m_labeled_segmented_image[index];
        if(topleftRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //top
    if(y > 0)
    {
        index = (x) + (y-1)*m_image_width;
        topRegionLabel = m_labeled_segmented_image[index];
        if(topRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //topright
    if(x < (m_image_width-1) && y > 0)
    {
        index = (x+1) + (y-1)*m_image_width;
        toprightRegionLabel = m_labeled_segmented_image[index];
        if(toprightRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //left
    if(x > 0)
    {
        index = (x-1) + (y)*m_image_width;
        leftRegionLabel = m_labeled_segmented_image[index];
        if(leftRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //right
    if(x < (m_image_width-1) )
    {
        index = (x+1) + (y)*m_image_width;
        rightRegionLabel = m_labeled_segmented_image[index];
        if(rightRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //bottomleft
    if(x > 0 && y < (m_image_height-1) )
    {
        index = (x-1) + (y+1)*m_image_width;
        bottomleftRegionLabel = m_labeled_segmented_image[index];
        if(bottomleftRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //bottom
    if(y < (m_image_height-1))
    {
        index = (x) + (y+1)*m_image_width;
        bottomRegionLabel = m_labeled_segmented_image[index];
        if(bottomRegionLabel != pixelLabel) isBoundryPixel = true;
    }

    //bottomright
    if(x < (m_image_width-1) && y < (m_image_height-1))
    {
        index = (x+1) + (y+1)*m_image_width;
        bottomrightRegionLabel = m_labeled_segmented_image[index];
        if(bottomrightRegionLabel != pixelLabel) isBoundryPixel = true;
    }
    
    return isBoundryPixel;
    */
    throw std::exception("Not Implemented yet");
}



//void HSWO::merge_regions(int region1Label, int region2Label)
HSWO::Region* HSWO::merge_regions(Region* r1, int region1Label, Region* r2, int region2Label)
{	
    /*
    //we want to delete any one of the regions, and the deleted regions pixels will get the label of the remaining region
    //, so we will choose to keep region 1, and delete region 2
    Region* r2 = m_regions[region2Label];  //wiill be deallocated from memory at the end of the function
    m_regions.erase(region2Label);

    map<int, Region*>::iterator it = m_regions.begin();

    // NOW, VERY IMPORTANT, make every neighbor to the deleted region a neighbor for the new region (except the new region itself)
    while(it != m_regions.end())
    {
        it->second->adjacentRegionsLabels.erase(region2Label);
        int regionID = it->first;
        
        if(regionID != region1Label)
        {
            it->second->adjacentRegionsLabels.insert(region1Label);
        }		
        it++;
    }
    */

    //we want to delete any one of the regions, and the deleted regions pixels will get the label of the remaining region
    //, so we will choose to keep region 1, and delete region 2
    //Region* r1 = m_regions[region1Label];
    //Region* r2 = m_regions[region2Label];  //will be deallocated from memory at the end of the function
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
    for each(int currentRegion in r2->adjacentRegionsLabels)
    {
        if(currentRegion != region1Label)
        {
            m_regions[region1Label]->adjacentRegionsLabels.insert(currentRegion);
        }		
    }
    */

    //DON'T forget to add the deleted region pixels to the new region
  //  for each(const int& itPixels in r2->pixels_IDs)
  //  {		
  //      r1->pixels_IDs.push_back(itPixels);
		//r1->_pixelCount += 1;
  //  }

	r1->pixels_IDs.merge(r2->pixels_IDs);
	r1->_pixelCount += r2->_pixelCount;

    //update new region band sums
    r1->updateRegionBandsSumsByMergingWith(r2, this->m_nBands);
    //r1->updateRegionBandsSumsByMergingWith(m_pRegionsSums ,region1Label, region2Label);

    //deallocate the deleted region from the memory
	//TODO very important, delete this Region* object if m_regions is heap allocated
    //delete r2;	

	r2->region_deleted = true;

    r1->resetAllComputationFlags();

    return r1;
}

void HSWO::doStep()
{
	{
		/*for each(auto reg in m_regions)
		{
			if(reg.first != reg.second->_regionLabel)
			{
				printf(" **************** lalalallalalalala **************\n");
			}
		}*/

		/*
		printf("dumping the cluster image\n" );
		auto spClustersLabelsImage = getClustersLabelsImage();

		for (int y=0; y<m_image_height; y++)			
		{			
			string str;
			for (int x=0; x<m_image_width; x++)
			{				
				str += ts(spClustersLabelsImage.get()[x  + y * m_image_width]) + "\t";
			}
			printf(str.c_str());
			printf("\n");
		}
		printf("\n");
		*/
	}

    //for every region: 
    // compute dissimilarity with all its neighbors
    // select the pair of regions with the minimum dissimilarity overall the computed pairs
    //
	{
		int minDissimRegion1Label, minDissimRegion2Label;
		float minDissim = MY_FLT_MAX;
		Region *ptrR1=nullptr, *ptrR2=nullptr;

		// region growing step
		// compute the dissim among all pairs of adjacent regions

		for each(const auto& itCurrentRegion in m_regions)
		{
			float bestAdjDissim = MY_FLT_MAX;
			int bestAdjLabel=-1;
			Region *bestAdj=nullptr;		

			if(itCurrentRegion.second->_needsAdjRecomputation)
			{
				for each(const auto& adjacentRegion in itCurrentRegion.second->adjacentRegions)
				{						
					//if((itCurrentRegion.first) < adjacentRegion->_regionLabel)
					{
						float dissim = measureDissimilarity(itCurrentRegion.second ,adjacentRegion);
						if(dissim < minDissim)
						{
							minDissim = dissim;
							minDissimRegion1Label = itCurrentRegion.first;
							minDissimRegion2Label = adjacentRegion->_regionLabel ;
							ptrR1 = itCurrentRegion.second;
							ptrR2 = adjacentRegion;
						}				

						if(dissim < bestAdjDissim)
						{
							bestAdjDissim = dissim;
							bestAdjLabel = adjacentRegion->_regionLabel;
							bestAdj = adjacentRegion;
						}
					}
				}			

				itCurrentRegion.second->_needsAdjRecomputation = false; //minor potential optimization, make this bool flag update in a separate loop outside
				itCurrentRegion.second->_bestAdjDissim = bestAdjDissim;
				itCurrentRegion.second->_bestAdjLabel = bestAdjLabel;
				itCurrentRegion.second->_bestAdj = bestAdj;
			}
			else
			{
				if(itCurrentRegion.second->_bestAdjDissim < minDissim)
				{
					minDissim = itCurrentRegion.second->_bestAdjDissim;
					minDissimRegion1Label = itCurrentRegion.first;
					minDissimRegion2Label = itCurrentRegion.second->_bestAdjLabel ;
					ptrR1 = itCurrentRegion.second;
					ptrR2 = itCurrentRegion.second->_bestAdj;
				}				
			}
		}

		_last_merge_dissim_value = minDissim;

		{
			//now we can merge the minimum dissimilarity pair	
			Region* spMergedRegion = merge_regions(ptrR1, minDissimRegion1Label, ptrR2, minDissimRegion2Label);

			//make all those who are adjacents to the old deleted regions (r1 and r2), need adj recomputation
			for each(const auto& currentAdjRegion in spMergedRegion->adjacentRegions)
			{		
				currentAdjRegion->resetAdjComputationFlag();
			}

#if SPECTRAL_DYNAMIC_PROGRAMMING
			// make all regions -that had the old deleted regions (r1 and r2) as their best regions- needs region recomputation
			for each(const auto& itCurrentRegion in m_regions)
			{			
				int regionLabel = itCurrentRegion.first;
				if (itCurrentRegion.second->_bestRegionLabel == minDissimRegion1Label ||
					itCurrentRegion.second->_bestRegionLabel == minDissimRegion2Label)
				{					
					itCurrentRegion.second->resetRegionComputationFlag();
				}
			}
#endif
		}
	}
	// spectral clustering step
    // compute the dissim among all pairs of non-adjacent regions
    // them merge all pairs that has dissim value less than or equal to _last_merge_dissim_value
	vector<Region>& regionsObjects = (*m_regionsObjects);
    
	//if(m_regions.size() > m_min_no_of_clusters && (m_regions.size() <= 378))	
	if(m_regions.size() > m_min_no_of_clusters)	
    {
        int minDissimRegion1Label, minDissimRegion2Label;
        float minDissim = _last_merge_dissim_value*_spclustWeight;
        Region *ptrR1=nullptr, *ptrR2=nullptr;

        //for all non-adjacent regions
        //std::list<int> merge_ids;
        //std::list<int> merge_labels1;
        //std::list<int> merge_labels2;
        //std::list<Region*> merge_Regions1;
        //std::list<Region*> merge_Regions2;
		//int iterations = 0, nocomp = 0;

		//memory blocking to enhance cache performance and reduce the impact of stride-n 
		// access of the first dimension of regions
		
		//int BlockWidth = 256;
		int BlockWidth = m_image_height*m_image_width;		
		int nBlocks = m_image_height*m_image_width/BlockWidth;
		
		//	for (int B_ind = 0; B_ind < nBlocks; ++B_ind)
		//int start_offset = B_ind*BlockWidth;
		//int end_offset = start_offset + BlockWidth;
		//for (int band = start_offset; band < (end_offset); ++band)

		for (int B_dim1_ind = 0; B_dim1_ind < nBlocks; ++B_dim1_ind)
		{
			int start_offset_dim1 = B_dim1_ind*BlockWidth;
			int end_offset_dim1 = start_offset_dim1 + BlockWidth;

			for (int B_dim2_ind = 0; B_dim2_ind < nBlocks; ++B_dim2_ind)
			{
				int start_offset_dim2 = B_dim2_ind*BlockWidth;
				int end_offset_dim2 = start_offset_dim2 + BlockWidth;


				//for each(const auto& itCurrentRegion1 in m_regions)												
				//for (auto itCurrentRegion1__ = regionsObjects.begin() ; itCurrentRegion1__ != regionsObjects.end();itCurrentRegion1__++)
				for (auto itCurrentRegion1__ = regionsObjects.begin()+start_offset_dim1
					; itCurrentRegion1__ != regionsObjects.begin()+end_offset_dim1;itCurrentRegion1__++)
				{
					//iterations++;      						
					Region* itCurrentRegion1 = &(*itCurrentRegion1__);
					if(itCurrentRegion1->region_deleted) continue;

					//int region1Label = itCurrentRegion1.first;
					int region1Label = itCurrentRegion1->_regionLabel;
					//cout<<region1Label << " : ";

					float bestRegionDissim = MY_FLT_MAX;
					int bestRegionLabel=-1;
					Region *bestRegion=nullptr;	
					//const auto& region1End = itCurrentRegion1.second->adjacentRegions.end();
					const auto& region1End = itCurrentRegion1->adjacentRegions.end();

					//set<int> tempR1Adjs; 
					//std::copy(itCurrentRegion1.second->adjacentRegionsLabels.begin(), itCurrentRegion1.second->adjacentRegionsLabels.end(), tempR1Adjs.begin());
#if SPECTRAL_DYNAMIC_PROGRAMMING
					if(itCurrentRegion1->_needsRegionRecomputation)
#endif
					{				
						//for each(const auto& itCurrentRegion2 in m_regions)
						//for (auto itCurrentRegion2__ = regionsObjects.begin() ; itCurrentRegion2__ != regionsObjects.end();itCurrentRegion2__++)
						for (auto itCurrentRegion2__ = regionsObjects.begin()+start_offset_dim2
							; itCurrentRegion2__ != regionsObjects.begin()+end_offset_dim2;itCurrentRegion2__++)
						{
							Region* itCurrentRegion2 = &(*itCurrentRegion2__);
							if(itCurrentRegion2->region_deleted) continue;
							//int region2Label = itCurrentRegion2.first;
							int region2Label = itCurrentRegion2->_regionLabel;
							//cout<<", " << region2Label;
							//if (region1Label < region2Label)
							{
								if(region1Label == region2Label) continue;
								//if (itCurrentRegion1.second->adjacentRegions.find(itCurrentRegion2.second) == region1End)							
								if (itCurrentRegion1->adjacentRegions.find(itCurrentRegion2) == region1End)							
								{
									//we have a non-adjacent region
									//mark the pair for merge if dissim is less than or equal the _last_merge_dissim_value from
									// the region growing step
									//float dissim = measureDissimilarity(itCurrentRegion1.second ,itCurrentRegion2.second);
									float dissim = measureDissimilarity(itCurrentRegion1,itCurrentRegion2);
									if(dissim < minDissim)
									{
										//mark the pair for merging, merge step is done later
										//merge_ids.push_back((region1Label-1)+ this->m_image_width*(region2Label-1));
										//merge_labels1.push_back(region1Label);
										//merge_labels2.push_back(region2Label);
										//merge_Regions1.push_back(itCurrentRegion1.second);
										//merge_Regions2.push_back(itCurrentRegion2.second);

										minDissim = dissim;
										minDissimRegion1Label = region1Label;
										minDissimRegion2Label = region2Label;
										//ptrR1 = itCurrentRegion1.second;
										ptrR1 = itCurrentRegion1;
										//ptrR2 = itCurrentRegion2.second;							
										ptrR2 = itCurrentRegion2;							
									}

									if(dissim < bestRegionDissim)
									{
										bestRegionDissim = dissim;
										bestRegionLabel = itCurrentRegion2->_regionLabel;
										bestRegion = itCurrentRegion2;
									}
								}
								else
								{
									//we find that region2 is an adjacent
									//remove it from the temp search vector (OPT)
									//itCurrentRegion1.second->adjacentRegionsLabels.erase(region2Label);						
								}
							}					
						}
						//cout<<endl;

						//TODO: VIP, this block is invalid when using memory blocking, mustn't be used !, fixed in the next block
						//itCurrentRegion1->_needsRegionRecomputation = false; 
						//itCurrentRegion1->_bestRegionDissim = bestRegionDissim;
						//itCurrentRegion1->_bestRegionLabel = bestRegionLabel;
						//itCurrentRegion1->_bestRegion = bestRegion;

						if((B_dim2_ind == nBlocks - 1))
						{
							itCurrentRegion1->_needsRegionRecomputation = false; 
						}
						if(bestRegionDissim < itCurrentRegion1->_bestRegionDissim)
						{
							itCurrentRegion1->_bestRegionDissim = bestRegionDissim;
							itCurrentRegion1->_bestRegionLabel = bestRegionLabel;
							itCurrentRegion1->_bestRegion = bestRegion;
						}
					}
#if SPECTRAL_DYNAMIC_PROGRAMMING
					else
					{	
						//nocomp++;
						if(itCurrentRegion1->_bestRegionDissim < minDissim)
						{
							minDissim = itCurrentRegion1->_bestRegionDissim;
							minDissimRegion1Label = itCurrentRegion1->_regionLabel;
							minDissimRegion2Label = itCurrentRegion1->_bestRegionLabel ;
							ptrR1 = itCurrentRegion1;
							ptrR2 = itCurrentRegion1->_bestRegion;

							//printf("reg %d, from label = %d, best reg = %d\n", minDissimRegion2Label, ptrR2->_regionLabel, ptrR2->_bestRegionLabel );
						}
					}			
#endif
				}
			}
		}

        //cout << "\n=================="<<endl;
        // merge all pairs with dissim less than or equal the _last_merge_dissim_value from
        // the region growing step
        //for each (int merge_id in merge_ids)
        //for each (int merge_id in merge_labels1)
        //{
        //	//int region1Label = merge_id % this->m_image_width +1;
        //	//int region2Label = merge_id / this->m_image_width +1;
        //	int region1Label = merge_id;
        //	int region2Label = merge_labels2.front();

        //	merge_regions(merge_Regions1.front(), region1Label, merge_Regions2.front(), region2Label);
        //	//merge_labels1.pop_front();
        //	merge_labels2.pop_front();
        //	merge_Regions1.pop_front();
        //	merge_Regions2.pop_front();
        //}		
        
		//printf("iterations = %d, nocomp = %d, percent = %f\n",iterations, nocomp, (float)nocomp/(float)iterations*100.0f);        
		//MERGE!
        if(ptrR1 != nullptr)
        {			
			//printf("doing spclust merge\n");
		    Region* spMergedRegion = merge_regions(ptrR1, minDissimRegion1Label, ptrR2, minDissimRegion2Label);

			//make all those who are adjacents to the old deleted regions (r1 and r2), need adj recomputation
			for each(const auto& currentAdjRegion in spMergedRegion->adjacentRegions)
			{		
				currentAdjRegion->resetAdjComputationFlag();
			}

#if SPECTRAL_DYNAMIC_PROGRAMMING
			// make all regions -that had the old deleted regions (r1 and r2) as their best regions- needs region recomputation
            for each(const auto& itCurrentRegion in m_regions)
            {			
                int regionLabel = itCurrentRegion.first;
                if (itCurrentRegion.second->_bestRegionLabel == minDissimRegion1Label ||
                    itCurrentRegion.second->_bestRegionLabel == minDissimRegion2Label)
                {
                    itCurrentRegion.second->resetRegionComputationFlag();
                }
            }
            
			//// recompute dissim of every other region again to this new region, of course store best dissims info for each
   //         // region
			//for each(auto& itCurrentRegion in regionsObjects)
			for (auto itCurrentRegion = regionsObjects.begin()
				; itCurrentRegion != regionsObjects.end();itCurrentRegion++)

			{			
				int regionLabel = itCurrentRegion->_regionLabel;

				if(itCurrentRegion->_needsRegionRecomputation == false)
				{
					if(itCurrentRegion->_regionLabel == spMergedRegion->_regionLabel) continue;
					float dissim = measureDissimilarity(&(*itCurrentRegion) ,spMergedRegion);
					if(dissim < itCurrentRegion->_bestRegionDissim)
					{
						itCurrentRegion->_bestRegionDissim = dissim;
						itCurrentRegion->_bestRegionLabel = spMergedRegion->_regionLabel;
						itCurrentRegion->_bestRegion = spMergedRegion;
					}					
				}				
			}			
#endif
        }	
    }	
}

void HSWO::InitializeHost()
{
    /// Build initial set of regions (each pixel is a region)
	//m_regions.clear();
	vector<Region>& regionsObjects = *(m_regionsObjects);

    int nPixels = m_image_height * m_image_width;
    for(int y = 0;y < m_image_height;y++){
        for(int x = 0;x < m_image_width;x++){

            auto createNewRegion = [&](int newRegionLabel)->Region* {
                regionsObjects[newRegionLabel-1].Construct(newRegionLabel, m_nBands);

				//Region* ptrNewRegion = new Region(newRegionLabel, m_nBands);
				Region* ptrNewRegion = &(regionsObjects[newRegionLabel-1]);
                this->m_regions[newRegionLabel] = ptrNewRegion ;
				
                
				ptrNewRegion->pixels_IDs.push_back(newRegionLabel);
				ptrNewRegion->_pixelCount += 1;
                //m_labeled_segmented_image[RegionLabel] = RegionLabel;

                //calculate initial regions means
                for (int band = 0; band < m_nBands; band++)
                {
                    //(ptrNewRegion->sumOfPixels.get())[band] = (m_pDataCube.get())[(newRegionLabel-1)+band*m_image_height*m_image_width];
					(ptrNewRegion->sumOfPixels)[band] = (m_pDataCube.get())[(newRegionLabel-1)+band*m_image_height*m_image_width];
                    //this->m_pRegionsSums[(RegionLabel-1)*m_nBands + band] = (m_pDataCube.get())[(RegionLabel-1)+band*m_image_height*m_image_width];
                }
                return ptrNewRegion;
            };
            // new logic, I decided to make the regions and pixels labels "1 based" to easily match 
            // the logic of RHSEG author 
            int RegionLabel = (x + y * m_image_width) + 1;
            
            Region* ptrCurrentRegion = nullptr;
            if(m_regions.find(RegionLabel) == m_regions.end())
            {
                ptrCurrentRegion = createNewRegion(RegionLabel);
            }
            else
            {
                ptrCurrentRegion = m_regions[RegionLabel];
            }
            


            // add boundary regions (which are the MAX_REGION_ADJ surrounding pixels) for this current region(x,y)
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
}

void HSWO::startRun()
{
    InitializeHost();
    while(m_regions.size() > m_min_no_of_clusters)
    {
        doStep();		
    }
}

void HSWO::continueRun()
{
    while(m_regions.size() > m_min_no_of_clusters)
    {
        doStep();		
    }
}

void HSWO::start_CPU_Run_Then_GPU()
{
	cout << "Started Mixed CPU / GPU thread # "<< myThreadID <<endl;
	
	InitializeHost();
	while(m_regions.size() > m_min_no_of_clusters && canWork)
	{
		doStep();		
	}
	
	cout << "Paused in the middle of Mixed CPU / GPU thread # "<< myThreadID << " with ("<<m_regions.size()<<") regions\n";

	canWork = true;

	if(m_regions.size() > m_min_no_of_clusters)
	{	
		TO_GPU_HSWO(this->m_deviceData.m_threadsPerBlock);

		while(m_regions.size() > m_min_no_of_clusters)
		{		
			Device_doStep();
		}
		DeviceExit(this);
	}	
	
	cout << "Finished Mixed CPU / GPU thread # "<< myThreadID ;

	if(nextToTryGPU != nullptr)
	{
		cout<< ", trying next thread # "<< nextToTryGPU->myThreadID <<"...\n";	
		nextToTryGPU->canWork = false;
	}
	cout<<"\n";
}

shared_ptr<int> HSWO::getClustersLabelsImage()
{
    shared_ptr<int> image(new int[m_image_width*m_image_height]);

	
    for(int i=0; i<m_image_height*m_image_width;i++) 
    {
        image.get()[i]= -1;
		 //image.get()[i]= 1;
    }
    set<int> allpixels;

	
    for each(auto& itCurrentRegion in m_regions)
    {
        for each(int it in itCurrentRegion.second->pixels_IDs)
        {
            int pixelLabel = it;
            if(allpixels.find(pixelLabel-1) != allpixels.end()) throw exception("HSWO: Must not happen");
            allpixels.insert(pixelLabel-1);
            image.get()[pixelLabel-1]=itCurrentRegion.first;			
        }		
    }
	
    return image;
}

int HSWO::getRegionsCount()
{
    return this->m_regions.size();
}


void HSWO::InitializeHost( HSWO& pre_segmentation )
{

}


