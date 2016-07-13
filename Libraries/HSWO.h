#pragma once
#include <memory>
#include <hash_map>
#include <set>
#include <list>
#include "./definitions.h"
#include <iostream>

#if USE_AMP
#include <amp.h>
#endif

#if USE_CUDA

#ifndef CUDA_SAFE_CALL
#define CUDA_SAFE_CALL(code) {                                              \
    cudaError err = code;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }
#endif

#endif

namespace HyperSpectralToolbox
{
	class HSWO
	{
	public:		
		//static float all_sums_of_pixels[1024][220];
		static int time_count;
		static bool silentMode;
		class Region
		{
		public:
			float sumOfPixels[HARD_CODED_N_BANDS];
			//std::shared_ptr<float> sumOfPixels;
			//float* sumOfPixels;
			std::list<int> pixels_IDs;  //potenial optimization: try vector instead, may be faster like herb sutter mentioned in : https://channel9.msdn.com/Events/Build/2014/2-661 (min 50)
										//update: I tried to remove the list merging in merge_regions() [the only refernce to pixels_IDs] and only changed run time from 508 to 506 !
										//seems no tuning is needed
			int _pixelCount;
			//std::set<int> adjacentRegionsLabels;			
			std::set<Region*> adjacentRegions;

			Region()
			{
			}

			Region(int regionLabel, int nBands)
			{
				//_nBands = nBands;
				//sumOfPixels.reset(new float [nBands]); //minor potential optimization: move this bands sums arrays outside the regions, and 
														// and do the dissim calculations on that outside matrix instead of accessing them indivdullay 
														// form inside the regions collection(may not be benefetial due to ununoforma nature of data access in RHSEG)

				//sumOfPixels = &(HSWO::all_sums_of_pixels[regionLabel -1][0]);
				//sumOfPixels = nullptr;
				//sumOfPixels = new float[_nBands];
				_regionLabel = regionLabel;
				_pixelCount = 0;

				resetAllComputationFlags();
				//std::cout << " created new Region " << _regionLabel << std::endl;

				region_deleted = false;
			}

			void Construct(int regionLabel, int nBands)
			{
				//_nBands = nBands;
				//sumOfPixels.reset(new float [nBands]); //minor potential optimization: move this bands sums arrays outside the regions, and 
														// and do the dissim calculations on that outside matrix instead of accessing them indivdullay 
														// form inside the regions collection(may not be benefetial due to ununoforma nature of data access in RHSEG)
				//sumOfPixels = &(HSWO::all_sums_of_pixels[regionLabel -1][0]);
				//sumOfPixels = nullptr;
				//sumOfPixels = new float[_nBands];
				_regionLabel = regionLabel;
				_pixelCount = 0;

				resetAllComputationFlags();
				//std::cout << " created new Region " << _regionLabel << std::endl;

				region_deleted = false;
			}


			//inline void HSWO::Region::updateRegionBandsSumsByMergingWith(float* ptrSums, int region1Label, int region2Label);
			inline void updateRegionBandsSumsByMergingWith(Region* region2, int nBands);

			void resetAllComputationFlags()
			{
				resetAdjComputationFlag();
				resetRegionComputationFlag();			
			}

			void resetAdjComputationFlag();

			void resetRegionComputationFlag();

			int _regionLabel;			
			
			int _bestAdjLabel;
			Region* _bestAdj;
			float _bestAdjDissim;
			
			int _bestRegionLabel;
			Region* _bestRegion;
			float _bestRegionDissim;

			bool _needsAdjRecomputation;
			bool _needsRegionRecomputation;
			//float* sumOfPixels;
			//int _nBands;

			bool region_deleted;
		};

	public:

		HSWO(int nCols, int nRows, int nBands, float spclustWeight, int min_no_of_clusters, std::shared_ptr<float> spDataCube, unsigned int threadsPerBlock = 32);		

		~HSWO(void);

		void TO_GPU_HSWO(unsigned int threadsPerBlock = 128);

		void InitializeHost();
		void InitializeHost(HSWO& pre_segmentation);
		void Device_Initialize();

		void doStep();
		void Device_doStep();
		
		void startRun();
		void continueRun();		
		void Device_run();
		void start_CPU_Run_Then_GPU();

		int getRegionsCount();

		//return an image of int 2D array (stored as 1D array of width*height) representing the
		//clustered image
		std::shared_ptr<int> getClustersLabelsImage();
		
		float _last_merge_dissim_value;
		unsigned int _last_block_count;
		int _max_nRegions;

	private:

		void computeDataStatistics();
		inline float measureDissimilarity(Region* region1, Region* region2);
		//inline float HSWO::measureDissimilarity(Region* region1, int region1Label, Region* region2, int region2Label);
		//float Device_measureDissimilarity(Region& r1, Region& r2);

		inline void measureMeanVector_for_merging_candidates(Region& r1, Region& r2, float*& mean_vector_for_merging_candids);

		Region* merge_regions(Region* region1, int region1Label, Region* region2, int region2Label);
		Region* Device_merge_regions(int region1Label, Region* r1, int region2Label, Region* r2);

		inline bool getPixelBoundryRegions(int x, int y, int pixelLabel, int& topleftRegionLabel, int& topRegionLabel, int& toprightRegionLabel
			, int& leftRegionLabel, int& rightRegionLabel, int& bottomleftRegionLabel, int& bottomRegionLabel, int& bottomrightRegionLabel);

	public:
		int m_min_no_of_clusters;
		const int m_nBands;
		int m_image_width, m_image_height;
		float _spclustWeight;
		bool canWork;
		int myThreadID;
		HSWO* nextToTryGPU;

		//float* m_pRegionsSums; // (width*height  *  nbands) 2D sparse matrix of regions bands pixels sums		

		//std::shared_ptr<float> _spDataCube;
		std::shared_ptr<float> m_pDataCube;

		struct DeviceData
		{
			float* m_dev_pDataCube;

			float* m_dev_pRegionMeans;
			//float* m_pRegionsMeans;

			int* m_dev_pRegionAdjancencies;
			int* m_pRegionAdjancencies;

			float m_dev_nRegions;

			//std::hash_map<int ,Region*> m_regions;

			int* m_dev_RegionKeys;
			int* m_pRegionKeys;

			int* m_dev_RegionPixelsCount;
			int* m_pRegionPixelsCount;

			int* m_labeled_segmented_image;

			//int* ;
			float* m_dev_regionsAdjMinimums;
			int* m_dev_regionsAdjMinimumsLabels;
			int* m_dev_needsAdjRecomputation;

			float* m_dev_regionsMinimums;
			int* m_dev_regionsMinimumsLabels;
			int* m_dev_needsRecomputation;

			int currentStepRegionID_to_merge_1;
			int currentStepRegionID_to_merge_2;

			float* _temp_all_dissims;
			int _temp_all_dissims_size;

			unsigned int m_threadsPerBlock;

			#if USE_AMP
			//TODO change all these plain pointers to smart pointers, to avoid the existing mem leak 
			//as i was lazy to put the delete statments

			Concurrency::array<int, 1>* amp_keysFilter;
			Concurrency::array<int, 1>* amp_adj;
			Concurrency::array<float, 1>* amp_regionsSums;
			Concurrency::array<int, 1>* amp_regionsPixelsCount;
			Concurrency::array<float, 1>* amp_regionsAdjMinimums;
			Concurrency::array<int, 1>* amp_regionsAdjMinimumsLabels;
			Concurrency::array<float, 1>* amp_regionsMinimums;
			Concurrency::array<int, 1>* amp_regionsMinimumsLabels;
			
			Concurrency::array<int, 1>* amp_needsAdjRecomputation;
			Concurrency::array<int, 1>* amp_needsRecomputation;
			#endif

			//std::shared_ptr<Region*> allRegions;
			Region** allRegions;
		};
		DeviceData m_deviceData;

		float* m_tempPixelVector3;
		float* m_minDataValue;
		float* m_meanDataValue;
		float* m_varDataValue;

		std::hash_map<int ,Region*> m_regions; //potenial optimization: try vector instead, may be faster like herb sutter mentioned in : https://channel9.msdn.com/Events/Build/2014/2-661 (min 50)
													//the hash_map is currently only needed in RHSEG stitching, and can be used only for that while using vector for main operations
		//std::vector<Region> m_regions;
		int* m_labeled_segmented_image;

		std::vector<Region>* m_regionsObjects;

	private:
	};
}