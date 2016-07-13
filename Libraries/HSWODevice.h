#pragma once

#if USE_CUDA
#include <thrust/host_vector.h>

int do_kernel();

namespace HyperSpectralToolbox
{
	void DeviceInitializeData(HSWO* object);
	void DeviceExit(HSWO* object);
	
	void DeviceInitStep(HSWO* object);
	void DeviceExitStep(HSWO* object);

	void DeviceUpdateRegionMean(HSWO* object, int lastChangedRegionID, HSWO::Region* ptrLastChangedRegionID);
	void DeviceUploadInitialMeans(HSWO* object);
	void DeviceUploadInitialAdjacents(HSWO* object);


	void DeviceCalcAllDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim);
	void DeviceCalcRegionDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim);
	void DeviceResetOthersBestRegionComputationFlagsFrom(HSWO* object, int region1Label, int region2Label ) ;
	void DeviceRecomputeOthersBestRegionsDissimilarityTo(HSWO* object, int regionLabel);

	void DeviceCalcAllRegionsMeanVectors(HSWO* object);
	
	//template <typename T> 
	void DeviceThrustFill(int* pdeviceArray, int size, int newValue);
	void DeviceThrustFillFloat(float* pdeviceArray, int size, float newValue);
	void DeviceThrustRemove(int* pdeviceArray, int size, int valueToRemove);

	void DeviceThrustResetAdjacencyOfMergedRegion(int* m_deviceDatam_dev_pRegionAdjancencies, int lastMergedRegion_OffsetInMatrix);
	void DeviceThrustReplace(int* m_deviceData_DOT_m_dev_pRegionAdjancencies_PLUS_currentAdjacentRegion_OffsetInMatrix, int nAdjacentsOfAdjacent, int adjRegionLabel, int regionLabel );

	float computeMax(HSWO* object);
	float computeMin(HSWO* object);

	void runThrustExpr(thrust::host_vector<int>& h_vec, int N);
	void run_expr();
}

#elif USE_AMP

namespace HyperSpectralToolbox
{
	void DeviceInitializeData(HSWO* object);
	void DeviceExit(HSWO* object);
	
	void DeviceInitStep(HSWO* object);
	void DeviceExitStep(HSWO* object);

	void DeviceUpdateRegionMean(HSWO* object, int lastChangedRegionID, HSWO::Region* ptrLastChangedRegionID);
	void DeviceUploadInitialMeans(HSWO* object);
	void DeviceUploadInitialAdjacents(HSWO* object);


	void DeviceCalcAllDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim);
	void DeviceCalcRegionDissims(HSWO* object, unsigned int threadsPerBlock, int& label1, int& label2, float& minDissim);
	void DeviceResetOthersBestRegionComputationFlagsFrom(HSWO* object, int region1Label, int region2Label );
	void DeviceRecomputeOthersBestRegionsDissimilarityTo(HSWO* object, int regionLabel);

	void DeviceCalcAllRegionsMeanVectors(HSWO* object);
	
	//template <typename T> 
	void DeviceThrustFill(int* pdeviceArray, int size, int newValue);
	void DeviceThrustFillFloat(float* pdeviceArray, int size, float newValue);
	void DeviceThrustRemove(int* pdeviceArray, int size, int valueToRemove);

	void DeviceThrustResetAdjacencyOfMergedRegion(int* m_deviceDatam_dev_pRegionAdjancencies, int lastMergedRegion_OffsetInMatrix);
	void DeviceThrustReplace(int* m_deviceData_DOT_m_dev_pRegionAdjancencies_PLUS_currentAdjacentRegion_OffsetInMatrix, int nAdjacentsOfAdjacent, int adjRegionLabel, int regionLabel );

	//float reduction_simple_1_min(const std::vector<float>& source, int& minIndex);
	float reduction_simple_1_min_INPLACE(Concurrency::array<float,1>& source, int& minIndex);
}

#endif

