#include "RHSEG.h"
#include <cmath>
#include <climits>
#include <string>
#include <QFutureSynchronizer>

#include "./Libraries/definitions.h"
#include "./Libraries/Logger.h"
#include "./Libraries/DebuggingUtilities.h"

#include "myglobals.h"

using namespace HyperSpectralToolbox;
using namespace std;

HyperSpectralToolbox::RHSEG::RHSEG(int nCols, int nRows, int nBands, float spclustWeight, int minNoOfClusters, int nClustersIntermediate, shared_ptr<float> spDataCube, int nRecursiveLevels):
_nBands(nBands)
{
	_n_recur_levels = nRecursiveLevels;
	this->pDataCube = spDataCube;
	_min_no_of_clusters = minNoOfClusters;

	_final_result.reset();

	_imageWidth = nCols;
	_imageHeight = nRows;

	_spclustWeight = spclustWeight;
	_nClustersIntermediate = nClustersIntermediate;
}

RHSEG::~RHSEG(void)
{
}

void HyperSpectralToolbox::RHSEG::InitializeHost()
{

}

shared_ptr<HSWO> HyperSpectralToolbox::RHSEG::merge_splits( vector<shared_ptr<HSWO>>& quarters)
{
	stream_to_log("started merging a split ...");
	//! WARNING, draft code designed to be working with NxN images, logic needs revising when working
	// with rectangular images

	// MERGE DATA CUBES HERE ************************************************************
	int recur_ncols = quarters[0]->m_image_width;
	int recur_nrows = quarters[0]->m_image_height;
	int merged_ncols = recur_ncols*2;
	int merged_nrows = recur_nrows*2;

	shared_ptr<float> sp_merged_DataCube(new float[(merged_nrows)*(merged_ncols)*_nBands]);			
	float* merged_DataCube=sp_merged_DataCube.get();

	int recur_section  = 0;
	for each(auto quarter in quarters)
	{	
		int row_offset;
		int col_offset;
		switch (recur_section)
		{
		case 0:
			col_offset = 0;
			row_offset = 0;
			break;
		case 1:
			col_offset = recur_ncols;
			row_offset = 0;
			break;
		case 2:
			col_offset = 0;
			row_offset = recur_nrows;
			break;
		case 3:
			col_offset = recur_ncols;
			row_offset = recur_nrows;
			break;
		}
		for (int row = 0; row < recur_nrows; ++row)
		{
			for (int col = 0; col < recur_ncols; ++col)
			{
				for (int b=0; b<_nBands;b++)
				{
					int pixel_index = (col + col_offset) + (row + row_offset) * merged_ncols + b*(merged_nrows*merged_ncols);
					unsigned int recur_pixel_index = col + row * recur_ncols + b*(recur_ncols*recur_nrows);					
					merged_DataCube[pixel_index] = quarter->m_pDataCube.get()[recur_pixel_index];
				}
			}
		}
		recur_section++;
	}
	
	/*
	auto dumpImageCube = [&](shared_ptr<float> cube)->void{
		stream_to_log("dumping the data cube in merge step");
		for (int b =0; b<3;b++)
		{		
			for (int y=0; y<merged_nrows; y++)
			{			
				string str;
				for (int x=0; x<merged_ncols; x++)
				{				
					str += ts(cube.get()[x  + y * merged_ncols + b * merged_nrows * merged_ncols]) + "\t\t";
				}
				stream_to_log(str.c_str());
			}
			stream_to_log("");
		}
	};
	dumpImageCube(sp_merged_DataCube);
	*/

	//***********************************************************************************	

	HSWO* hswoTopLeft =  (quarters[0].get());
	HSWO* hswoTopRight = (quarters[1].get());
	HSWO* hswoBotLeft = (quarters[2].get());
	HSWO* hswoBotRight = (quarters[3].get());
	
	int nBands = hswoTopLeft->m_nBands;
	shared_ptr<HSWO> sp_merged_hswo(new HSWO(hswoTopLeft->m_image_width*2, hswoTopLeft->m_image_height*2, nBands , _spclustWeight
		, _nClustersIntermediate, sp_merged_DataCube));

	HSWO* merged_hswo = sp_merged_hswo.get();

	int merged_width = merged_hswo->m_image_width;
	int section_width = hswoTopLeft->m_image_width;
	int section_height = hswoTopLeft->m_image_height;

	int merged_region_counter = 0;
	set<int> allpixels;

	map<int, int> section_to_merged_regionID;
	map<int, int> merged_pixelID_to_merged_regionID;

	auto copyRegionsToMerged = [&](HSWO* ptrHswoToCopy, int quarterNo)->void
	{	
		section_to_merged_regionID.clear();		

		for each(auto currentRegionID_cell in ptrHswoToCopy->m_regions )
		{
			int currentRegionID = currentRegionID_cell.first;

			merged_region_counter = merged_region_counter+1;
			section_to_merged_regionID[currentRegionID] = merged_region_counter;
			merged_hswo->m_regions[merged_region_counter] = new HSWO::Region(merged_region_counter, _nBands);
		}

		//adjacencies
		for each( auto currentRegionID_cell in ptrHswoToCopy->m_regions )
		{
			int currentRegionID = currentRegionID_cell.first;
			HSWO::Region& currentRegion = *(currentRegionID_cell.second);

			int currentMergedHSWORegionID = section_to_merged_regionID[currentRegionID];
			HSWO::Region& currentMergedHSWORegion = *(merged_hswo->m_regions[currentMergedHSWORegionID]);

			for each(auto currentAdjacentNeighbor in currentRegion.adjacentRegions )
			{
				// I have a section adjacent region ID

				// get its corresponding merged_ID
				int currentMergedHSWOAdjacentNeighborID = section_to_merged_regionID[currentAdjacentNeighbor->_regionLabel];

				// add the new adjacent region merged ID to its corresponding merged hswo
				// region
				currentMergedHSWORegion.adjacentRegions.insert(merged_hswo->m_regions[currentMergedHSWOAdjacentNeighborID]);
			}

			// also copy the region means from hswoTopLeft
			for(int b=0; b<nBands; b++)
			{
				//currentMergedHSWORegion.sumOfPixels[b] = currentRegion.sumOfPixels[b];
				(currentMergedHSWORegion.sumOfPixels)[b] = (currentRegion.sumOfPixels)[b];
			}
		}

		//============================ PIXELS ====================================
		for each(auto currentRegionID_cell in ptrHswoToCopy->m_regions )
		{
			int currentRegionID = currentRegionID_cell.first;
			HSWO::Region* ptrCurrentRegion = currentRegionID_cell.second;
			
			int currentMergedHSWORegionID = section_to_merged_regionID[currentRegionID];
			HSWO::Region& currentMergedHSWORegion = *(merged_hswo->m_regions[currentMergedHSWORegionID]);

			for each( int index in ptrCurrentRegion->pixels_IDs )
			{
				//matlab only======
				int zero_based_x = (index-1) % section_width;
				int zero_based_y = (floor((float)(index-1)/section_width));
				//=================
				int x = zero_based_x+1;
				int y = zero_based_y+1;
				int new_merged_PixelID = -1;

				switch (quarterNo)
				{
				case 0:				
					// ********************* FOR TOP LEFT *****************************
					new_merged_PixelID = (y-1)*merged_width+ x;
					// ****************************************************************
					break;

				case 1:
					// ********************* FOR TOP Right *****************************
					new_merged_PixelID = (y-1)*merged_width+ x+section_width;
					// ****************************************************************
					break;

				case 2:
					// ********************* FOR BOT LEFT *****************************
					new_merged_PixelID = (y-1)*merged_width+ x + (section_width*merged_width);
					// ****************************************************************
					break;

				case 3:
					// ********************* FOR BOT RIGHT *****************************
					new_merged_PixelID = (y-1)*merged_width+ x + (section_width*merged_width) + section_width;
					// ****************************************************************
					break;
				}

				currentMergedHSWORegion.pixels_IDs.push_back(new_merged_PixelID);
				currentMergedHSWORegion._pixelCount += 1;
				merged_pixelID_to_merged_regionID[new_merged_PixelID] = currentMergedHSWORegionID;

				/*
				if (allpixels.find(new_merged_PixelID) != allpixels.end())
				{
					throw std::exception("HSWO: Must not happen,one pixel in more than 1 region");
				}

				allpixels.insert(new_merged_PixelID);
				*/
			}
		}
	};

	copyRegionsToMerged(hswoTopLeft, 0);
	copyRegionsToMerged(hswoTopRight, 1);
	copyRegionsToMerged(hswoBotLeft, 2);
	copyRegionsToMerged(hswoBotRight, 3);

	///  ****************************************************************


	///  ****************************************************************
	/// Finally, merge the boundary regions that became adjacent to each other as a result of image merging
	///  ****************************************************************

	// first, regions on the vertical boundary (at x = section_width)
	for (int y = 1; y<=merged_hswo->m_image_height; y++)
	{
		int x = section_width;

		int zerobasedX=x-1, zerobasedY= y-1;
		int currentPixelIndex = (zerobasedX) + (zerobasedY) * merged_hswo->m_image_width;
		int currentPixelLabel = currentPixelIndex +1; // only in matlab code because of 1 based structures

		int current_region_ID_of_current_pixel = merged_pixelID_to_merged_regionID[currentPixelLabel];
		HSWO::Region* current_region_of_current_pixel  = merged_hswo->m_regions[current_region_ID_of_current_pixel];

		auto region_adj_logic = [&](int pixelLabel)->void{
			int adj_region_ID_of_adj_pixel = merged_pixelID_to_merged_regionID[pixelLabel]; 
			if (adj_region_ID_of_adj_pixel != current_region_ID_of_current_pixel )
			{
				HSWO::Region* adj_region_of_adj_pixel  = merged_hswo->m_regions[adj_region_ID_of_adj_pixel];
				
				current_region_of_current_pixel->adjacentRegions.insert(adj_region_of_adj_pixel);
				adj_region_of_adj_pixel->adjacentRegions.insert(current_region_of_current_pixel);
			}
		};

		// add boundary regions (which are the MAX_REGION_ADJ surrounding pixels) for this current pixel(x,y)
		//topleft
		if (zerobasedX > 0 && zerobasedY > 0 )
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}
		//top
		if (zerobasedY > 0)
		{
			int pixelIndex = (zerobasedX) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//topright
		if ( zerobasedX < (merged_hswo->m_image_width - 1) && zerobasedY > 0 )
		{
			int pixelIndex = (zerobasedX + 1) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//left
		if (zerobasedX > 0)
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//right
		if (zerobasedX < (merged_hswo->m_image_width - 1))
		{

			int pixelIndex = (zerobasedX + 1) + (zerobasedY) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottomleft
		if (zerobasedX > 0 && zerobasedY < (merged_hswo->m_image_height - 1) )
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottom
		if (zerobasedY < (merged_hswo->m_image_height - 1))
		{
			int pixelIndex = (zerobasedX) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottomright
		if (zerobasedX < (merged_hswo->m_image_width - 1)  && zerobasedY < (merged_hswo->m_image_height - 1) )
		{
			int pixelIndex = (zerobasedX + 1) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}	

	}

	///**************************************************************************
	/// second, regions on the horizontal boundary (at y=section_height)
	///**************************************************************************

	for (int x = 1; x<=merged_hswo->m_image_width;x++)
	{
		int y = section_height;

		int zerobasedX=x-1, zerobasedY= y-1;
		int currentPixelIndex = (zerobasedX) + (zerobasedY) * merged_hswo->m_image_width;
		int currentPixelLabel = currentPixelIndex +1; // only in matlab code because of 1 based structures

		int current_region_ID_of_current_pixel = merged_pixelID_to_merged_regionID[currentPixelLabel];
		HSWO::Region* current_region_of_current_pixel  = merged_hswo->m_regions[current_region_ID_of_current_pixel];

		auto region_adj_logic = [&](int pixelLabel)->void{
			int adj_region_ID_of_adj_pixel = merged_pixelID_to_merged_regionID[pixelLabel]; 
			if (adj_region_ID_of_adj_pixel != current_region_ID_of_current_pixel )
			{
				HSWO::Region* adj_region_of_adj_pixel  = merged_hswo->m_regions[adj_region_ID_of_adj_pixel];

				current_region_of_current_pixel->adjacentRegions.insert(adj_region_of_adj_pixel);
				adj_region_of_adj_pixel->adjacentRegions.insert(current_region_of_current_pixel);
			}
		};

		// add boundary regions (which are the MAX_REGION_ADJ surrounding pixels) for this current pixel(x,y)
		//topleft
		if (zerobasedX > 0 && zerobasedY > 0 )
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//top
		if (zerobasedY > 0)
		{
			int pixelIndex = (zerobasedX) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//topright
		if (zerobasedX < (merged_hswo->m_image_width - 1) && zerobasedY > 0)
		{
			int pixelIndex = (zerobasedX + 1) + (zerobasedY - 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//left
		if ( zerobasedX > 0)
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//right
		if (zerobasedX < (merged_hswo->m_image_width - 1))
		{
			int pixelIndex = (zerobasedX + 1) + (zerobasedY) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottomleft
		if( zerobasedX > 0 && zerobasedY < (merged_hswo->m_image_height - 1))
		{
			int pixelIndex = (zerobasedX - 1) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottom
		if (zerobasedY < (merged_hswo->m_image_height - 1))
		{
			int pixelIndex = (zerobasedX) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}

		//bottomright
		if (zerobasedX < (merged_hswo->m_image_width - 1)  && zerobasedY < (merged_hswo->m_image_height - 1))
		{ 
			int pixelIndex = (zerobasedX + 1) + (zerobasedY + 1) * merged_hswo->m_image_width;
			pixelIndex = pixelIndex+1;// only in matlab code because of 1 based structures
			region_adj_logic(pixelIndex);
		}
	}
	
	stream_to_log("finished merging the split, emitting signal to GUI...");
	//emit hswoReadyForDrawing(sp_merged_hswo);
	stream_to_log("signal emitted ...");

	return sp_merged_hswo;
}

//shared_ptr<HSWO> HyperSpectralToolbox::RHSEG::intermediateStep( int region_label_offset, int recur_level
MyFutureResult HyperSpectralToolbox::RHSEG::intermediateStep( int region_label_offset, int recur_level
	, int section_no, const int ncols, const int nrows, shared_ptr<float> sp_dataCube)
{	
	//double threshold = 0.0;

	// Declare arrays which will hold processing window seam information

	//unsigned int index_data_size = 1;
	//index_data_size = ncols * params.seam_size;
	//vector<Index> row_seam_index_data(index_data_size);
	//index_data_size = nrows * params.seam_size;
	//vector<Index> col_seam_index_data(index_data_size);
	
#pragma region Dumping data to debug output
	/*
	auto dumpImageCube = [&](shared_ptr<float> cube)->void{
		stream_to_log("dumping the data cube in intermediate step");
		for (int b =0; b<1;b++)
		{		
			for (int y=0; y<nrows; y++)
			{			
				string str;
				for (int x=0; x<ncols; x++)
				{				
					str += ts(cube.get()[x  + y * ncols + b * nrows * ncols]) + "\t";
				}
				stream_to_log(str.c_str());
			}
			stream_to_log("");
		}
	};
	//dumpImageCube(sp_dataCube);
	
	auto dumpHswoImage = [](shared_ptr<HSWO> theHswo)->void{
		stream_to_log("dumping the hswo image" );
		for (int b =0; b<1;b++)
		{		
			for (int y=0; y<theHswo->m_image_height; y++)			
			{			
				string str;
				for (int x=0; x<theHswo->m_image_width; x++)
				{				
					str += ts(theHswo->m_pDataCube.get()[x  + y * theHswo->m_image_width + b * theHswo->m_image_height * theHswo->m_image_width]) + "\t";
				}
				stream_to_log(str.c_str());
			}
			stream_to_log("");
		}		
	};
	*/
#pragma endregion Dumping data to debug output

	stream_to_log("started recurisve internediate step ...");

	if(recur_level == _n_recur_levels)
	{
		stream_to_log("recurisve internediate step in deepest level, running hswo ...");
		//the deepest level of recursion
		shared_ptr<HSWO> hswo(new HSWO(ncols, nrows, _nBands, _spclustWeight, _nClustersIntermediate, sp_dataCube));
		
		clock_t startClock2 = clock();

		QFuture<void> future = QtConcurrent::run(hswo.get(), &HSWO::startRun); // hswo->run()
		stream_to_log("finished hswo.");
		//dumpHswoImage(hswo);
		future.waitForFinished();

		printf("Finished CPU Quarter in RHSEG\n");
		clock_t endClock2 = clock();
		float timeInSeconds = (float)(endClock2 - startClock2)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		printf("CPU Quarter Time = %d seconds , %d milliseconds", (int)timeInSeconds, milliSeconds);
		
		stream_to_log("emitting the signal");
		//emit hswoReadyForDrawing(hswo);
		stream_to_log("emitted");

		MyFutureResult res;
		res.future = future;
		res.spHswo = hswo;

		return res;
	}
	else if (recur_level < _n_recur_levels)
	{
		// Need to set up and make a call to a deeper level of recursion
		stream_to_log("indermediate step in middle level ...");

		unsigned short int recur_section = 0;
		unsigned int recur_region_index, recur_nregions;
		unsigned int recur_region_classes_size, class_label_offset = 0;
		unsigned int recur_npixels;
		int col_offset = 0, row_offset = 0;

		// TODO revise these variables
		//int recur_pixel_ncols = pixel_ncols / 2;
		//int recur_pixel_nrows = 1;
		//recur_pixel_nrows = pixel_nrows / 2;

		// TODO revise this if statement
		//if (pixel_ncols != ncols) 
		//recur_pixel_ncols = pixel_ncols;

		// TODO revise this if statement
		//if (pixel_nrows != nrows)
		//recur_pixel_nrows = pixel_nrows;

		//recur_npixels = recur_pixel_ncols * recur_pixel_nrows;

		//	vector<Pixel> recur_pixel_data(recur_npixels, Pixel());
		//auto recur_region_classes = shnew(vector<shptr(HSWO::Region)>,());

		//nregions = 0;
		int recur_ncols = ncols / 2;
		int recur_nrows = nrows / 2;

		unsigned short int proc_section, next_recur_level = recur_level + 1;
		//unsigned short int next_stride = stride / params.nb_strides;
		//unsigned short int next_nb_sections = stride;
		vector<shared_ptr<HSWO>> quartersResults;
		QFutureSynchronizer<void> synchronizer;

		for (proc_section = 0; proc_section < 4; ++proc_section)
		{
			//TODO revise the stride variable logic
			switch (proc_section)
			{
			case 0:
				col_offset = 0;
				row_offset = 0;
				//recur_section = section;
				break;
			case 1:
				col_offset = recur_ncols;
				row_offset = 0;
				//recur_section = section + stride;
				break;
			case 2:
				col_offset = 0;
				row_offset = recur_nrows;
				//recur_section = section + 2 * stride;
				break;
			case 3:
				col_offset = recur_ncols;
				row_offset = recur_nrows;
				//recur_section = section + 3 * stride;
				break;
			}
			//class_label_offset = nregions;

			shared_ptr<float> recur_pixel_data(new float[recur_ncols*recur_nrows*_nBands]);
			float* p_recur_pixel_data = recur_pixel_data.get();
			float* p_DataCube = sp_dataCube.get();
			
			//TODO this if statement was added here in original Tilton RHSEG
			//because of the intermediate step output parameter AFAIK
			//but I use it now for creation of each segment pixels map (image)
			//if (recur_level >= params.ionb_levels)
			//{
			//qDebug() << "dumping the quarter image ";
			for (int row = 0; row < recur_nrows; ++row)
			{
				//stringstream str;				
				for (int col = 0; col < recur_ncols; ++col)
				{					
					{
						int pixel_index = (col + col_offset) + (row + row_offset) * ncols + 0*(ncols*nrows);
						unsigned int recur_pixel_index = col + row * recur_ncols + 0*(recur_ncols*recur_nrows);
						//str << p_DataCube[pixel_index] << "\t";
					}

					unsigned int recur_pixel_index;
					for (int b=0; b<_nBands;b++)
					{
						int pixel_index = (col + col_offset) + (row + row_offset) * ncols + b*(ncols*nrows);
						recur_pixel_index = col + row * recur_ncols + b*(recur_ncols*recur_nrows);
						p_recur_pixel_data[recur_pixel_index] = p_DataCube[pixel_index];						
					}
					//str << p_recur_pixel_data[recur_pixel_index] << "\t";
				}
				//qDebug()<<str.str().c_str();
			}
			//}

			stream_to_log("calling a new rhseg for a quarter ...");

			//call intermediateStep for every quarter	
			auto quarterResult = intermediateStep(class_label_offset, next_recur_level, recur_section
				, recur_ncols, recur_nrows, recur_pixel_data);

			stream_to_log("returned from rhseg call of rhseg in intermediate step.");

			quartersResults.push_back(quarterResult.spHswo);	
			synchronizer.addFuture(quarterResult.future);

			//if (recur_level >= params.ionb_levels)
			//{
			//for (row = 0; row < recur_pixel_nrows; ++row)
			//{
			//for (col = 0; col < recur_pixel_ncols; ++col)
			//{
			//pixel_index = (col + col_offset) + (row + row_offset) * pixel_ncols;
			//recur_pixel_index = col + row * recur_pixel_ncols;
			//pixel_data[pixel_index] = recur_pixel_data[recur_pixel_index];
			//}
			//}
			//}
			//else
			//{
			//if ((next_recur_level == params.ionb_levels) && (params.nb_sections > 1))
			//{
			//spatial_data.save_region_class_label_map(recur_section);
			//save_pixel_data(recur_section, recur_pixel_data, temp_data);
			//}
			//}

			//TODO add threshold logic
			//if (max_threshold < threshold)
			//max_threshold = threshold;

			//if (((unsigned int) (region_classes.size())) < (nregions + recur_nregions))
			//region_classes.resize(nregions + recur_nregions);

			//nregions += recur_nregions;	

		} // for (proc_section = 0; proc_section < params.nb_strides; ++proc_section)	

		//combine the last 4 results in 1 new region set
		//vector<HSWO*> veccy;
		//for each(auto quarter in quartersResults)
		{
			//veccy.push_back(quarter.get());
		}
		
		//synchronize all quarters futures
		printf("waiting all futures to finish\n");
		synchronizer.waitForFinished();
		printf("finished all futures\n");

		shared_ptr<HSWO> sp_merged_hswo = merge_splits(quartersResults);	

		stream_to_log("splits merged ...");

		emit stepCompleted(recur_level);	
		emit hswoReadyForDrawing(sp_merged_hswo);

		stream_to_log("running hswo on the merged HSWO result ...");
		//call hswo for the result
		if(recur_level == 1)
		{
			sp_merged_hswo->m_min_no_of_clusters = _min_no_of_clusters;
			sp_merged_hswo->continueRun();
		}
		else
		{
			sp_merged_hswo->continueRun();
		}

		stream_to_log("finished merged HSWO run.");
		
		//return the final result
		printf("Level finished\n");

		MyFutureResult res;
		res.spHswo = sp_merged_hswo;
		return res;
	}
}

void HyperSpectralToolbox::RHSEG::run()
{
	//the only 1 call
	auto result = intermediateStep(0, 1, -1, _imageWidth, _imageHeight, pDataCube);
	_final_result = result.spHswo;
	stream_to_log("finished all RHSEG.");
}

int HyperSpectralToolbox::RHSEG::getRegionsCount()
{
	if((_final_result).get() != nullptr)
	{
		return _final_result->m_regions.size();
	}
	return 0;
}

shared_ptr<int> HyperSpectralToolbox::RHSEG::getClustersLabelsImage()
{
	shared_ptr<int> image(new int[_imageWidth*_imageHeight]);
	if((_final_result).get() != nullptr)
	{
		for(int i=0; i<_imageHeight*_imageWidth;i++) 
		{
			image.get()[i]= -1;
		}
		set<int> allpixels;

		for each(auto itCurrentRegion in _final_result->m_regions)
		{
			for each(int pixelLabel in itCurrentRegion.second->pixels_IDs)
			{				
				if(allpixels.find(pixelLabel) != allpixels.end()) throw exception("HSWO: Must not happen");
				allpixels.insert(pixelLabel);
				image.get()[(pixelLabel-1)]=itCurrentRegion.first;				
			}		
		}
	}
	return image;
}


