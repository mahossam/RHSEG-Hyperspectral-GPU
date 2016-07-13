//ONLY if using visual studio 2012 with amp
#if USE_AMP

#include "CPUAlgorithm.h"
#include "ui_CPUAlgorithm.h"
#include <QTGui>
#include <stdexcept> // stdexcept header file contains runtime_error

#include "Libraries/definitions.h"
#include "Libraries/Logger.h"
#include "Libraries/HSWO.h"
#include "Libraries/RHSEG.h"
#include "Libraries/RHSEGDevice.h"

#include "Libraries/HSWODevice.h"

#include "myglobals.h"

using namespace HyperSpectralToolbox;

void CPUAlgorithm::ApplyGPUHSWO()
{
	try
	{
		if (_imageHeight <=0 || _imageWidth <= 0)
		{
			QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
			return;
		}			
		int min_co_of_clusters = ui->spinDesiredRegions->value();
		
		HSWO hswo(_imageWidth, _imageHeight, _nBands , ui->spinSpclustWeight->value(), min_co_of_clusters, spDataCube, ui->spinThreadsPerBlock->value());

		clock_t startClock = clock();  

		hswo.Device_Initialize();
		//int nsteps = 0;
		while(hswo.getRegionsCount() > min_co_of_clusters)
		{
			/// Calculating HSWO (computation is offloaded to a worker thread)==============================================
			qtawait_local_statments(hswo.Device_doStep(););
			//hswo.Device_doStep();
			//nsteps++;
			//printf("GPU step %d done\n",nsteps);

			//setWindowTitle((ts(hswo.getRegionsCount()) + " regions, last dissim value = " + ts(hswo._last_merge_dissim_value )
			//	+ ", r1 = " + ts(hswo.m_deviceData.currentStepRegionID_to_merge_1)
			//	+ ", r2 = " + ts(hswo.m_deviceData.currentStepRegionID_to_merge_2)
			//	+ ", TBP = "+ ts(hswo.m_deviceData.m_threadsPerBlock)
			//	+ ", last #Blocks = " + ts(hswo._last_block_count)).c_str());				
			
			//std::cout<<
			//	ts(hswo.getRegionsCount()) + " regions, last dissim value = " + ts(hswo._last_merge_dissim_value )
			//	+ ", r1 = " + ts(hswo.m_deviceData.currentStepRegionID_to_merge_1)
			//	+ ", r2 = " + ts(hswo.m_deviceData.currentStepRegionID_to_merge_2)
			//	+ ", TBP = "+ ts(hswo.m_deviceData.m_threadsPerBlock)
			//	+ ", last #Blocks = " + ts(hswo._last_block_count) 
			//	<<endl;

			/*
			int* clustersLabelsImage = new int[_imageWidth*_imageHeight];
			hswo.Device_getClustersLabelsImage(clustersLabelsImage);

			/// Drawing the image on the label box ===============================================================================

			QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);

			QPainter painter(&resultImage);
			painter.fillRect(resultImage.rect(), Qt::white);

			int lastColorIndex = 0;
			map<int, int> clusterIndex_to_ColorIndex;

			for (int irow = 0; irow < _imageHeight; irow++)
			{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
			painter.setPen(QPen(Qt::darkGray));
			/////painter.setPen(QPen(QColor::fromRgb(pixel_color, pixel_color, pixel_color,255)));

			//if(clusterIndex_to_ColorIndex.find(cluster_Index) == clusterIndex_to_ColorIndex.end())
			//{
			//	clusterIndex_to_ColorIndex[cluster_Index] = lastColorIndex;
			//	lastColorIndex++;
			//}
			//painter.setPen(QPen(_pixelColors[ clusterIndex_to_ColorIndex[cluster_Index] ]));

			painter.drawPoint(QPoint(icol, irow));
			}
			}


			_resultImage = resultImage;

			_zoomFactor = 27;
			QImage zoomedImg = _resultImage.scaled(_imageWidth * _zoomFactor, _imageHeight * _zoomFactor,Qt::KeepAspectRatio);


			QPainter overlayPainter(&zoomedImg);
			for (int irow = 0; irow < _imageHeight; irow++)
			{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
			int cluster_Index =  clustersLabelsImage[icol + irow*_imageWidth];

			overlayPainter.setPen(QPen(Qt::lightGray));

			overlayPainter.drawRect(QRectF((_zoomFactor) * icol 
			, (_zoomFactor) * irow
			, _zoomFactor, _zoomFactor));

			overlayPainter.setPen(QPen(Qt::white));

			overlayPainter.drawText(QRectF((_zoomFactor) * icol +3
			, (_zoomFactor) * irow +3
			, _zoomFactor, _zoomFactor),tr("%1").arg(cluster_Index));
			}
			}

			_resultImage = zoomedImg;

			UpdateSolutionImage();

			delete[] clustersLabelsImage;
			clustersLabelsImage = NULL;

			//loop till released by UI		
			//_threadWork = false;
			//while(_threadWork!= true) _pAlgorithmThread->MySleep(1);
			//_pAlgorithmThread->MySleep(10);
			*/
		}		
		
		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));


		shared_ptr<int> spClustersLabelsImage = hswo.getClustersLabelsImage();

		QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);

		QPainter painter(&resultImage);
		painter.fillRect(resultImage.rect(), Qt::white);

		int lastColorIndex = 1;
		map<int, int> clusterIndex_to_ColorIndex;

		for (int irow = 0; irow < _imageHeight; irow++)
		{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
				int cluster_index = spClustersLabelsImage.get()[icol + irow*_imageWidth];
				if(clusterIndex_to_ColorIndex.find(cluster_index) == clusterIndex_to_ColorIndex.end())
				{
					clusterIndex_to_ColorIndex[cluster_index] = lastColorIndex;
					lastColorIndex++;
				}
				painter.setPen(QPen(_pixelColors[ clusterIndex_to_ColorIndex[cluster_index] ]));

				painter.drawPoint(QPoint(icol, irow));
			}
		}

		_resultImage = resultImage;
		UpdateSolutionImage();
	}
	catch(std::runtime_error)
	{
		QMessageBox::critical(this, tr("General runtime exception in AMP GPU HSWO"), "");
		std::cout<<"General runtime exception in AMP GPU HSWO"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in AMP GPU HSWO"), e.what());
		std::cout<<"Exception in AMP GPU HSWO : "<<e.what() <<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in AMP GPU HSWO"), "");
		std::cout<<"Unknown exception in AMP GPU HSWO"<<endl;
	}
}
#endif