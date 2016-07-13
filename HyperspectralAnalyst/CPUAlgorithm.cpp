#include "CPUAlgorithm.h"
#include "ui_CPUAlgorithm.h"
#include <QTGui>
//#include <QTGui/QMessageBox>
//#include <QTGui/QPainter>
//#include <QTGui/QFileDialog>
//#include <QTGui/QImageReader>
#include <stdexcept> // stdexcept header file contains runtime_error

#include "Libraries/definitions.h"
#include "Libraries/Logger.h"
#include "Libraries/HSWO.h"
#include "Libraries/RHSEG.h"
#include "Libraries/RHSEGDevice.h"

//#include <thrust/host_vector.h>
//#include <thrust/sort.h>

// Required to include CUDA vector types
//#include <vector_types.h>
//#include <cutil_inline.h>

#include "Libraries/HSWODevice.h"

#include "myglobals.h"

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
//extern "C" void runTest(const int argc, const char** argv, 
//                        char* data, int2* data_int2, unsigned int len);
//extern "C" void runTest2(float * h_A, float * h_B, float * h_C);
//extern "C" void runThrustExpr(thrust::host_vector<int>& h_vec, int N);

using namespace HyperSpectralToolbox;

CPUAlgorithm::CPUAlgorithm(QWidget *parent) : QDialog(parent), ui(new Ui::CPUAlgorithm)
{
	spDataCube.reset();

	useGPU = false;

	ui->setupUi(this);

	InitializeApp();
	InitializeForm();	
}

CPUAlgorithm::~CPUAlgorithm()
{
	delete ui;	
}

void CPUAlgorithm::changeEvent(QEvent *e)
{
	QDialog::changeEvent(e);
	switch (e->type()) {
	case QEvent::LanguageChange:
		ui->retranslateUi(this);
		break;
	default:
		break;
	}
}

void CPUAlgorithm::InitializeApp()
{
	Logger::InitializeLogger();

	_zoomFactor = 1;
	_nEndmembers = 16;
	_bDrawEndmembers = true;
	_endmembersColor = Qt::black;

	_imageHeight = 0;
	_imageWidth = 0;	

/*
	unsigned char colormap_values[] = {0,127,200,255};
	short unsigned int index, red_value, green_value, blue_value;

	int red_index, green_index, blue_index;
	index = 0;
	for (red_index = 0; red_index < 4; red_index++)
	{
		for (green_index = 3; green_index >= 0; green_index--)
		{			
			for (blue_index = 0; blue_index < 4; blue_index++)
			{
				red_value = colormap_values[red_index];
				green_value = colormap_values[green_index];
				blue_value = colormap_values[blue_index];
				
				//if (red_value > 0)
					//red_value = ((red_value+1)*256) - 1;
				//else
					//red_value = 0;
				//if (blue_value > 0)
					//blue_value = ((blue_value+1)*256) - 1;
				///else
				//	blue_value = 0;
				//if (green_value > 0)
					//green_value = ((green_value+1)*256) - 1;
				///else
					//green_value = 0;
				
				_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);

				index++;
			}
		}
	}
	index = 21;
	red_value = 64; green_value = 127; blue_value = 191;
	//red_value = ((red_value+1)*256) - 1;
	//blue_value = ((blue_value+1)*256) - 1;
	//green_value = ((green_value+1)*256) - 1;
	_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);

	index = 63;
	red_value = 127; green_value = 63; blue_value = 191;
	//red_value = ((red_value+1)*256) - 1;
	//blue_value = ((blue_value+1)*256) - 1;
	//green_value = ((green_value+1)*256) - 1;
	_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);

	index = 64;
	red_value = 63; green_value = 191; blue_value = 127;
	//red_value = ((red_value+1)*256) - 1;
	//blue_value = ((blue_value+1)*256) - 1;
	//green_value = ((green_value+1)*256) - 1;
	_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);

	index = 65;
	red_value = 127; green_value = 191; blue_value = 63;
	//red_value = ((red_value+1)*256) - 1;
	//blue_value = ((blue_value+1)*256) - 1;
	//green_value = ((green_value+1)*256) - 1;
	_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);

	index = 66;
	red_value = 191; green_value = 63; blue_value = 127;
	//red_value = ((red_value+1)*256) - 1;
	//blue_value = ((blue_value+1)*256) - 1;
	//green_value = ((green_value+1)*256) - 1;
	_pixelColors[index] = QColor::fromRgb(red_value,green_value,blue_value,255);
*/

	/*
	_pixelColors[0] = Qt::black;// QColor::fromRgb(255, 254, 137,255);	
	_pixelColors[1] = Qt::cyan;// QColor::fromRgb(255, 254, 137,255);
	_pixelColors[2] = Qt::gray;
	_pixelColors[3] = Qt::red;//QColor::fromRgb(255, 89, 1,255);
	_pixelColors[4] = Qt::blue;//QColor::fromRgb(160, 78, 158,255);
	_pixelColors[5] = Qt::magenta;//QColor::fromRgb(255, 254, 137,255);
	_pixelColors[6] = Qt::yellow;//QColor::fromRgb(12, 255, 7,255);
	_pixelColors[7] = Qt::green;//QColor::fromRgb(89, 1, 255,255);
	_pixelColors[8] = Qt::darkCyan;//QColor::fromRgb(3, 28, 241,255);
	_pixelColors[9] = Qt::darkGray;
	_pixelColors[10] = Qt::darkRed;//QColor::fromRgb(5, 255, 133,255);
	_pixelColors[11] = Qt::darkBlue;//QColor::fromRgb(101, 173, 255,255);
	_pixelColors[12] = Qt::darkMagenta;//QColor::fromRgb(255, 2, 251,255);
	_pixelColors[13] = Qt::darkYellow;//QColor::fromRgb(172, 175, 84,255);
	_pixelColors[14] = Qt::darkGreen;//QColor::fromRgb(3, 171, 255,255);

	auto blendRGB = [](QColor col1, QColor col2){return QColor(col1.red()+col2.red()/2, col1.green()+col2.green()/2, col1.blue()+col2.blue()/2 );};

	_pixelColors[15] = blendRGB(Qt::red, Qt::cyan);
	_pixelColors[16] = QColor::fromRgb(110, 200, 255,255);//blendRGB(Qt::red, Qt::yellow);
	_pixelColors[17] = blendRGB(Qt::green, Qt::magenta);
	_pixelColors[18] = QColor::fromRgb(160, 78, 158,255);//blendRGB(Qt::green, Qt::yellow);
	_pixelColors[19] = QColor::fromRgb(5, 255, 133,255);//blendRGB(Qt::blue, Qt::magenta);
	_pixelColors[20] = blendRGB(Qt::blue, Qt::yellow);
	*/
	
	QColor mix = QColor::fromRgb(255,255,255);

	auto generateRandomColor = [](QColor mix)->QColor{
		int random = rand();
		double red = random % 256;
		random = rand();
		double green = random % 256;
		random = rand();
		double blue = random % 256;

		// mix the color
		red = (red + mix.red()) / 2;
		green = (green + mix.green()) / 2;
		blue = (blue + mix.blue()) / 2;		

		auto fRand = [](double fMin, double fMax)->double {
			double f = (double)rand() / RAND_MAX;
			return fMin + f * (fMax - fMin);
		};

		/*
		double constrast = fRand(.5, 1.3);
		red = red * constrast; if(red > 255.0) red = 255;
		green =green * constrast; if(green > 255.0) green = 255;
		blue = blue * constrast; if(blue > 255.0) blue = 255;
		*/

		QColor color = QColor::fromRgb(red, green, blue);
		return color;
	};

	for (int i=0; i< 200; i++)
	{
		_pixelColors[i] = generateRandomColor(mix);
	}
	
	/*
	//TODO, added this only for auto prifiling , remove for normal operation ==================
	qtawait_local_statments(

		ifstream f("F:\\mahmoud\\COLLEGE\\Research\\Thesis\\Data\\aviris\\f920612t01p02_r02c\\for_results\\middle_32_float.txt");

	f>>_imageWidth>>_imageHeight>>_nBands;

	spDataCube.reset(new float [_nBands*_imageHeight*_imageWidth]);

	float lastVal=0;

	///The new format ...
	for (int i = 0; i < _nBands; i++)
	{
		for (int irow = 0; irow < _imageHeight; irow++)
		{
			for (int icol = 0; icol < _imageWidth; icol++)			
			{		
				f>>lastVal ;
				(spDataCube.get())[icol  + irow * _imageWidth + i * _imageHeight * _imageWidth] = lastVal;
			}
		}
	}
	f.close();
	);	

	setWindowTitle(" All channels loaded.");	
	on_btnRunCUDA_clicked();

	QCoreApplication::exit();
	exit(0);
	*/
}

void CPUAlgorithm::InitializeForm()
{
	//threadDrawEndmembers = new QThread(gcnew ThreadStart(this, &nativeHyperSpectral::MainForm::DrawEndmembers));
	//threadProgress = new Thread(gcnew ThreadStart(this, &nativeHyperSpectral::MainForm::ShowProgress));
}

void CPUAlgorithm::LoadSampleImageFromFile()
{	
	QString fileName;
	QFileDialog::Options options;
	QString selectedFilter;
	fileName = QFileDialog::getOpenFileName(this,
		tr("Select Data File ..."),
		"",
		tr("All Files (*);;Text Files (*.txt)"),
		&selectedFilter,
		options);

	if (fileName.isEmpty())
	{
		QMessageBox::warning(this,tr("No File is Selected"),tr("please select a data file to open !"));
		return;
	}


	qtawait_local_statments(
		/// ///////////////////////////////////////////////////////////////////////
		/// Loading the new image file in a different thread using qtawait_local_statments
		/// ///////////////////////////////////////////////////////////////////////

		ifstream f(fileName.toUtf8().data());

		f>>_imageWidth>>_imageHeight>>_nBands;

		spDataCube.reset(new float [_nBands*_imageHeight*_imageWidth]);

		float lastVal=0;

		///The new format ...
		for (int i = 0; i < _nBands; i++)
		{
			for (int irow = 0; irow < _imageHeight; irow++)
			{
				for (int icol = 0; icol < _imageWidth; icol++)			
				{		
					f>>lastVal ;
					(spDataCube.get())[icol  + irow * _imageWidth + i * _imageHeight * _imageWidth] = lastVal;
				}
			}
		}
		f.close();
	);	

	setWindowTitle(" All channels loaded.");	
}

void CPUAlgorithm::StartComuptation()
{
	//btnStop->Enabled = true;
	//progressBar1->Visible = true;
	//progressBar1->Value = 0;
	//progressBar1->Maximum = _imageWidth * _imageHeight;
	//lblProgress->Text = "";
	//btnCluster->Enabled = false;
	//btnClassify->Enabled = false;

	//threadProgress->Start();	

//	_pAlgorithmThread->start();
	
	
	if(useGPU)
	{
		ApplyGPUHSWO();
	}
	else
	{
		ApplyHSWO();
	}	
}

void CPUAlgorithm::StopComputation()
{
	//threadProgress->Stop();	
	//_pAlgorithmThread->stop();
}
void CPUAlgorithm::DrawSolution()
{
}

#if USE_CUDA
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
		cout<<"inside GPU HSEG"<<endl;

		clock_t startClock = clock();  

		hswo.Device_Initialize();
		//int nsteps = 0;
		qtawait_local_statments(
		while(hswo.getRegionsCount() > min_co_of_clusters)
		{
			/// Calculating HSWO (computation is offloaded to a worker thread)==============================================
			hswo.Device_doStep();

			//clock_t midClock = clock();  
			//if(HSWO::silentMode)
			//{
			//	if((midClock - startClock) > 1000*30)
			//	{
			//		break;
			//	}
			//}


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
		);

		DeviceExit(&hswo);		

		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));
		cout<< "Time = " << timeInSeconds << ", " << milliSeconds <<endl;
		if(HSWO::silentMode) 
		{
			std::exit(0);
			return;
		}

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

		cout<< "comp time = " << HSWO::time_count <<endl;
		HSWO::time_count = 0;

		//cout<< "speedup over CPU HSEG = " << float(490000)/float(endClock - startClock)  <<endl;
		//HSWO::time_count = 0;
	}
	catch(std::runtime_error)
	{
		QMessageBox::critical(this, tr("General runtime exception in GPU HSWO"), "");
		std::cout<<"General runtime exception in GPU HSWO"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in GPU HSWO"), e.what());
		std::cout<<"Exception in GPU HSWO : "<<e.what()<<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in GPU HSWO"), "");
		std::cout<<"Unknown exception in GPU HSWO"<<endl;
	}
}
#endif

void CPUAlgorithm::ApplyHSWO()
{
	try
	{
		if (_imageHeight <=0 || _imageWidth <= 0)
		{
			QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
			return;
		}

		/// Calculating HSWO ==============================================
		int min_co_of_clusters = ui->spinDesiredRegions->value();

		HSWO hswo(_imageWidth, _imageHeight, _nBands, ui->spinSpclustWeight->value(), min_co_of_clusters, spDataCube);

		//	hswo.runHSWO();

		clock_t startClock = clock();  

		hswo.InitializeHost();  

		while(hswo.getRegionsCount() > min_co_of_clusters)
		{
			qtawait_local_statments(hswo.doStep(););					

			//setWindowTitle(tqs(hswo.getRegionsCount())+" regions of "+tqs(min_co_of_clusters)+" desired , last dissim value = "+tqs(hswo._last_merge_dissim_value));

			/*		
			auto spClustersLabelsImage = hswo.getClustersLabelsImage();
			int* clustersLabelsImage = spClustersLabelsImage.get();

			int zoomFactor = 27;
			QImage zoomedImg = _resultImage.scaled(_imageWidth * zoomFactor, _imageHeight * zoomFactor,Qt::KeepAspectRatio);

			QPainter overlayPainter(&zoomedImg);
			for (int irow = 0; irow < _imageHeight; irow++)
			{
				for (int icol = 0; icol < _imageWidth; icol++)
				{
					int cluster_Index =  clustersLabelsImage[icol + irow*_imageWidth];

					overlayPainter.setPen(QPen(Qt::lightGray));
				
					overlayPainter.drawRect(QRectF((zoomFactor) * icol 
						, (zoomFactor) * irow
						, zoomFactor, zoomFactor));

					overlayPainter.setPen(QPen(Qt::white));

					overlayPainter.drawText(QRectF((zoomFactor) * icol +3
						, (zoomFactor) * irow +3
						, zoomFactor, zoomFactor),tr("%1").arg(cluster_Index));
				}
			}

			_resultImage = zoomedImg;
			*/

			//UpdateSolutionImage();
			//cout<<"regions = "<<hswo.getRegionsCount()<<endl;
		}	

		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));
		
		//cout<< "comp time = " << HSWO::time_count <<endl;
		//HSWO::time_count = 0;

				
		/// Drawing the image on the label box ===============================================================================
		auto spClustersLabelsImage = hswo.getClustersLabelsImage();		
		QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);
		QPainter painter(&resultImage);
		painter.fillRect(resultImage.rect(), Qt::white);

		int lastColorIndex = 1;
		map<int, int> clusterIndex_to_ColorIndex;

		for (int irow = 0; irow < _imageHeight; irow++)
		{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
				/////painter.setPen(QPen(Qt::darkGray));
				int cluster_index = spClustersLabelsImage.get()[icol + irow*_imageWidth];
				//painter.setPen(QPen(QColor::fromRgb(pixel_color, pixel_color, pixel_color,255)));

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

		/*
		int* clustersLabelsImage = new int[_imageWidth*_imageHeight];
		hswo.getClustersLabelsImage(clustersLabelsImage);

		QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);

		QPainter painter(&resultImage);
		painter.fillRect(resultImage.rect(), Qt::white);

		int lastColorIndex = 1;
		map<int, int> clusterIndex_to_ColorIndex;

		for (int irow = 0; irow < _imageHeight; irow++)
		{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
				int cluster_index = clustersLabelsImage[icol + irow*_imageWidth];
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

		delete[] clustersLabelsImage;
		clustersLabelsImage = NULL;
		*/

		/*
		while(hswo.getRegionsCount() > min_co_of_clusters)
		{
			hswo.doStep();	

			setWindowTitle(tr("%1").arg(hswo.getRegionsCount())+" regions, last dissim value = "+tr("%1").arg(hswo._last_merge_dissim_value));


			int* clustersLabelsImage = new int[_imageWidth*_imageHeight];
			hswo.getClustersLabelsImage(clustersLabelsImage);

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
					//float pixel_color =  pDataCube[icol + irow*_imageWidth + 0]*(256.0/16.5);				
					painter.setPen(QPen(Qt::darkGray));
					//painter.setPen(QPen(QColor::fromRgb(pixel_color, pixel_color, pixel_color,255)));

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
		}
		*/
	}
	catch(std::runtime_error)
	{
		QMessageBox::critical(this, tr("General runtime exception in HSWO"), "");
		std::cout<<"General runtime exception in HSWO"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in HSWO"), e.what());
		std::cout<<"Exception in HSWO : "<<e.what()<<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in HSWO"), "");
		std::cout<<"Unknown exception in HSWO"<<endl;
	}	
}

void CPUAlgorithm::ApplyRHSEG()
{
	try
	{
		if (_imageHeight <=0 || _imageWidth <= 0)
		{
			QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
			return;
		}

		/// Calculating RHSEG ==============================================
		int min_co_of_clusters = ui->spinDesiredRegions->value();
		
		clock_t startClock = clock();

		shared_ptr<RHSEG> spRhseg(new RHSEG(_imageWidth, _imageHeight, _nBands, ui->spinSpclustWeight->value(), min_co_of_clusters, ui->spinNClustersIntermediate->value(), spDataCube, ui->spinLevels->value()));				
		
		{
			qRegisterMetaType<std::shared_ptr<HyperSpectralToolbox::HSWO>>("std::shared_ptr<HyperSpectralToolbox::HSWO>");
			bool succ = QObject::connect(spRhseg.get(), SIGNAL(hswoReadyForDrawing(std::shared_ptr<HyperSpectralToolbox::HSWO>))
				, this, SLOT(on_DrawHSWO(std::shared_ptr<HyperSpectralToolbox::HSWO>)), Qt::BlockingQueuedConnection);

			QProgressBar progress(this);
			progress.setMinimum(0);
			progress.setMaximum(spRhseg->_n_recur_levels);
			progress.setFormat("%v");			
			QObject::connect(spRhseg.get(), SIGNAL(stepCompleted(int)), &progress, SLOT(setValue(int)));		

			ui->LayoutDialog->addWidget(&progress);			

			qtawait_local_statments(spRhseg->run());

			ui->LayoutDialog->removeWidget(&progress);
		}
		
		stream_to_log("return to ApplyRHSEG after RHSEG run on diff thread");
		
		auto spClustersLabelsImage = spRhseg->getClustersLabelsImage();		

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

		/*
		int zoomFactor = 27;
		QImage zoomedImg = _resultImage.scaled(_imageWidth * zoomFactor, _imageHeight * zoomFactor,Qt::KeepAspectRatio);


		QPainter overlayPainter(&zoomedImg);
		for (int irow = 0; irow < _imageHeight; irow++)
		{
		for (int icol = 0; icol < _imageWidth; icol++)
		{
		int cluster_Index =  spClustersLabelsImage.get()[icol + irow*_imageWidth];

		overlayPainter.setPen(QPen(Qt::lightGray));

		overlayPainter.drawRect(QRectF((zoomFactor) * icol 
		, (zoomFactor) * irow
		, zoomFactor, zoomFactor));

		overlayPainter.setPen(QPen(Qt::white));

		overlayPainter.drawText(QRectF((zoomFactor) * icol +3
		, (zoomFactor) * irow +3
		, zoomFactor, zoomFactor),tr("%1").arg(cluster_Index));
		}
		}

		_resultImage = zoomedImg;
		UpdateSolutionImage();
		*/

		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));

		stream_to_log("finished ApplyRHSEG");


		ifstream ifstr;
		ifstr.open("f:\\mahmoud\\Dropbox\\personal\\Thesis\\results\\paper 4\\accuracy_assesment\\files.txt");
		int n = 0;
		ifstr >> n;
		n++;
		ifstr.close();

		ofstream ofstr("f:\\mahmoud\\Dropbox\\personal\\Thesis\\results\\paper 4\\accuracy_assesment\\files.txt");
		ofstr << n;
		ofstr.close();

		char buff[200];
		itoa(n, buff, 10);
		_resultImage.save((string("f:\\mahmoud\\Dropbox\\personal\\Thesis\\results\\paper 4\\accuracy_assesment\\") + string(buff) + ".bmp").c_str(), "bmp");

		ofstream ofstr2("f:\\mahmoud\\Dropbox\\personal\\Thesis\\results\\paper 4\\accuracy_assesment\\"+string(buff)+".txt");
		ofstr2 << "classes = " << min_co_of_clusters << endl;
		ofstr2 << "levels = " << spRhseg->_n_recur_levels<< endl;
		ofstr2 << "spectral weight = " << spRhseg->_spclustWeight<< endl;
		ofstr2 << "interm clusters = " << spRhseg->_nClustersIntermediate<< endl;
		ofstr2 << "computation duration = " << (int)timeInSeconds << " seconds, " << milliSeconds <<" millisec"<< endl;
		ofstr2.close();
	}
	catch(std::runtime_error)
	{
		QMessageBox::critical(this, tr("General runtime exception in RHSEG"), "");
		std::cout<<"General runtime exception in RHSEG"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in RHSEG"), e.what());
		std::cout<<"Exception in RHSEG : "<<e.what()<<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in RHSEG"), "");
		std::cout<<"Unknown exception in RHSEG"<<endl;
	}
}

void CPUAlgorithm::ApplyGPURHSEG()
{
	try
	{
		if (_imageHeight <=0 || _imageWidth <= 0)
		{
			QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
			return;
		}

		/// Calculating RHSEG ==============================================
		int min_co_of_clusters = ui->spinDesiredRegions->value();
		
		clock_t startClock = clock();

		shared_ptr<RHSEGDevice> spRhseg(new RHSEGDevice(_imageWidth, _imageHeight, _nBands, ui->spinSpclustWeight->value(), min_co_of_clusters, ui->spinNClustersIntermediate->value(), spDataCube, ui->spinLevels->value(), ui->spinThreadsPerBlock->value()));				
		
		{
			qRegisterMetaType<std::shared_ptr<HyperSpectralToolbox::HSWO>>("std::shared_ptr<HyperSpectralToolbox::HSWO>");
			bool succ = QObject::connect(spRhseg.get(), SIGNAL(Device_hswoReadyForDrawing(std::shared_ptr<HyperSpectralToolbox::HSWO>))
				, this, SLOT(on_DrawGPUHSWO(std::shared_ptr<HyperSpectralToolbox::HSWO>)), Qt::BlockingQueuedConnection);

			QProgressBar progress(this);
			progress.setMinimum(0);
			progress.setMaximum(spRhseg->_n_recur_levels);
			progress.setFormat("%v");			
			QObject::connect(spRhseg.get(), SIGNAL(stepCompleted(int)), &progress, SLOT(setValue(int)));		

			ui->LayoutDialog->addWidget(&progress);			

			qtawait_local_statments(spRhseg->run());			

			ui->LayoutDialog->removeWidget(&progress);
		}
		
		stream_to_log("return to ApplyGPURHSEG after RHSEG run on diff thread");
		
		auto spClustersLabelsImage = spRhseg->getClustersLabelsImage();		

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

		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));

		stream_to_log("finished ApplyGPURHSEG");
	}
	catch(std::runtime_error)
	{
		QMessageBox::critical(this, tr("General runtime exception in GPU RHSEG"), "");
		std::cout<<"General runtime exception in GPU RHSEG"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in GPU RHSEG"), e.what());
		std::cout<<"Exception in GPU RHSEG : "<<e.what()<<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in GPU RHSEG"), "");
		std::cout<<"Unknown exception in GPU RHSEG"<<endl;
	}
}

void CPUAlgorithm::UpdateSolutionImage()
{
	QImage zoomedImg = _resultImage.scaled(_imageWidth * _zoomFactor, _imageHeight * _zoomFactor,Qt::KeepAspectRatio);

	/*
	if (_AMC.m_endmembers != NULL && _AMC.m_endmembers > 0 && _bDrawEndmembers)
	{
	QPainter painter(&zoomedImg);
	//		g->SmoothingMode = System::Drawing::Drawing2D::SmoothingMode::AntiAlias;	

	for (int ie = 0; ie < _nEndmembers; ie++)
	{
	int length = _zoomFactor;
	if (_zoomFactor < 3) length = 4;

	painter.setPen(QPen(QBrush(_endmembersColor), 1));
	painter.drawRect(QRectF((_zoomFactor) * _AMC.m_endmembers[ie].X 
	, (_zoomFactor) * _AMC.m_endmembers[ie].Y
	, length, length));

	painter.setPen(QPen(QBrush(_endmembersColor), 2));
	painter.drawEllipse(QRectF((_zoomFactor) * _AMC.m_endmembers[ie].X - 1.5f
	, (_zoomFactor) * _AMC.m_endmembers[ie].Y - 1.5f
	, length + 1.5f, length + 1.5f));
	}
	}
	*/
	ui->imgBox->setPixmap(QPixmap::fromImage(zoomedImg));	
}

bool CPUAlgorithm::CheckCorrectness(float *A, float *B, float* C)
{
	bool success = true;
	for(int i=0; i< CUDA_N; i++)
	{
		if(C[i] != A[i] + B[i])
		{
			success = false;
			break;
		}
	}
	return success;
}


#pragma region old AMC code

void CPUAlgorithm::ApplyAMC()
{
	if (_imageHeight <=0 || _imageWidth <= 0)
	{
		QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
		return;
	}

	time_t startTime = time(NULL);

	{
		//
		//setWindowTitle("Extracting endmembers ... "); //
		//
	}

	/// Extracting endmembers

	_AMC.ExtractEndmembers((double*)spDataCube.get(), _nBands, _imageWidth, _imageHeight, _nEndmembers);
	_nEndmembers = _AMC.m_nEndmembers;

	{
		//
		setWindowTitle(tr("%1").arg(_nEndmembers)+" endmembers Found, Unmixing using Constraint Least Squares ... "); //
		//
	}

	/// Calculating abundances

	//bool succeeded = _AMC.CalculateAbundanceFractions(pDataCube, _imageWidth, _imageHeight, N_CHANNELS);
	bool succeeded =true;
	double durationInSec = difftime(time(NULL), startTime);

	{   ////////////////////////////////////////////////////////////////////////////////////
		//---------------------------------------------------------------------------------

		setWindowTitle(tr("%1").arg(_nEndmembers)+tr(" endmembers Found, Last Computation Duration = %1 sec").arg(durationInSec)); //
		if(!succeeded)
		{
			Logger::stringLine("unmix failed!");
			QMessageBox::warning(this,tr("%1").arg(_nEndmembers) + tr(" endmembers Found, Unmixing FAILED, for more info see "),tr(Logger::getLogFileName().c_str()));
			return;
		}
		else
		{
			Logger::stringLine("unmix done!");
		}
		//---------------------------------------------------------------------------------
		////////////////////////////////////////////////////////////////////////////////////
	}


	/// Drawing the image on the label box

	QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);

	QPainter painter(&resultImage);
	painter.fillRect(resultImage.rect(), Qt::white);

	//get the maxes of each pixel
	for (int irow = 0; irow < _imageHeight; irow++)
	{
		for (int icol = 0; icol < _imageWidth; icol++)
		{
			double maxAbundValue = -INT_MAX;
			int maxIndex = 0;
			bool anEndmember = false;

			for (int e = 0; e < _nEndmembers; e++)
			{
				double v = _AMC.m_abundanceFractions[e + icol*_nEndmembers+irow*_imageWidth*_nEndmembers];
				if (maxAbundValue < v)
				{
					maxAbundValue = v;
					maxIndex = e;
				}
			}
			{
				painter.setPen(QPen(_pixelColors[maxIndex]));
				painter.drawPoint(QPoint(icol, irow));
			}
		}
	}

	/// Output endmember locations
	QString str="";
	for (int i=0; i< _nEndmembers; i++)
	{
		str += QString().sprintf("/ %d , %d ", _AMC.m_endmembers[i].X, _AMC.m_endmembers[i].Y);
	}
	this->ui->lblEndmemberLocations->setText(str);


	_resultImage = resultImage;
	this->ui->imgBox->setPixmap(QPixmap::fromImage(_resultImage));

	//if (threadDrawEndmembers != nullptr && !(threadDrawEndmembers->ThreadState == ThreadState::Running))
	//{
	//threadDrawEndmembers->Start();
	//}
}
void CPUAlgorithm::DrawEndmembers()
{
	//while(true)
	//{
	//	UpdateProcessedImage();
	//	_bDrawEndmembers = !_bDrawEndmembers;
	//	Thread::Sleep(500);
	//}
}

void CPUAlgorithm::ShowComputationProgress()
{
	//while (threadProgress->ThreadState == ThreadState::Running)
	//{
	//	//if (_AMC.m_abundanceFractions != NULL)
	//	{
	//		Int32 i = ((int)(((double)progressBar1->Value / (double)progressBar1->Maximum) * 100.0));
	//		lblProgress->Text = i.ToString() + " %";
	//		progressBar1->Value = _AMC.m_currentX + _AMC.m_currentY * _imageWidth;
	//	}
	//	Thread::Sleep(1000);
	//}
}


#pragma endregion old AMC code

#pragma region event handlers
//////////////////////////////////////////////////////////////////////////
/// Events
//////////////////////////////////////////////////////////////////////////

void CPUAlgorithm::on_btnLoadImage_clicked()
{
	LoadSampleImageFromFile();
}

void CPUAlgorithm::on_btnRun_clicked()
{
	//delete the data after drawing to imagebox
	useGPU = false;
	StartComuptation();
}

void CPUAlgorithm::on_btnStop_clicked()
{
	StopComputation();
}

void CPUAlgorithm::on_btnZoomIn_clicked()
{
	_zoomFactor ++;	
	UpdateSolutionImage();
}

void CPUAlgorithm::on_btnZoomOut_clicked()
{
	_zoomFactor --;
	if(_zoomFactor < 1) _zoomFactor = 1;	
	UpdateSolutionImage();
}

void CPUAlgorithm::on_btnRunCUDA_clicked()
{
	/* 
	// input data
 	int len = 16;
    // the data has some zero padding at the end so that the size is a multiple of
    // four, this simplifies the processing as each thread can process four
    // elements (which is necessary to avoid bank conflicts) but no branching is
    // necessary to avoid out of bounds reads
    char str[] = { 82, 111, 118, 118, 121, 42, 97, 121, 124, 118, 110, 56,
                   10, 10, 10, 10};

    // Use int2 showing that CUDA vector types can be used in cpp code
    int2 i2[16];
    for( int i = 0; i < len; i++ )
    {
        i2[i].x = str[i];
        i2[i].y = 10;
    }

    // run the device part of the program
    //runTest(NULL, NULL, str, i2, len);	

	QMessageBox::information(this, "CUDA Output ...", tr(str));
	QString strOutput = "";
    for( int i = 0; i < len; i++ )
    {        
		strOutput += (char)(i2[i].x);
	}
	QMessageBox::information(this, "CUDA Output ...", strOutput);
	*/


	//float* host_A = new float[CUDA_N];
	//float* host_B = new float[CUDA_N];
	//float* host_C = new float[CUDA_N];

	//for(int i=0; i< CUDA_N; i++)
	//{
	//	host_A[i] = 1;
	//	host_B[i] = 3;
	//	host_C[i] = 0;
	//}


	//runTest2(host_A, host_B, host_C);
	//if(CheckCorrectness(host_A, host_B, host_C))
	{
		//QMessageBox::information(this, "Correct", "correct");
	}

	//delete host_A;
	//delete host_B;
	//delete host_C;
	
	/*
	int N = 1<<20;

	thrust::host_vector<int> h_vec(N), h_vec_2(N);	

	// generate 16M random numbers on the host	
    thrust::generate(h_vec.begin(), h_vec.end(), rand);
	h_vec_2=h_vec;

	ui->lblOutputLog->setText(ui->lblOutputLog->text()+tr("Number of random elements = %1\n").arg(h_vec.size()));

	
	time_t t = time(0);
	ui->lblOutputLog->setText(ui->lblOutputLog->text()+tr("\nStarted CPU Sorting vector ..."));
	thrust::sort(h_vec.begin(), h_vec.end());
	ui->lblOutputLog->setText(ui->lblOutputLog->text()+tr("\nDone CPU Sorting vector after ")+tr("%1 milliSeconds\n").arg(time(0)-t));
	

	
	time_t tGPU = time(0);
	ui->lblOutputLog->setText(ui->lblOutputLog->text()+tr("\nStarted GPU Sorting vector ...\n"));	
	runThrustExpr(h_vec_2,N);
	ui->lblOutputLog->setText(ui->lblOutputLog->text()+tr("\nDone GPU Sorting vector after ")+tr("%1 milliSeconds").arg(time(0)-tGPU));
	*/
	
	useGPU = true;
	StartComuptation();
}

void CPUAlgorithm::on_btnHSWOStep_clicked()
{	
}

void CPUAlgorithm::onClicked_btnRHSEG()
{	
	ApplyRHSEG();	
	UpdateSolutionImage();
}

void CPUAlgorithm::on_btnGPURHSEG_clicked()
{		
	ApplyGPURHSEG();
	UpdateSolutionImage();
	
	//qtawait_local_statments(do_kernel());
}


void CPUAlgorithm::on_DrawHSWO(shared_ptr<HSWO> spHSWO)
{
	stream_to_log("inside on_DrawHSWO slot");	
	
	auto spClustersLabelsImage = spHSWO->getClustersLabelsImage();

	QImage resultImage(QSize(spHSWO->m_image_width, spHSWO->m_image_height), QImage::Format_RGB32);

	QPainter painter(&resultImage);
	painter.fillRect(resultImage.rect(), Qt::white);

	int lastColorIndex = 1;
	map<int, int> clusterIndex_to_ColorIndex;

	for (int irow = 0; irow < spHSWO->m_image_height; irow++)
	{
		for (int icol = 0; icol < spHSWO->m_image_width; icol++)
		{
			int cluster_index = spClustersLabelsImage.get()[icol + irow*spHSWO->m_image_width];

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
	stream_to_log("finished slot");	
}

void CPUAlgorithm::on_DrawGPUHSWO(shared_ptr<HSWO> spHSWO)
{
	stream_to_log("inside on_DrawGPUHSWO slot");	

	auto spClustersLabelsImage = spHSWO->getClustersLabelsImage();

	QImage resultImage(QSize(spHSWO->m_image_width, spHSWO->m_image_height), QImage::Format_RGB32);

	QPainter painter(&resultImage);
	painter.fillRect(resultImage.rect(), Qt::white);

	int lastColorIndex = 1;
	map<int, int> clusterIndex_to_ColorIndex;

	for (int irow = 0; irow < spHSWO->m_image_height; irow++)
	{
		for (int icol = 0; icol < spHSWO->m_image_width; icol++)
		{
			int cluster_index = spClustersLabelsImage.get()[icol + irow*spHSWO->m_image_width];

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
	stream_to_log("finished slot");	
}

void CPUAlgorithm::on_btnTest1_clicked()
{
	try
	{
		if (_imageHeight <=0 || _imageWidth <= 0)
		{
			QMessageBox::warning(this, tr("No image is loaded !"), tr("please press load image button"));
			return;
		}

		/// Calculating HSWO ==============================================
		int min_co_of_clusters = ui->spinDesiredRegions->value();

		HSWO hswo(_imageWidth, _imageHeight, _nBands, ui->spinSpclustWeight->value(), min_co_of_clusters, spDataCube);

		clock_t startClock = clock();  

		hswo.InitializeHost();  

		while(hswo.getRegionsCount() > 512)
		{
			qtawait_local_statments(hswo.doStep(););					

			setWindowTitle(tqs(hswo.getRegionsCount())+" regions of "+tqs(min_co_of_clusters)+" desired , last dissim value = "+tqs(hswo._last_merge_dissim_value));			
		}	

		hswo.TO_GPU_HSWO(ui->spinThreadsPerBlock->value());

		while(hswo.getRegionsCount() > min_co_of_clusters)
		{
			qtawait_local_statments(hswo.Device_doStep(););					

			setWindowTitle(tr("%1").arg(hswo.getRegionsCount()) +" regions"
				+". last dissim value = " +tr("%1").arg(hswo._last_merge_dissim_value)
				+". TBP = "+tr("%1").arg(hswo.m_deviceData.m_threadsPerBlock)
				+". last #Blocks = "+tr("%1").arg(hswo._last_block_count));				
		}	
		DeviceExit(&hswo);

		clock_t endClock = clock();
		float timeInSeconds = (float)(endClock - startClock)/CLOCKS_PER_SEC;
		int milliSeconds = ( float((float)timeInSeconds - (int)(timeInSeconds))  )*1000.0f;

		ui->lblOutputLog->setText(tr("Time = %1 seconds , %2 milliseconds").arg((int)timeInSeconds).arg(milliSeconds));
				
		/// Drawing the image on the label box ===============================================================================
		auto spClustersLabelsImage = hswo.getClustersLabelsImage();		
		QImage resultImage(QSize(_imageWidth, _imageHeight), QImage::Format_RGB32);
		QPainter painter(&resultImage);
		painter.fillRect(resultImage.rect(), Qt::white);

		int lastColorIndex = 1;
		map<int, int> clusterIndex_to_ColorIndex;

		for (int irow = 0; irow < _imageHeight; irow++)
		{
			for (int icol = 0; icol < _imageWidth; icol++)
			{
				/////painter.setPen(QPen(Qt::darkGray));
				int cluster_index = spClustersLabelsImage.get()[icol + irow*_imageWidth];
				//painter.setPen(QPen(QColor::fromRgb(pixel_color, pixel_color, pixel_color,255)));

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
		QMessageBox::critical(this, tr("General runtime exception in HSWO"), "");
		std::cout<<"General runtime exception in HSWO"<<endl;
	}
	catch(std::exception e)
	{
		QMessageBox::critical(this, tr("Exception in HSWO"), e.what());
		std::cout<<"Exception in HSWO : "<<e.what()<<endl;
	}
	catch(...)
	{
		QMessageBox::critical(this, tr("Unknown exception in HSWO"), "");
		std::cout<<"Unknown exception in HSWO"<<endl;
	}	
}

void CPUAlgorithm::on_btnTest2_clicked()
{
#if USE_CUDA
	HyperSpectralToolbox::run_expr();
#endif
}

#pragma endregion event handlers