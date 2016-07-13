#ifndef CPUALGORITHM_H
#define CPUALGORITHM_H

#include <QtGui/QDialog>

#include "Libraries/AMC.h" 
#include "Libraries/HSWO.h" 
#include "Libraries/DebuggingUtilities.h"
#include <QtCore/QThread>

namespace Ui 
{
	class CPUAlgorithm ;
}

class CPUAlgorithm : public QDialog {
	Q_OBJECT

public:

	//QThread threadLoadChannels;
	//QAlgorithmThread* _pAlgorithmThread;
	bool useGPU;
	//QThread* threadProgress;
	//QThread* threadDrawEndmembers;

protected:
	int _zoomFactor;
	int _imageWidth;
	int _imageHeight;
	int _nBands;
	int _nEndmembers;
	bool _bDrawEndmembers;
	shared_ptr<float> spDataCube;

	QColor _endmembersColor;
	QColor _pixelColors[200];
	QImage _resultImage;

	HyperSpectralToolbox::AMC _AMC;	

private:
	void LoadSampleImageFromFile();
	void StartComuptation();
	void DrawSolution();

	void UpdateSolutionImage();
	void DrawEndmembers();
	void ShowComputationProgress();

	void StopComputation();
public:
	CPUAlgorithm(QWidget *parent = 0);

	void InitializeForm();
	void InitializeApp();
	void ApplyAMC();	
	void ApplyHSWO();
	void ApplyGPUHSWO();
	void ApplyRHSEG();
	void ApplyGPURHSEG();

	~CPUAlgorithm();

	public slots:
		void on_DrawHSWO(std::shared_ptr<HyperSpectralToolbox::HSWO> spHSWO);		
		void on_DrawGPUHSWO(std::shared_ptr<HyperSpectralToolbox::HSWO> spHSWO);		
		//template<class Typy>
		//void DefaultSlot();		

protected:
	void changeEvent(QEvent *e);

private:
	Ui::CPUAlgorithm *ui;
	bool CheckCorrectness(float * A, float* B, float* C);

	private slots:
		void on_btnStop_clicked();
		void on_btnLoadImage_clicked();
		void on_btnRun_clicked();
		void on_btnZoomIn_clicked();
		void on_btnZoomOut_clicked();	
		void on_btnRunCUDA_clicked();
		void on_btnHSWOStep_clicked();		
		void onClicked_btnRHSEG();		
		void on_btnGPURHSEG_clicked();	
		void on_btnTest1_clicked();
		void on_btnTest2_clicked();
};

#endif // CPUALGORITHM_H
