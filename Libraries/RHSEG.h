#pragma once
#include <QObject>
#include <map>
//#include <list>
#include <set>
#include <string>
#include <vector>
#include <QFuture>

#include "./definitions.h"
#include "myglobals.h"

#include "HSWO.h"

//struct MyFutureResult{
//	std::shared_ptr<HyperSpectralToolbox::HSWO> spHswo;
//	QFuture<void> future;
//};

namespace HyperSpectralToolbox
{

	class RHSEG: public QObject
	{	
		Q_OBJECT
	public:
		RHSEG(int nCols, int nRows, int nBands, float spclustWeight, int minNoOfClusters, int nClustersIntermediate, std::shared_ptr<float> spDataCube, int nRecursiveLevels);
		~RHSEG(void);		

		void InitializeHost();
		void run();		

		int getRegionsCount();

		std::shared_ptr<int> getClustersLabelsImage();
		float _last_merge_dissim_value;
		int _last_image_number;		
		
		float _spclustWeight;
		int _nClustersIntermediate;
		int _n_recur_levels;


		int _min_no_of_clusters;
		const int _nBands;
		int _imageWidth, _imageHeight;

		signals:
		void hswoReadyForDrawing(std::shared_ptr<HyperSpectralToolbox::HSWO>);
		void stepCompleted(int);

	
	private:
		//std::shared_ptr<HSWO> intermediateStep(int region_label_offset, int recur_level
		MyFutureResult intermediateStep(int region_label_offset, int recur_level
			, int section_no, const int ncols, const int nrows, std::shared_ptr<float> up_dataCube);
		std::shared_ptr<HSWO> merge_splits(std::vector<std::shared_ptr<HSWO>>& quarters);

		std::shared_ptr<float> pDataCube;			

		std::shared_ptr<HSWO> _final_result;		
	};

}
