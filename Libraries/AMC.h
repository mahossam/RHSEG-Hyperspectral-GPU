#pragma once
#include <map>
#include "./matrix/newmat.h"
#include "./matrix/newmatrc.h"
#include "./matrix/newmatio.h"
#include "./definitions.h"

namespace HyperSpectralToolbox
{
	//////////////////////////////////////////////////////////////////////////
	/// AMC algorithm , based on Erosion and Dilation operator from Mathematical Morphology for endmember extraction
	/// and LeastSquaresMethod for Spectral Linear Unmixing
	//////////////////////////////////////////////////////////////////////////
	class Utilities
	{
	public:
		double static getSpectralAngle(MatrixCol& a, MatrixCol& b);
	};

	class AMC 
	{
	public:
		class Point2D
		{
		public: Point2D(int x, int y)
				{
					X = x; Y = y;
				}
		public: Point2D()
				{
					X = 0; Y = 0;
				}
		public: int X, Y;
		};

		Point2D* m_endmembers;

		double* m_abundanceFractions;
		int m_currentX , m_currentY, m_nEndmembers;  //for progress monitoring
		std::map<double,Point2D> all_endmembers;

	protected: short* se;
			   int seSize;
	protected: 
		double* pSummations;  //height*  width
		int m_width, m_height;


	private: int srcStride ;
			 int m_nBands;
	private: double* m_pDataCube;

	public: AMC() ;
	public: ~AMC() ;

	public: AMC(short* se);

	public:
		//////////////////////////////////////////////////////////////////////////
		/// Extracts Endmembers from given image at dataCube
		/// dataCube is nBands * width * height
		//////////////////////////////////////////////////////////////////////////			
		void ExtractEndmembers(double* dataCube, int nBands, int width, int height, int max_endmembers);

	public:
		//////////////////////////////////////////////////////////////////////////
		/// The function is used to calculates the abundance images for the given "datacube"
		/// the output abundance array (m_abundanceFractions) dimensions are (height(or nRows) * width(or nCols) * nEndmembers)
		/// Note: "dataCube" must be (bands * height * width )
		//////////////////////////////////////////////////////////////////////////
		bool CalculateAbundanceFractions(double* dataCube, int width, int height, int nBands);

	private: inline double SAD(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY);
	private: inline double SID(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY);
	private: inline double mixed_SID_SAD(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY);


	protected: 	
		class LeastSquares
		{

		public:
			LeastSquares();				 
			
		private:
			// last row value in Matrix and last b vector element
			// for constraint Sum xi = 1 (GAMMA=weight)
			int GAMMA ;  


		private: double find_max(double x, double y);
		private: double round(double x);

		public:
			//////////////////////////////////////////////////////////////////////////
			/// The function calculates the abundance images for the given datacube and endmembers of the image
			/// the abundance array dimensions are nRows * nCols * nEndmembers // (height * width * nEndmembers)
			/// dataCube must be (bands * height * width )
			//////////////////////////////////////////////////////////////////////////
			bool LeastSquaresMinimization(double* dataCube, int width, int height, int nBands, int nEndmembers, Point2D* endmembers, double* abundanceFractions, int& currentX, int& currentY);
		
		private: bool doLeastSquares(Matrix& A, double* dataCube, int width, int height, int nBands, int nEndmembers, double* abundanceFractions, int& currentX, int& currentY);

		};
	};
}


