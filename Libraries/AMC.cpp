#include "AMC.h"
#include <cmath>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "./Libraries/definitions.h"
#include "./Libraries/Logger.h"
#include "./Libraries/DebuggingUtilities.h"

#include "myglobals.h"

using namespace HyperSpectralToolbox;
using namespace std;
//short global_se[9*9] = {  1,1,1,1,1,1, 1,1, 1 ,  1,1,1,1,1,1, 1,1, 1 ,  1,1,1,1,1,1, 1,1, 1 \
//		,1,1,1,1,1,1, 1,1, 1    ,1,1,1,1,1,1, 1,1, 1    ,1,1,1,1,1,1, 1,1, 1   \
//		,1,1,1,1,1,1, 1,1, 1    ,1,1,1,1,1,1, 1,1, 1    ,1,1,1,1,1,1, 1,1, 1    };

short global_se[3*3]={  1, 1, 1 ,  1, 1, 1 ,  1, 1, 1  };
#define PI 3.1415926535897932384626433832795
#define LEAST_SQUARES_TOLERANCE .001
#define SAD_TOLERANCE .005
#define NCANDIDATE_ENDMEMBERS(nClasses) (nClasses+12)

double Utilities::getSpectralAngle(MatrixCol& a, MatrixCol& b)
{
	return acos(DotProd(a, b) / (a.Norm2()*b.Norm2()));
}


AMC::AMC()
{
	se = &global_se[0];
	seSize=3;
	m_pDataCube = NULL;
	m_nBands = 0;
	srcStride = -1;
	m_currentX = m_currentY = 0;
	m_abundanceFractions = NULL;
	m_endmembers =  NULL;
	m_nEndmembers= 0;	
	pSummations = NULL;
}


AMC::AMC(short* se)
{
	this->se = se;
}

void AMC::ExtractEndmembers(double* dataCube, int nBands, int width, int height, int max_endmembers)
{
	//for debugging in log file
	IntegerDebuggingUtilities::getInstance().AddArray("array_bad_channels");

	//====================== processing
	// process the filter
	m_nBands = nBands;
	m_pDataCube = dataCube;
	m_width = width;
	m_height = height;

	if(m_endmembers != NULL)
	{
		delete[] m_endmembers;
		m_endmembers = NULL;
	}

	int pixelSize = 1 ;

	// processing start and stop X,Y positions
	int startX = 0;
	int startY = 0;
	int stopX = startX + width;
	int stopY = startY + height;

	// loop and array indexes
	int t, ir, jr, i, j, ci, cir, cj, cjr;
	// structuring element's radius
	int r = seSize >> 1;

	/// do the job

	//initialize Summations needed for SID
	//
	//double *pSummations=&summations[0];
	pSummations=new double[width*height];
	double *pMEI= new double[width*height];

	for (int y = startY; y < stopY; y++)
	{
		for (int x = startX; x < stopX; x++)
		{
			pSummations[x+y*width] = 0;
			pMEI[x+y*width]=0;
		}
	}

	double* src = (double*)dataCube;
	for (int l = 0; l < nBands; l++)
	{
		for (int y = startY; y < stopY; y++)
		{
			for (int x = startX; x < stopX; x++)
			{
				double* Value = src + (x  + y * width + l * height * width);
				pSummations[x+y*width] += *Value;
			}
		}
	}


	double maxComulativeDist = -INT_MAX, minComulativeDist = INT_MAX, comulativeDist = 0.0;
	int minVectorX=0, minVectorY=0, maxVectorX=0, maxVectorY=0;

	// for each line
	for (int y = startY; y < stopY; y++)
	{
		// for each pixel
		for (int x = startX; x < stopX; x++)
		{
			m_currentX = x;
			m_currentY = y;

			// for each structuring element's row
			for (i = 0; i < seSize; i++)
			{
				ir = i - r;
				t = y + ir;

				// skip row
				if (t < startY)
					continue;
				// break
				if (t >= stopY)
					break;

				// for each structuring element's column
				for (j = 0; j < seSize; j++)
				{
					jr = j - r;
					t = x + jr;

					// skip column
					if (t < startX)
						continue;
					if (t < stopX)
					{
						if (se[i+j*seSize] == 1)
						{

							comulativeDist = 0;
							// for each structuring element's row (to compute the cumulative Distance for dataCube(ir, jr)
							for (ci = 0; ci < seSize; ci++)
							{
								cir = ci - r;
								t = y + cir;

								// skip row
								if (t < startY)
									continue;
								// break
								if (t >= stopY)
									break;

								// for each structuring element's (to compute the cumulative Distance for dataCube(ir, jr)
								for (cj = 0; cj < seSize; cj++)
								{
									cjr = cj - r;
									t = x + cjr;

									// skip column
									if (t < startX)
										continue;
									if (t < stopX)
									{
										if (cir == ir && cjr == jr) continue;
										//comulativeDist += mixed_SID_SAD(ref srcData, x + cjr, y + cir, x + jr, y + ir); //src[ir * stride + jr];  //TODO
										comulativeDist += SID(x + cjr, y + cir, x + jr, y + ir); //src[ir * stride + jr];  //TODO
										//comulativeDist += SAD(x + cjr, y + cir, x + jr, y + ir); //src[ir * stride + jr];  //TODO
									}
								}
							}
							if (comulativeDist > maxComulativeDist)
							{
								maxComulativeDist = comulativeDist;
								maxVectorY = y + ir;
								maxVectorX = x+jr;
							}
							if (comulativeDist < minComulativeDist)
							{
								minComulativeDist = comulativeDist;
								minVectorY = y+ir;
								minVectorX = x+jr;
							}
						}
					}
				}
			}

			//double value = SAD(maxVectorX, maxVectorY, minVectorX, minVectorY);
			double value = SID(maxVectorX, maxVectorY, minVectorX, minVectorY);

			//updating endmembers
			pMEI[maxVectorX+maxVectorY*width] += value;
			//pMEI[x+y*width] += value;

			// reinitialize distances
			maxComulativeDist = -INT_MAX;
			minComulativeDist = INT_MAX;
		}

	}


	for (int y = startY; y < stopY; y++)
	{
		for (int x = startX; x < stopX; x++)
		{
			all_endmembers.insert(std::make_pair(pMEI[x+y*width],Point2D(x,y)));
		}
	}
	
	Point2D* candidate_endmembers = new Point2D[NCANDIDATE_ENDMEMBERS(max_endmembers)];

	map<double,Point2D>::reverse_iterator it = all_endmembers.rbegin();
	for (int i=0; it != all_endmembers.rend() && i< NCANDIDATE_ENDMEMBERS(max_endmembers); it++,i++)
	{
		candidate_endmembers[i] = it->second;
	}

	int nFoundEndmembers = NCANDIDATE_ENDMEMBERS(max_endmembers);
	for(int iEnd1 =0; iEnd1 <NCANDIDATE_ENDMEMBERS(max_endmembers); iEnd1++)
	{
		if(candidate_endmembers[iEnd1].X == -1 ) continue;
		//choosing only unique endmembers by measuring the SAD between every vector pair
		for(int iEnd2 =iEnd1+1; iEnd2 < NCANDIDATE_ENDMEMBERS(max_endmembers); iEnd2++)
		{
			if(iEnd1 == iEnd2) continue;
			if(candidate_endmembers[iEnd2].X == -1) continue;
			if(SID(candidate_endmembers[iEnd1].X,candidate_endmembers[iEnd1].Y, candidate_endmembers[iEnd2].X, candidate_endmembers[iEnd2].Y) < SAD_TOLERANCE)
			{
				candidate_endmembers[iEnd2].X = -1;
				candidate_endmembers[iEnd2].Y = -1;
				nFoundEndmembers--;
			}
		}
	}

	m_endmembers = new AMC::Point2D[nFoundEndmembers];
	m_nEndmembers = nFoundEndmembers;

	//copy the endmembers to the final array:
	for(int ice =0, e=0; ice <NCANDIDATE_ENDMEMBERS(max_endmembers); ice++)
	{
		if(candidate_endmembers[ice].X == -1 ) continue;
		m_endmembers[e].X = candidate_endmembers[ice].X;
		m_endmembers[e].Y = candidate_endmembers[ice].Y;
		e++;
	}


	delete[] pMEI;
	delete[] pSummations;
	m_currentX = 0;
	m_currentY = 0;

	{
		//
		//Logger::stringLine("=======");
		//Logger::stringLine("bad channels : ");
		//std::set<int>* pArr = IntegerDebuggingUtilities::getInstance().GetArray("array_bad_channels");
		//for(std::set<int>::iterator it=pArr->begin(); it != pArr->end(); it++)
		//{
		//Logger::decimal(*it);Logger::sstring(", ");
		//}
		//Logger::endl();
		//Logger::stringLine("=======");		
	}
}

bool AMC::CalculateAbundanceFractions(double* dataCube, int width, int height, int nBands)
{
	Logger::stringLine("CalculateAbundanceFractions ...");
	AMC::LeastSquares ls;
	if(m_abundanceFractions != NULL)
	{
		delete[] m_abundanceFractions;
		m_abundanceFractions = NULL;
	}
	if(m_abundanceFractions == NULL) m_abundanceFractions = new double[height * width * m_nEndmembers];
	bool succeeded = ls.LeastSquaresMinimization(dataCube, width, height,nBands,  m_nEndmembers, m_endmembers, m_abundanceFractions, m_currentX, m_currentY);

	if (succeeded)
	{	
		Logger::stringLine("CalculateAbundanceFractions finished successfully!");
	}
	else
	{
		Logger::stringLine("CalculateAbundanceFractions FAILED !");
	}

	return succeeded;
}

inline double AMC::SAD(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY)
{
	double dotProduct = 0, vec1Square = 0, vec2Square = 0;
	double* src = m_pDataCube;
	for (int i = 0; i < m_nBands; i++)
	{
		double maxValue = *(src + (maxVectorX  + maxVectorY * m_width + i * m_height * m_width));

		double minValue = *(src + (minVectorX  + minVectorY * m_width + i * m_height * m_width));

		dotProduct += maxValue * minValue;
		vec1Square += maxValue * maxValue;
		vec2Square += minValue * minValue;
	}

	return acos((double) dotProduct/(sqrtf(vec1Square)*sqrtf(vec2Square)));
}


inline double AMC::SID(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY)
{
	double maxVectorSum = pSummations[maxVectorX+ maxVectorY*m_width], minVectorSum = pSummations[minVectorX+ minVectorY*m_width];
	double* src = m_pDataCube;
	double sid = 0.0;

	for (int i = 0; i < m_nBands; i++)
	{
		double maxValue = *(src + (maxVectorX  + maxVectorY * m_width + i * m_height * m_width));
		double minValue = *(src + (minVectorX  + minVectorY * m_width + i * m_height * m_width));

		double pl = (double)maxValue / maxVectorSum;
		double ql = (double)minValue / minVectorSum;

		if(pl <= 0.0 || ql <= 0.0)
		{
			LOG_DISABLE
				Logger::stringLine("----");
				Logger::sstring("c: ");Logger::decimal(i);
				Logger::endl();Logger::sstring("!! pl = ");Logger::decimal(pl);Logger::endl();
				Logger::endl();Logger::sstring("!! ql = ");Logger::decimal(ql);Logger::endl();
				Logger::stringLine("----");
			LOG_END
				IntegerDebuggingUtilities::getInstance().AddArrayValue("array_bad_channels", i+1);
			continue;
		}
		sid += (  (pl * log10f((float)(pl / ql)))        + (ql * log10f((float)(ql / pl)))  );
	}
	return sid;
}

inline double AMC::mixed_SID_SAD(int maxVectorX, int maxVectorY, int minVectorX, int minVectorY)
{
	return SID(maxVectorX, maxVectorY, minVectorX, minVectorY) * tan(SAD(maxVectorX, maxVectorY, minVectorX, minVectorY));
}


AMC::~AMC(void)
{
	if(m_abundanceFractions != NULL)
	{
		delete[] m_abundanceFractions;
		m_abundanceFractions = NULL;
	}

	if(m_endmembers != NULL)
	{
		delete[] m_endmembers;
		m_endmembers = NULL;
	}

	if(pSummations != NULL)
	{
		delete[] pSummations;
		pSummations = NULL;
	}

}


//////////////////////////////////////////////////////////////////////////
/// LeastSqaures Method
//////////////////////////////////////////////////////////////////////////


AMC::LeastSquares::LeastSquares()
{
	GAMMA = 10;//.0;
}

double AMC::LeastSquares::find_max(double x, double y)
{
	return (x > y ? x : y);
}

double AMC::LeastSquares::round(double x)
{
	double n;

	if (x >= 0.0)
		n = x + .5;
	else
	{
		n = -x + .5;
		n = -n;
	}
	return n;
}

bool AMC::LeastSquares::LeastSquaresMinimization(double* dataCube, int width, int height, int nBands, int nEndmembers, AMC::Point2D* endmembers, double* abundanceFractions, int& currentX, int& currentY)
{
	Logger::stringLine("LeastSquaresMinimization ...");

	Matrix A(nBands, nEndmembers);

	Logger::stringLine("\nThe Matrix A = ");
	char buf[10]="";
	for (int e = 0; e < nEndmembers; e++)
	{
		for (int b = 0; b < nBands; b++)
		{
			Point2D point = endmembers[e];
			A.element(b, e) = dataCube[point.X+ point.Y*width+b*height*width];
			gcvt(A.element(b, e), 4, buf);
			Logger::sstring(buf);
			Logger::sstring(" ");
		}
		Logger::stringLine("");			
	}
	
	

	bool succeeded = doLeastSquares(A, dataCube, width, height, nBands, nEndmembers, abundanceFractions, currentX, currentY);

	if (succeeded)
	{	
		Logger::stringLine("LeastSquaresMinimization successful !");
	}
	else
	{
		Logger::stringLine("LeastSquaresMinimization FAILED !");
	}
	return succeeded;
}

bool AMC::LeastSquares::doLeastSquares(Matrix& A, double* dataCube, int width, int height, int nBands, int nEndmembers, double* abundanceFractions, int& currentX, int& currentY)
{
	Logger::stringLine("doLeastSquares ...");
	currentX = 0;
	currentY = 0;

	int row, col;
	int nrows = height, ncols = width;
	int band;
	int i, j, k, iterations;
	Matrix* A_tilde;

	int index = -1;
	double max1, max2, max_total = 0.0;
	double change, mu, deviation;

	double anglefield[255*255];

	bool bError=false;

	// here we go... 

	// ATTENTION: Internally we work here with col-oriented matrixfile,
	// but the user has to enter the spectra row-wise for his/her's
	// convenience...  That means: Don't mix row- and col-orientation
	// in source code and modules messages output!
	//
	//Spectral Matrix is stored in A now (diagonally flipped to input
	// file) Generally:    n: cols ....  for matrix A
	//             m: rows
	//                 |
	//                 |
	//

	// 1. Check matrix orthogonality:
	//    Ref: Youngsinn Sohn, Roger M. McCoy 1997: Mapping desert shrub
	//    rangeland using spectral unmixing and modeling spectral
	//    mixtrues with TM data. Photogrammetric Engineering &
	//    Remote Sensing,  Vol.63,  No6.
	//
	//
	// 2. Beside checking matrix orthogonality we find out the maximum
	//    entry of the matrix for configuring stepsize mu later.  

	double curr_angle = 0.0;

	for (i = 0; i < A.Ncols(); i++) // go columnwise through matrix
	{
		MatrixCol Avector1(&A,LoadOnEntry,i);
		max1 = Avector1.Maximum1(INT_MIN, index);  // get the max. element of this vector 
		for (j = 0; j < A.Ncols(); j++)
		{
			if (j != i)
			{
				MatrixCol Avector2(&A,LoadOnEntry,j);  // get next col in A 
				max2 = Avector2.Maximum1(INT_MIN, index); // get the max. element of this vector 
				max_total = find_max(max1, max2); // find max of matrix A 

				curr_angle = Utilities::getSpectralAngle(Avector1, Avector2);                 // check vector angle 
				anglefield[i+ j*255] = ((float)curr_angle / PI) * 180;     // save angle in degree 
			}
		}
	}

	Logger::stringLine("  -Checking linear dependencies (orthogonality check) of Matrix A");

	// print out the result 
	for (i = 0; i < A.Ncols(); i++)
	{
		for (j = 0; j < A.Ncols(); j++)
		{
			if (j != i)
			{
				char buff[500];
				// internally this is col and not row certainly 
				sprintf(buff, "		-Angle between row %d and row %d: %f degree",(i + 1), (j + 1), anglefield[i+ j*255]);
				Logger::stringLine(buff);
			}

		}
	}

	// check it 
	//Logger::stringLine("Checking Orthogonality");
	//bError = false;
	//
	//for (i = 0; i < A.Ncols(); i++)
	//{
	//	for (j = 0; j < A.Ncols(); j++)
	//	{
	//		if (j != i)
	//		{
	//			if (anglefield[i+ j*255] < 8.0) //was 8.0 in the original algorithm of GRASS
	//			{
	//				char buff[500];
	//				// internally this is col and not row certainly
	//				sprintf(buff, "		-ERROR: Spectral entries row %d: and row %d: in your matrix are linear dependent!",(i+1), (j+1));

	//				Logger::stringLine(buff);
	//				Logger::stringLine("		-You have to revise your reference spectra");
	//				bError = true;
	//			}
	//		}
	//	}
	//}
	//


	if (!bError)
	{
		Logger::stringLine("Spectral matrix is o.k. Proceeding...");
	}
	else
	{
		Logger::stringLine("doLeastSquares FAILED !");
		return false;
	}
	//if (!flag.quiet->answer)
	//Console.WriteLine("Spectra matrix is o.k. Proceeding...\n\n");



	/// Begin calculations  //////////////////////////////////////////////////////////////////////////
	// 1. Constraint SUM xi = 1
	//   add last row "1" elements to Matrix A, store in A_tilde
	//   A_tilde is one row-dimension more than A 
	A_tilde = new Matrix(A.Nrows() + 1, A.Ncols());  // memory allocation 

	for (i = 0; i < A.Nrows(); i++)   // copy row wise 
		for (j = 0; j < A.Ncols(); j++)  // copy col wise 
			A_tilde->element(i, j) = A.element(i, j);

	// fill last row with 1 elements 
	for (j = 0; j < A.Ncols(); j++)
		(*A_tilde).element(A.Nrows(), j) = (double)(GAMMA);

	/// ---- now we have an overdetermined (non-square) system 
	// We have a least square problem here: error minimization
	//                             T          -1         T
	// unknown fraction = [A_tilde * A_tilde]  * A_tilde * b
	//
	// A_tilde is the non-square matrix with first constraint in last row.
	// b is pixel vector from satellite image
	//
	// Solve this by deriving above equation and searching the
	// minimum of this error function in an iterative loop within
	// both constraints.


	// calculate the transpose of A_tilde
	Matrix A_tilde_trans = A_tilde->t();

	// initialize some values 
	// step size must be small enough for convergence  of iteration:
	//  mu=0.000001;      step size for spectra in range of W/m^2/um
	//  mu=0.000000001;   step size for spectra in range of mW/m^2/um
	//  mu=0.000001;      step size for spectra in range of reflectance
	////
	// check  max_total for number of digits to configure mu size


	ColumnVector* fraction = NULL;

	//mu = 0.000001 * pow(10, -1 * ceil(log10f((float)max_total)));
	mu = 0.000000001;

	ColumnVector startvector (A.Ncols());                // length: no. of spectra 

	ColumnVector A_times_startvector (A_tilde->Nrows());  // length: no. of bands   
	//MultipliedMatrix* A_times_startvector ;				// length: no. of bands   

	ColumnVector errorvector (A_tilde->Nrows());          // length: no. of bands   
	//SubtractedMatrix* errorvector ;

	ColumnVector temp (A_tilde->Ncols());                 // length: no. of spectra 
	//MultipliedMatrix* temp;

	Matrix A_tilde_trans_mu(A_tilde->Ncols(), A_tilde->Nrows());
	//ScaledMatrix* A_tilde_trans_mu;

	// Now we can calculate the fractions pixelwise 
	//nrows // get geographical region 
	//ncols

	Logger::sstring("mu = ");Logger::decimal(mu);Logger::endl();

	for (row = 0; row < nrows; row++)             // rows loop in images 
	{
		//if (!flag2.veryquiet->answer)
		//G_percent(row, nrows, 1);
		//
		//for (band = 0; band < nBands; band++) //get one row for all bands
		//{
		//if (G_get_map_row (cellfd[band], cell[band], row) < 0)
		//Application.Exit();
		//}

		for (col = 0; col < ncols; col++) // cols loop, work pixelwise for all bands 
		{
			currentX = col;
			currentY = row;

			// get pixel values of each band and store in b vector:
			ColumnVector b_gamma(A_tilde->Nrows());              // length: no. of bands + 1 (GAMMA)
			for (band = 0; band < nBands; band++)
				b_gamma.element(band) = (double)(dataCube[col+ row*width+band*height*width]);
			// add GAMMA for 1. constraint as last element
			b_gamma.element(nBands) = (double)(GAMMA);

			{
				//
				Logger::endl();
				Logger::sstring("  row,  col\n");
				Logger::sstring("( ");Logger::decimal(row+1);Logger::sstring(", ");Logger::decimal(col+1);Logger::sstring(") ,   pixelVector = ");//Logger::openStream() << b_gamma;Logger::closeStream();
				Logger::endl();
				Logger::stringLine("-----------------------------------------");				
			}

			// calculate fraction vector for current pixel
			// Result is stored in fractions vector
			// with second constraint: Sum x_i = 1

			change = 1000;  // initialize 
			deviation = 1000;
			iterations = 0;
			for (k = 0; k < (A_tilde->Ncols()); k++)  // no. of spectra times 
				startvector.element(k) = (double)((1.0 / A_tilde->Ncols()));

			// get start vector and initialize it with equal fractions:
			// using the neighbor pixelvector as startvector

			// solve with iterative solution: 
			while (fabs(change) > LEAST_SQUARES_TOLERANCE) //0.0001
			{
				// go a small step into direction of negative gradient
				A_times_startvector = ((*A_tilde) * startvector);
				errorvector =  *(A_times_startvector.Evaluate()) - b_gamma;
				A_tilde_trans_mu =  mu * *(A_tilde_trans.Evaluate());
				temp = *(A_tilde_trans_mu.Evaluate()) * *(errorvector.Evaluate());

				{
					LOG_DISABLE
						//Logger::getSingletonLoggerStream().width(5);
						//Logger::getSingletonLoggerStream()<<iterations<<":   startV = "<<startvector<<endl;
						//Logger::getSingletonLoggerStream()<<"     :   A*stv = "<<A_times_startvector<<endl;
						//Logger::getSingletonLoggerStream()<<"     :   errorV = "<<errorvector<<endl;
						//Logger::getSingletonLoggerStream()<<"     :   temp = "<<temp<<endl;
						LOG_END
				}
				startvector = startvector - *(temp.Evaluate()); /// update startvector 

				//if one element gets negative, set it to zero
				for (k = 0; k < (A_tilde->Ncols()); k++)  // no. of spectra times
				{
					if (startvector.element(k) < 0.0)
						startvector.element(k) = 0.0;
				}

				// Check the deviation
				MatrixCol errorCol(errorvector.Evaluate(),LoadOnEntry, 0);
				double errorVectorNorm = errorCol.Norm2();
				change = deviation - errorVectorNorm;
				deviation = errorVectorNorm;

				//LOG_START
				//Logger::getSingletonLoggerStream()<<"     :   change = "<<change<<endl;
				//Logger::getSingletonLoggerStream()<<"     :   deviation = "<<deviation<<endl;
				//Logger::LogStringLine("-----------------------------------------");
				//LOG_END

				iterations++;
			} // while

			// length: no. of spectra
			double error = deviation / MatrixCol(&b_gamma,LoadOnEntry).Norm2();
			fraction = &startvector;

			{
				Logger::decimal(iterations,5);Logger::sstring(":");
				Logger::sstring("     :   fraction = ");//Logger::openStream() << *fraction;Logger::closeStream();
				Logger::endl();
				Logger::stringLine("=====================================================");
			}

			///  end of second constraint 
			// store fractions in resulting rows of resulting files
			// (number of bands = vector dimension)

			/// write result percent 
			for (i = 0; i < A.Ncols(); i++)  // no. of spectra 
				abundanceFractions[i + col*nEndmembers+row*width*nEndmembers] = fraction->element(i);

			// save error and iterations
			//error_cell[col] = (CELL) (100 * error);
			//iter_cell[col] = iterations;

		} // columns loop 

		/// write the resulting rows into output files:
		//for (i = 0; i < A.ncols; i++)   // no. of spectra 
		//  G_put_map_row(resultfd[i], result_cell[i]);
		//if (error_fd > 0)
		// G_put_map_row(error_fd, error_cell);
		//if (iter_fd > 0)
		// G_put_map_row(iter_fd, iter_cell);

	} // rows loop 

	Logger::stringLine("doLeastSquares finished successfully !");

	currentX = 0;
	currentY = 0;
	//delete A_tilde_trans_mu;
	delete A_tilde;

	return true;
}




