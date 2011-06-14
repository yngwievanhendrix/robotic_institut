/**

 \file wavelets.h
 (c) 2006 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        A trous wavelet decomposition and reconstruction<br>
 COMPONENT:      wavelets.cpp, wavelets.h
 \brief Wavelet decomposition
 
 **************************************************************************/

/***************************************************************************
 *
 *
 *  Provided classes:
 *
 *		CWaveletSkeleton
 *
 *  Provided structures:
 *
 *		none
 *
 *	Defined constants:
 *
 *		none
 *
 *  Defined types:
 *
 *		none
 *
 *	Provided methods and operators:
 *
 *		none
 *
 ****************************************************************************
 *
 *  These files require: 
 *		std::vector
 *		std::list
 *		std::string
 *		math.h
 *
 ***************************************************************************/

/****************************************************************************
 *
 * Required includes 
 *
 ***************************************************************************/

#ifndef __WAVELETS
#define __WAVELETS

#include "math/include/vectorops.h"

namespace Math {

	/**
	A trous wavelet decompositon and reconstruction class
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
	*/
	class CWaveletSkeleton {
	private:
		//! The input data: time
		DoubleVector m_x;
		//! The input data: signal
		DoubleVector m_y;
		//! The detail coefficients
		DoubleMatrix m_W;
		//! The approximation coefficients
		DoubleVector m_c;
		//! The approximation coefficients (all of them)
		DoubleMatrix m_AllCoeffsC;
		//! The number of scales to use
		int m_nNumScales;
		//! Is the class initialised
		bool m_fIsInitialised;
	public:
		/**
		Initialise variables.
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		CWaveletSkeleton(void);
		/**
		Initialise variables and compute the a trous wavelet decomposition of x at times y, at NumScales different scales.
		@param x	time
		@param y	signal
		@param NumScales	maximal number of scales
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		CWaveletSkeleton(DoubleVector x, DoubleVector y, int NumScales);
		/**
		Initialise variables and compute the a trous wavelet decomposition of x at times y, at NumScales different scales.
		@param x	time
		@param y	signal
		@param NumScales	maximal number of scales
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		bool Init(DoubleVector x, DoubleVector y, int NumScales);
		/**
		Get the wavelet reconstruction and trend
		@param FirstScale	the first scale to use in reconstruction
		@param LastScale	the last scale to use in reconstruction
		@param[out] appx	the reconstructed signal 
		@param[out] trend	the trend part of the signal
		@param[out] Delay	the shift induced by filtering
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		bool GetAppx(int FirstScale, int LastScale, DoubleVector &appx, DoubleVector &trend, int &Delay);
		/**
		Get the last available value of the wavelet reconstruction and trend 
		@param FirstScale	the first scale to use in reconstruction
		@param LastScale	the last scale to use in reconstruction
		@param[out] appx	the reconstructed signal
		@param[out] trend	the trend part of the signal
		@param[out] Delay	the shift induced by filtering
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		bool GetAppx(int FirstScale, int LastScale, double &appx, double &trend, int &Delay);
		/**
		Compute and get a vector of relative scale energies
		@param energies	the relative energies (return value)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
		*/
		void GetEnergies(DoubleVector &energies);
		/**
		Compute and get a vector of relative scale energies from start to end
		@param energies	the relative energies (return value)
		@param start	first index to consider
		@param end		last index to consider
		@returns boolean value indicating success
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-22
		*/
		bool GetEnergies(DoubleVector &energies, int start, int end);
		void GetW(DoubleMatrix &W);
		void GetC(DoubleVector &c);
		/** 
		Push a new measurement and compute the decomposition
		@param x	time stamp of the new measurement
		@param y	position of the new measurement
		@returns boolean value indicating success
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-04-10
		*/
		bool Push(double x, double y);
		/**
		Get the length of the data
		@returns	the length of the data stored in the class
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-04-10
		*/
		int GetSize(void);
	};
}

#endif
