/**

 \file vectorops.h
 (c) 2006 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        Vector operations<br>
 COMPONENT:      vectorops.cpp, vectorops.h
 \brief Vector and data manipulation methods
 
 **************************************************************************/

/***************************************************************************
 *
 *
 *  Provided classes:
 *
 *		none
 *
 *  Provided structures:
 *
 *		none
 *
 *	Defined constants:
 *
 *		GEOM_EPS, MIN, MAX, SGN
 *
 *  Defined types:
 *
 *		DoubleVector, LongVector, BoolVector
 *		DoubleList, DoubleMatrix, LongMatrix
 *		Cube
 *
 *	Provided methods and operators:
 *
 *		- addition, subtraction and scalar multiplication of vectors
 *		- dot product
 *		- computation of maxima, minima, mean, sum and variance of vectors or 
 *		  parts of vectors
 *		- computation of index of the minimum and maximum of a vector
 *		- rounding
 *		- RMS error of two vectors
 *		- clipping of vectors against a cube
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

#ifndef __VECTOROPS
#define __VECTOROPS

#define NOMINMAX
#include <algorithm>

#define MIN std::min
#define MAX std::max

#include <vector>
#include <list>
#include <string>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#ifndef M_PI
	//! define PI
    #define M_PI 3.1415926535897932384626433832795
#endif

#ifndef SGN
	//! geometric accuracy
	#define GEOM_EPS 0.000001
	//! sign function
	#define SGN(X) (((X) > GEOM_EPS) ? 1 : (((X) < -GEOM_EPS) ? -1 : 0))
#endif

namespace Math {
	//! vector of doubles
	typedef std::vector<double> DoubleVector;
	//! vector of integers
	typedef std::vector<int> LongVector;
	//! vector of bools
	typedef std::vector<bool> BoolVector;
	//! matrix of doubles
	typedef std::vector<DoubleVector> DoubleMatrix;
	//! vector of double vectors
	typedef std::vector<DoubleVector> DoubleVectorVector;
	//! matrix of integers
	typedef std::vector<LongVector> LongMatrix;
	//! list of doubles
	typedef std::list<double> DoubleList;

	//! Generate a subview of a vector (MATLAB colon operator)
	DoubleVector subview(const DoubleVector& v, int i, int j);
	//! Compute the standard deviation of a vector
	double stddev(const DoubleVector& v, size_t start = 0, size_t end = -1);
	//! compute the norm of a vector
	double norm(const DoubleVector& v, int norm = 2, bool doRoot = true);
	//! compute the signs of a vector
	DoubleVector sign(const DoubleVector&v);
	//! compute element-wise power of a vector
	DoubleVector Pow(const DoubleVector &v, double e);
	//! compute absolute values of a vector
	DoubleVector Abs(const DoubleVector &v);

	//! addition of two double vectors
	DoubleVector operator + (const DoubleVector& x, const DoubleVector& y);
	//! addition of a double vector and a scalar
	DoubleVector operator + (const DoubleVector& x, const double& s);
	//! addition of a double vector and a scalar
	DoubleVector operator + (const double& s, const DoubleVector& x);
	//! addition of a double vector and a scalar
	DoubleVector operator - (const DoubleVector& x, const double& s);
	//! subtraction of two double vectors
	DoubleVector operator - (const DoubleVector& x, const DoubleVector& y);
	//! scalar multiplication of a double vector and a double
	DoubleVector operator * (const DoubleVector& x, const double& r);
	//! scalar multiplication of a double and a double vector
	DoubleVector operator * (const double& r, const DoubleVector& x);
	//! dot product of two double vectors
	double operator * (const DoubleVector& x, const DoubleVector& y);
	//! element-wise product of two double vectors
	DoubleVector ElementProduct(const DoubleVector& x, const DoubleVector& y);

	//! addition of two int vectors
	LongVector operator + (const LongVector& x, const LongVector& y);
	//! subtraction of two int vectors
	LongVector operator - (const LongVector& x, const LongVector& y);
	//! scalar multiplication of a int vector and a int
	LongVector operator * (const LongVector& x, const int& r);
	//! scalar multiplication of a int and a int vector
	LongVector operator * (const int& r, const LongVector& x);
	//! dot product of two int vectors
	int operator * (const LongVector& x, const LongVector& y);

	//! cross product of two 3D vectors
	DoubleVector cross(const DoubleVector x, const DoubleVector y);

	/**
		A cube in 3D space, represented by two points.
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-06-21
	*/
	class CCube {
	private:
		DoubleVector m_p, m_q;
	public:
		CCube(void);
		CCube(const DoubleVector &x, const DoubleVector &y);
		void Init(const DoubleVector &x, const DoubleVector &y);
		void Clip(DoubleVector &x);
		void Clip(double &x, double &y, double &z);
		void Clip(double* x);
	};

	/**
		Round a double to a certain accuracy

		@param x	the value to round
		@param acc	the rounding accuracy (default: 1)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Round(double x, int acc = 1);

	/**
		Find the maximum of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we look (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Max(DoubleVector &x, int start = 0);
	/**
		Find the index of the maximum of a double vector

		@param x	the vector
		@param absolute	take absolute values first
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-21
	*/
	int MaxInd(DoubleVector &x, bool absolute = false, double *max = NULL);
	/**
		Find the minimum of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we look (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Min(DoubleVector &x, int start = 0);
	/**
		Find the index of the minimum of a double vector

		@param x	the vector
		@param absolute	take absolute values first
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-29
	*/
	int MinInd(DoubleVector &x, bool absolute = false);
	/**
		Find the maximum of a double vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we look
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Max(DoubleVector &x, int start, int end);
	/**
		Find the minimum of a double vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we look
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Min(DoubleVector &x, int start, int end);
	/**
		Find the maximum of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we look (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Max(LongVector &x, int start = 0);
	/**
		Find the index of the maximum of an int vector

		@param x	the vector
		@param absolute	take absolute values first
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-21
	*/
	int MaxInd(LongVector &x, bool absolute = false);
	/**
		Find the minimum of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we look (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Min(LongVector &x, int start = 0);
	/**
		Find the index of the minimum of an int vector

		@param x	the vector
		@param absolute	take absolute values first
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-29
	*/
	int MinInd(LongVector &x, bool absolute = false);
	/**
		Find the maximum of an int vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we look
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Max(LongVector &x, int start, int end);
	/**
		Find the minimum of an int vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we look
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Min(LongVector &x, int start, int end);

	/**
		Compute the sum of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we sum (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Sum(DoubleVector &x, int start = 0);
	/**
		Compute the sum of a double vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we sum
		@param end	up to which index do we sum
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Sum(DoubleVector &x, int start, int end);
	/**
		Compute the sum of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from which index on do we sum (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Sum(LongVector &x, int start = 0);
	/**
		Compute the sum of an int vector, starting at <i>start</i> and ending at <i>end</i>

		@param x	the vector
		@param start	from which index on do we sum
		@param end	up to which index do we sum
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	int Sum(LongVector &x, int start, int end);

	/**
		Compute the mean of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Mean(const DoubleVector &x, int start = 0);
	/**
		Compute the mean of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Mean(const DoubleVector &x, int start, int end);
	/**
		Compute the mean of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Mean(const LongVector &x, int start = 0);
	/**
		Compute the mean of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Mean(const LongVector &x, int start, int end);

	/**
		Compute the variance of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Var(LongVector &x, int start = 0);
	/**
		Compute the variance of a double vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Var(LongVector &x, int start, int end);
	/**
		Compute the variance of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start (default: 0)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Var(DoubleVector &x, int start = 0);
	/**
		Compute the variance of an int vector, starting at <i>start</i>

		@param x	the vector
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-17
	*/
	double Var(DoubleVector &x, int start, int end);
	/**
		Compute the RMS difference between x(start:end) and y(start:end)

		@param x	the first vector
		@param y	the second vector (a list)
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-21
	*/
	double RMSError(DoubleVector &x, DoubleList &y, int start, int end);
	/**
		Compute the RMS difference between x(start:end) and y(start:end)

		@param x	the first vector
		@param y	the second vector
		@param start	from where we start
		@param end	up to which index do we look
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2006-11-21
	*/
	double RMSError(DoubleVector &x, DoubleVector &y, int start, int end);
	/** 
		Compute the jitter of the input signal

		@param x	the input signal
		@param dt	the sampling size
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-11-24
	*/
	double Jitter(DoubleVector &x, double time);
	/** 
	Compute the jitter of the 3D input signal

	@param x	the input signal
	@param y	the input signal
	@param z	the input signal
	@param dt	the sampling size
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-11-24
	*/
	double Jitter(DoubleVector &x, DoubleVector &y, DoubleVector &z, double time);
	/**
		Find local maxima and minima of a vector.
		@param v	the vector to look at
		@param t	threshold (minimum difference between a maximum and a minimum)
		@param max	the indices of the maxima (return value)
		@param maxVals	the values corresponding to the indices in max (return value)
		@param min	the indices of the minima (return value)
		@param minVals	the values corresponding to the indices in min (return value)
		@returns	boolean value indicating success
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-01-10
	*/
	bool LMaxMin(DoubleVector &v, double t, LongVector &max, DoubleVector &maxVals, LongVector &min, DoubleVector &minVals);
}

#endif
