#ifndef __MATRIX_H
#define __MATRIX_H

#include <iostream>
#include <math.h>
#include <vector>

// #include "general/include/TransformMatrix.h"

#ifndef M_PI
	#define M_PI 3.1415926535897932384626433832795
#endif

#ifdef __USE_LAPACK
	// LaPack++
	#include "lapackpp/include/gmd.h"
	#include "lapackpp/include/lavd.h"
	#include "lapackpp/include/laindex.h"
	#include "lapackpp/include/laslv.h"
	#include "lapackpp/include/blas1pp.h"
	#include "lapackpp/include/blas2pp.h"
	#include "lapackpp/include/blas3pp.h"
#else
	// GSL
	#include "gsl/gsl_linalg.h"
	#include "gsl/gsl_eigen.h"
	#include "gsl/gsl_blas.h"
	#include "gsl/gsl_math.h"
	#include "gsl/gsl_statistics_double.h"
#endif

namespace Math {
	const double m_dEps = 1e-8;

	class Matrix {
	public:
	#ifndef __USE_LAPACK
		double *_values;
		int _size[2];
	#else
		LaGenMatDouble _matrix;
	#endif

		Matrix();
		Matrix(int m, int n = 0);
		Matrix(const Matrix& m);
	#ifdef _TRANSFORMMATRIX_H_
		inline Matrix(const TransformMatrix& m) {
			#ifndef __USE_LAPACK
				_values = NULL;
				_size[0] = _size[1] = 0;
			#endif
			resize(4,4);
			for(int i=0; i<3; i++)
				for(int j=0; j<4; j++)
					set(m(i+1,j+1), i, j);
			set(0.0,3,0);
			set(0.0,3,1);
			set(0.0,3,2);
			set(1,3,3);
		}
	#endif
	#ifndef __USE_LAPACK
		~Matrix(void);
	#endif

		inline double& operator()(int i, int j) { 
			#ifndef __USE_LAPACK
				return _values[_size[1]*i + j];
			#else
				return _matrix(i,j); 
			#endif
		}
		inline const double& operator()(int i, int j) const {
			#ifndef __USE_LAPACK
				return _values[_size[1]*i + j];
			#else
				return _matrix(i,j); 
			#endif
		}
		inline bool operator==(const Matrix& s) const {
			if(_size[0] != s._size[0] || _size[1] != s._size[1])
				return false;
			for(int i=0; i<_size[0]; i++) {
				for(int j=0; j<_size[1]; j++) {
					if(fabs((*this)(i,j) - s(i,j)) > m_dEps)
						return false;
				}
			}
			return true;
		}
		inline bool operator!=(const Matrix& s) const {
			return !((*this)==s);
		}
		inline Matrix& operator=(const Matrix& s) { 
			#ifndef __USE_LAPACK
				resize(s._size[0], s._size[1]);
				_size[0] = s._size[0]; _size[1] = s._size[1];
				for(int i=0; i<_size[0]*_size[1];i++)
					_values[i] = s._values[i];
				return (*this);
			#else
				_matrix = s._matrix; 
				return (*this); 
			#endif
		}

		Matrix copy(int fromRow, int fromCol, int numRows, int numCols);
		int size(int dim) const;
		void inject(const Matrix& m);
		void scale(double s);
		void resize(int m, int n = 0);
		void set(double *vals, int m, int n);
		inline void set(double val, int i, int j) {
			#ifndef __USE_LAPACK
				_values[_size[1]*i + j] = val;
			#else
				_matrix(i,j) = val; 
			#endif
		}

		inline void add(const Matrix& mat) 
		{
			#ifndef __USE_LAPACK
				for(int i=0; i<_size[0]*_size[1]; i++)
					_values[i] += mat._values[i];
			#else
				for(int i=0; i<_size[0]; i++)
					for(int j=0; j<_size[1]; j++)
						_matrix(i,j) += _mat(i,j); 
			#endif
		}
		inline void sub(const Matrix& mat) 
		{
			#ifndef __USE_LAPACK
				for(int i=0; i<_size[0]*_size[1]; i++)
					_values[i] -= mat._values[i];
			#else
				for(int i=0; i<_size[0]; i++)
					for(int j=0; j<_size[1]; j++)
						_matrix(i,j) -= _mat(i,j); 
			#endif
		}

		void leftMultiply(const Matrix& left, Matrix &result, bool transpose_A = false, bool transpose_B = false, double alpha = 1, double beta = 0) const;
		void rightMultiply(const Matrix& right, Matrix &result, bool transpose_A = false, bool transpose_B = false, double alpha = 1, double beta = 0) const;

		void linearSolve(Matrix &X, const Matrix& B, bool overwrite = false) const;
		void QRSolve(Matrix &X, const Matrix& B, bool overwrite = false) const;
#ifndef __USE_LAPACK
		void SVD(Matrix &U, Matrix &S, Matrix &V) const;
		bool orthonormalise(void);
#endif
		double trace() const;
		double norm2() const;

	#ifdef _TRANSFORMMATRIX_H_
		bool ToTransformMatrix(TransformMatrix &matrix);
	#endif

		bool GetEigenvalues(std::vector<double> &eigs_r, std::vector<double> &eigs_i);
		bool GetEigenvalues(std::vector<double> &eigs_r, std::vector<double> &eigs_i, std::vector<double> &eigs_norm);
		bool GetEigenvectors(Matrix &V, bool sort = false);

		static void MatrixMultiply(const Matrix& A, const Matrix& B, Matrix &C, bool transpose_A = false, bool transpose_B = false, double alpha = 1, double beta = 0);
		
		static void invert(const Matrix &X, Matrix &Xinv);
		static double Max(const std::vector<double> &v);
		static double Max(const double *v, int num);
		static double Norm(const std::vector<double> &v);
		static double Norm(const double *v, int num);
		
		bool TransformErrors(double &escale, double &erot, double &etrans);

		bool covariance(Matrix &r);
		bool PCA(Matrix &p, Matrix &v);

		inline static Matrix zeros (int m, int n=0) {
			Matrix mat(m,n);
			#ifndef __USE_LAPACK
				for(int i=0; i<n*m; i++)
					mat._values[i] = 0;
			#else
				mat._matrix = LaGenMatDouble::zeros(m,n);
			#endif
			return mat;
		}
		inline static Matrix ones (int m, int n=0) {
			Matrix mat(m,n);
			#ifndef __USE_LAPACK
				for(int i=0; i<n*m; i++)
					mat._values[i] = 1;
			#else
				mat._matrix = LaGenMatDouble::ones(m,n);
			#endif
			return mat;
		}
		inline static Matrix eye (int m) {
			Matrix mat(m,m);
			#ifndef __USE_LAPACK
				for(int i=0; i<m; i++) {
					for(int j=0; j<m; j++) {
						mat(i,j) = 0;
						if(i==j)
							mat(i,j) = 1;
					}
				}
			#else
				mat._matrix = LaGenMatDouble::eye(m,m);
			#endif
			return mat;
		}
	};

	inline std::ostream& operator<<(std::ostream& s, const Matrix& M) { 
	#ifndef __USE_LAPACK
		for(int i=0; i<M._size[0]; i++) {
			for(int j=0; j<M._size[1]; j++) {
				s << M(i,j);
				if(j < M._size[1]-1) {
					s << " ";
				}
			} 
			s << std::endl;
		}
	#else
		s << M._matrix; 
	#endif
		return s;
	}

	/**
	Generate a random matrix
	@param[in]	rotmax	maximal rotation around x, y, z
	@param[in]	posmax	maximal position in x, y, z
	@param[out]	M		resulting Matrix
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-01-17
	*/
	bool RandMatrix(Matrix &M, double rotmax = 5, double posmax = 10);
}

#endif
