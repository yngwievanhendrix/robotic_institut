#include "math/include/matrix.h"
#include <assert.h>

using namespace Math;


Math::Matrix::Matrix() { 
#ifndef __USE_LAPACK
	_values = NULL;
	_size[0] = _size[1] = 0;
#endif	
}

Math::Matrix::Matrix(int m, int n /* = 0 */) {
#ifndef __USE_LAPACK
	_values = NULL;
	_size[0] = _size[1] = 0;
#endif	
	resize(m,n);
}

Math::Matrix::Matrix(const Matrix& m) {
#ifndef __USE_LAPACK
	_values = NULL;
	_size[0] = _size[1] = 0;
#endif	
	inject(m);
}

#ifndef __USE_LAPACK
Math::Matrix::~Matrix(void) {
	if(_values != NULL)
		free(_values);
	_values = NULL;
	_size[0] = _size[1] = 0;
}
#endif

Math::Matrix Math::Matrix::copy(int fromRow, int fromCol, int numRows, int numCols)
{
	Matrix copy(numRows, numCols);
	for (int i=fromRow; i<fromRow+numRows; i++)
		for (int j=fromCol; j<fromCol+numCols; j++)
		{
			double d = (*this)(i,j);
			copy(i-fromRow, j-fromCol) = d;
		}
	return copy;
}


int Math::Matrix::size(int dim) const {
	if(dim < 0 || dim > 1)
		return -1;
#ifndef __USE_LAPACK
	return _size[dim];
#else
	return _matrix.size(dim);
#endif
}

void Math::Matrix::inject(const Matrix& m) {
	resize(m.size(0), m.size(1));

	for(int i=0; i<m.size(0); i++)
		for(int j=0; j<m.size(1); j++)
			set(m(i,j), i,j);
}

void Math::Matrix::scale(double s) {
#ifndef __USE_LAPACK
	for(int i=0; i<_size[0]*_size[1]; i++)
		_values[i] *= s;
#else
	_matrix.scale(s);
#endif
}

void Math::Matrix::resize(int m, int n /* = 0 */) {
#ifndef __USE_LAPACK
	if(_values != NULL)
		free(_values);

	n = n==0?m:n;
	_values = (double *)malloc(m*n*sizeof(double));
	_size[0] = m; _size[1] = n;
#else
	_matrix.resize(m, n==0?m:n);
#endif
}

void Math::Matrix::leftMultiply(const Matrix& left, Matrix &result, bool transpose_A, bool transpose_B, double alpha, double beta) const {
#ifndef __USE_LAPACK
	gsl_matrix_view _left = gsl_matrix_view_array(left._values, left._size[0], left._size[1]);
	gsl_matrix_view _right = gsl_matrix_view_array(_values, _size[0], _size[1]);
	gsl_matrix_view _output = gsl_matrix_view_array(result._values, result._size[0], result._size[1]);
	gsl_blas_dgemm(transpose_A?CblasTrans:CblasNoTrans, transpose_B?CblasTrans:CblasNoTrans, alpha, &_left.matrix, &_right.matrix, beta, &_output.matrix);
//	for(int i=0; i<result.size(0); i++)
//		for(int j=0; j<result.size(1); j++)
//			result._values[i*result.size(1) + j] = gsl_matrix_get(&_output.matrix, i, j);
#else
	Blas_Mat_Mat_Mult(left._matrix, this->_matrix, result._matrix, transpose_A, transpose_B, alpha, beta);
#endif
}

void Math::Matrix::rightMultiply(const Matrix& right, Matrix &result, bool transpose_A, bool transpose_B, double alpha, double beta) const {
#ifndef __USE_LAPACK
	gsl_matrix_view _left = gsl_matrix_view_array(_values, _size[0], _size[1]);
	gsl_matrix_view _right = gsl_matrix_view_array(right._values, right._size[0], right._size[1]);
	gsl_matrix_view _output = gsl_matrix_view_array(result._values, result._size[0], result._size[1]);
	gsl_blas_dgemm(transpose_A?CblasTrans:CblasNoTrans, transpose_B?CblasTrans:CblasNoTrans, alpha, &_left.matrix, &_right.matrix, beta, &_output.matrix);
//	for(int i=0; i<result.size(0); i++)
//		for(int j=0; j<result.size(1); j++)
//			result._values[i*result.size(1) + j] = gsl_matrix_get(&_output.matrix, i, j);
#else
	Blas_Mat_Mat_Mult(this->_matrix, right._matrix, result._matrix, transpose_A, transpose_B, alpha, beta);
#endif
}

void Math::Matrix::linearSolve(Matrix &X, const Matrix& B, bool overwrite) const {
	X.resize(_size[1], B.size(1));
#ifndef __USE_LAPACK
	double *v = NULL;
	gsl_matrix_view _A;
	if(!overwrite) {
		v = (double *)malloc(_size[0]*_size[1]*sizeof(double));
		for(int i=0; i<_size[0]*_size[1]; i++)
			v[i] = _values[i];
		_A = gsl_matrix_view_array(v, _size[0], _size[1]);
	} else {
		_A = gsl_matrix_view_array(_values, _size[0], _size[1]);
	}

	gsl_matrix_view x;
	x  = gsl_matrix_view_array(X._values, X._size[0], X._size[1]);
	gsl_matrix_view b;
	b  = gsl_matrix_view_array(B._values, B._size[0], B._size[1]);

	gsl_permutation *p = gsl_permutation_alloc(_size[0]);
	gsl_vector* tau = gsl_vector_alloc(_size[1]);
	gsl_vector* residual = gsl_vector_alloc(B._size[0]);
	gsl_linalg_QR_decomp (&_A.matrix, tau);
	
	for(int i=0; i<B.size(1); i++)
		gsl_linalg_QR_lssolve(&_A.matrix, tau, &gsl_matrix_column(&b.matrix,i).vector, &gsl_matrix_column(&x.matrix,i).vector, residual);

	gsl_permutation_free(p);
	gsl_vector_free(residual);
	gsl_vector_free(tau);
	if(v != NULL)
		free(v);
#else
	LaLinearSolve(_matrix, X._matrix, B._matrix);
#endif
}

void Math::Matrix::QRSolve(Matrix &X, const Matrix& B, bool overwrite) const {
#ifndef __USE_LAPACK
	double *v = NULL;
	gsl_matrix_view _A;
	if(!overwrite) {
		v = (double *)malloc(_size[0]*_size[1]*sizeof(double));
		for(int i=0; i<_size[0]*_size[1]; i++)
			v[i] = _values[i];
		_A = gsl_matrix_view_array(v, _size[0], _size[1]);
	} else {
		_A = gsl_matrix_view_array(_values, _size[0], _size[1]);
	}
	
	gsl_matrix_view x = gsl_matrix_view_array(X._values, X._size[0], X._size[1]);
	gsl_matrix_view b = gsl_matrix_view_array(B._values, B._size[0], B._size[1]);

	gsl_vector *tau = gsl_vector_alloc((_size[0]<=_size[1])?_size[0]:_size[1]);
	gsl_vector *res = gsl_vector_alloc(B.size(0));

	gsl_linalg_QR_decomp(&_A.matrix, tau);
	for(int i=0; i<B.size(1); i++)
		gsl_linalg_QR_lssolve(&_A.matrix, tau, &gsl_matrix_column(&b.matrix,i).vector, &gsl_matrix_column(&x.matrix,i).vector, res);

	gsl_vector_free(tau);
	gsl_vector_free(res);
	if(v != NULL)
		free(v);
#else
	LaQRLinearSolve(_matrix, X._matrix, B._matrix);
#endif
}

/*
@author Volker Martens
@param U has same size as this matrix (m x n)
@param S is square with size (n x n) and singular values in main diagonale
@param V is square with size (n x n)
*/
#ifndef __USE_LAPACK
void Math::Matrix::SVD(Matrix &U, Matrix &S, Matrix &V) const 
{
	int m = _size[0];
	int n = _size[1];
	gsl_vector *workspace = gsl_vector_alloc(n);
	gsl_vector *s = gsl_vector_alloc(n);

	U = *this;
	gsl_matrix_view u = gsl_matrix_view_array(U._values, m, n);
	gsl_matrix_view v = gsl_matrix_view_array(V._values, n, n);
	gsl_linalg_SV_decomp(&u.matrix, &v.matrix, s, workspace);
	for (int i=0; i<n; i++)
		S(i,i)=s->data[i];

	gsl_vector_free(workspace);
	gsl_vector_free(s);
}
#endif

double Math::Matrix::norm2() const 
{
	int m = _size[0];
	int n = _size[1];
	gsl_vector *workspace = gsl_vector_alloc(n);
	gsl_vector *s = gsl_vector_alloc(n);

	Matrix U = *this;
	Matrix V(n,n);
	gsl_matrix_view u = gsl_matrix_view_array(U._values, m, n);
	gsl_matrix_view v = gsl_matrix_view_array(V._values, n, n);
	gsl_linalg_SV_decomp(&u.matrix, &v.matrix, s, workspace);

	double maxSingVal = -10000000;
	for (int i=0; i<n; i++)
		if (s->data[i] > maxSingVal)
			maxSingVal = s->data[i];
	gsl_vector_free(workspace);
	gsl_vector_free(s);
	return maxSingVal;
}

double Math::Matrix::trace() const 
{
	int dim = GSL_MIN(_size[0], _size[1]);
	double trace = 0;
	for (int i=0; i<dim; i++)
		trace += (*this)(i,i);
	return trace;
}


void Math::Matrix::set(double *vals, int m, int n) {
	resize(m,n);
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			set(vals[n*i+j],i,j);
}

void Math::Matrix::MatrixMultiply(const Matrix& A, const Matrix& B, Matrix &C, bool transpose_A, bool transpose_B, double alpha, double beta) {
#ifndef __USE_LAPACK
	gsl_matrix_view _left = gsl_matrix_view_array(A._values, A._size[0], A._size[1]);
	gsl_matrix_view _right = gsl_matrix_view_array(B._values, B._size[0], B._size[1]);
	gsl_matrix_view _output = gsl_matrix_view_array(C._values, C._size[0], C._size[1]);
	gsl_blas_dgemm(transpose_A?CblasTrans:CblasNoTrans, transpose_B?CblasTrans:CblasNoTrans, alpha, &_left.matrix, &_right.matrix, beta, &_output.matrix);
//	for(int i=0; i<C.size(0); i++)
//		for(int j=0; j<C.size(1); j++)
//			C._values[i*C.size(1) + j] = gsl_matrix_get(&_output.matrix, i, j);
#else
	Blas_Mat_Mat_Mult(A._matrix, B._matrix, C._matrix, transpose_A, transpose_B, alpha, beta);
#endif
}

#ifdef _TRANSFORMMATRIX_H_
bool Math::Matrix::ToTransformMatrix(TransformMatrix &matrix) {
	if (size(0) != 4 || size(1) != 4)
		return false;

	for(int i=0; i<3; i++) {
		for(int j=0; j<4; j++) {
			matrix(i+1,j+1) = (*this)(i,j);
		}
	}
	return true;
}
#endif

bool Math::Matrix::GetEigenvalues(std::vector<double> &eigs_r, std::vector<double> &eigs_i, std::vector<double> &eigs_norm) {
	if(!GetEigenvalues(eigs_r, eigs_i))
		return false;

	eigs_norm.resize(eigs_i.size());
	for(size_t i=0; i<eigs_norm.size(); i++)
		eigs_norm[i] = sqrt(eigs_i[i] * eigs_i[i] + eigs_r[i] * eigs_r[i]);

	return true;
}

bool Math::Matrix::GetEigenvalues(std::vector<double> &eigs_r, std::vector<double> &eigs_i) {
#ifndef __USE_LAPACK
	if(_size[0] != _size[1])
		return false;
	
	gsl_matrix_view tmp = gsl_matrix_view_array(_values, _size[0], _size[1]);
	gsl_matrix *m = gsl_matrix_alloc(_size[0], _size[1]);
	gsl_matrix_memcpy(m, &tmp.matrix);

	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(_size[0]);
	gsl_vector_complex *eval = gsl_vector_complex_alloc(_size[0]);
	if(gsl_eigen_nonsymm(m, eval, w) != 0)
		return false;

	eigs_i.resize(_size[0]);
	eigs_r.resize(_size[0]);
	for(int i=0; i<_size[0]; i++) {
		eigs_r[i] = gsl_vector_complex_get(eval, i).dat[0];
		eigs_i[i] = gsl_vector_complex_get(eval, i).dat[1];
	}
	gsl_eigen_nonsymm_free(w);
	gsl_vector_complex_free(eval);

	gsl_matrix_free(m);

	return true;
#else
	return false;
#endif
}

double Math::Matrix::Max(const std::vector<double> &v) {
	double max = v[0];
	for(size_t i=1; i<v.size(); i++)
		if(v[i] > max)
			max = v[i];

	return max;
}

double Math::Matrix::Max(const double *v, int num) {
	double max = v[0];
	for(int i=1; i<num; i++)
		if(v[i] > max)
			max = v[i];

	return max;
}

double Math::Matrix::Norm(const std::vector<double> &v) {
	double norm = 0;
	for(size_t i=0; i<v.size(); i++)
		norm += v[i] * v[i];

	return sqrt(norm);
}

double Math::Matrix::Norm(const double *v, int num) {
	double norm = 0;
	for(int i=0; i<num; i++)
		norm += v[i] * v[i];

	return sqrt(norm);
}


// @author Volker Martens
void Math::Matrix::invert(const Matrix &X, Matrix &Xinv)
{
	int dim = X.size(0);
	assert(dim == X.size(1));
	assert(dim == Xinv.size(0));
	assert(dim == Xinv.size(1));

	int s;
    gsl_matrix* ludecomp = gsl_matrix_calloc(dim,dim);
	gsl_permutation* perm = gsl_permutation_alloc(dim);
	gsl_matrix_view x = gsl_matrix_view_array(X._values, X._size[0], X._size[1]);
	gsl_matrix_view xinv = gsl_matrix_view_array(Xinv._values, Xinv._size[0], Xinv._size[1]);
	gsl_matrix_memcpy(ludecomp, &x.matrix);
    gsl_linalg_LU_decomp(ludecomp, perm, &s);
	gsl_linalg_LU_invert(ludecomp, perm, &xinv.matrix);

	gsl_permutation_free(perm);
}


bool Math::Matrix::TransformErrors(double &escale, double &erot, double &etrans) {
	if(size(0) != 4 || size(1) != 4)
		return false;

	std::vector<double> er(4,0), ei(4,0), en(4,0);
	std::vector<double> trans(3,0);

	if(!GetEigenvalues(er, ei, en))
		return false;

	// scale error
	size_t num = en.size();
	for(size_t i=0; i<num; i++)
		en.push_back(1.0 / en[i]);

	escale = fabs(Matrix::Max(en) - 1);

	// translation error
	trans[0] = (*this)(0,3);
	trans[1] = (*this)(1,3);
	trans[2] = (*this)(2,3);
	etrans = Matrix::Norm(trans);

	// rotation error
	if(Matrix::Norm(ei) < 1e-10)
		erot = 0;
	else {
		// find complex eigenvalues
		double e[2][2];
		int k=0;
		for(size_t i=0; i<ei.size() && k < 2; i++) {
			if(ei[i] > 1e-10)
				e[0][0] = er[i], e[0][1] = ei[i], k++;
			else if(ei[i] < -1e-10)
				e[1][0] = er[i], e[1][1] = ei[i], k++;
		}
		// determine rotation angle
		erot = fabs(atan2(e[0][1] - e[1][1],e[0][0] + e[1][0])) / M_PI * 180.0;
	}	

	return true;
}

bool Math::Matrix::covariance( Matrix &r )
{
	#ifdef __USE_LAPACK
		return false;
	#else
		gsl_vector_view a, b;
		gsl_matrix_view m = gsl_matrix_view_array(_values, _size[0], _size[1]);
		size_t i, j;
		r.resize(m.matrix.size2);

		for (i = 0; i < m.matrix.size2; i++) {
			for (j = 0; j < m.matrix.size2; j++) {
				a = gsl_matrix_column (&(m.matrix), i);
				b = gsl_matrix_column (&(m.matrix), j);
				r(i,j) = gsl_stats_covariance (a.vector.data, a.vector.stride,
					b.vector.data, b.vector.stride, a.vector.size);
			}
		}
		return true;
	#endif
}

bool Math::Matrix::GetEigenvectors( Matrix &V, bool sort )
{
	#ifdef __USE_LAPACK
		return false;
	#else
		if(_size[0] != _size[1])
			return false;
		V.resize(_size[0]);

		gsl_matrix_view tmp = gsl_matrix_view_array(_values, _size[0], _size[1]);
		gsl_matrix *m = gsl_matrix_alloc(_size[0], _size[1]);
		gsl_matrix_memcpy(m, &tmp.matrix);

		gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(_size[0]);
		gsl_vector_complex *eval = gsl_vector_complex_alloc(_size[0]);
		gsl_matrix_complex *evec = gsl_matrix_complex_alloc(_size[0],_size[0]);
		if(gsl_eigen_nonsymmv(m, eval, evec, w) != 0)
			return false;
		if(sort) {
			if(gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC) != 0)
				return false;
		} 

		for(int i=0; i<_size[0]; i++) {
			for(int j=0; j<_size[0]; j++) {
				V(i,j) = gsl_matrix_complex_get(evec, i, j).dat[0];
			}
		}

		gsl_eigen_nonsymmv_free(w);
		gsl_vector_complex_free(eval);
		gsl_matrix_complex_free(evec);

		gsl_matrix_free(m);

		return true;
	#endif
}

bool Math::Matrix::PCA( Matrix &p, Matrix &v )
{
	// get covariance matrix
	Matrix cov;
	if(!this->covariance(cov))
		return false;

	// get eigenvector matrix
	if(!cov.GetEigenvectors(v, true))
		return false;

	// transform point cloud
	p.resize(_size[0], _size[1]);
	this->rightMultiply(v, p);

	return true;
}

#ifndef __USE_LAPACK
bool Math::Matrix::orthonormalise( void )
{
	int m = _size[0];
	int n = _size[1];
	if(m != n || m != 4)
		return false;

	gsl_vector *workspace = gsl_vector_alloc(3);
	gsl_vector *s = gsl_vector_alloc(3);

	gsl_matrix *r = gsl_matrix_alloc(3,3);
	gsl_matrix *u = gsl_matrix_alloc(3,3);
	gsl_matrix *v = gsl_matrix_alloc(3,3);
	for(m=0; m<3; m++) {
		for(n=0; n<3; n++) {
			gsl_matrix_set(u, m, n, (*this)(m,n));
		}
	}
	gsl_linalg_SV_decomp(u, v, s, workspace);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, u, v, 0, r);
	for(m=0; m<3; m++) {
		for(n=0; n<3; n++) {
			(*this)(m,n) = gsl_matrix_get(r, m, n);
		}
	}

	gsl_vector_free(workspace);
	gsl_vector_free(s);
	gsl_matrix_free(v);
	gsl_matrix_free(u);
	gsl_matrix_free(r);

	return true;
}
#endif

bool Math::RandMatrix(Matrix &M, double rotmax, double posmax) {
	//! Degrees to radians
	double CAL_DEG2RAD = M_PI/180.0;

	//! Degrees to radians
	double CAL_RAD2DEG = 180.0/M_PI;

	// prepare output matrix
	M.resize(4,4);

	// rotation matrices
	Matrix rotX = Matrix::eye(4);
	Matrix rotY = Matrix::eye(4);
	Matrix rotZ = Matrix::eye(4);
	Matrix temp = Matrix::eye(4);

	// generate random rotation
	double rot[3] = { 
		CAL_DEG2RAD*2*rotmax*((double)rand()/(double)RAND_MAX - 0.5), 
		CAL_DEG2RAD*2*rotmax*((double)rand()/(double)RAND_MAX - 0.5),
		CAL_DEG2RAD*2*rotmax*((double)rand()/(double)RAND_MAX - 0.5) 
	};

	rotX(1,1) = cos(rot[0]);
	rotX(1,2) = sin(rot[0]);
	rotX(2,1) = -sin(rot[0]);
	rotX(2,2) = cos(rot[0]);

	rotY(0,0) = cos(rot[1]);
	rotY(0,2) = sin(rot[1]);
	rotY(2,0) = -sin(rot[1]);
	rotY(2,2) = cos(rot[1]);

	rotZ(0,0) = cos(rot[2]);
	rotZ(0,1) = sin(rot[2]);
	rotZ(1,0) = -sin(rot[2]);
	rotZ(1,1) = cos(rot[2]);

	rotX.rightMultiply(rotY, temp, false, false, 1);
	temp.rightMultiply(rotZ, M, false, false, 1);

	// generate random position offset
	double pos[3] = {
		2*posmax*((double)rand()/(double)RAND_MAX - 0.5),
		CAL_DEG2RAD*180*((double)rand()/(double)RAND_MAX),
		CAL_DEG2RAD*360*((double)rand()/(double)RAND_MAX) 
	};

	M(0,3) = pos[0]*sin(pos[1])*cos(pos[2]);
	M(1,3) = pos[0]*sin(pos[1])*sin(pos[2]);
	M(2,3) = pos[0]*cos(pos[1]);

	return true;
}

