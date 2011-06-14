#include "math/include/Horn.h"

using namespace std;
using namespace Math;

Math::Horn::Horn(void) {
	Alloc();
	Init(3, NULL, NULL);
}

Math::Horn::Horn(int n, gsl_matrix* dataPointsLeft, gsl_matrix* dataPointsRight) {
	Alloc();
	Init(n, dataPointsLeft, dataPointsRight);
}

Math::Horn::Horn(std::vector<Math::tDataObject> &left, std::vector<Math::tDataObject> &right) {
	if(left.size() != right.size())
		return;
	
	// allocate memory
	Alloc();
	int n = left.size();
	_dataPointsLeft = gsl_matrix_calloc(3, n);
	_dataPointsRight = gsl_matrix_calloc(3, n);
	_needToFreePoints = true;
	for(int i=0; i<n; i++) {
		for(int j=0; j<3; j++) {
			if(left[i].RotationType == Math::rotMATRIX) {
				gsl_matrix_set(_dataPointsLeft, j, i, left[i].matrix(j,3));
				gsl_matrix_set(_dataPointsRight, j, i, right[i].matrix(j,3));
			} else {
				gsl_matrix_set(_dataPointsLeft, j, i, left[i].trans[j]);
				gsl_matrix_set(_dataPointsRight, j, i, right[i].trans[j]);
			}
		}
	}
	Init(n, _dataPointsLeft, _dataPointsRight);
}

Math::Horn::Horn(std::vector<Math::Matrix> &left, std::vector<Math::Matrix> &right) {
	if(left.size() != right.size())
		return;

	// allocate memory
	Alloc();
	int n = left.size();
	_dataPointsLeft = gsl_matrix_calloc(3, n);
	_dataPointsRight = gsl_matrix_calloc(3, n);
	_needToFreePoints = true;
	for(int i=0; i<n; i++) {
		for(int j=0; j<3; j++) {
			if (left[i].size(0) != 4 || left[i].size(1) != 4 || right[i].size(0) != 4 || right[i].size(1) != 4) {
				gsl_matrix_free(_dataPointsLeft); _dataPointsLeft = NULL;
				gsl_matrix_free(_dataPointsRight); _dataPointsRight = NULL;
				_needToFreePoints = false;
				return;
			}
			
			gsl_matrix_set(_dataPointsLeft, j, i, left[i](j,3));
			gsl_matrix_set(_dataPointsRight, j, i, right[i](j,3));
		}
	}
	Init(n, _dataPointsLeft, _dataPointsRight);
}
Math::Horn::Horn(std::vector<Math::DoubleVector> &left, std::vector<Math::Matrix> &right) {
	if(left.size() != right.size())
		return;

	// allocate memory
	Alloc();
	int n = left.size();
	_dataPointsLeft = gsl_matrix_calloc(3, n);
	_dataPointsRight = gsl_matrix_calloc(3, n);
	_needToFreePoints = true;
	for(int i=0; i<n; i++) {
		for(int j=0; j<3; j++) {
			if (left[i].size() != 3 || right[i].size(0) != 4 || right[i].size(1) != 4) {
				gsl_matrix_free(_dataPointsLeft); _dataPointsLeft = NULL;
				gsl_matrix_free(_dataPointsRight); _dataPointsRight = NULL;
				_needToFreePoints = false;
				return;
			}

			gsl_matrix_set(_dataPointsLeft, j, i, left[i][j]);
			gsl_matrix_set(_dataPointsRight, j, i, right[i](j,3));
		}
	}
	Init(n, _dataPointsLeft, _dataPointsRight);
}
Math::Horn::Horn(std::vector<Math::Matrix> &left, std::vector<Math::DoubleVector> &right) {
	if(left.size() != right.size())
		return;

	// allocate memory
	Alloc();
	int n = left.size();
	_dataPointsLeft = gsl_matrix_calloc(3, n);
	_dataPointsRight = gsl_matrix_calloc(3, n);
	_needToFreePoints = true;
	for(int i=0; i<n; i++) {
		for(int j=0; j<3; j++) {
			if (right[i].size() != 3 || left[i].size(0) != 4 || left[i].size(1) != 4) {
				gsl_matrix_free(_dataPointsLeft); _dataPointsLeft = NULL;
				gsl_matrix_free(_dataPointsRight); _dataPointsRight = NULL;
				_needToFreePoints = false;
				return;
			}

			gsl_matrix_set(_dataPointsLeft, j, i, left[i](j,3));
			gsl_matrix_set(_dataPointsRight, j, i, right[i][j]);
		}
	}
	Init(n, _dataPointsLeft, _dataPointsRight);
}
Math::Horn::Horn(std::vector<Math::DoubleVector> &left, std::vector<Math::DoubleVector> &right) {
	if(left.size() != right.size())
		return;

	// allocate memory
	Alloc();
	int n = left.size();
	_dataPointsLeft = gsl_matrix_calloc(3, n);
	_dataPointsRight = gsl_matrix_calloc(3, n);
	_needToFreePoints = true;
	for(int i=0; i<n; i++) {
		for(int j=0; j<3; j++) {
			if (left[i].size() != 3 || right[i].size() != 3) {
				gsl_matrix_free(_dataPointsLeft); _dataPointsLeft = NULL;
				gsl_matrix_free(_dataPointsRight); _dataPointsRight = NULL;
				_needToFreePoints = false;
				return;
			}

			gsl_matrix_set(_dataPointsLeft, j, i, left[i][j]);
			gsl_matrix_set(_dataPointsRight, j, i, right[i][j]);
		}
	}
	Init(n, _dataPointsLeft, _dataPointsRight);
}

void Math::Horn::Alloc(void) {
	_matrixLeftToRight	= gsl_matrix_calloc(3,3);
	_eigenvector		= gsl_vector_calloc(4);
	_translationLtoR	= gsl_vector_calloc(3);
	
	_leftCentroid = gsl_vector_calloc(3);
	_rightCentroid = gsl_vector_calloc(3);
	_matrixM = gsl_matrix_calloc(3,3);
	_matrixN = gsl_matrix_calloc(4,4);
	_eigenvalues = gsl_vector_calloc(4);
	_eigenvectors = gsl_matrix_calloc(4,4);
	_work = gsl_eigen_symmv_alloc(4);
	_rotatedLeft	= gsl_vector_calloc(3);
	_rotatedRight = gsl_vector_calloc(3);
	_errorVector = gsl_vector_calloc(3);
	
	_dataPointsLeftOld	= gsl_matrix_calloc(3,3); 
	_dataPointsRightOld	= gsl_matrix_calloc(3,3); 

	_needToFreePoints = false;

	_n = 3;
}

void Math::Horn::Init(int n, gsl_matrix* dataPointsLeft, gsl_matrix* dataPointsRight)
{
	_dataPointsLeft = dataPointsLeft;
	_dataPointsRight = dataPointsRight;
	_maxPosEigenvalue = 0.0;
	_scale = 1.0;
	gsl_vector_set_all(_rightCentroid, 0);
	gsl_vector_set_all(_leftCentroid, 0);

	if(_dataPointsLeft != NULL && _dataPointsRight != NULL) {
		if(_n != n) {
			gsl_matrix_free(_dataPointsLeftOld);
			gsl_matrix_free(_dataPointsRightOld);

			_dataPointsLeftOld = gsl_matrix_calloc(3,n);
			_dataPointsRightOld = gsl_matrix_calloc(3,n);
		}
	}

	_n = n;
	_isCalculated = false;
}

Math::Horn::~Horn() {
	gsl_matrix_free(_matrixLeftToRight);
	gsl_vector_free(_eigenvector);
	gsl_vector_free(_translationLtoR);

	gsl_vector_free(_leftCentroid);
	gsl_vector_free(_rightCentroid);
	gsl_matrix_free(_matrixM);
	gsl_matrix_free(_matrixN);

	gsl_vector_free(_eigenvalues);
	gsl_matrix_free(_eigenvectors);
	gsl_eigen_symmv_free(_work);

	gsl_vector_free(_rotatedLeft);
	gsl_vector_free(_rotatedRight);
	gsl_vector_free(_errorVector);


	if(_dataPointsLeftOld != NULL)
		gsl_matrix_free(_dataPointsLeftOld);
	if(_dataPointsRightOld != NULL)
		gsl_matrix_free(_dataPointsRightOld);

	if(_needToFreePoints) {
		if(_dataPointsRight != NULL)
			gsl_matrix_free(_dataPointsRight);
		if(_dataPointsLeft != NULL)
			gsl_matrix_free(_dataPointsLeft);
	}
}

void Math::Horn::calculateTransformation(void) 
{
		findCentroids(_rightCentroid, _leftCentroid);
		relativizeDataPoints(_rightCentroid, _leftCentroid);
		construct3x3Matrix(_matrixM);
		construct4x4Matrix(_matrixN, _matrixM);
		findEigenvector(_matrixN);

		buildRotationMatrix();
		calculateScale();

		calculateTranslation(_rightCentroid, _leftCentroid);

		_isCalculated = true;
}

void Math::Horn::findCentroids(gsl_vector* rightCentroid, gsl_vector* leftCentroid) {
	gsl_vector_view actRight;
	gsl_vector_view actLeft;
	for (int i=0; i<_n; i++) {
		actRight = gsl_matrix_column(_dataPointsRight,i);
		gsl_vector_add(rightCentroid,&(actRight.vector));
		actLeft = gsl_matrix_column(_dataPointsLeft,i);
		gsl_vector_add(leftCentroid,&(actLeft.vector));
	}
	gsl_vector_scale(rightCentroid, (1.0/_n));
	gsl_vector_scale(leftCentroid, (1.0/_n));
}

void Math::Horn::relativizeDataPoints(gsl_vector* rightCentroid, gsl_vector* leftCentroid) {
	gsl_vector_view actRight;
	gsl_vector_view actLeft;
	gsl_matrix_memcpy(_dataPointsLeftOld, _dataPointsLeft);
	gsl_matrix_memcpy(_dataPointsRightOld,_dataPointsRight);
	for (int i=0; i<_n; i++) {
		actRight = gsl_matrix_column(_dataPointsRight,i);
		gsl_vector_sub(&(actRight.vector),rightCentroid);
		actLeft = gsl_matrix_column(_dataPointsLeft,i);
		gsl_vector_sub(&(actLeft.vector),leftCentroid);
	}
}

void Math::Horn::construct3x3Matrix(gsl_matrix* matrixM) {
	double S = 0.0;
	// calculate Sxx, Sxy, ... ,Szz
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			S = 0.0;
			for (int k=0; k<_n; k++) {
				S = S + gsl_matrix_get(_dataPointsLeft,i,k) * gsl_matrix_get(_dataPointsRight,j,k);
			}
			gsl_matrix_set(matrixM, i, j, S);
		}
	}
}

void Math::Horn::construct4x4Matrix(gsl_matrix* matrixN, gsl_matrix* matrixM) {
	double Sxx = gsl_matrix_get(matrixM,0,0);
	double Sxy = gsl_matrix_get(matrixM,0,1);
	double Sxz = gsl_matrix_get(matrixM,0,2);
	double Syx = gsl_matrix_get(matrixM,1,0);
	double Syy = gsl_matrix_get(matrixM,1,1);
	double Syz = gsl_matrix_get(matrixM,1,2);
	double Szx = gsl_matrix_get(matrixM,2,0);
	double Szy = gsl_matrix_get(matrixM,2,1);
	double Szz = gsl_matrix_get(matrixM,2,2);
	gsl_matrix_set(matrixN,0,0,(Sxx+Syy+Szz));
	gsl_matrix_set(matrixN,0,1,(Syz-Szy));
	gsl_matrix_set(matrixN,1,0,(Syz-Szy));
	gsl_matrix_set(matrixN,2,0,(Szx-Sxz));
	gsl_matrix_set(matrixN,0,2,(Szx-Sxz));
	gsl_matrix_set(matrixN,0,3,(Sxy-Syx));
	gsl_matrix_set(matrixN,3,0,(Sxy-Syx));
	gsl_matrix_set(matrixN,1,1,(Sxx-Syy-Szz));
	gsl_matrix_set(matrixN,1,2,(Sxy+Syx));
	gsl_matrix_set(matrixN,2,1,(Sxy+Syx));
	gsl_matrix_set(matrixN,1,3,(Szx+Sxz));
	gsl_matrix_set(matrixN,3,1,(Szx+Sxz));
	gsl_matrix_set(matrixN,2,2,(-Sxx+Syy-Szz));
	gsl_matrix_set(matrixN,2,3,(Syz+Szy));
	gsl_matrix_set(matrixN,3,2,(Syz+Szy));
	gsl_matrix_set(matrixN,3,3,(-Sxx-Syy+Szz));
}

void Math::Horn::findEigenvector(gsl_matrix* matrixN) {
	gsl_eigen_symmv (matrixN, _eigenvalues, _eigenvectors, _work);
	// sort in descending order in numerical value
	gsl_eigen_symmv_sort(_eigenvalues, _eigenvectors, GSL_EIGEN_SORT_VAL_DESC);
	_maxPosEigenvalue = gsl_vector_get(_eigenvalues,0);	
	gsl_matrix_get_col(_eigenvector, _eigenvectors, 0);
	// derive unit vector
	gsl_vector_scale(_eigenvector,(1.0/gsl_blas_dnrm2(_eigenvector)));
}

void Math::Horn::buildRotationMatrix() {
	double q0 = gsl_vector_get(_eigenvector, 0);
	double qx = gsl_vector_get(_eigenvector, 1);
	double qy = gsl_vector_get(_eigenvector, 2);
	double qz = gsl_vector_get(_eigenvector, 3);
	
	double q0_2 = q0*q0;
	double qx_2 = qx*qx;
	double qy_2 = qy*qy;
	double qz_2 = qz*qz;

	gsl_matrix_set(_matrixLeftToRight,0,0,(q0_2+qx_2-qy_2-qz_2));
	gsl_matrix_set(_matrixLeftToRight,0,1,(2*(qx*qy-q0*qz)));
	gsl_matrix_set(_matrixLeftToRight,0,2,(2*(qx*qz+q0*qy)));
	gsl_matrix_set(_matrixLeftToRight,1,0,(2*(qx*qy+q0*qz)));
	gsl_matrix_set(_matrixLeftToRight,1,1,(q0_2-qx_2+qy_2-qz_2));
	gsl_matrix_set(_matrixLeftToRight,1,2,(2*(qy*qz-q0*qx)));
	gsl_matrix_set(_matrixLeftToRight,2,0,(2*(qz*qx-q0*qy)));
	gsl_matrix_set(_matrixLeftToRight,2,1,(2*(qz*qy+q0*qx)));
	gsl_matrix_set(_matrixLeftToRight,2,2,(q0_2-qx_2-qy_2+qz_2));
}

void Math::Horn::calculateScale() {
	// symmetrische Skalierung
	double sumNormRr = 0.0;
	double sumNormRl = 0.0;
	gsl_vector_view RrView;
	gsl_vector_view RlView;
	double scalarProdRr;
	double scalarProdRl;
	for (int i=0; i<_n; i++) {
		RrView = gsl_matrix_column(_dataPointsRight,i);
		RlView = gsl_matrix_column(_dataPointsLeft,i);
		gsl_blas_ddot(&RrView.vector,&RrView.vector, &scalarProdRr);
		gsl_blas_ddot(&RlView.vector,&RlView.vector, &scalarProdRl);
		sumNormRr = sumNormRr + scalarProdRr;
		sumNormRl = sumNormRl + scalarProdRl;
	}
	_scale = sqrt(sumNormRr/sumNormRl); 
}

void Math::Horn::calculateTranslation(gsl_vector* rightCentroid, gsl_vector* leftCentroid) {

	gsl_vector *rotatedLeft = gsl_vector_calloc(3);
	//_translation = rightCentroid - R(leftCentroid)
	gsl_blas_dgemv(CblasNoTrans,1.0,_matrixLeftToRight,leftCentroid,1.0,rotatedLeft);
	gsl_vector_memcpy(_translationLtoR,rightCentroid);
	gsl_vector_sub(_translationLtoR,rotatedLeft);
	gsl_vector_free(rotatedLeft);
}

double Math::Horn::errorSumOfSquares(std::vector<std::vector<double> > *errors) {
	double error = 0.0;
	double skalarProd = 0.0;
	gsl_vector_view leftVector;
	gsl_vector_view rightVector;
	
	std::vector<double> v(_n,0);
	if(errors != NULL)
		errors->resize(4, v);
	for (int i=0; i<_n; i++) {
		gsl_vector_set_all(_rotatedLeft,0.0);
		leftVector = gsl_matrix_column(_dataPointsLeftOld,i);
		rightVector = gsl_matrix_column(_dataPointsRightOld,i);
		// R(l)
		gsl_blas_dgemv(CblasNoTrans,1.0,_matrixLeftToRight,&leftVector.vector,1.0,_rotatedLeft);
//		gsl_vector_scale(rotatedLeft,_scale);
		// R(l)+t
		gsl_vector_add(_rotatedLeft,_translationLtoR);
		gsl_vector_memcpy(_errorVector,&rightVector.vector);
		// r-[R(l)+t]
		gsl_vector_sub(_errorVector,_rotatedLeft);
		// ||r-[R(l)+t]||_2^2 = (r-[R(l)+t])^T*(r-[R(l)+t])
		gsl_blas_ddot(_errorVector, _errorVector, &skalarProd);
		if(errors != NULL) {
			errors->at(0)[i] = gsl_vector_get(_errorVector, 0);
			errors->at(1)[i] = gsl_vector_get(_errorVector, 1);
			errors->at(2)[i] = gsl_vector_get(_errorVector, 2);
			errors->at(3)[i] = sqrt(skalarProd);
		}
		error = error + skalarProd;
	}

	error = sqrt(error);
	return error;
}

double Math::Horn::maxError(void) {
	double error = 0.0;
	double skalarProd = 0.0;
	gsl_vector_view leftVector;
	gsl_vector_view rightVector;

	for (int i=0; i<_n; i++) {
		gsl_vector_set_all(_rotatedLeft,0.0);
		leftVector = gsl_matrix_column(_dataPointsLeftOld,i);
		rightVector = gsl_matrix_column(_dataPointsRightOld,i);
		// R(l)
		gsl_blas_dgemv(CblasNoTrans,1.0,_matrixLeftToRight,&leftVector.vector,1.0,_rotatedLeft);
	//	gsl_vector_scale(rotatedLeft,_scale);
		// R(l)+t
		gsl_vector_add(_rotatedLeft,_translationLtoR);
		gsl_vector_memcpy(_errorVector,&rightVector.vector);
		// r-[R(l)+t]
		gsl_vector_sub(_errorVector,_rotatedLeft);
		// ||r-[R(l)+t]||_2^2 = (r-[R(l)+t])^T*(r-[R(l)+t])
		gsl_blas_ddot(_errorVector, _errorVector, &skalarProd);
		error = dmax(error,skalarProd);
	}

	error = sqrt(error);
	return error;
}

double Math::Horn::dmax(double a, double b) {
	if (a>=b) {
		return a;
	}
	return b;
}

bool Math::Horn::getQuaternion(double quat[4]) {
	if(_isCalculated) {
		quat[0] = gsl_vector_get(_eigenvector, 0);
		quat[1] = gsl_vector_get(_eigenvector, 1);
		quat[2] = gsl_vector_get(_eigenvector, 2);
		quat[3] = gsl_vector_get(_eigenvector, 3);
		return true;
	} else 
		return false;
}

bool Math::Horn::getTranslation(double trans[3]) {
	if(_isCalculated) {
		trans[0] = gsl_vector_get(_translationLtoR, 0);
		trans[1] = gsl_vector_get(_translationLtoR, 1);
		trans[2] = gsl_vector_get(_translationLtoR, 2);
		return true;
	} else 
		return false;
}

bool Math::Horn::getScale(double &scale) {
	if(_isCalculated) {
		scale = _scale;
		return true;
	} else 
		return false;
}

bool Math::Horn::getTransform(double quat[4], double trans[3], double &scale) {
	if(_isCalculated) {
		quat[0] = gsl_vector_get(_eigenvector, 0);
		quat[1] = gsl_vector_get(_eigenvector, 1);
		quat[2] = gsl_vector_get(_eigenvector, 2);
		quat[3] = gsl_vector_get(_eigenvector, 3);
		trans[0] = gsl_vector_get(_translationLtoR, 0);
		trans[1] = gsl_vector_get(_translationLtoR, 1);
		trans[2] = gsl_vector_get(_translationLtoR, 2);
		scale = _scale;
		return true;
	} else 
		return false;
}

bool Math::Horn::getTransform(double matrixRowWise[12]) {
	if(_isCalculated) {
		matrixRowWise[0] = gsl_matrix_get(_matrixLeftToRight, 0, 0);
		matrixRowWise[1] = gsl_matrix_get(_matrixLeftToRight, 0, 1);
		matrixRowWise[2] = gsl_matrix_get(_matrixLeftToRight, 0, 2);
		matrixRowWise[3] = gsl_vector_get(_translationLtoR, 0);
		matrixRowWise[4] = gsl_matrix_get(_matrixLeftToRight, 1, 0);
		matrixRowWise[5] = gsl_matrix_get(_matrixLeftToRight, 1, 1);
		matrixRowWise[6] = gsl_matrix_get(_matrixLeftToRight, 1, 2);
		matrixRowWise[7] = gsl_vector_get(_translationLtoR, 1);
		matrixRowWise[8] = gsl_matrix_get(_matrixLeftToRight, 2, 0);
		matrixRowWise[9] = gsl_matrix_get(_matrixLeftToRight, 2, 1);
		matrixRowWise[10] = gsl_matrix_get(_matrixLeftToRight, 2, 2);
		matrixRowWise[11] = gsl_vector_get(_translationLtoR, 2);
		return true;
	} else 
		return false;
}

bool Math::Horn::getTransform(Matrix &matrix) {
	if(_isCalculated) {
		matrix = Matrix::eye(4);
		for(int i=0; i<3; i++)
			for(int j=0; j<3; j++)
				matrix(i,j) = gsl_matrix_get(_matrixLeftToRight, i, j);
		matrix(0,3) = gsl_vector_get(_translationLtoR, 0);
		matrix(1,3) = gsl_vector_get(_translationLtoR, 1);
		matrix(2,3) = gsl_vector_get(_translationLtoR, 2);
		return true;
	} else {
		return false;
	}
}

bool Math::Horn::getMeanError(double &error) {
	if(_isCalculated) {
		error = errorSumOfSquares() / sqrt((double)_n);
		return true;
	} else 
		return false;
}

bool Math::Horn::getMaxError(double &error) {
	if(_isCalculated) {
		error = maxError();
		return true;
	} else 
		return false;
}

