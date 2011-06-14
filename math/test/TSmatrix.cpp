#include "cpptest/include/cpptest.h"
#include "TSmatrix.h"
#include "../include/matrix.h"

using namespace std;
using namespace Math;

TSmatrix::TSmatrix()
{
	prepare();
	TEST_ADD(TSmatrix::testInversion);
	TEST_ADD(TSmatrix::testMatrixMultiply);
}

TSmatrix::~TSmatrix()
{
	close();
}
    
void TSmatrix::prepare()
{
};

void TSmatrix::close()
{

};

void TSmatrix::testMatrixMultiply()
{
	Matrix test= Matrix::ones(3,3);
	Matrix AAA = Matrix::ones(3,1);
	Matrix BBB = Matrix::zeros(3,1);
	Matrix::MatrixMultiply(test, AAA, BBB);
}

void TSmatrix::testInversion()
{
	Math::Matrix m = Math::Matrix::zeros(4,4);
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			if (j<=i)
				m(i,j) = 1;
	m(3,3) = 1;
	Math::Matrix minv(4,4);
	Math::Matrix::invert(m, minv);
	

	Math::Matrix identity = Math::Matrix::eye(4);
	Math::Matrix minv2(4,4);
	m.linearSolve(minv2, identity);

	minv.sub(minv2);
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			TEST_ASSERT(fabs(minv(i,j)) < 0.0001); 
}

