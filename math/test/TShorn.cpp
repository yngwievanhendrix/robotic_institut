#include "cpptest/include/cpptest.h"
#include "TShorn.h"
#include "../include/Horn.h"

using namespace std;
using namespace Math;

TShorn::TShorn()
{
	prepare();
	TEST_ADD(TShorn::testHorn);
}

TShorn::~TShorn()
{
	close();
}

void TShorn::prepare()
{
};

void TShorn::close()
{
};

void TShorn::testHorn()
{
	using namespace Math;
	
	// random system matrix
	Matrix system(4,4), system_, tmp(4,4);
	RandMatrix(system, 180, 1000);

	// random points
	int num = 200;
	std::vector<Matrix> left(num, Matrix::eye(4));
	std::vector<Matrix> right(num, Matrix::eye(4));
	for(int i=0; i<num; i++) {
		RandMatrix(left[i], 180, 1000);
		Matrix::MatrixMultiply(system, left[i], right[i]);
	}
	Horn horn(left, right);
	horn.calculateTransformation();
	horn.getTransform(system_);

	Matrix::invert(system_, tmp);
	Matrix::MatrixMultiply(system, tmp, system_);

	TEST_ASSERT_MSG(horn.errorSumOfSquares() < 0.01, "Large error detected");

	TEST_ASSERT_MSG(fabs(system_.norm2() - 1) < 0.01, "Horn returned non-matching matrix" );
}

