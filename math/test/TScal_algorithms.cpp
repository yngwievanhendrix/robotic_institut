#include "cpptest/include/cpptest.h"
#include "TScal_algorithms.h"

using namespace std;
using namespace Math;

TScal_algorithms::TScal_algorithms()
{
	prepare();
	TEST_ADD(TScal_algorithms::test_Tensorproduct);
	TEST_ADD(TScal_algorithms::test_CalibAlgs);
}

TScal_algorithms::~TScal_algorithms()
{
	close();
}
    
void TScal_algorithms::prepare()
{
	m_system(4,4);
	RandMatrix(m_system, 180, 1000);
	m_offset(4,4);
	RandMatrix(m_offset, 180, 100);

	int num = 50;
	m_leftlist.resize(num);
	m_rightlist.resize(num);
	Matrix temp(4,4);
	for (int i=0; i<num; i++)
	{
		m_leftlist[i].resize(4,4);
		m_rightlist[i].resize(4,4);
		RandMatrix(m_leftlist[i], 180, 1000);
		// right = system * left * offset;
		Matrix::MatrixMultiply(m_system, m_leftlist[i], temp);
		Matrix::MatrixMultiply(temp, m_offset, m_rightlist[i]);
	}
};

void TScal_algorithms::close()
{

};


void TScal_algorithms::test_Tensorproduct()
{
	DoubleVector tempPos(3);
	DoubleVector tempQuat(4);
	DoubleVector q(4);
	gsl_matrix* tempTensor = gsl_matrix_alloc(4,4);
	Matrix temp(4,4);
	Matrix temp2(4,4);
	double scale, rot, trans;
	for (int i = 0; i < m_leftlist.size(); i++) //we have n-1 linearly independent relations between the views
	{
		Math::MatToQuat(m_leftlist[i], tempPos, tempQuat);
		Math::quaternion2Tensor(tempTensor, tempQuat);
		Math::tensor2Quaternion(tempTensor, tempQuat);
		Math::QuatToMat(tempPos, tempQuat, temp);
		Math::Matrix::invert(temp, temp2);
		Math::Matrix::MatrixMultiply(temp2, m_leftlist[i], temp);
		temp.TransformErrors(scale, rot, trans);
		TEST_ASSERT_MSG(scale < 0.01,"Scale wrong for quaternion and tensor product");
		TEST_ASSERT_MSG(rot < 0.01,"Rotation wrong for quaternion and tensor product");
		TEST_ASSERT_MSG(trans < 0.01,"Translation wrong for quaternion and tensor product");
	}
}



void TScal_algorithms::test_CalibAlgs()
{
	Calibration::CCalibAlgorithm* alg;
	int num = 0;
	DoubleVector errors(9);
	Matrix system(4,4);
	Matrix offset(4,4);
	Matrix offsetinv(4,4);
	Matrix systeminv(4,4);

	double scale, rot, trans;

	char msg[1024];

	for(Calibration::kCalibrationAlgorithms k = Calibration::kCALIB_ALG_MIN; k < Calibration::kCALIB_ALG_COUNT; k++) {

		alg = m_bk.Get(k);
		alg->SetPoints(m_leftlist, m_rightlist);
		
		alg->ComputeMatrices(num, errors);
		alg->GetMatrices(system, offset);

		sprintf(msg, "%s", (alg->GetName() + ": ").c_str());

		Matrix::invert(offset, offsetinv);
		Matrix::invert(system, systeminv);

		systeminv.rightMultiply(m_system, system);
		offset.rightMultiply(m_offset, offsetinv);

		system.TransformErrors(scale, rot, trans);

		TEST_ASSERT_MSG(scale < 0.01,strcat(msg, "Scale wrong for system matrix"));
		TEST_ASSERT_MSG(rot < 0.01,strcat(msg, "Rotation wrong for system matrix"));
		TEST_ASSERT_MSG(trans < 0.01,strcat(msg, "Translation wrong for system matrix"));

		offsetinv.TransformErrors(scale, rot, trans);
		if(alg->GetType() != Calibration::kCALIB_QR15) {
			TEST_ASSERT_MSG(scale < 0.01,strcat(msg, "Scale wrong for offset matrix"));
			TEST_ASSERT_MSG(rot < 0.01,strcat(msg, "Rotation wrong for offset matrix"));
		}
		TEST_ASSERT_MSG(trans < 0.01,strcat(msg, "Translation wrong for offset matrix"));
	}
}
