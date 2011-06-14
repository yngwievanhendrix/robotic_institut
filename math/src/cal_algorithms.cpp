#include "math/include/cal_algorithms.h"
#include <fstream>

using namespace Calibration;

Calibration::kCalibrationDataTypeAccepted Calibration::ConvertDataType( kROTTYPE type )
{
	switch(type) {
		case rotEULER:
		case rotMATRIX:
		case rotQUATERNION:
		case rotAXISANGLE:
		case rotYPR:
			return kCALIB_DTYPE_POSROT;
		case rotFORCETORQUE:
			return kCALIB_DTYPE_FORCETORQUE;
		case rotNONE:
		case rotUNKNOWN:
		default:
			return kCALIB_DTYPE_UNKNOWN;
	}
}

Calibration::CCalibAlgorithm::CCalibAlgorithm(kCalibrationAlgorithms type, std::string name, kCalibrationDataTypeAccepted typeAcceptedLeft, kCalibrationDataTypeAccepted typeAcceptedRight) :
	m_Type(type),
	m_szName(name),
	m_fPointsSet(false),
	m_fIsCalibrated(false),
	m_fDefaultOrtho(false)

{ 	
	m_TypeAccepted[0] = typeAcceptedLeft;
	m_TypeAccepted[1] = typeAcceptedRight;
	m_Matrices[0].resize(4,4); m_Matrices[0] = Matrix::eye(4);
	m_Matrices[1].resize(4,4); m_Matrices[1] = Matrix::eye(4);
}

CAL_ECODE Calibration::CCalibAlgorithm::SetPoints(const std::vector<tDataObject>& system_1, const std::vector<tDataObject>& system_2) {
	// check dimensions
	if(system_1.size() != system_2.size())
		return TRC_VEC_INVALID_DIM;

	// clear
	m_Data[0].clear();
	m_Data[1].clear();
	m_Points[0].clear();
	m_Points[1].clear();

	// check data type
	for(size_t i=0; i<system_1.size(); i++)
		if(!(m_TypeAccepted[0] & ConvertDataType(system_1[i].RotationType)) || !(m_TypeAccepted[1] & ConvertDataType(system_2[i].RotationType)) )
			return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	if(m_TypeAccepted[0] == kCALIB_DTYPE_POSROT) {
		// check rotation type
		for(size_t i=0; i<system_1.size(); i++)
			if((system_1[i].RotationType != rotMATRIX && system_1[i].RotationType != rotQUATERNION))
				return CAL_UNSUPPORTED_ROTTYPE;

		// store
		m_Points[0].resize(system_1.size());

		for(size_t i=0; i<system_1.size(); i++) {
			if(system_1[i].RotationType == rotMATRIX)
				m_Points[0][i].inject(system_1[i].matrix);
			else if(system_1[i].RotationType == rotQUATERNION)
				QuatToMat(system_1[i].trans, system_1[i].quat, m_Points[0][i]);
		}
	} else if(m_TypeAccepted[0] == kCALIB_DTYPE_FORCETORQUE) {
		// check size
		for(size_t i=0; i<system_1.size(); i++) {
			if(system_1[i].forcetorque.size() != 6)
				return TRC_VEC_INVALID_DIM;
		}

		// store
		m_Data[0].resize(system_1.size());

		for(size_t i=0; i<system_1.size(); i++) {
			m_Data[0][i] = system_1[i].forcetorque;
		}
	} else {
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
	}

	if(m_TypeAccepted[1] == kCALIB_DTYPE_POSROT) {
		// check rotation type
		for(size_t i=0; i<system_2.size(); i++)
			if((system_2[i].RotationType != rotMATRIX && system_2[i].RotationType != rotQUATERNION))
				return CAL_UNSUPPORTED_ROTTYPE;

		// store
		m_Points[1].resize(system_1.size());

		for(size_t i=0; i<system_2.size(); i++) {
			if(system_2[i].RotationType == rotMATRIX)
				m_Points[1][i].inject(system_2[i].matrix);
			else if(system_2[i].RotationType == rotQUATERNION)
				QuatToMat(system_2[i].trans, system_2[i].quat, m_Points[1][i]);
		}
	} else if(m_TypeAccepted[0] == kCALIB_DTYPE_FORCETORQUE) {
		// check size
		for(size_t i=0; i<system_2.size(); i++) {
			if(system_2[i].forcetorque.size() != 6)
				return TRC_VEC_INVALID_DIM;
		}

		// store
		m_Data[1].resize(system_2.size());

		for(size_t i=0; i<system_2.size(); i++) {
			m_Data[1][i] = system_2[i].forcetorque;
		}
	} else {
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
	}

	m_fPointsSet = true;

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::SetPoints(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_1_pos.size() != system_1_rot.size() || system_1_pos.size() != system_2_pos.size() || system_2_pos.size() != system_2_rot.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_1_pos.size();

	m_Points[0].resize(num);
	m_Points[1].resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_1_pos[i], system_1_rot[i], m_Points[0][i]);
		QuatToMat(system_2_pos[i], system_2_rot[i], m_Points[1][i]);
	}

	m_fPointsSet = true;
	m_fIsCalibrated = false;

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::SetPoints(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<Matrix>& system_2) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_1_pos.size() != system_1_rot.size() || system_1_pos.size() != system_2.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_1_pos.size();

	m_Points[0].resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_1_pos[i], system_1_rot[i], m_Points[0][i]);
		if(system_2[i].size(0) != 4 || system_2[i].size(1) != 4)
			return TRC_MAT_INVALID_DIM;
	}

	m_Points[1] = system_2;

	m_fPointsSet = true;
	m_fIsCalibrated = false;

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::SetPoints(const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_2_pos.size() != system_2_rot.size() || system_2_pos.size() != system_1.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_2_pos.size();

	m_Points[1].resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_2_pos[i], system_2_rot[i], m_Points[1][i]);
		if(system_1[i].size(0) != 4 || system_1[i].size(1) != 4)
			return TRC_MAT_INVALID_DIM;
	}

	m_Points[0] = system_1;

	m_fPointsSet = true;
	m_fIsCalibrated = false;

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::SetPoints(const std::vector<Matrix>& system_1, const std::vector<Matrix>& system_2) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_2.size() != system_1.size())
		return TRC_VEC_INVALID_DIM;

	for(size_t i=0; i<system_1.size(); i++)
		if(system_1[i].size(0) != 4 || system_1[i].size(1) != 4 || system_2[i].size(0) != 4 || system_2[i].size(1) != 4)
			return TRC_MAT_INVALID_DIM;

	m_Points[0] = system_1;
	m_Points[1] = system_2;

	m_fPointsSet = true;
	m_fIsCalibrated = false;

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::GetMatrices(Matrix &system, Matrix &offset) {
	if(!m_fIsCalibrated)
		return CAL_NOT_CALIBRATED;

	offset.resize(4,4);
	offset.inject(m_Matrices[0]);
	system.resize(4,4);
	system.inject(m_Matrices[1]);

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(DoubleVector& err, bool ortho) {
	if(m_Points[0].size() > 0 && m_Points[1].size() > 0)
		return ComputeDeviations(m_Points[0], m_Points[1], err, ortho);
	else if(m_Points[0].size() > 0 && m_Data[1].size() > 0)
		return ComputeDeviations(m_Points[0], m_Data[1], err, ortho);
	else if(m_Data[0].size() > 0 && m_Points[1].size() > 0)
		return ComputeDeviations(m_Data[0], m_Points[1], err, ortho);
	else if(m_Data[0].size() > 0 && m_Data[1].size() > 0)
		return ComputeDeviations(m_Data[0], m_Data[1], err, ortho);
	else
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(const std::vector<tDataObject>& system_1, const std::vector<tDataObject>& system_2, DoubleVector &err, bool ortho) {
	// check dimensions
	if(system_1.size() != system_2.size())
		return TRC_VEC_INVALID_DIM;

	std::vector<DoubleVector> data[2];
	std::vector<Matrix> points[2];

	// check data type
	for(size_t i=0; i<system_1.size(); i++)
		if(!(m_TypeAccepted[0] & ConvertDataType(system_1[i].RotationType)) || !(m_TypeAccepted[1] & ConvertDataType(system_2[i].RotationType)) )
			return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	if(m_TypeAccepted[0] == kCALIB_DTYPE_POSROT) {
		// check rotation type
		for(size_t i=0; i<system_1.size(); i++)
			if((system_1[i].RotationType != rotMATRIX && system_1[i].RotationType != rotQUATERNION))
				return CAL_UNSUPPORTED_ROTTYPE;

		// store
		points[0].resize(system_1.size());

		for(size_t i=0; i<system_1.size(); i++) {
			if(system_1[i].RotationType == rotMATRIX)
				points[0][i].inject(system_1[i].matrix);
			else if(system_1[i].RotationType == rotQUATERNION)
				QuatToMat(system_1[i].trans, system_1[i].quat, points[0][i]);
		}
	} else if(m_TypeAccepted[0] == kCALIB_DTYPE_FORCETORQUE) {
		// check size
		for(size_t i=0; i<system_1.size(); i++) {
			if(system_1[i].forcetorque.size() != 6)
				return TRC_VEC_INVALID_DIM;
		}

		// store
		data[0].resize(system_1.size());

		for(size_t i=0; i<system_1.size(); i++) {
			data[0][i] = system_1[i].forcetorque;
		}
	} else {
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
	}

	if(m_TypeAccepted[1] == kCALIB_DTYPE_POSROT) {
		// check rotation type
		for(size_t i=0; i<system_2.size(); i++)
			if((system_2[i].RotationType != rotMATRIX && system_2[i].RotationType != rotQUATERNION))
				return CAL_UNSUPPORTED_ROTTYPE;

		// store
		points[1].resize(system_1.size());

		for(size_t i=0; i<system_2.size(); i++) {
			if(system_2[i].RotationType == rotMATRIX)
				points[1][i].inject(system_2[i].matrix);
			else if(system_2[i].RotationType == rotQUATERNION)
				QuatToMat(system_2[i].trans, system_2[i].quat, points[1][i]);
		}
	} else if(m_TypeAccepted[0] == kCALIB_DTYPE_FORCETORQUE) {
		// check size
		for(size_t i=0; i<system_2.size(); i++) {
			if(system_2[i].forcetorque.size() != 6)
				return TRC_VEC_INVALID_DIM;
		}

		// store
		data[1].resize(system_2.size());

		for(size_t i=0; i<system_2.size(); i++) {
			data[1][i] = system_2[i].forcetorque;
		}
	} else {
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
	}

	if(points[0].size() > 0 && points[1].size() > 0)
		return ComputeDeviations(points[0], points[1], err, ortho);
	else if(points[0].size() > 0 && data[1].size() > 0)
		return ComputeDeviations(points[0], data[1], err, ortho);
	else if(data[0].size() > 0 && points[1].size() > 0)
		return ComputeDeviations(data[0], points[1], err, ortho);
	else if(data[0].size() > 0 && data[1].size() > 0)
		return ComputeDeviations(data[0], data[1], err, ortho);
	else
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(const std::vector<Matrix>& system_1, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	if(!m_fIsCalibrated)
		return CAL_NOT_CALIBRATED;

	if(system_1.size() != system_2.size())
		return TRC_VEC_INVALID_DIM;

	Matrix tmp1(4,4), tmp2(4,4), tmp3(4,4), e(4,4), P1(4,4), P2(4,4);
	double es, er, et;

	// posDevMean, posDevRMS, posDevMax, rotDevMean, rotDevRMS, rotDevMax
	err.resize(9, 0);
	err = 0 * err;
	int num = (int)system_1.size();

	DoubleVector pos(3,0), quat(4,0);

	for(int i=0; i<num; i++) {
		if(system_1[i].size(0) != 4 || system_2[i].size(0) != 4)
			return TRC_MAT_INVALID_DIM;

		// generate rotation matrices
		P1 = system_1[i];
		P2 = system_2[i];

		m_Matrices[1].linearSolve(tmp1, Matrix::eye(4));

		Matrix::MatrixMultiply(tmp1, P2, tmp2);
		Matrix::MatrixMultiply(tmp2, m_Matrices[0], tmp1);

		P1.linearSolve(tmp2, Matrix::eye(4));

		if(ortho)
			tmp1.orthonormalise();

		Matrix::MatrixMultiply(tmp2, tmp1, e);

		/* m_Matrices[0].linearSolve(tmp1, Matrix::eye(4));
		Matrix::MatrixMultiply(P1, tmp1, tmp2);
		Matrix::MatrixMultiply(m_Matrices[1], tmp2, tmp1);

		P2.linearSolve(tmp2, Matrix::eye(4));

		if(ortho)
			tmp1.orthonormalise();

		Matrix::MatrixMultiply(tmp2, tmp1, e); */

		e.TransformErrors(es, er, et);

		err[0] += et;
		err[1] += et*et;
		err[2] = (et > err[2])?et:err[2];
		err[3] += er;
		err[4] += er*er;
		err[5] = (er > err[5])?er:err[5];
		err[6] += es;	
		err[7] += es*es;
		err[8] = (es > err[8])?es:err[8];
	}

	err[0] /= (double)num;
	err[3] /= (double)num;
	err[6] /= (double)num;

	err[1] = sqrt(err[1] / (double)num);
	err[4] = sqrt(err[4] / (double)num);
	err[7] = sqrt(err[7] / (double)num);

	// check if system supports rotation and scaling
	if(m_Type == kCALIB_QR15) {
		err[3] = err[4] = err[5] = err[6] = err[7] = err[8] = -1;
	}

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos, DoubleVector &err, bool ortho) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_1_pos.size() != system_1_rot.size() || system_1_pos.size() != system_2_pos.size() || system_2_pos.size() != system_2_rot.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_1_pos.size();

	std::vector<Matrix> s1, s2;
	s1.resize(num);
	s2.resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_1_pos[i], system_1_rot[i], s1[i]);
		QuatToMat(system_2_pos[i], system_2_rot[i], s2[i]);
	}

	return ComputeDeviations(s1, s2, err, ortho);
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_1_pos.size() != system_1_rot.size() || system_1_pos.size() != system_2.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_1_pos.size();

	std::vector<Matrix> s1;
	s1.resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_1_pos[i], system_1_rot[i], s1[i]);
		if(system_2[i].size(0) != 4 || system_2[i].size(1) != 4)
			return TRC_MAT_INVALID_DIM;
	}

	return ComputeDeviations(s1, system_2, err, ortho);
}

CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations(const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos, DoubleVector &err, bool ortho) {
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// check dimensions
	if(system_2_pos.size() != system_2_rot.size() || system_2_pos.size() != system_1.size())
		return TRC_VEC_INVALID_DIM;

	size_t num = system_2_pos.size();

	std::vector<Matrix> s2;
	s2.resize(num);

	for(size_t i=0; i<num; i++) {
		QuatToMat(system_2_pos[i], system_2_rot[i], s2[i]);
		if(system_1[i].size(0) != 4 || system_1[i].size(1) != 4)
			return TRC_MAT_INVALID_DIM;
	}

	return ComputeDeviations(system_1, s2, err, ortho);
}

Calibration::CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations( const std::vector<DoubleVector>& system_1, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho )
{
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_FORCETORQUE) || !(m_TypeAccepted[1] & kCALIB_DTYPE_POSROT))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// TODO: compute deviations!

	// currently not supported
	return CAL_ALG_DOESNT_SUPPORT_DTYPE;
}

Calibration::CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations( const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2, DoubleVector &err, bool ortho )
{
	// check data type
	if(!(m_TypeAccepted[0] & kCALIB_DTYPE_POSROT) || !(m_TypeAccepted[1] & kCALIB_DTYPE_FORCETORQUE))
		return CAL_ALG_DOESNT_SUPPORT_DTYPE;

	// TODO: compute deviations!

	// currently not supported
	return CAL_ALG_DOESNT_SUPPORT_DTYPE;
}

Calibration::CAL_ECODE Calibration::CCalibAlgorithm::ComputeDeviations( const std::vector<DoubleVector>& system_1, const std::vector<DoubleVector>& system_2, DoubleVector &err, bool ortho )
{
	// currently not supported
	return CAL_ALG_DOESNT_SUPPORT_DTYPE;
}

CAL_ECODE Calibration::CCalibAlgorithmQR15::ComputeMatrices(int &numPoints, DoubleVector &err) {
	/* if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;
	
	numPoints = m_Points[0].size();

	Matrix A(3*numPoints, 15), b(3*numPoints, 1), x(15,1);
	A.inject(Matrix::zeros(3*numPoints, 15));
	b.inject(Matrix::zeros(3*numPoints, 1));

	// Ax=b im 2D-Fall
	// Zeile von A: a11-a21  a12-a22  -y1 -y2 -1 y1 y2 1
	// x: x1 x2 b11 b12 b13 b21 b22 b23
	// zugehöriger Eintrag in b: a23-a13
	for (int i = 0; i < numPoints; i++)
	{
		A(3*i, 0) = m_Points[0][i](0,0);
		A(3*i, 1) = m_Points[0][i](0,1);
		A(3*i, 2) = m_Points[0][i](0,2);
		A(3*i, 3) = -m_Points[1][i](0,3);
		A(3*i, 4) = -m_Points[1][i](1,3);
		A(3*i, 5) = -m_Points[1][i](2,3);
		A(3*i, 6) = -1;

		A(3*i+1, 0) = m_Points[0][i](1,0);
		A(3*i+1, 1) = m_Points[0][i](1,1);
		A(3*i+1, 2) = m_Points[0][i](1,2);
		A(3*i+1, 7) = -m_Points[1][i](0,3);
		A(3*i+1, 8) = -m_Points[1][i](1,3);
		A(3*i+1, 9) = -m_Points[1][i](2,3);
		A(3*i+1, 10) = -1;

		A(3*i+2, 0) = m_Points[0][i](2,0);
		A(3*i+2, 1) = m_Points[0][i](2,1);
		A(3*i+2, 2) = m_Points[0][i](2,2);
		A(3*i+2, 11) = -m_Points[1][i](0,3);
		A(3*i+2, 12) = -m_Points[1][i](1,3);
		A(3*i+2, 13) = -m_Points[1][i](2,3);
		A(3*i+2, 14) = -1;

		b(3*i, 0) = -m_Points[0][i](0,3);
		b(3*i+1, 0) = -m_Points[0][i](1,3);
		b(3*i+2, 0) = -m_Points[0][i](2,3);
	}   

	A.QRSolve(x, b, false);
	m_Matrices[0].inject(Matrix::eye(4));

	Matrix tmp(4);
	tmp.inject(Matrix::eye(4));

	m_Matrices[0](0,3) = -x(0,0);
	m_Matrices[0](1,3) = -x(1,0);
	m_Matrices[0](2,3) = -x(2,0);

	tmp(0,0) = x(3,0);
	tmp(0,1) = x(4,0);
	tmp(0,2) = x(5,0);
	tmp(0,3) = x(6,0);

	tmp(1,0) = x(7,0);
	tmp(1,1) = x(8,0);
	tmp(1,2) = x(9,0);
	tmp(1,3) = x(10,0);

	tmp(2,0) = x(11,0);
	tmp(2,1) = x(12,0);
	tmp(2,2) = x(13,0);
	tmp(2,3) = x(14,0);	

	tmp.linearSolve(m_Matrices[1], Matrix::eye(4));

	m_fIsCalibrated = true;

	ComputeDeviations(err, true);

	return CAL_SUCCESS; */
	
	if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;

	numPoints = m_Points[0].size();

	Matrix A(3*numPoints, 15), b(3*numPoints, 1), x(15,1);
	A.inject(Matrix::zeros(3*numPoints, 15));
	b.inject(Matrix::zeros(3*numPoints, 1));

	// Ax=b im 2D-Fall
	// Zeile von A: a11-a21  a12-a22  -y1 -y2 -1 y1 y2 1
	// x: x1 x2 b11 b12 b13 b21 b22 b23
	// zugehöriger Eintrag in b: a23-a13
	for (int i = 0; i < numPoints; i++)
	{
		A(3*i, 0) = m_Points[1][i](0,0);
		A(3*i, 1) = m_Points[1][i](0,1);
		A(3*i, 2) = m_Points[1][i](0,2);
		A(3*i, 3) = -m_Points[0][i](0,3);
		A(3*i, 4) = -m_Points[0][i](1,3);
		A(3*i, 5) = -m_Points[0][i](2,3);
		A(3*i, 6) = -1;

		A(3*i+1, 0) = m_Points[1][i](1,0);
		A(3*i+1, 1) = m_Points[1][i](1,1);
		A(3*i+1, 2) = m_Points[1][i](1,2);
		A(3*i+1, 7) = -m_Points[0][i](0,3);
		A(3*i+1, 8) = -m_Points[0][i](1,3);
		A(3*i+1, 9) = -m_Points[0][i](2,3);
		A(3*i+1, 10) = -1;

		A(3*i+2, 0) = m_Points[1][i](2,0);
		A(3*i+2, 1) = m_Points[1][i](2,1);
		A(3*i+2, 2) = m_Points[1][i](2,2);
		A(3*i+2, 11) = -m_Points[0][i](0,3);
		A(3*i+2, 12) = -m_Points[0][i](1,3);
		A(3*i+2, 13) = -m_Points[0][i](2,3);
		A(3*i+2, 14) = -1;

		b(3*i, 0) = -m_Points[1][i](0,3);
		b(3*i+1, 0) = -m_Points[1][i](1,3);
		b(3*i+2, 0) = -m_Points[1][i](2,3);
	}   

	A.QRSolve(x, b, false);
	m_Matrices[0].inject(Matrix::eye(4));

	m_Matrices[0](0,3) = x(0,0);
	m_Matrices[0](1,3) = x(1,0);
	m_Matrices[0](2,3) = x(2,0);

	m_Matrices[1](0,0) = x(3,0);
	m_Matrices[1](0,1) = x(4,0);
	m_Matrices[1](0,2) = x(5,0);
	m_Matrices[1](0,3) = x(6,0);

	m_Matrices[1](1,0) = x(7,0);
	m_Matrices[1](1,1) = x(8,0);
	m_Matrices[1](1,2) = x(9,0);
	m_Matrices[1](1,3) = x(10,0);

	m_Matrices[1](2,0) = x(11,0);
	m_Matrices[1](2,1) = x(12,0);
	m_Matrices[1](2,2) = x(13,0);
	m_Matrices[1](2,3) = x(14,0);	

	m_fIsCalibrated = true;

	ComputeDeviations(err, true);

	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithmForceTorque::ComputeMatrices(int &numPoints, DoubleVector &err) {
	if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;

	//Linear force bias:
	const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmder;
    gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc (T, numPoints, 4);

	//initial guess
	double cfx = 1;
	double cfy = 2;
	double cfz = 3;
	double fmagnitude = 4;

	double null_i = 0;
	gsl_vector *f;

	for(int i = 0; i < numPoints; i++){
		null_i = sqrt((m_Data[0][i].at(0) - cfx)*(m_Data[0][i].at(0) - cfx) + (m_Data[0][i].at(1) - cfy)*(m_Data[0][i].at(1) - cfy) + (m_Data[0][i].at(2) - cfz)*(m_Data[0][i].at(2) - cfz)) - fmagnitude;
		gsl_vector_set(f,i,null_i);
	}

	gsl_multifit_function_fdf fdf;
//	fdf.f = f;
	

	m_bias[0] = 1;
	m_bias[1] = 2;
	m_bias[2] = 3;
	m_bias[3] = 4;
	m_bias[4] = 5;
	m_bias[5] = 6;

	m_rot_sensor_effector = m_rot_sensor_effector.zeros(3,3);

	m_centroid_sensor[0] = 8;
	m_centroid_sensor[1] = 9;
	m_centroid_sensor[2] = 10;

	m_magnitude = 14.763;
	m_rot_tool_effector = m_rot_tool_effector.zeros(4,4);

	m_centroid_tool[0] = 1;
	m_centroid_tool[1] = 2;
	m_centroid_tool[2] = 3;

	m_fIsCalibrated = true;
	
	return CAL_SUCCESS;
}

CAL_ECODE Calibration::CCalibAlgorithmQR24::ComputeMatrices(int &numPoints, DoubleVector &err) {
	if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;

	int a = m_Points[0].size();

	Matrix ETMK(12*a, 24); ETMK = Matrix::zeros(12*a, 24);
	Matrix b(12*a, 1); b = Matrix::zeros(12*a, 1);
	Matrix lMatrix(4,4), rMatrix(4,4);
	lMatrix = Matrix::zeros(4);
	rMatrix = Matrix::zeros(4);

	// durch alle Versuche v
	for (int v=0; v<a; v++) {
		rMatrix.inject(m_Points[0][v]);
		lMatrix.inject(m_Points[1][v]);
		int i=0, j=0, z=0;
		for (int m=0; m<12; m++) {
			i = m/4;
			j = m%4;
			z = v*12 + m;
			for (int k=0; k<4; k++) {
				// ET
				if (k<3) {
					ETMK(z, 4*k + j) = lMatrix(i, k);
				} else {
					// k=3 (4. Zeile) gilt (0 0 0 1)
					if (j==3) {
						b(z, 0) = b(z, 0) - lMatrix(i,3) / (m_fCalcInMetres ? 1000.0 : 1.0);
					}
				}
				// MK
				ETMK(v*12+m, 12 + i * 4 + k) = - rMatrix(k, j) / ((m_fCalcInMetres && j==3 && k < 3) ? 1000.0 : 1.0);
			}
		}
	}

	m_Matrices[0].resize(4,4); m_Matrices[0] = Matrix::zeros(4);
	m_Matrices[1].resize(4,4); m_Matrices[1] = Matrix::zeros(4);
	
	Matrix x(24,1);
	
	ETMK.QRSolve(x, b, true);
	for (int m=0; m<12; m++) {
		m_Matrices[1](m/4, m%4) = x(m+12, 0) * ((m_fCalcInMetres && m%4==3) ? 1000.0 : 1.0);
		m_Matrices[0](m/4, m%4) = x(m, 0) * ((m_fCalcInMetres && m%4==3) ? 1000.0 : 1.0);
	}
	m_Matrices[1](3,0) = m_Matrices[1](3,1) = m_Matrices[1](3,2) = 0;
	m_Matrices[0](3,0) = m_Matrices[0](3,1) = m_Matrices[0](3,2) = 0;

	m_Matrices[1](3,3) = 1.0f;
	m_Matrices[0](3,3) = 1.0f;

	m_fIsCalibrated = true;

	ComputeDeviations(err, true);

	return CAL_SUCCESS;
}



Matrix Math::crossprod (double x, double y, double z)
{
	Matrix m = Matrix::zeros(3,3);
	m(0, 1) = -z;
	m(0, 2) = y;
	m(1, 0) = z;
	m(1, 2) = -x;
	m(2, 0) = -y;
	m(2, 1) = x;
	return m;
}


void Math::rodrigues(const Matrix &mat, double* rod)  // rod of length 3
{
	int m = mat.size(0);
	int n = mat.size(1);
	assert(m == n);
	assert(m == 3);
	
//	double* rod = new double[3];

	// project the rotation matrix to SO(3);
	Matrix U(3,3);
	Matrix S(3,3);
	Matrix V(3,3);
	mat.SVD(U, S, V);

	Matrix UxVt(3,3);
	U.rightMultiply(V, UxVt, false, true);
	double trace = (UxVt.trace()-1)*0.5;  // max (3-1)*0.5 = 1 => sqrt(trace) is real
	double theta = acos(trace);  // only the real part, no imaginary, please
    
	if (sin(theta) >= 0.00001)
	{
//		double diagElem = -0.5/sqrt(1-trace*trace);
		double vth = 1.0/(2*sin(theta));
//		double dvthdtheta = -vth * cos(theta) / sin(theta);

		rod[0] = UxVt(2,1)-UxVt(1,2);
		rod[1] = UxVt(0,2)-UxVt(2,0);
		rod[2] = UxVt(1,0)-UxVt(0,1);

		for (int i=0; i<3; i++)
			rod[i] *= vth*theta;
	}
	else
	{
		if (trace > 0)  // norm(om) == 0
		{
			rod[0] = rod[1] = rod[2] = 0;
		}
		else // norm(om) == PI
		{
			// out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
			for (int i=0; i<3; i++)
			{
				rod[i] = sqrt((UxVt(i,i)+1)*0.5);
			}
			if (UxVt(0,1) < 0)
				rod[1] = -rod[1];
			if (UxVt(0,2) < 0)
				rod[2] = -rod[2];
			for (int i=0; i<3; i++)
			{
				rod[i] *= theta;
			}
		}
	}
}


CAL_ECODE Calibration::CCalibAlgorithmTsaiLenz::ComputeMatrices(int &numPoints, DoubleVector &err) 
{
	if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;

	std::ofstream f;
	f.open("/home/ernst/test.m");

	numPoints = m_Points[0].size();
	int K = (numPoints * (numPoints - 1)) / 2, k = 0;

	Matrix A(3*K,3), B(3*K,1), Tinv(4,4), Rinv(4,4), r(3,1), x(4,1);
	std::vector<Matrix> R(K, Matrix::eye(4)), T(K, Matrix::eye(4));

	DoubleVector q1(4), q2(4), q(3);

	for(int i=0; i<numPoints; i++) {
		for(int j = i+1; j<numPoints; j++) {
			Matrix::invert(m_Points[1][i], Rinv);
			Rinv.leftMultiply(m_Points[1][j], R[k]);     // Transformation from i-th to j-th gripper pose
			Math::MatToQuat(R[k], q1);                    // ... and the corresponding quaternion
			q1 = q1 * (-2);

			Matrix::invert(m_Points[0][i], Tinv);
			Tinv.leftMultiply(m_Points[0][j], T[k]);      // Transformation from i-th to j-th camera pose
			Math::MatToQuat(T[k], q2);                    // ... and the corresponding quaternion
			q2 = q2 * (-2);
			
			A(3*k, 0) = 0;                                // A(3*k+(0:2), 0:2) = skew(q1(1:3) + q2(1:3));
			A(3*k + 1, 1) = 0;
			A(3*k + 2, 2) = 0;
			A(3*k, 1) = -q1[3] - q2[3];
			A(3*k, 2) = q1[2] + q2[2];
			A(3*k + 1, 0) = q1[3] + q2[3];
			A(3*k + 1, 2) = -q1[1] - q2[1];
			A(3*k + 2, 0) = -q1[2] - q2[2];
			A(3*k + 2, 1) = q1[1] + q2[1];

			B(3*k, 0) = q2[1] - q1[1];                    // B(3*k+(0:2)) = q2(1:3) - q1(1:3);
			B(3*k + 1, 0) = q2[2] - q1[2];
			B(3*k + 2, 0) = q2[3] - q1[3];

			k++;
		}
	}

	A.QRSolve(r, B); // solve A*r = B

	// compute quaternion from this solution
	double s = sqrt(1 + r(0,0)*r(0,0) + r(1,0)*r(1,0) + r(2,0)*r(2,0));
	q1[1] = r(0,0) / s;
	q1[2] = r(1,0) / s;
	q1[3] = r(2,0) / s;
	q1[0] = sqrt(1 - q1[1]*q1[1] - q1[2]*q1[2] - q1[3]*q1[3]); //TODO: CHECK!!!
	// compute rotation matrix
	Math::QuatToMat(DoubleVector(3,0), q1, m_Matrices[0]);
	
	// calculate translational component
	k = 0;
	r.resize(4,1);
	for(int i=0; i<numPoints; i++) {
		for(int j = i+1; j<numPoints; j++) {
			// A(3*k+(0:2), 0:2) = R[k] - eye(3);
			for(int l=0; l<3; l++) {
				for(int m=0; m<3; m++) {
					A(3*k + l, m) = R[k](l,m) - ((m == l) ? 1 : 0);
				}
			}

			// B(3*k+(0:2)) = tmp(0:2,0:2)*T[k](0:2,3) - R[k](0:2,3);
			for(int m=0; m<4; m++)
				x(m,0) = T[k](m,3);
			m_Matrices[0].rightMultiply(x, r);

			for(int m=0; m<3; m++)
				B(3*k+m, 0) = r(m,0) - R[k](m,3);

			k++;
		}
	}
	r.resize(3,1);
	A.QRSolve(r, B);
	
	for(int m=0; m<3; m++)
		m_Matrices[0](m,3) = r(m,0);

	Matrix::invert(m_Matrices[0], m_Matrices[0]);

	f.close();

	// compute second matrix by averaging
	DoubleVector tempPos(3);
	DoubleVector tempQuat(4);
	DoubleVector avgPos(3);
	gsl_matrix* tempTensor = gsl_matrix_alloc(4,4);
	gsl_matrix* avgTensor = gsl_matrix_alloc(4,4);
	gsl_matrix_set_zero(avgTensor);
	Matrix fourthMat(4,4);
	Matrix temp(4,4);
	Matrix leftinv(4,4);
	for (int i = 0; i < numPoints; i++) //we have n-1 linearly independent relations between the views
	{
		Matrix::invert(m_Points[0][i], leftinv);
		Matrix::MatrixMultiply(leftinv, m_Matrices[0], temp);
		Matrix::MatrixMultiply(temp, m_Points[1][i], fourthMat);

		Math::MatToQuat(fourthMat, tempPos, tempQuat);
		Math::quaternion2Tensor(tempTensor, tempQuat);
		for (int k=0; k<3; k++)
			avgPos[k] += tempPos[k];
		gsl_matrix_add(avgTensor, tempTensor);
	}
	for (int i=0; i<3; i++)
		avgPos[i] /= numPoints;
	Math::tensor2Quaternion(avgTensor, tempQuat);
	Math::QuatToMat(avgPos, tempQuat, m_Matrices[1]);

	gsl_matrix_free(avgTensor);
	gsl_matrix_free(tempTensor);

	Matrix::invert(m_Matrices[0], m_Matrices[0]);
	Matrix::invert(m_Matrices[1], m_Matrices[1]);

	Matrix tmp = m_Matrices[0];
	m_Matrices[0] = m_Matrices[1];
	m_Matrices[1] = tmp;

	m_fIsCalibrated = true;
	this->ComputeDeviations(err);

/*  VOLKER'S ORIGINAL CODE: BIAS ON FIRST MEASUREMENT!!


	numPoints = m_Points[0].size();

	// get first matrices of both lists and invert them
	Matrix P1inv(4,4);
	Matrix::invert(m_Points[0][0], P1inv);
	Matrix Q1inv(4,4);
	Matrix::invert(m_Points[1][0], Q1inv);

	int num = m_Points[0].size();
	std::vector<Math::Matrix> pvec(num-1);
	std::vector<Math::Matrix> qvec(num-1);
	Matrix A(3*(num-1),3);
	Matrix b(3*(num-1),1);
	double rodP[3];
	double rodQ[3];
	for (int i = 1; i < num; i++) //we have n-1 linearly independent relations between the views
	{
		Matrix *Pi = &m_Points[0][i];
		Matrix *Qi = &m_Points[1][i];

		Matrix P(4,4);
		Matrix Q(4,4);
		Matrix::MatrixMultiply(*Pi, P1inv, P);
		Matrix::MatrixMultiply(*Qi, Q1inv, Q);
		pvec[i-1] = P;
		qvec[i-1] = Q;

		// Rodrigues formula
		Matrix Prot = P.copy(0,0,3,3);
		Matrix Qrot = Q.copy(0,0,3,3);
		rodrigues(Prot, rodP);
		rodrigues(Qrot, rodQ);

		b((i-1)*3  , 0) = rodQ[0]-rodP[0];
		b((i-1)*3+1, 0) = rodQ[1]-rodP[1];
		b((i-1)*3+2, 0) = rodQ[2]-rodP[2];

		for (int k=0; k<3; k++)
			rodP[k] += rodQ[k];

		Matrix crossproduct = crossprod(rodP[0], rodP[1], rodP[2]);
		for (int k=0; k<3; k++)
			for (int m=0; m<3; m++)
				A((i-1)*3+k, m) = crossproduct(k,m);
	}
	// compute rotation
	Matrix xr(3,1);
	A.linearSolve(xr, b);

	double norm = xr.norm2();
	double scaleFactor = 2/sqrt(1+norm*norm);
	double xx = xr(0,0) * scaleFactor;
	double yy = xr(1,0) * scaleFactor;
	double zz = xr(2,0) * scaleFactor;
	Matrix rcg2 = crossprod(xx, yy ,zz);
	rcg2.scale(sqrt(4-(xx*xx + yy*yy + zz*zz)));

	Matrix rcg3(3, 3);
	rcg3(0,0) = xx*xx;
	rcg3(0,1) = xx*yy;
	rcg3(0,2) = xx*zz;

	rcg3(1,0) = yy*xx;
	rcg3(1,1) = yy*yy;
	rcg3(1,2) = yy*zz;

	rcg3(2,0) = zz*xx;
	rcg3(2,1) = zz*yy;
	rcg3(2,2) = zz*zz;

	rcg2.add(rcg3);
	rcg2.scale(0.5);
	double rcgMainDiagonale = 1-(xx*xx + yy*yy + zz*zz)/2;
	for (int i=0; i<3; i++)
		rcg2(i,i) += rcgMainDiagonale;
	
	Matrix Qresult(3,1);
	Matrix Pt(3,1);
	Matrix Qt(3,1);

	// Translation		
	for (int i = 1; i < num; i++) //we have n-1 linearly independent relations between the views
	{			
		Matrix Pi = pvec[i-1];
		for (int k=0; k<3; k++)
			Pt(k,0) = Pi(k, 3);

		Matrix Qi = qvec[i-1];
		for (int k=0; k<3; k++)
			Qt(k,0) = Qi(k, 3);

		Matrix::MatrixMultiply(rcg2, Qt, Qresult);
		Qresult.sub(Pt);

		for (int k=0; k<3; k++)
		{
			b((i-1)*3+k,0) = Qresult(k,0);
			for (int m=0; m<3; m++)
			{
				if (m==k)
					A((i-1)*3+k, m) = Pi(k,m)-1;
				else
					A((i-1)*3+k, m) = Pi(k,m);
			}
		}
	}
	// compute translation
	Matrix xt(3,1);
	A.linearSolve(xt, b);

	for (int i=0; i<3; i++)
	{
		m_Matrices[0](i,3) = xt(i,0);
		for (int j=0; j<3; j++)
		{
			m_Matrices[0](i,j) = rcg2(i,j);
		}
	}
	
	// compute second matrix by averaging
	DoubleVector tempPos(3);
	DoubleVector tempQuat(4);
	DoubleVector avgPos(3);
	gsl_matrix* tempTensor = gsl_matrix_alloc(4,4);
	gsl_matrix* avgTensor = gsl_matrix_alloc(4,4);
	gsl_matrix_set_zero(avgTensor);
	Matrix fourthMat(4,4);
	Matrix temp(4,4);
	Matrix leftinv(4,4);
	for (int i = 0; i < num; i++) //we have n-1 linearly independent relations between the views
	{
		Matrix::invert(m_Points[0][i], leftinv);
		Matrix::MatrixMultiply(leftinv, m_Matrices[0], temp);
		Matrix::MatrixMultiply(temp, m_Points[1][i], fourthMat);

		Math::MatToQuat(fourthMat, tempPos, tempQuat);
		Math::quaternion2Tensor(tempTensor, tempQuat);
		for (int k=0; k<3; k++)
			avgPos[k] += tempPos[k];
		gsl_matrix_add(avgTensor, tempTensor);
	}
	for (int i=0; i<3; i++)
		avgPos[i] /= num;
	Math::tensor2Quaternion(avgTensor, tempQuat);
	Math::QuatToMat(avgPos, tempQuat, m_Matrices[1]);

	gsl_matrix_free(avgTensor);
	gsl_matrix_free(tempTensor);

	Matrix::invert(m_Matrices[0], m_Matrices[0]);
	Matrix::invert(m_Matrices[1], m_Matrices[1]);

	Matrix tmp = m_Matrices[0];
	m_Matrices[0] = m_Matrices[1];
	m_Matrices[1] = tmp;

	m_fIsCalibrated = true;
	this->ComputeDeviations(err); */

	return CAL_SUCCESS;
}


CAL_ECODE Calibration::CCalibAlgorithmDualQuaternion::ComputeMatrices(int &numPoints, DoubleVector &err) {
	if(!m_fPointsSet)
		return CAL_NOT_INITIALISED;

	// TODO!!

	return CAL_ALG_DOESNT_EXIST;
}


void Math::QuatToMat(const DoubleVector &pos, const DoubleVector &quat, Matrix &M) {
	if(M.size(0) != 4 || M.size(1) != 4)
		M.resize(4,4);

	// assign translation offset
	M(3,3) = 1; M(0,3) = pos[0]; M(1,3) = pos[1]; M(2,3) = pos[2];

	// compute rotational part
	M(0,0) = quat[0]*quat[0]+quat[1]*quat[1]-quat[2]*quat[2]-quat[3]*quat[3];
	M(1,0) = 2*(quat[1]*quat[2]+quat[0]*quat[3]);
	M(2,0) = 2*(quat[1]*quat[3]-quat[0]*quat[2]);
	M(3,0) = 0;
	//
	M(0,1) = 2*(quat[1]*quat[2]-quat[0]*quat[3]);
	M(1,1) = quat[0]*quat[0]-quat[1]*quat[1]+quat[2]*quat[2]-quat[3]*quat[3];
	M(2,1) = 2*(quat[2]*quat[3]+quat[0]*quat[1]);
	M(3,1) = 0;
	//
	M(0,2) = 2*(quat[1]*quat[3]+quat[0]*quat[2]);
	M(1,2) = 2*(quat[2]*quat[3]-quat[0]*quat[1]);
	M(2,2) = quat[0]*quat[0]-quat[1]*quat[1]-quat[2]*quat[2]+quat[3]*quat[3];
	M(3,2) = 0;
}

void Math::MatToQuat(const Matrix &M, DoubleVector &pos, DoubleVector &quat) {
	// extract position vector
	pos.resize(3);
	pos[0] = M(0,3); pos[1] = M(1,3); pos[2] = M(2,3);

	// extract quaternion
	quat.resize(4);

	// determine matrix trace
	double T = M(0,0) + M(1,1) + M(2,2) + 1;
	double s;

	if(T > 1e-6) {
		s = sqrt(T) * 2;
		quat[1] = ( M(1,2) - M(2,1) ) / s;
		quat[2] = ( M(2,0) - M(0,2) ) / s;
		quat[3] = ( M(0,1) - M(1,0) ) / s;
		quat[0] = -0.25 * s;
	} else {
		if ( M(0,0) > M(1,1) && M(0,0) > M(2,2) )  {	// Column 0: 
			s  = sqrt( 1.0 + M(0,0) - M(1,1) - M(2,2) ) * 2;
			quat[1] = 0.25 * s;
			quat[2] = (M(0,1) + M(1,0) ) / s;
			quat[3] = (M(2,0) + M(0,2) ) / s;
			quat[0] = -(M(1,2) - M(2,1) ) / s;
		} else if ( M(1,1) > M(2,2) ) {			// Column 1: 
			s  = sqrt( 1.0 + M(1,1) - M(0,0) - M(2,2) ) * 2;
			quat[1] = (M(0,1) + M(1,0) ) / s;
			quat[2] = 0.25 * s;
			quat[3] = (M(1,2) + M(2,1) ) / s;
			quat[0] = -(M(2,0) - M(0,2) ) / s;
		} else {						// Column 2:
			s  = sqrt( 1.0 + M(2,2) - M(0,0) - M(1,1) ) * 2;
			quat[1] = (M(2,0) + M(0,2) ) / s;
			quat[2] = (M(1,2) + M(2,1) ) / s;
			quat[3] = 0.25 * s;
			quat[0] = -(M(0,1) - M(1,0) ) / s;
		}
	}

	return;

	/* // find squares of the entries of the quaternion
	quat[0] = 1 + M(0,0) + M(1,1) + M(2,2);
	quat[1] = 1 + M(0,0) - M(1,1) - M(2,2);
	quat[2] = 1 - M(0,0) + M(1,1) - M(2,2);
	quat[3] = 1 - M(0,0) - M(1,1) + M(2,2);

	// take maximum
	int ind = 0;
	for(int i=1; i<4; i++) {
		if(quat[i] > quat[ind]) {
			ind = i;
		}
	}

	// extract first component
	quat[ind] = sqrt(quat[ind]) / 2;

	// extract remaining components
	switch(ind) {
		case 0:
			quat[1] = 1/(4*quat[ind]) * (M(2,1) - M(1,2));
			quat[2] = 1/(4*quat[ind]) * (M(0,2) - M(2,0));
			quat[3] = 1/(4*quat[ind]) * (M(1,0) - M(0,1));
			break;
		case 1:
			quat[0] = 1/(4*quat[ind]) * (M(2,1) - M(1,2));
			quat[2] = 1/(4*quat[ind]) * (M(1,0) + M(0,1));
			quat[3] = 1/(4*quat[ind]) * (M(0,2) + M(2,0));
			break;
		case 2:
			quat[0] = 1/(4*quat[ind]) * (M(0,2) - M(2,0));
			quat[1] = 1/(4*quat[ind]) * (M(1,0) + M(0,1));
			quat[3] = 1/(4*quat[ind]) * (M(2,1) + M(1,2));
			break;
		case 3:
			quat[0] = 1/(4*quat[ind]) * (M(1,0) - M(0,1));
			quat[1] = 1/(4*quat[ind]) * (M(0,2) + M(2,0));
			quat[2] = 1/(4*quat[ind]) * (M(2,1) + M(1,2));
	} */
}

void Math::MatToQuat( const Matrix &M, DoubleVector &quat )
{
	DoubleVector pos(3,0);
	MatToQuat(M, pos, quat);
}

CAL_ECODE Math::ConvertRotation(Calibration::tDataObject& obj, Calibration::kROTTYPE dest) {
	if(obj.RotationType == dest)
		return Calibration::CAL_SUCCESS;

	// currently supported conversions: Matrix <-> Quaternion
	if(obj.RotationType == Calibration::rotMATRIX && dest == Calibration::rotQUATERNION) {
		MatToQuat(obj.matrix, obj.trans, obj.quat);
	} else if(obj.RotationType == Calibration::rotQUATERNION && dest == Calibration::rotMATRIX) {
		QuatToMat(obj.trans, obj.quat, obj.matrix);
	} else {
		return Calibration::CAL_CONVERSION_NOT_SUPPORTED;
	}

	obj.RotationType = dest;
	return Calibration::CAL_SUCCESS;
}

Calibration::CCalibAlgorithm *Calibration::CCalibAlgorithmBookKeeper::Get(kCalibrationAlgorithms which) {
	// check if algorithm exists
	for(size_t i=0; i<m_pAlgorithms.size(); i++)
		if(m_pAlgorithms[i]->GetType() == which)
			return m_pAlgorithms[i];

	// if not, create new algorithm if possible
	switch (which)
	{
	case kCALIB_QR15:
		m_pAlgorithms.push_back(new CCalibAlgorithmQR15);
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	case kCALIB_QR24:
		m_pAlgorithms.push_back(new CCalibAlgorithmQR24);
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	case kCALIB_QR24M:
		m_pAlgorithms.push_back(new CCalibAlgorithmQR24(true));
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	case kCALIB_DUALQUATERNION:
		m_pAlgorithms.push_back(new CCalibAlgorithmDualQuaternion);
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	case kCALIB_TSAILENZ:
		m_pAlgorithms.push_back(new CCalibAlgorithmTsaiLenz);
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	case kCALIB_FORCETORQUE:
		m_pAlgorithms.push_back(new CCalibAlgorithmForceTorque);		
		return m_pAlgorithms[m_pAlgorithms.size()-1];
		break;
	default:
		return NULL;
	}
}

Calibration::CCalibAlgorithmBookKeeper::~CCalibAlgorithmBookKeeper() {
	for (size_t i=0; i<m_pAlgorithms.size(); i++)
	{
		if (m_pAlgorithms[i] != NULL)
		{
			delete m_pAlgorithms[i];
			m_pAlgorithms[i] = NULL;
		}
	}
}



void Math::quaternion2Tensor(gsl_matrix* mat, DoubleVector &q)
{
	for (int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			gsl_matrix_set(mat, i, j, q[i]*q[j]);
		}
	}
}


void Math::tensor2Quaternion(gsl_matrix* matrixN, DoubleVector &q)
{
	gsl_vector* eigenvalues = gsl_vector_calloc(4);
	gsl_matrix* eigenvectors = gsl_matrix_calloc(4,4);
	gsl_eigen_symmv_workspace* work = gsl_eigen_symmv_alloc(4);
	gsl_vector*	eigenvector = gsl_vector_calloc(4);
	gsl_eigen_symmv (matrixN, eigenvalues, eigenvectors, work);
	// sort in descending order in numerical value
	gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_VAL_DESC);
	gsl_matrix_get_col(eigenvector, eigenvectors, 0);
	// derive unit vector
	gsl_vector_scale(eigenvector,(1.0/gsl_blas_dnrm2(eigenvector)));

   	gsl_vector_free(eigenvalues);
  	gsl_matrix_free(eigenvectors);
	gsl_eigen_symmv_free(work);

	q[0] = gsl_vector_get(eigenvector, 0);  // 0
	q[1] = gsl_vector_get(eigenvector, 1);  // x
	q[2] = gsl_vector_get(eigenvector, 2);  // y
	q[3] = gsl_vector_get(eigenvector, 3);  // z
	gsl_vector_free(eigenvector);
}

