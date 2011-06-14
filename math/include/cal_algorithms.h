#ifndef __cal_algorithms_h
#define __cal_algorithms_h

#include "math/include/vectorops.h"
#include "math/include/matrix.h"
#include "math/include/dataobject.h"

#include "math/include/cal_reply_codes.h"
#include "math/include/trc_reply_codes.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

namespace Calibration {

	using namespace Math; 
	using namespace Tracking;
	
	//! Available algorithms
	typedef enum {
		kCALIB_QR24,
		kCALIB_QR24M,
		kCALIB_QR15,
		kCALIB_TSAILENZ,
		kCALIB_DUALQUATERNION,
		kCALIB_FORCETORQUE,
		kCALIB_ALG_COUNT,
		kCALIB_ALG_MAX = kCALIB_TSAILENZ,
		// kCALIB_ALG_MAX = kCALIB_DUALQUATERNION
		kCALIB_ALG_MIN = kCALIB_QR24
		// kCALIB_ALG_MIN = kCALIB_QR24
	} kCalibrationAlgorithms;

	//! Data type accepted by algorithm
	typedef enum {
		kCALIB_DTYPE_POSROT = 1,
		kCALIB_DTYPE_FORCETORQUE = 2,
		kCALIB_DTYPE_ALL = 3,
		kCALIB_DTYPE_UNKNOWN = 2048
	} kCalibrationDataTypeAccepted;

	kCalibrationDataTypeAccepted ConvertDataType(kROTTYPE type);

	//! Iterate through the algorithms (prefix increment)
	inline kCalibrationAlgorithms &operator++(kCalibrationAlgorithms &a) {
		return a = kCalibrationAlgorithms(a + 1);
	}

	//! Iterate through the algorithms (postfix increment)
	inline kCalibrationAlgorithms operator++(kCalibrationAlgorithms &a, int) {
		kCalibrationAlgorithms r(a);
		a = kCalibrationAlgorithms(a + 1);
		return r;
	}

	/** 
	Abstract calibration algorithm class
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-19
	*/
	class CCalibAlgorithm {
	protected:
		//! calibration matrices
		Matrix m_Matrices[2];
		//! type of algorithm
		kCalibrationAlgorithms m_Type;
		//! name of the algorithm
		std::string m_szName;
		//! points are set
		bool m_fPointsSet;
		//! calibration state
		bool m_fIsCalibrated;
		//! Calibration errors
		DoubleVector m_Errors;
		//! points
		std::vector<Matrix> m_Points[2];
		//! other data types
		std::vector<DoubleVector> m_Data[2];
		//! data type accepted
		kCalibrationDataTypeAccepted m_TypeAccepted[2];
		//! default value for orthonormalisation
		bool m_fDefaultOrtho;

		

	public:
//**********************************************//
		//			!only for force-torque:				//
		//**********************************************//
		//!sensor specific:
		//!sensor bias values
		double m_bias[6];
		//!rotation matrix from sensor to endeffector
		Matrix m_rot_sensor_effector;
		//
		//!tool specific:
		//distance to the centroid in sensor coodinates
		double m_centroid_sensor[3];
		//! force magnitude
		double m_magnitude;
		//!tool matrix: matrix from endeffector to tool
		Matrix m_rot_tool_effector;
		//distance to the centroid in tool coodinates
		double m_centroid_tool[3];
		//**********************************************//

		CCalibAlgorithm(kCalibrationAlgorithms type, std::string name, kCalibrationDataTypeAccepted typeAcceptedLeft = kCALIB_DTYPE_POSROT, kCalibrationDataTypeAccepted typeAcceptedRight = kCALIB_DTYPE_POSROT);

		CAL_ECODE SetPoints(const std::vector<tDataObject>& system_1, const std::vector<tDataObject>& system_2);
		CAL_ECODE SetPoints(const std::vector<Matrix>& system_1, const std::vector<Matrix>& system_2);
		CAL_ECODE SetPoints(const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos);
		CAL_ECODE SetPoints(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<Matrix>& system_2);
		CAL_ECODE SetPoints(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos);

		//! Compute the calibration errors
		CAL_ECODE ComputeDeviations(DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<tDataObject>& system_1, const std::vector<tDataObject>& system_2, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<DoubleVector>& system_1_rot, const std::vector<DoubleVector>& system_1_pos, const std::vector<DoubleVector>& system_2_rot, const std::vector<DoubleVector>& system_2_pos, DoubleVector &err, bool ortho = false);

		CAL_ECODE ComputeDeviations(const std::vector<Matrix>& system_1, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<DoubleVector>& system_1, const std::vector<Matrix>& system_2, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<Matrix>& system_1, const std::vector<DoubleVector>& system_2, DoubleVector &err, bool ortho = false);
		CAL_ECODE ComputeDeviations(const std::vector<DoubleVector>& system_1, const std::vector<DoubleVector>& system_2, DoubleVector &err, bool ortho = false);

		//! Get the calibration matrices
		CAL_ECODE GetMatrices(Matrix &system, Matrix &offset);

		//! Compute the calibration matrices
		virtual CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err) = 0;

		//! get the algorithm's type
		inline kCalibrationAlgorithms GetType(void) {
			return m_Type;
		}

		//! get default value for orthonormalisation
		inline bool GetOrtho(void) {
			return m_fDefaultOrtho;
		}

		//! get the algorithm's accepted data type
		inline void GetAcceptedDataTypes(kCalibrationDataTypeAccepted &left, kCalibrationDataTypeAccepted &right) {
			left = m_TypeAccepted[0];
			right = m_TypeAccepted[1];
		}

		//! get the algorithm's name
		inline std::string GetName(void) {
			return m_szName;
		}
	};

	/**
	QR 15 calibration algorithm
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-19
	*/
	class CCalibAlgorithmQR15 : public CCalibAlgorithm {
	public:
		inline CCalibAlgorithmQR15(void) : CCalibAlgorithm(kCALIB_QR15, "QR15") { m_fDefaultOrtho = true; }
		CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err);
	};

	/**
	QR 24 calibration algorithm
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-19
	*/
	class CCalibAlgorithmQR24 : public CCalibAlgorithm {
	private:
		bool m_fCalcInMetres;
	public:
		inline CCalibAlgorithmQR24(bool metres = false) : CCalibAlgorithm(kCALIB_QR24, metres ? "QR24M" : "QR24") { m_fCalcInMetres = metres; m_fDefaultOrtho = true; }
		inline void SetCalcMode(bool metres) { m_fCalcInMetres = metres; }
		CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err);
	};

	/**
	Tsai-Lenz calibration algorithm
	@author Volker Martens (martens@rob.uni-luebeck.de) @date 2008-08-20
	*/
	class CCalibAlgorithmTsaiLenz : public CCalibAlgorithm {
	public:
		inline CCalibAlgorithmTsaiLenz(void) : CCalibAlgorithm(kCALIB_TSAILENZ, "TsaiLenz") { }
		CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err);
	};

	/**
	Dual Quaternion calibration algorithm
	@author Volker Martens (martens@rob.uni-luebeck.de) @date 2008-08-20
	*/
	class CCalibAlgorithmDualQuaternion : public CCalibAlgorithm {
	public:
		inline CCalibAlgorithmDualQuaternion(void) : CCalibAlgorithm(kCALIB_DUALQUATERNION, "DUALQUATERNION") { }
		CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err);
	};

	/**
	Force-Torque sensor calibration algorithm
	@author Lars Richter (richter@rob.uni-luebeck.de) @date 2010-02-02
	*/
	class CCalibAlgorithmForceTorque : public CCalibAlgorithm {
	public:
		inline CCalibAlgorithmForceTorque(void) : CCalibAlgorithm(kCALIB_FORCETORQUE, "FORCETORQUE", kCALIB_DTYPE_FORCETORQUE, kCALIB_DTYPE_POSROT) { }
		CAL_ECODE ComputeMatrices(int &numPoints, DoubleVector &err);
	};

	/**
	Book keeping for calibration algorithms
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-26
	*/
	class CCalibAlgorithmBookKeeper {
	private:
		//! Make copy constructor private to forbid copies
		inline CCalibAlgorithmBookKeeper(const CCalibAlgorithmBookKeeper &c) {};
		//! Make assignment operator private to forbid copies
		inline CCalibAlgorithmBookKeeper& operator=(const CCalibAlgorithmBookKeeper &c) { return *this; };
		//! Storage for created algorithms (will be deleted on destruction)
		std::vector<CCalibAlgorithm *> m_pAlgorithms;
	public:
		/**
		Get algorithm of specified type
		@param	which	type of the algorithm
		@returns	pointer to the new algorithm (will be deleted on destruction of this class!)
		@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-26
		*/
		CCalibAlgorithm *Get(kCalibrationAlgorithms which);
		inline CCalibAlgorithmBookKeeper(void) {}
		~CCalibAlgorithmBookKeeper(void);
	};
}

namespace Math 
{
	/**
	Generate transformation matrix from translation and quaternion
	@param[in]	pos	a DoubleVector with the position part of the matrix
	@param[in]	quat	a DoubleVector with the quaternion of the matrix
	@param[out]	M	a matrix containing the transformation described by pos and quat
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-03-14
	*/
	void QuatToMat(const DoubleVector &pos, const DoubleVector &quat, Matrix &M);
	/**
	Generate translation and quaternion from transformation matrix
	@param[in]	M	a matrix containing the transformation described by pos and quat
	@param[out]	pos	a DoubleVector with the position part of the matrix
	@param[out]	quat	a DoubleVector with the quaternion of the matrix
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2007-03-14
	*/
	void MatToQuat(const Matrix &M, DoubleVector &pos, DoubleVector &quat);
	void MatToQuat(const Matrix &M, DoubleVector &quat);
	/**
	Convert rotation part of a data object
	@param[out]	obj	the data object to convert
	@param[in]	dest	the destination rotation type
	@author Floris Ernst (ernst@rob.uni-luebeck.de) @date 2008-08-22
	*/
	Calibration::CAL_ECODE ConvertRotation(Calibration::tDataObject& obj, Calibration::kROTTYPE dest);

	Matrix crossprod (double x, double y, double z);

	void rodrigues(const Matrix &mat, double* rod);

	void quaternion2Tensor(gsl_matrix* mat, DoubleVector &q);

	void tensor2Quaternion(gsl_matrix* matrixN, DoubleVector &q);

}

#endif
