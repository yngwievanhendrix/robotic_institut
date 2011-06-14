#ifndef __dataobject_h
#define __dataobject_h

#include "math/include/vectorops.h"
#include "math/include/matrix.h"

namespace Math {
	typedef enum {
		rotMATRIX = 0,
		rotQUATERNION,
		rotAXISANGLE,
		rotYPR,
		rotEULER,
		rotFORCETORQUE,
		rotNONE,
		rotUNKNOWN
	} kROTTYPE;

	typedef struct _tDataObject {
		inline _tDataObject(void) {
			RotationType = rotMATRIX;
			axis.resize(3,0);
			angles.resize(3,0);
			angle = 0;
			trans.resize(3,0);
			quat.resize(4,0); quat[0] = 1;
			matrix = Matrix::eye(4);
			forcetorque.resize(6,0);
			vis = true;
			quality = -1;
			ts = -1;
		}
		inline _tDataObject(kROTTYPE rottype) {
			RotationType = rottype;
			axis.resize(3,0);
			angles.resize(3,0);
			angle = 0;
			trans.resize(3,0);
			quat.resize(4,0); quat[0] = 1;
			matrix = Matrix::eye(4);
			forcetorque.resize(6,0);
			vis = true;
			quality = -1;
			ts = -1;
		}
		//! Rotation type of the stored data
		kROTTYPE RotationType;
		//! Axis
		DoubleVector axis;
		//! Angle (axis-angle)
		double angle;
		//! Angles (Euler, YPR)
		DoubleVector angles;
		//! Translation
		DoubleVector trans;
		//! Quaternion
		DoubleVector quat;
		//! Matrix
		Matrix matrix;
		//! Force-Torque vector
		DoubleVector forcetorque;
		//! visibility
		bool vis;
		//! quality
		double quality;
		//! timestamp
		double ts;
	} tDataObject;
}

#endif

