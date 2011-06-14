#ifndef __distortion_h
#define __distortion_h

#include "math/include/dataobject.h"

namespace Math {
	bool ComputeDistortion(std::vector<Matrix> left, std::vector<DoubleVector> right, std::vector<DoubleVector> &errVec, DoubleVector &errors, Matrix &matrix);
	DoubleVector ProjectToPlane(const DoubleVector n, const DoubleVector v, const DoubleVector p);
	DoubleVector InterpolateOnPlane(const DoubleVector p[4], const DoubleVector x, const DoubleVector v[4]);
	DoubleVector InterpolateOnLine(const DoubleVector p, const DoubleVector q, const DoubleVector x, DoubleVector v[2]);
	DoubleVector ProjectToLine(const DoubleVector p, const DoubleVector q, const DoubleVector x);

	class CDCube : public CCube {
	private:
		static int m_Faces[4][6];
		DoubleVector m_Points[8];
		DoubleVector m_Values[8];
	public:
		CDCube(DoubleVector p, DoubleVector q);
		CDCube(DoubleVector v[8]);
		CDCube(void);
		void Init(DoubleVector p, DoubleVector q);
		void SetValues(DoubleVector v[8]);
		void SetValue(DoubleVector v, int i);

		DoubleVector Trilinear(const DoubleVector p, double tolerance = 0);
		DoubleVector InterpolateOnFace(const int f, const DoubleVector p);
		DoubleVector ProjectToFace(const int f, const DoubleVector p);
		bool Inside(const DoubleVector p, double tolerance = 0);
	};
}

#endif
