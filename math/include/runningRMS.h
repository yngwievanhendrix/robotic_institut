#ifndef __runningRMS_h
#define __runningRMS_h

#include "math/include/vectorops.h"

namespace Math {
	class CRunningRMS {
	public:
		CRunningRMS(void);
		CRunningRMS(const DoubleVector &v);
		void Init(const DoubleVector &v);
		void Init(int n);
		double RMS(void);
		double Push(double val);
	private:
		int m_nSize;
		int m_nCurrent;
		DoubleVector m_vStorage;
		double m_dRMS_squared;
	};
}

#endif