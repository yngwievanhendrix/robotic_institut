/**

 \file wavelets.cpp
 (c) 2006 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        A trous wavelet decomposition and reconstruction<br>
 COMPONENT:      wavelets.cpp, wavelets.h
 \brief Wavelet decomposition
 
 **************************************************************************/

/***************************************************************************
 *
 * Required includes 
 *
 **************************************************************************/

#include "math/include/wavelets.h"
#include <iostream>

/***************************************************************************
 * 
 *  Function bodies
 *
 **************************************************************************/

using namespace Math;

Math::CWaveletSkeleton::CWaveletSkeleton(void) {
	m_fIsInitialised = false;
}

Math::CWaveletSkeleton::CWaveletSkeleton(DoubleVector x, DoubleVector y, int NumScales) {
	m_fIsInitialised = Init(x, y, NumScales);
}

void Math::CWaveletSkeleton::GetW(DoubleMatrix &W) {
	W = m_W;
}

void Math::CWaveletSkeleton::GetC(DoubleVector &c) {
	c = m_c;
}

bool Math::CWaveletSkeleton::Init(DoubleVector x, DoubleVector y, int NumScales) {
	//if(pow((double)2,(NumScales-1)) > y.size()) {
	//	return false;
	//}

	// prepare storage
	m_y = y;
	m_x = x;
	m_c = y;
	m_W.resize(NumScales, DoubleVector(y.size(), 0));
	m_AllCoeffsC.resize(0, DoubleVector(y.size(), 0));
	m_nNumScales = NumScales;
	DoubleVector cOld;

/*	// open output file
	FILE *file;
	file = fopen("wavskel.txt", "w");

	for(unsigned int t=0;t<m_y.size();t++) {
		fprintf(file, "%.12lf ", m_x[t]);
	}
	fprintf(file, "\n");

	for(unsigned int t=0;t<m_y.size();t++) {
		fprintf(file, "%.12lf ", m_y[t]);
	}
	fprintf(file, "\n"); */

	// compute coefficients
	int k = 0;
	for(int j=0;j<=NumScales - 1;j++) {
		m_AllCoeffsC.push_back(m_c);
		cOld = m_c;
		for(unsigned int t=0;t<m_y.size();t++) {
			k = t-(int)floor(pow((double)2, j));
			m_c[t] = 0.5*(cOld[MAX(0,k)] + cOld[t]);
			m_W[j][t] = cOld[t] - m_c[t];
//			fprintf(file, "%.12lf ", m_W[j][t]);
		}
//		fprintf(file, "\n");
	}
	m_AllCoeffsC.push_back(m_c);

//	for(unsigned int t=0;t<m_y.size();t++) {
//		fprintf(file, "%.12lf ", m_c[t]);
//	}

//	fclose(file);

	m_fIsInitialised = true;
	return true;
}

bool Math::CWaveletSkeleton::GetAppx(int FirstScale, int LastScale, double &appx, double &trend, int &Delay) {
	// check for validity
	if(!m_fIsInitialised)
		return false;
	if(LastScale == -1)
		LastScale = m_nNumScales - 1;
	if(FirstScale < 0 || LastScale >= m_nNumScales || LastScale < FirstScale)
		return false;

	// compute delay
	Delay = (int)floor(pow(2.0, FirstScale-1)) - 1;
	Delay = MAX(0, Delay);

	// compute trend
	trend = 0;
	size_t last = m_x.size() - 1;
	for(int i=LastScale+1; i<m_nNumScales; i++) 
		trend += m_W[i][last];
	trend += m_c[last];

	// compute approximation
	appx = 0;
	for(int i=FirstScale; i<=LastScale; i++) 
		appx += m_W[i][last];

	// done
	return true;
}

bool Math::CWaveletSkeleton::GetAppx(int FirstScale, int LastScale, DoubleVector &appx, DoubleVector &trend, int &Delay) {
	// check for validity
	if(!m_fIsInitialised)
		return false;
	if(LastScale == -1)
		LastScale = m_nNumScales - 1;
	if(FirstScale < 0 || LastScale >= m_nNumScales || LastScale < FirstScale)
		return false;

	// compute delay
	Delay = (int)floor(pow(2.0, FirstScale-1)) - 1;
	Delay = MAX(0, Delay);

	// compute trend
	trend = DoubleVector(m_y.size(), 0);
	for(int i=LastScale+1; i<m_nNumScales; i++) 
		trend = trend + m_W[i];
	trend = trend + m_c;

	// compute approximation
	appx = DoubleVector(m_y.size(), 0);
	for(int i=FirstScale; i<=LastScale; i++) 
		appx = appx + m_W[i];

	// done
	return true;
}

bool Math::CWaveletSkeleton::GetEnergies(DoubleVector &energies, int start, int end) {
	if(start < 0 || start >= end || end > (int)m_y.size())
		return false;

	int i,j;

	// create local copies of y, c and W
	DoubleVector y(end-start + 1), c(end-start + 1);
	DoubleMatrix W(m_nNumScales, DoubleVector(end-start + 1, 0));
	for(i=0; i<end-start+1; i++) {
		y[i] = m_y[start + i];
		c[i] = m_c[start + i];
		for(j = 0; j<m_nNumScales; j++) {
			W[j][i] = m_W[j][start + i];
		}
	}

	// compute energies
	energies = DoubleVector(m_nNumScales,0);
	double all = (y-c)*(y-c);
	for(unsigned int i=0; i<energies.size(); i++) {
		energies[i] = (W[i]*W[i])/all;
	}

	return true;
}

void Math::CWaveletSkeleton::GetEnergies(DoubleVector &energies) {
	energies = DoubleVector(m_nNumScales,0);
	double all = (m_y-m_c)*(m_y-m_c);
	for(unsigned int i=0; i<energies.size(); i++) {
		energies[i] = (m_W[i]*m_W[i])/all;
	}
}

bool Math::CWaveletSkeleton::Push(double x, double y) {
	if(!m_fIsInitialised)
		return false;
	
	// push new values
	m_x.push_back(x);
	m_y.push_back(y);

	// compute next coefficients
	int k = 0;
	int t = (int)m_y.size() - 1;

	m_AllCoeffsC[0].push_back(y);
	for(int j=1;j<=m_nNumScales;j++) {
		k = t - (int)floor(pow((double)2, j));
		m_AllCoeffsC[j].push_back(0.5*(m_AllCoeffsC[j-1][MAX(0,k)] + m_AllCoeffsC[j-1][t]));
		m_W[j-1].push_back(m_AllCoeffsC[j-1][t] - m_AllCoeffsC[j][t]);
	}
	m_c.push_back(m_AllCoeffsC[m_nNumScales][m_AllCoeffsC[m_nNumScales].size() - 1]);

	return true;
}

int Math::CWaveletSkeleton::GetSize(void) {
	return (int)m_x.size();
}
