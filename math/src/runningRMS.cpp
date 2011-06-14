#include "math/include/runningRMS.h"

Math::CRunningRMS::CRunningRMS( void )
{
	m_nSize = 0;
	m_vStorage.resize(0, 0);
}

Math::CRunningRMS::CRunningRMS( const DoubleVector &v )
{
	Init(v);
}

void Math::CRunningRMS::Init( const DoubleVector &v )
{
	m_nSize = v.size();
	m_nCurrent = v.size() - 1;
	m_vStorage.resize(v.size());
	m_dRMS_squared = 0;
	for(size_t i=0; i<v.size(); i++) {
		m_vStorage[i] = v[i] * v[i];
		m_dRMS_squared += v[i];
	}
}

void Math::CRunningRMS::Init( int n )
{
	Init(DoubleVector(n, 0));
}

double Math::CRunningRMS::RMS( void )
{
	return sqrt(m_dRMS_squared / (double)m_nSize);
}

double Math::CRunningRMS::Push( double val )
{
	m_nCurrent++;
	while(m_nCurrent >= m_nSize)
		m_nCurrent -= m_nSize;
	m_dRMS_squared -= m_vStorage[m_nCurrent];
	m_vStorage[m_nCurrent] = val * val;
	m_dRMS_squared += m_vStorage[m_nCurrent];
	return RMS();
}