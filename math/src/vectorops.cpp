/**

 \file vectorops.cpp
 (c) 2006 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        Vector operations<br>
 COMPONENT:      vectorops.cpp, vectorops.h
 \brief Vector and data manipulation methods
 
 **************************************************************************/

/***************************************************************************
 *
 * Required includes 
 *
 **************************************************************************/

#include "math/include/vectorops.h"

/***************************************************************************
 * 
 *  Function bodies
 *
 **************************************************************************/

using namespace Math;

Math::CCube::CCube(void) {
	Init(DoubleVector(3,0), DoubleVector(3,1));
}

Math::CCube::CCube(const DoubleVector &x, const DoubleVector &y) {
	Init(x,y);
}

void Math::CCube::Init(const DoubleVector &x, const DoubleVector &y) {
	// check dimensions
	if(x.size() != 3 || y.size() != 3) {
		Init(DoubleVector(3,0), DoubleVector(3,1));
		return;
	}

	// initialise
	m_p.resize(3); m_q.resize(3);
	for(int i=0; i<3; i++) {
		m_p[i] = MIN(x[i], y[i]);
		m_q[i] = MAX(x[i], y[i]);
	}
}

void Math::CCube::Clip(DoubleVector &x) {
	if(x.size() != 3)
		return;

	double v[3] = { x[0], x[1], x[2] };
	Clip(v);
	x[0] = v[0]; x[1] = v[1]; x[2] = v[2];
}

void Math::CCube::Clip(double &x, double &y, double &z) {
	double v[3] = { x, y, z };
	Clip(v);
	x = v[0]; y = v[1]; z = v[2];
}

void Math::CCube::Clip(double* x) {
	// perform clipping in three dimensions
	for(int i=0; i<3; i++) {
		x[i] = MIN(m_q[i], MAX(m_p[i], x[i]));
	}
}

DoubleVector Math::subview(const DoubleVector& v, int i, int j) { 
	if(j <= 0)
		j += v.size() - 1;
	
	if(i < 0) i = 0;
	if(j >= (int)v.size()) j = (int)v.size()-1;
	
	assert(i <= j);

	DoubleVector w;
	w.resize(j-i+1,0);
	for(int k=i; k<=j; k++)
		w[k-i] = v[k];
	return w;
}

double Math::stddev(const DoubleVector& v, size_t start, size_t end) { 
	end = end == -1 ? v.size() : end;
	double m = Mean(v, start, end);
	DoubleVector v_ = subview(v, start, end);
	return sqrt((1.0/((double)(end-start-1)))*((v_-m)*(v_-m)));
}
DoubleVector Math::operator + (const DoubleVector& x, const DoubleVector& y) {
	DoubleVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = x[n] + y[n];
	}
	return res;
}

DoubleVector Math::operator - (const DoubleVector& x, const DoubleVector& y) {
	DoubleVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = x[n] - y[n];
	}
	return res;
}

DoubleVector Math::operator + (const DoubleVector& x, const double& s) {
	DoubleVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = x[n] + s;
	}
	return res;
}

DoubleVector Math::operator + (const double& s, const DoubleVector& x) {
	return x + s;
}

DoubleVector Math::operator - (const DoubleVector& x, const double& s) {
	return x + (-s);
}

DoubleVector Math::operator * (const DoubleVector& x, const double& r) {
	DoubleVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = r*x[n];
	}
	return res;
}

DoubleVector Math::operator * (const double& r, const DoubleVector& x) {
	return x*r;
}

LongVector Math::operator + (const LongVector& x, const LongVector& y) {
	LongVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = x[n] + y[n];
	}
	return res;
}

LongVector Math::operator - (const LongVector& x, const LongVector& y) {
	LongVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = x[n] - y[n];
	}
	return res;
}

LongVector Math::operator * (const LongVector& x, const int& r) {
	LongVector res(x.size(), 0);
	for(unsigned int n=0;n<x.size();n++) {
		res[n] = r*x[n];
	}
	return res;
}

LongVector Math::operator * (const int& r, const LongVector& x) {
	return x*r;
}

double Math::operator * (const DoubleVector& x, const DoubleVector& y) {
	double res = 0;
	for(unsigned int n=0;n<x.size();n++) {
		res += x[n]*y[n];
	}
	return res;
}

DoubleVector Math::ElementProduct(const DoubleVector& x, const DoubleVector& y) {
	DoubleVector s = x;
	for(size_t n=0;n<x.size();n++) {
		s[n] *= y[n];
	}
	return s;
}

int Math::operator * (const LongVector& x, const LongVector& y) {
	int res = 0;
	for(unsigned int n=0;n<x.size();n++) {
		res += x[n]*y[n];
	}
	return res;
}

double Math::Round(double x, int acc) {
	double up = ceil(x*acc), down = floor(x*acc);
	if(up - x*acc < 0.5) {
		return up/acc;
	} else {
		return down/acc;
	}
}

double Math::Max(DoubleVector &x, int start) {
	return Max(x, start, (int)x.size() - 1);
}

double Math::Max(DoubleVector &x, int start, int end) {
	int n = start;
	double val = x[start];
	for(n=start;n<=end;n++)
		val = MAX(val, x[n]);
	return val;
}

int Math::MaxInd(DoubleVector &x, bool absolute, double *max) {
	int n = 0;
	int ind = n;
	for(n=0;n<(int)x.size();n++) {
		if(absolute) {
			if(fabs(x[n]) > fabs(x[ind])) {
				ind = n;
			}
		} else {
			if(x[n] > x[ind]) {
				ind = n;
			}
		}
	}
	if(max != NULL)
		*max = absolute ? fabs(x[ind]) : x[ind];
	return ind;
}

double Math::Min(DoubleVector &x, int start) {
	return Min(x, start, (int)x.size() - 1);
}

double Math::Min(DoubleVector &x, int start, int end) {
	int n = start;
	double val = x[start];
	for(n=start;n<=end;n++)
		val = MIN(val, x[n]);
	return val;
}

int Math::MinInd(DoubleVector &x, bool absolute) {
	int n = 0;
	int ind = n;
	for(n=0;n<(int)x.size();n++) {
		if(absolute) {
			if((fabs(x[n]) < fabs(x[ind]) && x[n] == x[n]) || x[ind] != x[ind]) {
				ind = n;
			}
		} else {
			if((x[n] < x[ind] && x[n] == x[n]) || x[ind] != x[ind]) {
				ind = n;
			}
		}
	}
	return ind;
}

int Math::Max(LongVector &x, int start) {
	return Max(x, start, (int)x.size() - 1);
}

int Math::Max(LongVector &x, int start, int end) {
	int n = start;
	int val = x[start];
	for(n=start;n<=end;n++)
		val = MAX(val, x[n]);
	return val;
}

int Math::MaxInd(LongVector &x, bool absolute) {
	int n = 0;
	int ind = n;
	for(n=0;n<=(int)x.size();n++) {
		if(absolute) {
			if(abs(x[n]) > abs(x[ind])) {
				ind = n;
			}
		} else {
			if(x[n] > x[ind]) {
				ind = n;
			}
		}
	}
	return ind;
}

int Math::Min(LongVector &x, int start) {
	return Min(x, start, (int)x.size() - 1);
}

int Math::Min(LongVector &x, int start, int end){
	int n = start;
	int val = x[start];
	for(n=start;n<=end;n++)
		val = MIN(val, x[n]);
	return val;
}

int Math::MinInd(LongVector &x, bool absolute) {
	int n = 0;
	int ind = n;
	for(n=0;n<=(int)x.size();n++) {
		if(absolute) {
			if(abs(x[n]) < abs(x[ind])) {
				ind = n;
			}
		} else {
			if(x[n] < x[ind]) {
				ind = n;
			}
		}
	}
	return ind;
}

double Math::Sum(DoubleVector &x, int start) {
	return Sum(x, start, (int)x.size() - 1);
}

double Math::Sum(DoubleVector &x, int start, int end) {
	double val = 0;
	for(int n=start;n<=end;n++)
		val += x[n];
	return val;
}

int Math::Sum(LongVector &x, int start) {
	return Sum(x, start, (int)x.size() - 1);
}
int Math::Sum(LongVector &x, int start, int end) {
	int val = 0;
	for(int n=start;n<=end;n++)
		val += x[n];
	return val;
}

double Math::Mean(const DoubleVector &x, int start) {
	return Mean(x, start, (int)x.size() - 1);
}

double Math::Mean(const DoubleVector &x, int start, int end) {
	double retVal = 0;
	for(int i=start;i<=end;i++) {
		retVal += x[i] / (end - start + 1);
	}
	return retVal;
}

double Math::Mean(const LongVector &x, int start) {
	return Mean(x, start, (int)x.size() - 1);
}

double Math::Mean(const LongVector &x, int start, int end) {
	double retVal = 0;
	for(int i=start;i<=end;i++) {
		retVal += (double)x[i] / (end - start + 1);
	}
	return retVal;
}

double Math::Var(LongVector &x, int start) {
	return Var(x, start, (int)x.size() - 1);
}

double Math::Var(LongVector &x, int start, int end) {
	double val = 0;
	double mean = Mean(x, start, end);
	for(int n=start;n<=end;n++) {
		val += (x[n] - mean)*(x[n] - mean);
	}
	return val/(end - start);
}
double Math::Var(DoubleVector &x, int start) {
	return Var(x, start, (int)x.size() - 1);
}
double Math::Var(DoubleVector &x, int start, int end) {
	double val = 0;
	double mean = Mean(x, start, end);
	for(int n=start;n<=end;n++) {
		val += (x[n] - mean)*(x[n] - mean);
	}
	return val/(end - start);
}

double Math::RMSError(DoubleVector &x, DoubleVector &y, int start, int end) {
	if(end == -1)
		end = (int)MIN(x.size(), y.size()) - 1;
	DoubleVector v(end - start + 1);
	for(int i=0;i<end-start+1;i++) {
		v[i] = (x[start + i] - y[start + i]) * (x[start + i] - y[start + i]);
	}
	return sqrt(Sum(v) / v.size());
}

double Math::RMSError(DoubleVector &x, DoubleList &y, int start, int end) {
	if(end == -1)
		end = (int)MIN(x.size(), y.size());
	DoubleVector v(end - start + 1);
	DoubleList::iterator it;
	it = y.begin();
	for(int i=0;i<start;i++)
		it++;
	for(int i=0;i<end-start+1;i++, it++) {
		v[i] = (x[start + i] - (*it)) * (x[start + i] - (*it));
	}
	return sqrt(Sum(v) / v.size());
}

bool Math::LMaxMin(DoubleVector &v, double t, LongVector &max, DoubleVector &maxVals, LongVector &min, DoubleVector &minVals) {
	min.clear();
	max.clear();
	minVals.clear();
	maxVals.clear();

	if(v.size() == 1) {
		max.push_back(0);
		min.push_back(0);
		maxVals.push_back(v[0]);
		minVals.push_back(v[0]);
		return true;
	}

	int k = 0;
	int dir = 0;
	if(v[1] >= v[2])
		dir = 1;
	
	int curmax = k;
	int curmin = k;

	for(int n=k+1; n<(int)v.size(); n++) {
		if(v[n] < v[curmax] - t && dir == 1) {
			if(curmax != 0) {
				max.push_back(curmax);
				maxVals.push_back(v[curmax]);
			}
			dir = 0;
			curmin = n;
		} else if(v[n] > v[curmin] + t && dir == 0) {
			if(curmin != 0) {
				min.push_back(curmin);
				minVals.push_back(v[curmin]);
			}
			dir = 1;
			curmax = n;
		}

		if(v[n] > v[curmax] && dir == 1)
			curmax = n;
		else if(v[n] < v[curmin] && dir == 0)
			curmin = n;
	}

	if(min.size() == 0) {
		min.push_back(curmin);
		minVals.push_back(v[curmin]);
	}

	if(max.size() == 0) {
		max.push_back(curmax);
		maxVals.push_back(v[curmax]);
	}

	return true;
}

double Math::Jitter( DoubleVector &x, double time ) {
	double jitter = 0;
	for(size_t i=0; i<x.size()-1; i++)
		jitter += fabs(x[i] - x[i+1]);
	
	jitter /= (time * (x.size() - 2));

	return jitter;
}

double Math::Jitter( DoubleVector &x, DoubleVector &y, DoubleVector &z, double time ) {
	double jitter = 0;
	size_t s = x.size() < y.size() ? x.size() : y.size();
	s = s < z.size() ? s : z.size();

	for(size_t i=0; i<s-1; i++) 
		jitter += sqrt( (x[i] - x[i+1]) * (x[i] - x[i+1]) + (y[i] - y[i+1]) * (y[i] - y[i+1]) + (z[i] - z[i+1]) * (z[i] - z[i+1]) );

	jitter /= (time * (x.size() - 2));

	return jitter;
}

Math::DoubleVector Math::cross( const DoubleVector x, const DoubleVector y )
{
	assert(x.size() == 3 && y.size() == 3);

	DoubleVector v(3,0);
	v[0] = x[1]*y[2] - x[2]*y[1];
	v[1] = x[2]*y[0] - x[0]*y[2];
	v[2] = x[0]*y[1] - y[1]*y[0];
	return v;
}

double Math::norm( const DoubleVector& v, int norm, bool doRoot )
{
	double n = 0;
	double tmp;
	if(norm >= 1) {
		for(size_t i=0; i<v.size(); i++) {
			tmp = 1;
			for(int j=1; j<=norm; j++)
				tmp *= fabs(v[i]);
			n += tmp;
		}
		if(doRoot)
			return pow(n, 1.0/norm);
		else
			return n;
	} else {
		n = fabs(v[0]);
		for(size_t i=1; i<v.size(); i++) {
			tmp = fabs(v[i]);
			if(tmp > n)
				n = tmp;
		}
		return n;
	}
}

DoubleVector Math::sign( const DoubleVector&v )
{
	DoubleVector s = v;
	for(size_t i=0; i<v.size(); i++) {
		s[i] = (v[i] > 0 ? 1 : (v[i] < 0 ? -1 : 0));
	}
	return s;
}

Math::DoubleVector Math::Pow( const DoubleVector &v, double e )
{
	if(e == 1)
		return v;

	DoubleVector s(v.size(), 1);
	if(e == 0)
		return s;

	for(size_t i=0; i<v.size(); i++) {
		s[i] = pow(v[i], e);
	}
	return s;
}

Math::DoubleVector Math::Abs( const DoubleVector &v )
{
	DoubleVector s(v.size(), 0);
	for(size_t i=0; i<v.size(); i++) {
		s[i] = fabs(v[i]);
	}
	return s;
}