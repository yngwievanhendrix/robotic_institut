#include "math/include/distortion.h"
#include "math/include/Horn.h"

// bottom, front, right, back, left, top
int Math::CDCube::m_Faces[4][6] = { 
	{0,0,1,2,3,5}, 
	{1,4,5,6,7,4}, 
	{2,5,6,7,4,7},
	{3,1,2,3,0,6} };

bool Math::ComputeDistortion(std::vector<Matrix> left, std::vector<DoubleVector> right, std::vector<DoubleVector> &errVec, DoubleVector &errors, Matrix &matrix) {
	matrix.resize(4,4);
	Horn horn(left, right);
	horn.calculateTransformation();
	horn.getTransform(matrix);

	errors.resize(3,0);
	errors[1] = horn.errorSumOfSquares(&errVec) / sqrt((double)left.size());
	errors[0] = Mean(errVec[3]);
	errors[2] = Max(errVec[3]);

	return true;
}

Math::DoubleVector Math::ProjectToPlane( const DoubleVector n, const DoubleVector v, const DoubleVector p )
{
	// normalize n
	DoubleVector _n = 1.0 / norm(n) * n;

	// compute distance
	double d = _n * (p - v);

	// compute projected point
	return p - d * _n;
}

Math::DoubleVector Math::ProjectToLine( const DoubleVector p, const DoubleVector q, const DoubleVector x )
{
	DoubleVector v = q - p;

	// compute lambda
	double lambda = (v * x - v * p ) / (v * v);

	// compute projected point
	return p + lambda * v;
}

Math::DoubleVector Math::InterpolateOnPlane( const DoubleVector p[4], const DoubleVector x, const DoubleVector v[4] )
{
	// compute projected points and interpolate on these lines
	DoubleVector _x[2], _v[2], V[2];	
	for(int i=0; i<2; i++) {
		V[0] = v[2*i]; V[1] = v[(2*i+1) % 4];
		_x[i] = ProjectToLine(p[2*i], p[(2*i+1) % 4], x);
		_v[i] = InterpolateOnLine(p[2*i], p[(2*i+1) % 4], x, V);
	}

	return InterpolateOnLine(_x[0], _x[1], x, _v);
}

Math::DoubleVector Math::InterpolateOnLine( const DoubleVector p, const DoubleVector q, const DoubleVector x, DoubleVector v[2] )
{
	// project point onto line
	DoubleVector _x = ProjectToLine(p, q, x);
	
	// ratios
	double r[2] = { norm(p-_x) / norm(p-q), norm(q-_x) / norm(p-q) };

	// interpolate
	return r[1] * v[0] + r[0] * v[1];
}

Math::CDCube::CDCube( DoubleVector p, DoubleVector q ) : Math::CCube(p, q)
{
	Init(p, q);
}

Math::CDCube::CDCube( DoubleVector v[8] )
{
	DoubleVector p(3,0), q(3,0);
	p = v[0]; q = v[0];
	for(int i=0; i<8; i++) {
		m_Points[i]=v[i];
		for(int j=0; j<3; j++) {
			p[j] = MIN(p[j], v[i][j]);
			q[j] = MAX(q[j], v[i][j]);
		}
	}
}

Math::CDCube::CDCube( void ) : CCube()
{

}

Math::DoubleVector Math::CDCube::Trilinear( DoubleVector p, double tolerance )
{
	assert(Inside(p, tolerance));

	// DoubleVector v[6];
	DoubleVector V[2];

	// v[0].resize(m_Values[0].size(), 0);
	
	/* for(int i=0; i<6; i++)
		v[i] = InterpolateOnFace(i, p);

	V[0] = v[0]; V[1] = v[5];
	v[0] = InterpolateOnLine(ProjectToFace(0, p), ProjectToFace(5, p), p, V);
	V[0] = v[1]; V[1] = v[3];
	v[1] = InterpolateOnLine(ProjectToFace(1, p), ProjectToFace(3, p), p, V);
	V[0] = v[2]; V[1] = v[4];
	v[2] = InterpolateOnLine(ProjectToFace(2, p), ProjectToFace(4, p), p, V);

	return 1.0 / 3.0 * (v[0] + v[1] + v[2]); */

	V[0] = InterpolateOnFace(0, p);
	V[1] = InterpolateOnFace(5, p);

	return InterpolateOnLine(ProjectToFace(0, p), ProjectToFace(5, p), p, V);
}

bool Math::CDCube::Inside( const DoubleVector p, double tolerance )
{
	assert(p.size() == 3);
	DoubleVector n;
	double v;
	tolerance = -fabs(tolerance);

	// check for each face
	for(int i=0; i<6; i++) {
		n = cross(m_Points[m_Faces[0][i]] - m_Points[m_Faces[1][i]], m_Points[m_Faces[0][i]] - m_Points[m_Faces[3][i]]);
		v = n * (p - m_Points[m_Faces[0][i]]);
		if(v < tolerance)
			return false;
	}

	return true;
}

Math::DoubleVector Math::CDCube::InterpolateOnFace( const int f, const DoubleVector p )
{
	DoubleVector q[4], v[4];
	for(int i=0; i<4; i++) {
		q[i] = m_Points[m_Faces[i][f]];
		v[i] = m_Values[m_Faces[i][f]];
	}

	return InterpolateOnPlane(q, p, v);
}

Math::DoubleVector Math::CDCube::ProjectToFace( const int f, const DoubleVector p )
{
	return ProjectToPlane(cross(m_Points[m_Faces[1][f]] - m_Points[m_Faces[0][f]], m_Points[m_Faces[3][f]] - m_Points[m_Faces[0][f]]), 
		m_Points[m_Faces[0][f]], p);
}

void Math::CDCube::SetValues( DoubleVector v[8] )
{
	for(size_t i=0; i<8; i++) {
		assert(v[0].size() == v[i].size());
		m_Values[i] = v[i];
	}
}

void Math::CDCube::SetValue( DoubleVector v, int i )
{
	m_Values[i] = v;
}

void Math::CDCube::Init( DoubleVector p, DoubleVector q )
{
	int i;
	for(i=0; i<8; i++)
		m_Points[i].resize(3,0);

	for(i=0; i<3; i++) {
		m_Points[0][i] = MIN(p[i], q[i]);
		m_Points[6][i] = MAX(p[i], q[i]);
	}

	m_Points[3][0] = m_Points[4][0] = m_Points[7][0] = m_Points[0][0];
	m_Points[1][0] = m_Points[2][0] = m_Points[5][0] = m_Points[6][0];
	m_Points[1][1] = m_Points[4][1] = m_Points[5][1] = m_Points[0][1];
	m_Points[2][1] = m_Points[3][1] = m_Points[7][1] = m_Points[6][1];
	m_Points[1][2] = m_Points[2][2] = m_Points[3][2] = m_Points[0][2];
	m_Points[4][2] = m_Points[5][2] = m_Points[7][2] = m_Points[6][2];
}