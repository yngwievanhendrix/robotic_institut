#include "../include/cal_algorithms.h"
#include "cpptest/include/cpptest.h"

using namespace std;

class TSmatrix : public Test::Suite
{
public:

	TSmatrix();
	~TSmatrix();
    
private:

	void prepare();

	void close();

	void testInversion();

	void testMatrixMultiply();
};
