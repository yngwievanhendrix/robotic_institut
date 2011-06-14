#include "../include/cal_algorithms.h"
#include "cpptest/include/cpptest.h"
#include "math/include/math.h"

using namespace std;

class TScal_algorithms : public Test::Suite
{
public:

	Math::Matrix m_system;
	Math::Matrix m_offset;
	std::vector<Math::Matrix> m_leftlist;
	std::vector<Math::Matrix> m_rightlist;

	TScal_algorithms();
	~TScal_algorithms();

	Calibration::CCalibAlgorithmBookKeeper m_bk;
    
private:

	void prepare();

	void close();

	void test_Tensorproduct();

	void test_CalibAlgs();
};
