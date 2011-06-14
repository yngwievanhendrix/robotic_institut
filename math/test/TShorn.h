#include "../include/cal_algorithms.h"
#include "cpptest/include/cpptest.h"

using namespace std;

class TShorn : public Test::Suite
{
public:

	TShorn();
	~TShorn();

private:

	void prepare();

	void close();

	void testHorn();

};
