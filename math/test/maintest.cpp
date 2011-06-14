#include "TScal_algorithms.h"
#include "TSmatrix.h"
#include "TShorn.h"
#include "cpptest/include/cpptest.h"
#include <fstream>
#include <time.h>

#include "math/include/distortion.h"

using namespace std;
using namespace Math;

// Main test program
//
int main(int argc, char* argv[])
{
	DoubleVector p(3,0), q(3,0), r(3,0), t;
	q[0] = 1; q[1] = 1; q[2] = 1;
	
	DoubleVector v[8];
	for(size_t i=0; i<8; i++)
		v[i].resize(3,0);
	
	v[0][0] = 1;
	v[0][1] = 1;
	v[0][2] = -1;

	v[6][0] = -1;
	v[6][1] = 2;
	v[6][2] = 1;

	CDCube testcube(p, q);
	testcube.SetValues(v);

	int num = 25;

	cout << "R = zeros(" << num+1 << "," << num+1 << "," << num+1 << ",3);\n";
	int i,j,k;
	for(i=0; i<num+1; i++) {
		for(j=0; j<num+1; j++) {
			for(k=0; k<num+1; k++) {
				r[0] = 1.0/num*i; r[1] = 1.0/num*j; r[2] = 1.0/num*k;
				t = testcube.Trilinear(r);
				std::cout << "R(" << i+1 << "," << j+1 << "," << k+1 << ",1) = " << t[0] << ";\n";
				std::cout << "R(" << i+1 << "," << j+1 << "," << k+1 << ",2) = " << t[1] << ";\n";
				std::cout << "R(" << i+1 << "," << j+1 << "," << k+1 << ",3) = " << t[2] << ";\n";
			}
		}
	}
	cout << "close all; clc; for j=1:11\n surf(R(:,:,j)); waitforbuttonpress; end\n ";

/*	srand(time(NULL));
	Test::Suite trackingTestSuite;
	trackingTestSuite.add(auto_ptr<Test::Suite>(new TScal_algorithms));
	trackingTestSuite.add(auto_ptr<Test::Suite>(new TSmatrix));
	trackingTestSuite.add(auto_ptr<Test::Suite>(new TShorn));
	
	// html output 
	// get the path
	// char path[256]
	char file[256];
	string filePath;
		
	filePath = "htmlTestOutput.html";
	strcpy(file, filePath.c_str());	//to char*

//
	ofstream datei(file);
	if( datei )
	{
		auto_ptr<Test::Output> output (new Test::HtmlOutput);
		trackingTestSuite.run(*output,true);
		Test::HtmlOutput* html = dynamic_cast<Test::HtmlOutput*>(output.get());
		html->generate(datei, true, "TransformMatrix");
	}
	datei.close();
//
	// text output
	Test::TextOutput output(Test::TextOutput::Verbose);
	int val = trackingTestSuite.run(output);
	cin.get(); */

	return 0;
}
