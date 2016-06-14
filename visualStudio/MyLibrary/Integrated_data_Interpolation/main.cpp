#include "EdgeSpline.h"
#include <vector> 
#include <IGOR/IGORitx.hpp>

int main(){
	//---	example data	---
	std::vector <double> testdata(16);
	testdata[ 0] = 0.00386029;
	testdata[ 1] = 0.00555841;
	testdata[ 2] = 0.00762782;
	testdata[ 3] = 0.0113996; 
	testdata[ 4] = 0.0174413; 
	testdata[ 5] = 0.0350182;
	testdata[ 6] = 0.118844;
	testdata[ 7] = 1.29932;
	testdata[ 8] = 0.287288;
	testdata[ 9] = 0.068146;
	testdata[10] = 0.025427;
	testdata[11] = 0.0140719;
	testdata[12] = 0.00872065;
	testdata[13] = 0.00609149;
	testdata[14] = 0.00434196;
	testdata[15] = 0.00322576;
	std::vector <double> testdata_x(16);
	for(size_t i=0; i<16; ++i) testdata_x[i] = 1.0*i;
	//------------------------
	EdgeLogSplineArea spline;
	spline.set(testdata_x, testdata, 0.8, 16, spline.Natural);
	std::vector<double> x(128), y(128);
	for(size_t i=0; i<x.size(); ++i){
		x[i] = 0.0 + 1.0*testdata.size()/(x.size()-1)*i;
		y[i] = spline.get(x[i]);
	}
	IGORdata::write_itx(x, y, "test.itx", "x_interp", "y_interp");
	IGORdata::write_itx(testdata, "testdata.itx", "testdata");
}

















