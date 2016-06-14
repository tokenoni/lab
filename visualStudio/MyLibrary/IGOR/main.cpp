#include <IGOR\IGORitx.hpp>
#include <IGOR\IGORibw.hpp>
#include <stdio.h>

int main(){

	std::vector<int> intdata(0);
	std::vector<float> floatdata(0);
	std::vector<std::vector <float>> floatmatrix(0);

	for(size_t i=0; i<100; i++){
		intdata.push_back(i);
		floatdata.push_back((float)10.0*i);
	}
	for(size_t j=0; j<55; j++){
		floatmatrix.push_back(floatdata);
	}

	IGORdata::write_itx<int>(intdata, "intdata.itx", "intdata");
	IGORdata::write_ibw_onedim_double(intdata, intdata.size(), "intdata.ibw", "intdata_ibw");
//	IGORdata::	

}