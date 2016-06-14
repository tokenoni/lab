#include "MyWavelet.h"
#include <IGOR/IGORitx.hpp>

int main(){
//------	example data	------
	size_t size=10000;
	std::vector<double> example(size, 0.0);
	double freq;
	for(size_t i=size/10; i<size/10*9; ++i){
//		example[i] = exp(-(1.0*i-size/2)*(1.0*i-size/2)/(1.0*size/10*size/10)) * sin(1.0*i/20);
		freq = 10.0+10.0*i/10000;
		example[i] = sin(1.0*i/freq);
	}
	IGORdata::write_itx(example, "example.itx", "example");

	MyWavelet wavelet;
	wavelet.makeWave(10.0, 10.0);
	IGORdata::write_itx(wavelet.getWaveletShape_r(), "WaveletShape_r.itx", "WaveletShape_r");
	IGORdata::write_itx(wavelet.getWaveletShape_i(), "WaveletShape_i.itx", "WaveletShape_i");
	
	MyWaveletTransform_1a wavelet_trans_1a;
	std::vector<double> a(40);
	std::vector<std::vector<double>> result_r(0);
	for(size_t i=0; i<a.size(); ++i){	a[i] = 2.0*i+1;}
	MyWaveletTransform wavelet_trans;

	wavelet_trans.set(example, size, a, 10.0);
//	IGORdata::write_itx(wavelet_trans.getWab_r(), "result_r.itx", "result_r");
	IGORdata::write_itx(wavelet_trans.getSmoothedWab(100), "result_amp.itx", "result_amp");

	IGORdata::write_edgeVector(a, "a_v.itx", "a_v");
//	std::vetor<double> 
}