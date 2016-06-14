#include "ContinuousWaveletTransform.h"
#include "FFT2.h"
#include <IGOR/IGORitx.hpp>
#include <cmath>

void wavelet_test(){
	std::vector<double> test(20000);
	for(size_t i=0; i<test.size(); ++i){
		double f = 0.0000001*i;
		test[i] = sin(2.0*3.1415*i)+0.5*sin(1.0*i*i);
	};
//	IGORdata::write_itx(test, "test.itx", "test");
}

void FFT2complex_test(){
	std::vector<double> real(1024);
	std::vector<double> imag(1024);
	for(size_t i=0; i<1024;++i){
		double f = 0.0001*i;
//		real[i] = cos(2.0*3.1415*i+0.)+0.5*sin(;
//		real[i] = cos(2.0*3.1415*i/20+0.1*sin(0.1*i))+0.5*sin(1.0*i*i);
		real[i] = cos(2.0*M_PI*i/20);
		imag[i] = sin(2.0*M_PI*i*i/1000);
//		real[i] = cos(2.0*M_PI*i/15)+0.5*sin(1.0*i*i/5);
//		imag[i] = -sin(2.0*M_PI*i/15)+0.5*cos(1.0*i*i/5);
	}
	IGORdata::write_itx(real, "test.itx", "test");

	FFT2 fft2;
	fft2.set(real, 1000, 1.0, FFT2::Rectangular);
	FFT2Complex fft2c;
	fft2c.set(real, imag, 1000, 1.0, FFT2::Rectangular);
	
	std::vector<double> result(fft2.size());
	for(size_t i=0; i<fft2.size(); ++i)	result[i] = fft2.real(i);
	IGORdata::write_itx(result, "result.itx", "result");

	std::vector<double> result_r(fft2.size());
	std::vector<double> result_i(fft2.size());
	for(size_t i=0; i<fft2.size(); ++i){
		result_r[i] = fft2c.real(i);
		result_i[i] = fft2c.imag(i);
	}
	IGORdata::write_itx(result_r, result_i, "result_c.itx", "result_r", "result_i");
	std::vector<double> power=fft2.Power();
	IGORdata::write_itx(power, "power.itx", "power");

	ContinuousWaveletTransform CWT;
	std::vector<double> a(100);
	for(size_t i=0; i<100; ++i) a[i] = (i+1.0);
	CWT.set(real, 1024, a, 5.0, CWT.Gabor);
	IGORdata::write_itx(a, "freq.itx", "freq");
//	IGORdata::write_itx(CWT.getPower(), "CWTtest.itx", "CWTtest");
	IGORdata::write_itx(CWT.getReal(), "CWTtest.itx", "CWTtest");
	IGORdata::write_itx(CWT.getImag(), "CWTtest.itx", "CWTtest");
}
int main(){
	//wavelet_test();
	FFT2complex_test();
}