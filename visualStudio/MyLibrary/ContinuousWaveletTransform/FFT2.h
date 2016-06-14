#pragma once
#include<vector>
#include <gsl/gsl_fft_real.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

class FFT2
{
public:
	double* f;
	double rate;
	bool Allocated;
	size_t num;
	//---	window functions	---
	enum FFTwindow{
		Rectangular,
		Hann,
		Hamming,
		Blackman
	};
	void RectangularWindow(){	return;	}
	void HannWindow(){		for(size_t i=0; i<num; ++i) f[i] = f[i]*(0.5-0.5*cos(2.0*M_PI*i/(num-1)));	}
	void HammingWindow(){	for(size_t i=0; i<num; ++i) f[i] = f[i]*(0.54-0.46*cos(2.0*M_PI*i/(num-1)));}
	void BlackmanWindow(){	for(size_t i=0; i<num; ++i) f[i] = f[i]*(0.42-0.5*cos(2.0*M_PI*i/(num-1))+0.08*cos(4.0*M_PI*i/(num-1)));}

//---	fft core functions	---
	template < class Tyv_ >
	void set(std::vector< Tyv_ >& stdv_raw_data, double rate, FFTwindow window);
	template < class Tyv_ >
	void get(std::vector< Tyv_ >& stdv_raw_data, double rate, FFTwindow window){set(stdv_raw_data, rate, window);};
	template < class Tyv_ >
	void set( Tyv_ & raw_data, size_t num_, double rate, FFTwindow window);
	
//---	practical functions	---
	double real(size_t i);
	double imag(size_t i);
	template< class Tyv_ >
	void Power(Tyv_& power);
	template< class Tyv_ >
	void CrossPower(Tyv_& crosspower, const FFT2& src);
	template< class Tyv_ >
	void Coherency(Tyv_& coherency, FFT2& src);
	template< class Tyv_ >
	void CrossPhase(Tyv_& crossphase, const FFT2& src);
	std::vector<double> Power();
	std::vector<double> Phase();
//---	preferred point number is power of 2	---
	static inline size_t getPower2Num(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0)+0.5);
		return (size_t)pow(2.0, 1.0*power2);	}
	static inline size_t getPower2NumFloor(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0) + 0.01);
		return (size_t)pow(2.0, 1.0*power2);	}
//---	check the memory allocation	---
	void div(double val);
	void div(const FFT2& src);
	void set_zero();
	size_t size(){return num;};

	inline void Allocate(size_t& num_tmp);
	FFT2(void);
	~FFT2(void);
};

//----------------------------------------------------------
//
//			class for Fourier transform for complex values
//
//----------------------------------------------------------
class FFT2Complex: public FFT2{
public:
	FFT2 fft2imag;
	template < class Tyv_ >
	void set(std::vector< Tyv_ >& stdv_real, std::vector< Tyv_ >& stdv_imag, double rate, FFTwindow window);
	template < class Tyv_ >
	void set( Tyv_ & raw_real, Tyv_ & raw_imag, size_t size_, double rate, FFTwindow window);
	double real(size_t i);
	double imag(size_t i);
};

class FFT2Cor: public FFT2
{
public:
	FFT2 fft0, fft1;
	
	template < class Tyv0_, class Tyv1_ >
	void get(std::vector< Tyv0_ >& src0, std::vector< Tyv1_>& src1, double rate, size_t avgnum, FFT2::FFTwindow);

	void add(FFT2& src0, FFT2& src1);
	void addCorrelation(FFT2& src0, FFT2& src1);
};

class FFT2Develop{
public:
	std::vector<FFT2> fft2set;
	void set();
	std::vector<std::vector<double>> Power();
	std::vector<std::vector<double>> Phase();
};

#include "FFT2.inl"

