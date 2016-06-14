#pragma once
#include<vector>
#ifndef __GSL_MATH_H__
#define _USE_MATH_DEFINES
#include <math.h>
#endif
#include <gsl/gsl_fft_real.h>

//---------------------------------------------------------------//
//																 //
//	base functions for sharing static variables and functions	 //
//																 //
//---------------------------------------------------------------//
class FFT2base{
public:
	double* f;
	size_t size;

	FFT2base(void);
	~FFT2base(void);
	FFT2base(const FFT2base& obj);
	FFT2base& operator = (const FFT2base& obj);

//---	Tip functions	---
	size_t Size() const {return size;}
	bool IsExist() const {return Allocated;}
//---	check the memory allocation	---
	void MemoryControl(size_t size_tmp);
	void Free();
	void Allocate(size_t size_tmp);

	void div(double val);
	void divReal(const FFT2base& src);	//	deviding real part of src.f
	void divAmp(const FFT2base& src);	//	deviding amplitude of src.f
	void add(const FFT2base& src);
	void addReal(const FFT2base& src);
	void mul(const double val);
	void mulAmp(const FFT2base& src);
	void addmul(const FFT2base& src0, const FFT2base& src1);	//	src.f(conjugate()) * f
	void addmuldagger(const FFT2base& src0, const FFT2base& src1);	//	src.f(conjugate()) * f
	void set_zero();
	void set(const FFT2base& src);
	void amp(const FFT2base& src);
	void addamp(const FFT2base& src);
	std::vector<double> power();
	std::vector<double> amplitude();
	std::vector<double> phase();
private:
	bool Allocated;
};

//-------------------------------------------//
//											 //
//	core class for Fast Fourier Transform	 //
//											 //
//-------------------------------------------//
class FFT2: public FFT2base
{
public:	
	double rate;
//---	window functions	---
	enum FFTwindow{
		Rectangular,
		Hann,
		Hamming,
		Blackman
	};
//---	preferred point number is power of 2	---
	static inline size_t getPower2Num(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0)+0.5);
		return (size_t)pow(2.0, 1.0*power2);	}
	static inline size_t getPower2NumFloor(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0) + 0.01);
		return (size_t)pow(2.0, 1.0*power2);	}

	void RectangularWindow(){	return;	}
	void HannWindow(){		for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.5-0.5*cos(2.0*M_PI*i/(  Size()-1)));	}
	void HammingWindow(){	for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.54-0.46*cos(2.0*M_PI*i/(Size()-1)));}
	void BlackmanWindow(){	for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.42-0.5*cos(2.0*M_PI*i/( Size()-1))+0.08*cos(4.0*M_PI*i/(Size()-1)));}

//---	Tip functions	---
	double operator [] (const size_t i) const {return f[i];}
	double Rate() const {return rate;}

//---	fft core functions	---
	template < class Tyv_ >
	void set(Tyv_& raw_data, size_t size_, double rate_, FFTwindow window);
	template < class Tyv_ >
	void set(Tyv_& raw_data, size_t size_, size_t offset_, double rate_, FFTwindow window);
	template < class Tyv_ >
	void get(std::vector< Tyv_ >& stdv_raw_data, double rate, FFTwindow window);
	template < class Tyv_ >
	void get(Tyv_& raw_data, size_t size_, double rate, FFTwindow window);
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
};

class FFT2Cor
{
public:
	template < class Tyv0_, class Tyv1_ >
	void get(std::vector< Tyv0_ >& src0, std::vector< Tyv1_>& src1, double rate, size_t avgnum, FFT2::FFTwindow);

	void add(FFT2& src0, FFT2& src1);
	void addCorrelation(FFT2& src0, FFT2& src1);
};

//----------------------------------------------------------------------//
//																		//
//	high level class for short time FFT, and its time development		//
//																		//
//----------------------------------------------------------------------//
class FFT2development{
public:
	FFT2development(void);
	~FFT2development(void);
	void Free();
	//---	tip functions	---
	size_t FFTsize()const{return fft_size;}
	size_t FFTnum()const{return fft_num;}

	//---	core functions	---
	template < class Tyv0_ >
	void set(Tyv0_& src0, size_t data_size, size_t fft_size, size_t fft_num, double rate, FFT2::FFTwindow window);
	std::vector<std::vector<double>> getPower(size_t avgnum);
	std::vector<std::vector<double>> getCrossPower(FFT2development& src, size_t avgnum);
	std::vector<std::vector<double>> getCrossPhase(FFT2development& src, size_t avgnum);
	std::vector<double> getFrequency();
	std::vector<double> getFFTtiming(double offset);
	

	std::vector<FFT2> fftv;
private:
	size_t avgnum;		//	number of average
	size_t fft_num;		//	number of time development 
	size_t fft_size;	//	number of fft point
	size_t start_p, increment;	//	point of start for fftdevelopment and increment
	double rate;
	FFT2base CrossPower, Amplitude0, Amplitude1, CrossPhase; 
};


#include "FFT2.inl"

