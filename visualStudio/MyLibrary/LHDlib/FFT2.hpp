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
	FFT2base(void);
	~FFT2base(void);
	FFT2base(const FFT2base& obj);
	FFT2base& operator = (const FFT2base& obj);
	FFT2base& copy(const FFT2base& obj);

//---	Tip functions	---
	size_t Size() const {return size;}
	bool IsExist() const {return Allocated;}

	std::vector<double> power() const;
	std::vector<double> amplitude() const;
	std::vector<double> phase() const;


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

protected:
	double* f;
	size_t size;

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
//---	window functions	---
	enum FFTwindow{
		Rectangular,
		Hann,
		Hamming,
		Blackman
	};
	FFT2(void):FFT2base(){};
	FFT2(const size_t size_, const double rate_, FFT2::FFTwindow window_);
	FFT2(const FFT2& obj);
	FFT2& operator = (const FFT2& obj);

//---	preferred point number is power of 2	---
	static inline size_t getPower2Num(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0)+0.5);
		return (size_t)pow(2.0, 1.0*power2);	}
	static inline size_t getPower2NumFloor(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0) + 0.01);
		return (size_t)pow(2.0, 1.0*power2);	}


//---	Tip functions	---
	double operator [] (const size_t i) const {return f[i];}
	double Rate() const {return rate;}

//---	fft core functions	---
	template < class Tyv_ >
	void set(const Tyv_& raw_data, size_t size_, double rate_, FFTwindow window);
	template < class Tyv_ >
	void set(const Tyv_& raw_data, size_t size_, size_t offset_, double rate_, FFTwindow window);
	
	void setEachpreparation(const size_t size_, double rate_, FFTwindow window);
	void setEach(const size_t index, const double data_i);
	void setEachFinish();

	template< class Tyv_ >
	void Power(Tyv_& power)const;
	template< class Tyv_ >
	void CrossPower(Tyv_& crosspower, const FFT2& src)const;
	template< class Tyv_ >
	void Coherency(Tyv_& coherency, FFT2& src)const;
	template< class Tyv_ >
	void CrossPhase(Tyv_& crossphase, const FFT2& src)const;
	std::vector<double> Power()const;
	std::vector<double> Phase()const;
	FFTwindow getWindow()const{return window;}

private:
	double rate;
	FFTwindow window;
	void RectangularWindow(){	return;	}
	void HannWindow(){		for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.5-0.5*cos(2.0*M_PI*i/(  Size()-1)));	}
	void HammingWindow(){	for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.54-0.46*cos(2.0*M_PI*i/(Size()-1)));}
	void BlackmanWindow(){	for(size_t i=0; i<Size(); ++i) f[i] = f[i]*(0.42-0.5*cos(2.0*M_PI*i/( Size()-1))+0.08*cos(4.0*M_PI*i/(Size()-1)));}
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
	std::vector<std::vector<double>> getCrossPower(const FFT2development& src, size_t avgnum);
	std::vector<std::vector<double>> getCrossPhase(const FFT2development& src, size_t avgnum);
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

