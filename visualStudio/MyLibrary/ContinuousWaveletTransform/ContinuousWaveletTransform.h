#pragma once
#include <vector>
#include "fft2.h"


class WaveletBase{
public:
	enum Shape {	Gabor,	MexicanHat,	Dummy	};
};
//-----------------------------------------------------------------//
//																   //
//	class for the wavelet functions with one "a" value			   //
//																   //
//-----------------------------------------------------------------//
class WaveletFunctions:public WaveletBase{
public:
	//---	constructor	---
	WaveletFunctions(void);
	~WaveletFunctions(void);
	size_t size(){return num;}
	void set(size_t num_, Shape shape_, double a_, double sigma_);
	void set(size_t num_, Shape shape_, double a_);
	double real_t(size_t i){return ft_real[i];}
	double imag_t(size_t i){return ft_imag[i];}
	double real(size_t i){if(IsComplex) return fft2c.real(i); else return fft2.real(i);}
	double imag(size_t i){return fft2c.imag(i);}
	double da();	//	return the frequency resolution da;
	double geta(){return a;}
private:
	Shape shape;		//	store the wavelet shape
	double a;			//	expansion factor (related to the frequency)
	double sigma;		//	shape factor defining the ratio between frequency resolution and time resolution
	double* ft_real;	double* ft_imag;	//	functions in time space
//	double* Fk_real;	double* Fk_imag;	//	functions in frequency space	(Fourier transformed of ft)
	size_t num;
	bool IsComplex;
	FFT2 fft2;
	FFT2Complex fft2c;
	
	void DefineWavelet(Shape shape_, double a_, double sigma_);
	
//---	functions for the memory control	---
	bool Allocated;
	void Allocation(size_t num_);
	void FreeMemory();
	void MemoryControl(size_t num_);
};

//-----------------------------------------------------------------//
//																   //
//	class for storing wavelet coefficients for 1 channel data	   //
//																   //
//-----------------------------------------------------------------//
class WaveletCoef{
private:
	std::vector<double> Wab;
	double a;
public:
	void resize(size_t anum){Wab.resize(anum);}
	double& operator[] (size_t i) {return Wab[i];}
};

//-----------------------------------------------------------------//
//																   //
//				class for the wavelet analysis					   //
//																   //
//-----------------------------------------------------------------//
class ContinuousWaveletTransform:public WaveletBase
{
public:
	ContinuousWaveletTransform(void){};
	~ContinuousWaveletTransform(void){};
	template< class TyW_ >
	void set( TyW_ & data, size_t num_, std::vector<double> a, Shape shape);
	template< class TyW_ >
	void set( TyW_ & data, size_t num_, std::vector<double> a, double sigma_, Shape shape);
	std::vector<std::vector<double>> getReal();
	std::vector<std::vector<double>> getImag();
	std::vector<std::vector<double>> getPower();
private:
	size_t num;
	FFT2 data_fft2;
	std::vector<WaveletFunctions> Wavelet;
	std::vector<WaveletCoef> Coef;
};


#include "ContinuousWaveletTransform.inl"