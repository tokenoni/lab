#ifndef MyWAVELET_H
#define MyWAVELET_H

#ifndef __GSL_MATH_H__
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include <vector>
//--------------------------------------------------------------------//
//																	  //
//			Class for storing the shape of  Gabor Wavelet			  //
//																	  //
//--------------------------------------------------------------------//
class MyWavelet
{
public:
	MyWavelet(void)	;
	~MyWavelet(void);
	MyWavelet(const MyWavelet& obj);
	MyWavelet& operator=(const MyWavelet& s);

	size_t size(){return num;}
	double geta(){return a;}
	double getsigma(){return sigma;}
	std::vector<double> getWaveletShape_r();
	std::vector<double> getWaveletShape_i();
	
	template<class TyW_ >
	double get_r(const TyW_& src, size_t size_src, double a_, size_t b_);
	template<class TyW_ >
	double get_i(const TyW_& src, size_t size_src, double a_, size_t b_);
	template<class TyW_ >
	double get_avg(const TyW_& src, size_t size_src, double a_, size_t b_);
	void makeWave(double sigma_, double a_);
private:
	double* wave_r;			//	real part
	double* wave_i;			//	imaginary part
	double* wave_amp;		//	a gaussian. used for average derivation
	double sigma;
	const double how_many_sigma;
	double a;
	size_t offset;
//	void make
//----	functions and variables for memory control	---	
	bool Allocated;
	size_t num;
	void MemoryControl(size_t num_);	
	void Allocation();
	void FreeMemory();
};

//--------------------------------------------------------------------//
//																	  //
//	Class for Continuous Wavelet Transform using Gabor Wavelet		  //
//			for a certain "a" value 	(Low Level Class)			  //
//																	  //
//--------------------------------------------------------------------//
class MyWaveletTransform_1a{
public:
	MyWaveletTransform_1a(void);
	~MyWaveletTransform_1a(void);
	MyWaveletTransform_1a(const MyWaveletTransform_1a& obj);
	MyWaveletTransform_1a& operator=(const MyWaveletTransform_1a& s);

	MyWavelet wavelet;
	size_t size(){return num;}
	template<class TyWT_ >
	void set(const TyWT_& src, size_t size_src, double a, double sigma_);
	std::vector<double> getResult_r();
	std::vector<double> getResult_i();
	std::vector<double> getResult_amp2(size_t avg_num=1);	//	avg_times average
	std::vector<double> getResult_amp_scaled(size_t avg_num=1);	//	avg_times average
	std::vector<double> getResult_real(size_t avg_num=1);	//	avg_times average
	std::vector<double> getCorrelation(MyWaveletTransform_1a& src);
	std::vector<double> getCorrelation(MyWaveletTransform_1a& src, size_t avg_num);
	std::vector<double> getCorrelationReal(MyWaveletTransform_1a& src);
	std::vector<double> getCorrelationReal(MyWaveletTransform_1a& src, size_t avg_num);
	std::vector<double> getPhase(MyWaveletTransform_1a& src);
	std::vector<double> getPhase(MyWaveletTransform_1a& src, size_t avg_num);
	//---	convenient functions	---
	double geta(){  return wavelet.geta();}
	double getSigma(){ return wavelet.getsigma();}
	//---	frequency, frequency_width	---
	double getf(double rate){ return rate/geta();};
	double getdf(double rate){return rate/(sqrt(2.0)*getSigma()*geta());};
	//---	time_width	---
	double getdt(double rate){return getSigma()*geta()/(sqrt(2.0)*rate);};
private:
	double* data_r;
	double* data_i;
	double* data_avg;
	bool Allocated;
	size_t num;
	void MemoryControl(size_t num_);	
	void Allocation();
public:
	void FreeMemory();
};

//--------------------------------------------------------------------//
//																	  //
//	Class for Continuous Wavelet Transform using Gabor Wavelet		  //
//										(High Level Class)			  //
//																	  //
//--------------------------------------------------------------------//
class MyWaveletTransform{
public:
	MyWaveletTransform(void){};
	~MyWaveletTransform(void){};
	MyWaveletTransform(const MyWaveletTransform& obj){};
	void FreeMemory();

	template<class TyWT_ >
	void set(const TyWT_& src, size_t size_src, std::vector<double> a_, double sigma_);

	std::vector<std::vector<double>> getCor(MyWaveletTransform& src);
	std::vector<std::vector<double>> getCor(MyWaveletTransform& src, size_t avg);
	std::vector<std::vector<double>> getCorReal(MyWaveletTransform& src);
	std::vector<std::vector<double>> getCorReal(MyWaveletTransform& src, size_t avg);
	std::vector<std::vector<double>> getPhase(MyWaveletTransform& src);
	std::vector<std::vector<double>> getPhase(MyWaveletTransform& src, size_t avg);
	std::vector<std::vector<double>> getAmp(size_t avg=1);
	std::vector<std::vector<double>> getAmpScaled(size_t avg=1);
	std::vector<std::vector<double>> getReal(size_t avg=1);
	size_t size(){return a.size();}
	size_t size2() {if(a.size()) return wavelet_trans[0].size(); else return 0;}
	std::vector<MyWaveletTransform_1a> wavelet_trans;

	std::vector<double> geta(){return a;}
	double getsigma(){ if(a.size()) return wavelet_trans[0].wavelet.getsigma(); else return 0;}

	//---	returning the frequency vector	---
	std::vector<double> getFrequency(double rate);
	std::vector<double> getFrequency(double rate, std::vector<double>& a_);
	std::vector<double> getFrequencyWidth(double rate);
	std::vector<double> getFrequencyWidth(double rate, std::vector<double>& a_, double sigma);
	std::vector<double> getTiming(double rate, double offset, size_t avgnum);
	static std::vector<double> getAvector(std::vector<double> frequency, double rate);

private:
	std::vector<double>	a;
};


#include "MyWavelet.inl"

#endif