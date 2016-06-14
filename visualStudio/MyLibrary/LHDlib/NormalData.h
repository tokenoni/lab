#ifndef NORMALDATA_H
#define NORMALDATA_H

#include <MyWavelet/MyWavelet.h>
#include <MyFunc/FFT2.h>

//----------------------------------------------------------//
//															//
//		class for data manipulating time-ordered data		//
//															//
//----------------------------------------------------------//
class NormalData{
public:
	NormalData(void);
	NormalData(const NormalData& obj);
	~NormalData(void);
	NormalData& operator=(const NormalData& obj);

	//---	basic functions	---
public:
	size_t Size(){return size;}
	double Rate(){return rate;}
	void setRate(double rate_){ rate = rate_;}
	
	inline double* get(){ return data;};
	inline double* get(size_t i){ return (data+i);}
	inline double operator [] (size_t i){return data[i];}

	//---	storing data with original size	---
	template< class TySrcy_, class TySrcx_>
	void setAsIs(TySrcy_& y_src, TySrcx_& x_src, size_t size);
	template< class TySrcy_>
	void setAsIs(TySrcy_& y_src, double offset_, double rate_, size_t size);

	//---	resizing	---
	template< class TySrcY_, class TySrcX_, class Tydest_ >
	void set(TySrcY_& y_src, TySrcX_& x_src, Tydest_& x_dest, size_t size_src, size_t size_dest);
	template< class TySrcY_, class TySrcX_, class Tydest_ >
	void set_averaged(TySrcY_& y_src, TySrcX_& x_src, Tydest_& x_dest, size_t size_src, size_t size_dest);
	//---	smoothing	---
	std::vector<double> getSmoothedVector(size_t avg_num);
	std::vector<double> getSmoothedVector_t(size_t avg_num);
//	std::vector<double> getSubVector(size_t avg_num);
//	std::vector<double> getSubVector_t(size_t avg_num);

	//---	FFT analysis	----
public:
	//FFT2 fft2;
	FFT2development fft2development;
	//std::vector<FFT2> fft_v;
	void setFFTdevelopment(double start_t, double end_t, size_t timing_num, size_t fft_point_num_, FFT2::FFTwindow fftwindow=FFT2::Blackman);
	std::vector<std::vector<double>> getFFTdevelopment(){return getFFTdevelopment(1);}
	std::vector<std::vector<double>> getFFTdevelopment(size_t avg_num);
	std::vector<std::vector<double>> getFFTCorrelation(size_t avg_num, NormalData& src);
	std::vector<std::vector<double>> getFFTPhase(size_t avg_num, NormalData& src);
	void FreeFFTdevelopment();
	std::vector<double> getFFTfrequency();
	std::vector<double> getFFTtiming(double start_t);


	//---	wavelet analysis	----
public:
//	void setWavelet(double a, double sigma);
//	void setWavelet(double a, double sigma, double start_t, double final_t);
	std::vector<std::vector<double>> getWavelet(size_t avg_num);
	std::vector<std::vector<double>> getWavelet(std::vector<double> a, double sigma, double start_t, double final_t, size_t avg_num);
	std::vector<std::vector<double>> getWaveletScaled(size_t avg_num);
	std::vector<std::vector<double>> getWaveletReal(size_t avg_num=1);
	
	void setWavelet(std::vector<double> a, double sigma, double start_t, double final_t);
	std::vector<std::vector<double>> getWaveletCorrelation(size_t avg_num, NormalData& src);
	std::vector<std::vector<double>> getWaveletCorrelationReal(size_t avg_num, NormalData& src);
	std::vector<std::vector<double>> getWaveletPhase(size_t avg_num, NormalData& src);
	void WaveletMemoryFree();
	std::vector<double> getWaveletFrequency();
	std::vector<double> getWaveletTiming(double start_t, size_t avgnum);

	MyWaveletTransform_1a wavelet1a;
	MyWaveletTransform wavelet;
	size_t start_p_wl, num_p_wl;
	
	//----	data strage	----
public:
	double *data;	//	corrected data
	double rate;
	double offset;	//	time = i/rate + offset	

	//----	Tip functions	-------------
	static void VectorAdd(std::vector<double>& dest, const std::vector<double>& src);
	static void VectorAdd(std::vector<double>& dest, const std::vector<double>& src, double scale);


	//----	memory control functions	---
	void MemoryControl(size_t size_);
	inline double getTiming(size_t t){return offset + 1.0*t/rate;}

private:
	void Allocation();
	void FreeMemory();
	bool Allocated;	//	true if the memory is allicated
	size_t size;
};


#include "NormalData.inl"
#endif