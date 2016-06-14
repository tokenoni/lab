#ifndef __FFT_REAL_HPP__
#define __FFT_REAL_HPP__
#include <gsl/gsl_fft_real.h>
//#include <gsl/gsl_fft_halfcomplex.h>
#include <vector>

static const double pi = 3.1415926535897932384626433832795028841971;

template < class Tyfft_ >
class MyFFT2{
public:
	std::vector<std::vector<double>> Power; //	Power spectrum
	std::vector<double> fft_timing;
	std::vector<double> frequency;
	double rate;

//---	window functions	---
	enum FFTwindow{
		Rectangular,
		Hann,
		Hamming,
		Blackman
	};

//---	summerize	---
	void FFT2vector(const std::vector< Tyfft_ >& rawdata, size_t freq_point_num, FFTwindow window, double rate_){
		rate = rate_;
		size_t size2Power = getPower2Num(freq_point_num);
		size_t time_point_num = size2Power/2;
		size_t time_slice_num = (size_t)floor(rawdata.size()/time_point_num);	//	number of time slice
		fft_timing.resize(time_slice_num);
		Power.resize(time_slice_num);
		double* src_x = new double[size2Power];
		double* amp_x = new double[size2Power];
		for(size_t t=0; t<rawdata.size(); t+=time_point_num){	//	time slice
			double t_avg=0.0;
			for(size_t i=0; i<size2Power; i++){
				src_x[i] = getVal(rawdata[t+i]);
				t_avg += getTime(i);
			}
			// fft
			FFT2(src_x, size2Power, window);
			FFT2Power(amp_x, src_x, size2Power);
			Power[t]=(Pointer2Stdv(amp_x, size2Power/2);
			//time
			fft_timing[t] = t_avg/size2Power;
		}
		//	frequency
		frequency.resize(size2Power/2);
		for(size_t i=0; i<size2Power/2; i++) frequency[i] = rate*i/size2Power;
		delete[] src_x;
		delete[] amp_x;
	}

private:
	virtual double getVal( const Tyfft_ x){return 1.0*x;}	//	correction (raw data -> physical value).
	virtual double getTime( const size_t i){return rate*i;}

	void RectangularWindow(double* x, size_t size){	return;	}
	void HannWindow(double* x, size_t size){
		for(size_t i=0; i<size; i++) x[i] = x[i]*(0.5-0.5*cos(2.0*pi*i/(size-1)));
	}
	void HammingWindow(double* x, size_t size){
		for(size_t i=0; i<size; i++) x[i] = x[i]*(0.54-0.46*cos(2.0*pi*i/(size-1)));
	}
	void BlackmanWindow(double* x, size_t size){
		for(size_t i=0; i<size; i++) 
			x[i] = x[i]*(0.42-0.5*cos(2.0*pi*i/(size-1))+0.08*cos(4.0*pi*i/(size-1));
	}

//---	preferred point number is power of 2	---
	inline size_t getPower2Num(size_t point_num){
		size_t power2 = (size_t)(log(1.0*point_num)/log(2.0) + 0.5);
		return (size_t)pow(2.0, 1.0*power2);
	}

//---	fft core functions	---
	inline void FFT2(double* x, size_t size2Power, FFTwindow window){
		//---	applying the window functions	---
		if(window == Rectangular)	RectangularWindow(x, size2Power);
		else if(window == Hann)		HannWindow(x, size2Power);
		else if(window == Hamming)	HammingWindow(x, size2Power);
		else if(window == Blackman)	BlackmanWindow(x, size2Power);
		//	fft
		gsl_fft_real_radix2_transform(x, 1, size2Power);
		for(size_t i=0; i<size2Power/2; i++) x[i] *= rate;
	}

	inline void FFT2Power(double* amp_x, const double* src_x, size_t size2Power){
		amp_x[0] = src_x[0]*src_x[0];
		for(size_t i=1; i<size2Power/2; i++)
			amp_x[i] = src_x[i]*src_x[i] + src_x[size2Power-i]*src_x[size2Power-i];
	}

	inline void FFT2Covariance(double* cov_xy, const double* src_x, const double* src_y, size_t size2Power){
		cov_xy[0] = src_x[0]*src_y[0];
		for(size_t i=1; i<size2Power/2; i++){
			cov_xy[i] = src_x[i]*src_y[i] - src_x[size2Power-i]*src_x[size2Power-i];
			cov_xy[size2Power-i] = - src_x[i]*src_y[size2Power-i] + src_x[i]*src_x[size2Power-i];
		}
	}

//---	Tip functions	---
public:
	std::vector<double> Pointer2Stdv(const double* src_x, size_t size){
		std::vector<double> tmp(size);
		for(size_t i=0; i<size; i++) tmp[i] = src_x[i];
		return tmp;
	}
	template< class Tyv_ >
	static inline std::vector<double> getEdgeVector(std::vector< Tyv_ >& src){
		std::vector<double> edge(src.size()+1);
		if(src.size()<=1) return edge;
		edge[0] = 1.0*src[0]-2.0*src[1];
		for(size_t i=1; i<src.size(); i++)	edge[i] = 0.5*src[i-1]+0.5*src[i];
		edge[src.size()] = 2.0*src[src.size()-1] - 1.0*src[src.size()-2];
		return edge;
	}

};


#endif