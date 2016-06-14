#ifndef __TIME_DOMAIN_DATA_FFT_HPP__
#define __TIME_DOMAIN_DATA_FFT_HPP__

#include "FFT2.hpp"
#include "TimeDomainData.hpp"
#include <vector>

class TimeDomainData_FFT{
public:
	bool set(const TimeDomainData& src, const size_t timeNum = 500, const size_t resampleScale=1, const double freqNum=1024, FFT2::FFTwindow window_=FFT2::Blackman);

	size_t f_size() const;	//	returns the frequency size
	size_t t_size() const;  //	returns the timing size
	std::vector<double> getTiming() const;
	std::vector<double> getFrequency() const;
	std::vector< std::vector <double> > getPower()const;
	std::vector< std::vector <double> > getFreqCenters(const size_t numPeak, const size_t t_avgnum, const size_t f_avgnum)const;

private:
	std::vector< FFT2 > fftDevelopment;	//	data of fft_development
	double rate, start_t, end_t;	//	parameters specifying short time FFT timng
	double delta_f;					//	for frequency
};

#include "TimeDomainData_FFT.inl"


#endif
