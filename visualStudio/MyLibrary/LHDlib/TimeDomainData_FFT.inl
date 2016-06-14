inline bool TimeDomainData_FFT::set
	(const TimeDomainData& src, const size_t timeNum, const size_t resampleNum, const double freqNum, FFT2::FFTwindow window_)
{
	if(src.size() <= timeNum) return false;
	size_t offset = floor((1.0*src.size()-freqNum*resampleNum)/(timeNum-1));
	fftDevelopment.resize(timeNum);
	rate = src.getRate()/offset;
	start_t = src.getStart_t() + 0.5/rate;
	end_t   = start_t + (timeNum-1) / rate;
	delta_f = src.getRate()/resampleNum/freqNum;

	if(resampleNum==1) {
		for(size_t t=0; t<timeNum; ++t){
			fftDevelopment[t].setEachpreparation(freqNum, src.getRate(), window_);
			for(size_t t2=0; t2<freqNum; ++t2)
				fftDevelopment[t].setEach(t2, src.get(t*offset + t2));
			fftDevelopment[t].setEachFinish();
		}
	}else{
		for(size_t t=0; t<timeNum; ++t){
			fftDevelopment[t].setEachpreparation(freqNum, src.getRate()/resampleNum, window_);
			for(size_t t2=0; t2<freqNum; ++t2)
				fftDevelopment[t].setEach(t2, src.getAvg(t*offset + t2*resampleNum, resampleNum));
			fftDevelopment[t].setEachFinish();
		}
	}
}

inline std::vector< std::vector <double> > TimeDomainData_FFT::getPower()const{
	std::vector<std::vector<double>> result(0);
	for(size_t t=0; t<fftDevelopment.size(); ++t)
		result.push_back(fftDevelopment[t].Power());
	return result;
}

inline std::vector<double> TimeDomainData_FFT::getFrequency()const{
	std::vector<double> freq(0);
	if(fftDevelopment.size()==0) return freq;
	freq.resize(fftDevelopment[0].Size()/2+1);
	for(size_t i=0; i<freq.size(); ++i)
		freq[i] = delta_f*i;
	return freq;
}

inline std::vector<double> TimeDomainData_FFT::getTiming()const{
	std::vector<double> timing(fftDevelopment.size());
	for(size_t t=0; t<fftDevelopment.size(); ++t)
		timing[t] = start_t + t/rate;
	return timing;
}

inline std::vector< std::vector <double> > TimeDomainData_FFT::getFreqCenters(const size_t numPeak, const size_t t_avgnum, const size_t f_avgnum)const{
	std::vector< std::vector<double>> result(fftDevelopment.size(), std::vector<double>(numPeak));
	std::vector< std::vector<double>> fftPower = getPower();
	std::vector< std::vector<double>> fftPower_smth(fftPower.size(), std::vector<double>(fftPower[0].size(), 0.0));

	for(size_t t=t_avgnum; t<fftPower.size()-t_avgnum-1; ++t){
		for(size_t f=f_avgnum; f<fftPower[0].size()-f_avgnum-1; ++f){
			for(int t2 = -t_avgnum; t2<=(int)t_avgnum; ++t2){
				for(int f2= -f_avgnum; f2<=(int)f_avgnum; ++f2){
	 				fftPower_smth[t][f] += fftPower[t+t2][f+f2];
				}
			}
		}
	}
	return fftPower_smth;
}
