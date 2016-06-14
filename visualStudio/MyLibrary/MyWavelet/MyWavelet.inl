//--------------------------------------------------------------------//
//																	  //
//			Class for storing the shape of  Gabor Wavelet			  //
//																	  //
//--------------------------------------------------------------------//
//-----		constructors and = operator		-----
inline MyWavelet::MyWavelet(void)	: how_many_sigma(3.0)
{
	Allocated = false;
	num = 0;
}
inline MyWavelet::~MyWavelet(void){FreeMemory();}
inline MyWavelet::MyWavelet(const MyWavelet& obj)	: how_many_sigma(3.0)
{
	Allocated = false;
	num = 0;
}
inline MyWavelet& MyWavelet::operator = (const MyWavelet& s){
	MemoryControl(s.num);
	for(size_t i=0; i<num; ++i){wave_r[i] = s.wave_r[i];	wave_i[i] = s.wave_i[i];	wave_amp[i] = s.wave_amp[i];}
	sigma = s.sigma; a=s.a;	offset = s.offset;
	return *this;
}

//------	core of the wavelet transform	---
template<class TyW_ >
inline double MyWavelet::get_r(const TyW_& src, size_t size_src, double a_, size_t b){
	if(b < num || b >= size_src-num) return 0.0;	//	return 0.0 if edge region
	double sum=0.0;
	for(size_t i=0; i<num; ++i) sum += src[i+b-offset]*wave_r[i];
	return sum;
}
template<class TyW_ >
inline double MyWavelet::get_i(const TyW_& src, size_t size_src, double a_, size_t b){
	if(b < num || b >= size_src-num) return 0.0;	//	return 0.0 if edge region
	double sum=0.0;
	for(size_t i=0; i<num; ++i) sum += src[i+b-offset]*wave_i[i];
	return sum;
}
template<class TyW_ >
inline double MyWavelet::get_avg(const TyW_& src, size_t size_src, double a_, size_t b){
	if(b < num || b >= size_src-num) return 0.0;	//	return 0.0 if edge region
	double sum=0.0;
	for(size_t i=0; i<num; ++i) sum += src[i+b-offset]*wave_amp[i];
	return sum;
}

inline void MyWavelet::makeWave(double sigma_, double a_){
	sigma = sigma_;	a = a_;
	double coef = 1.0/(2.0*sqrt(M_PI*a)*sigma);
	double delta_t = 1.0*a*sigma;
	offset = (size_t)floor(how_many_sigma*delta_t);
	size_t num_tmp = 2*offset-1;
	if(offset==0) num_tmp=0;

	//	memory allocation
	MemoryControl(num_tmp);
	for(size_t i=0; i<num; ++i){
		double coef_tmp = coef*exp(-(1.0*i-offset)*(1.0*i-offset)/(1.0*a*a*sigma*sigma));
		wave_r[i]   = coef_tmp * cos(2.0*M_PI*(1.0*i-offset)/a);
		wave_i[i]   =-coef_tmp * sin(2.0*M_PI*(1.0*i-offset)/a);
		wave_amp[i] = coef_tmp;
	}
	//	this part is nescesarry for the condition" sum of real part is 0"
	double plus_sum = 0.0;
	double mins_sum = 0.0;
	for(size_t i=0; i<num; ++i)	{
		if(wave_r[i] >= 0.0) plus_sum += wave_r[i];
		else				 mins_sum -= wave_r[i];
	}
	double coef_plus = mins_sum/plus_sum;
	for(size_t i=0; i<num; ++i)	{
		if(wave_r[i] >= 0.0) wave_r[i]*=coef_plus;
	}
/*
	
	double offset_r = 0.0;
	for(size_t i=0; i<num; ++i)	offset_r += wave_r[i];
	offset_r /= num;
	for(size_t i=0; i<num; ++i)	wave_r[i] -= offset_r;
*/
}

inline std::vector<double> MyWavelet::getWaveletShape_r(){
	std::vector<double> tmp(num);
	for(size_t i=0; i<num; ++i) tmp[i] = wave_r[i];
	return tmp;
}

inline std::vector<double> MyWavelet::getWaveletShape_i(){
	std::vector<double> tmp(num);
	for(size_t i=0; i<num; ++i) tmp[i] = wave_i[i];
	return tmp;
}

//----	functions and variables for memory control	---	
inline void MyWavelet::MemoryControl(size_t num_){
	if(Allocated && num != num_) FreeMemory();
	num = num_;
	if(!Allocated) Allocation();
}

inline void MyWavelet::Allocation(){
	if(!Allocated) {
		wave_r =   new double [num];
		wave_i =   new double [num];
		wave_amp = new double [num];
		Allocated = true;
	}
}

inline void MyWavelet::FreeMemory(){
	if(Allocated){
		delete [] wave_r;
		delete [] wave_i;
		delete [] wave_amp;
		Allocated = false;
	}
}
//--------------------------------------------------------------------//
//																	  //
//	Class for Continuous Wavelet Transform using Gabor Wavelet		  //
//																	  //
//--------------------------------------------------------------------//

inline MyWaveletTransform_1a::MyWaveletTransform_1a(void)
{
	Allocated = false;
	num = 0;
}
inline MyWaveletTransform_1a::~MyWaveletTransform_1a(void){FreeMemory();}

inline MyWaveletTransform_1a::MyWaveletTransform_1a(const MyWaveletTransform_1a& obj){
	Allocated = false;
	num=0;
};
inline MyWaveletTransform_1a& MyWaveletTransform_1a::operator=(const MyWaveletTransform_1a& s){
	MemoryControl(s.num);
	for(size_t i=0; i<num; ++i){data_r[i] = s.data_r[i];	data_i[i] = s.data_i[i];}
	wavelet = s.wavelet;
	return *this;
}


template<class TyWT_ >
inline void MyWaveletTransform_1a::set(const TyWT_& src, size_t size_src, double a, double sigma_){
	num = size_src;
	MemoryControl(num);
	wavelet.makeWave(sigma_, a);

	for(size_t i=0; i<num; ++i){
		data_r[i]   = wavelet.get_r(  src, num, a, i);
		data_i[i]   = wavelet.get_i(  src, num, a, i);
		data_avg[i] = wavelet.get_avg(src, num, a, i);
	}
}
inline std::vector<double> MyWaveletTransform_1a::getResult_r(){
	std::vector<double> result(num);
	for(size_t i=0; i<num; ++i) result[i] = data_r[i];
	return result;
}

inline std::vector<double> MyWaveletTransform_1a::getResult_i(){
	std::vector<double> result(num);
	for(size_t i=0; i<num; ++i) result[i] = data_i[i];
	return result;
}

//----	smoothing	----
inline std::vector<double> MyWaveletTransform_1a::getResult_amp2(size_t avg_num){
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size);
	for(size_t i=0; i<avg_size; ++i){
		double avg = 0.0;
		for(size_t j=0; j<avg_num; ++j)
			avg += data_r[i*avg_num+j]*data_r[i*avg_num+j] + data_i[i*avg_num+j]*data_i[i*avg_num+j];
		result[i] = avg/avg_num;
	}
	return result;
}

inline std::vector<double> MyWaveletTransform_1a::getResult_amp_scaled(size_t avg_num){
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size);
	for(size_t i=0; i<avg_size; ++i){
		double amp_avg = 0.0;
		double avg = 0.0;
		for(size_t j=0; j<avg_num; ++j){
			amp_avg += data_r[i*avg_num+j] * data_r[i*avg_num+j] 
					 + data_i[i*avg_num+j] * data_i[i*avg_num+j];
			avg += data_avg[i*avg_num+j]*data_avg[i*avg_num+j];
		}
		if(avg <= 0.0) result[i] = 0.0;
		else result[i] = sqrt(amp_avg/avg);
	}
	return result;
}

inline std::vector<double> MyWaveletTransform_1a::getResult_real(size_t avg_num){	//	avg_times average
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size);
	for(size_t i=0; i<avg_size; ++i){
		double avg = 0.0;
		for(size_t j=0; j<avg_num; ++j)
			avg += data_r[i*avg_num+j];
		result[i] = avg/avg_num;
	}
	return result;
}


//----	correlation	----
inline std::vector<double> MyWaveletTransform_1a::getCorrelation(MyWaveletTransform_1a& src){
	return getCorrelation(src, 1);
}
inline std::vector<double> MyWaveletTransform_1a::getCorrelation(MyWaveletTransform_1a& src, size_t avg_num){
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size, 0.0);
	if(num != src.size()) return result;

	for(size_t i=0; i<avg_size; ++i){
		double avg_r = 0.0;
		double avg_i = 0.0;
		double avg   = 0.0;
		for(size_t j=0; j<avg_num; ++j){
			avg_r += data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + data_i[i*avg_num+j]*src.data_i[i*avg_num+j];
			avg_i += data_i[i*avg_num+j]*src.data_r[i*avg_num+j] - data_r[i*avg_num+j]*src.data_i[i*avg_num+j];	
			avg   += sqrt( (    data_r[i*avg_num+j]*    data_r[i*avg_num+j] +     data_i[i*avg_num+j]*    data_i[i*avg_num+j])
				          *(src.data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + src.data_i[i*avg_num+j]*src.data_i[i*avg_num+j]));
////			avg   += sqrt(src.data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + src.data_i[i*avg_num+j]*src.data_i[i*avg_num+j]);
		}
		if(avg<=0.0) result[i] = 0.0;
		else result[i] = (avg_r*avg_r + avg_i*avg_i) /(avg*avg);
	}
	return result;
}
inline std::vector<double> MyWaveletTransform_1a::getCorrelationReal(MyWaveletTransform_1a& src){
	return getCorrelationReal(src, 1);
}
inline std::vector<double> MyWaveletTransform_1a::getCorrelationReal(MyWaveletTransform_1a& src, size_t avg_num){
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size, 0.0);
	if(num != src.size()) return result;

	for(size_t i=0; i<avg_size; ++i){
		double avg_r = 0.0;
		double avg_i = 0.0;
		double avg   = 0.0;
		for(size_t j=0; j<avg_num; ++j){
			avg_r += data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + data_i[i*avg_num+j]*src.data_i[i*avg_num+j];
			avg_i += data_i[i*avg_num+j]*src.data_r[i*avg_num+j] - data_r[i*avg_num+j]*src.data_i[i*avg_num+j];	
			avg   += sqrt( (    data_r[i*avg_num+j]*    data_r[i*avg_num+j] +     data_i[i*avg_num+j]*    data_i[i*avg_num+j])
				          *(src.data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + src.data_i[i*avg_num+j]*src.data_i[i*avg_num+j]));
		}
		if(avg<=0.0) result[i] = 0.0;
		else result[i] = avg_r/avg;
	}
	return result;
}
inline std::vector<double> MyWaveletTransform_1a::getPhase(MyWaveletTransform_1a& src){
	return getPhase(src, 1);
}
inline std::vector<double> MyWaveletTransform_1a::getPhase(MyWaveletTransform_1a& src, size_t avg_num){
	size_t avg_size = (size_t)floor(1.0*num/avg_num);
	std::vector<double> result(avg_size, 0.0);
	if(num != src.size()) return result;

	for(size_t i=0; i<avg_size; ++i){
		double avg_r = 0.0;
		double avg_i = 0.0;
		for(size_t j=0; j<avg_num; ++j){	//	src.conj * data
			avg_r += data_r[i*avg_num+j]*src.data_r[i*avg_num+j] + data_i[i*avg_num+j]*src.data_i[i*avg_num+j];
			avg_i += data_i[i*avg_num+j]*src.data_r[i*avg_num+j] - data_r[i*avg_num+j]*src.data_i[i*avg_num+j];	
		}
		if(avg_r == 0.0) result[i] = 0.0;
		else result[i] = atan(avg_i/avg_r);

		if(avg_r < 0.0 && avg_i > 0.0) result[i] += M_PI;
		if(avg_r < 0.0 && avg_i < 0.0) result[i] -= M_PI;
	}
	return result;
}

//----	functions and variables for memory control	---	
inline void MyWaveletTransform_1a::MemoryControl(size_t num_){
	if(Allocated && num != num_) FreeMemory();
	num = num_;
	if(!Allocated) Allocation();
}

inline void MyWaveletTransform_1a::Allocation(){
	if(!Allocated) {
		data_r =   new double [num];
		data_i =   new double [num];
		data_avg = new double [num];
		Allocated = true;
	}
}

inline void MyWaveletTransform_1a::FreeMemory(){
	if(Allocated){
		delete [] data_r;
		delete [] data_i;
		delete [] data_avg;
		Allocated = false;
	}
}

//--------------------------------------------------------------------//
//																	  //
//	Class for Continuous Wavelet Transform using Gabor Wavelet		  //
//										(High Level Class)			  //
//																	  //
//--------------------------------------------------------------------//

inline void MyWaveletTransform::FreeMemory(){
	for(size_t i=0; i<wavelet_trans.size(); ++i) 
		wavelet_trans[i].FreeMemory();
}

template<class TyWT_ >
inline void MyWaveletTransform::set
			(const TyWT_& src, size_t size_src, std::vector<double> a_, double sigma_){
	a = a_;
	wavelet_trans.resize(a.size());
	for(size_t i=0; i<a.size(); ++i) wavelet_trans[i].set( src, size_src, a[i], sigma_);
}

inline std::vector<std::vector<double>> MyWaveletTransform::getCor(MyWaveletTransform& src){
	return getCor(src, 1);
}

inline std::vector<std::vector<double>> MyWaveletTransform::getCor(MyWaveletTransform& src, size_t avg_num){
	std::vector<std::vector<double>> cor(0);
	if(size() != src.size() || size2() != src.size2()) return cor;
	cor.resize(size());
	for(size_t i=0; i<size(); ++i) cor[i] = wavelet_trans[i].getCorrelation(src.wavelet_trans[i], avg_num);
	return cor;
}

inline std::vector<std::vector<double>> MyWaveletTransform::getCorReal(MyWaveletTransform& src){
	return getCorReal(src, 1);
}
inline std::vector<std::vector<double>> MyWaveletTransform::getCorReal(MyWaveletTransform& src, size_t avg_num){
	std::vector<std::vector<double>> cor(0);
	if(size() != src.size() || size2() != src.size2()) return cor;
	cor.resize(size());
	for(size_t i=0; i<size(); ++i) cor[i] = wavelet_trans[i].getCorrelationReal(src.wavelet_trans[i], avg_num);
	return cor;
}

inline std::vector<std::vector<double>> MyWaveletTransform::getPhase(MyWaveletTransform& src){
	return getPhase(src, 1);
}
inline std::vector<std::vector<double>> MyWaveletTransform::getPhase(MyWaveletTransform& src, size_t avg_num){
	std::vector<std::vector<double>> phase(0);
	if(size() != src.size() || size2() != src.size2()) return phase;
	phase.resize(size());
	for(size_t i=0; i<size(); ++i) phase[i] = wavelet_trans[i].getPhase(src.wavelet_trans[i], avg_num);
	return phase;
}

inline std::vector<std::vector<double>> MyWaveletTransform::getAmp(size_t avg_num){
	std::vector<std::vector<double>> amp(size());
	if(size() ==0 ) return amp;
	for(size_t i=0; i<size(); ++i) amp[i] = wavelet_trans[i].getResult_amp2(avg_num);
	return amp;
}

inline std::vector<std::vector<double>> MyWaveletTransform::getAmpScaled(size_t avg_num){
	std::vector<std::vector<double>> amp(size());
	if(size() ==0 ) return amp;
	for(size_t i=0; i<size(); ++i) amp[i] = wavelet_trans[i].getResult_amp_scaled(avg_num);
	return amp;
}

inline std::vector<std::vector<double>> MyWaveletTransform::getReal(size_t avg_num){
	std::vector<std::vector<double>> amp(size());
	if(size() ==0 ) return amp;
	for(size_t i=0; i<size(); ++i) amp[i] = wavelet_trans[i].getResult_real(avg_num);
	return amp;
}


inline std::vector<double> MyWaveletTransform::getFrequency(double rate){
	return getFrequency(rate, a);
}
inline std::vector<double> MyWaveletTransform::getFrequency(double rate, std::vector<double>& a_){
	std::vector<double> freq(a_.size());
	for(size_t i=0; i<a_.size(); ++i)
		freq[i] =rate/a_[i];
	return freq;
}

inline std::vector<double> MyWaveletTransform::getAvector(std::vector<double> frequency, double rate){
	std::vector<double> a(frequency.size());
	for(size_t i=0; i<frequency.size(); ++i)
		a[i] = rate/frequency[i];
	return a;
}

inline std::vector<double> MyWaveletTransform::getFrequencyWidth(double rate){
	return getFrequencyWidth(rate, a, getsigma());
}
inline std::vector<double> MyWaveletTransform::getFrequencyWidth(double rate, std::vector<double>& a_, double sigma){
	std::vector<double> dfreq(a_.size());
	for(size_t i=0; i<a_.size(); ++i)
		dfreq[i] = rate/(sqrt(2.0)*a_[i]*sigma);
	return dfreq;
}

inline std::vector<double> MyWaveletTransform::getTiming(double rate, double offset, size_t avgnum){
	std::vector<double> tmp(0);
	if(a.size()==0) return tmp;
	size_t timing_num = (size_t)floor(1.0*wavelet_trans[0].size()/avgnum);
	tmp.resize(timing_num);
	for(size_t i=0; i<timing_num; ++i)
		tmp[i] = offset + 1.0/rate*avgnum*(i+0.5);
	return tmp;
}
