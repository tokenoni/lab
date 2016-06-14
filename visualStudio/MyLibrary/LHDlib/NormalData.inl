//----------------------------------------------------------//
//															//
//		class for data manipulating time-ordered data		//
//															//
//----------------------------------------------------------//
inline NormalData::NormalData(void){
	Allocated = false; 
	size = 0;};
inline NormalData::~NormalData(void){
	FreeMemory();
}
inline NormalData::NormalData(const NormalData& obj){
	if(obj.Allocated){
		MemoryControl(obj.size);
		for(size_t i=0; i<size; i++) data[i] = obj.data[i];
	}else{
		Allocated = false; size=0;
	}
}
inline NormalData& NormalData::operator=(const NormalData& obj){
	if(obj.Allocated){
		MemoryControl(obj.size);
		for(size_t i=0; i<size; i++) data[i] = obj.data[i];
	}else{
		Allocated = false; size=0;
	}
	return *this;
}

inline void NormalData::MemoryControl(size_t size_){
	//	if size is different, free the memory
	if(Allocated && size != size_)	FreeMemory();
	size = size_;
	if(!Allocated) Allocation();
}

inline void NormalData::Allocation(){
	//	if the memory is not allocated 
	if(!Allocated){
		data = new double[size];
		Allocated = true;
	}
}

inline void NormalData::FreeMemory(){
	if(Allocated){
		delete [] data;
		Allocated = false;
	}
}

//---	smoothing	---
inline std::vector<double> NormalData::getSmoothedVector(size_t avg_num){
	size_t num = (size_t)floor(1.0*size/avg_num);	//	number of data points of the smoothed data
	std::vector<double> smoothed_v(num);
	for(size_t i=0; i<num; ++i){
		double avg = 0.0;
		for(size_t j=0; j<avg_num; ++j){
			avg += data[i*avg_num+j];
		}
		smoothed_v[i] = avg/avg_num;
	}
	return smoothed_v;
}
//---	smoothing	---
inline std::vector<double> NormalData::getSmoothedVector_t(size_t avg_num){
	size_t num = (size_t)floor(1.0*size/avg_num);	//	number of data points of the smoothed data
	std::vector<double> smoothed_v(num);
	for(size_t i=0; i<num; ++i){
		double avg = 0.0;
		for(size_t j=0; j<avg_num; ++j){
			avg += getTiming(i*avg_num+j);
		}
		smoothed_v[i] = avg/avg_num;
	}
	return smoothed_v;
}

//---	FFT analysis	----
inline void NormalData::setFFTdevelopment(double start_t, double end_t, size_t timing_num, size_t fft_point_num_, FFT2::FFTwindow fftwindow){
	start_t -= offset;
	end_t   -= offset;
	//	deriving the index under consideration
	int start_p = (int)floor(start_t*rate);
	int end_p   = (int)ceil( end_t*rate);
	if(start_p<0) start_p=0;
	if(end_p  <0) end_p = 0;
	if(end_p >= size) end_p = size -1;
	if(start_p >= end_p) start_p = 0;
	size_t data_num = end_p- start_p +1;

	double* data_pointer = data+start_p;
	fft2development.set(data_pointer,data_num, fft_point_num_, timing_num, rate, fftwindow);
}

inline std::vector<std::vector<double>> NormalData::getFFTdevelopment(size_t avg_num){
	return fft2development.getPower(avg_num);
}

inline std::vector<std::vector<double>> NormalData::getFFTCorrelation(size_t avg_num, NormalData& src){
	return fft2development.getCrossPower(src.fft2development, avg_num);
}

inline std::vector<std::vector<double>> NormalData::getFFTPhase(size_t avg_num, NormalData& src){
	return fft2development.getCrossPhase(src.fft2development, avg_num);
}

inline void NormalData::FreeFFTdevelopment(){
	fft2development.Free();
}

inline std::vector<double> NormalData::getFFTfrequency(){
	return fft2development.getFrequency();
}

inline std::vector<double> NormalData::getFFTtiming(double start_t){
	return fft2development.getFFTtiming(offset+start_t);
}

//-------		wavelet analysis		--------
inline void NormalData::WaveletMemoryFree(){
	wavelet.FreeMemory();
}

inline std::vector<std::vector<double>> NormalData::getWavelet(size_t avg_num){
	return wavelet.getAmp(avg_num);
}

inline std::vector<std::vector<double>> NormalData::getWaveletReal(size_t avg_num){
	return wavelet.getReal(avg_num);
}

inline std::vector<std::vector<double>> NormalData::getWaveletScaled(size_t avg_num){
	return wavelet.getAmpScaled(avg_num);
}

inline void NormalData::setWavelet(std::vector<double> a, double sigma, double start_t, double final_t){
	size_t start_p = (size_t)floor(0.5+(start_t-offset)*rate);
	size_t final_p = (size_t)floor(0.5+(final_t-offset)*rate);
	if(final_p >= size) final_p = size -1;
	if(start_p >= final_p) start_p = 0;
	//	pointer offset is used
	wavelet.set(data+start_p, (final_p - start_p + 1), a, sigma);
	start_p_wl = start_p;
	num_p_wl   = final_p - start_p + 1;
}

inline std::vector<std::vector<double>> NormalData::getWaveletCorrelation(size_t avg_num, NormalData& src){
//	wavelet.set(data + src.start_p_wl, src.num_p_wl, src.wavelet.geta(), src.wavelet.getsigma());
	return wavelet.getCor(src.wavelet, avg_num);
}
inline std::vector<std::vector<double>> NormalData::getWaveletCorrelationReal(size_t avg_num, NormalData& src){
//	wavelet.set(data + src.start_p_wl, src.num_p_wl, src.wavelet.geta(), src.wavelet.getsigma());
	return wavelet.getCorReal(src.wavelet, avg_num);
}
inline std::vector<std::vector<double>> NormalData::getWaveletPhase(size_t avg_num, NormalData& src){
//	wavelet.set(data + src.start_p_wl, src.num_p_wl, src.wavelet.geta(), src.wavelet.getsigma());
	return wavelet.getPhase(src.wavelet, avg_num);
}

inline std::vector<double> NormalData::getWaveletFrequency(){
	return wavelet.getFrequency(rate);
}
inline std::vector<double> NormalData::getWaveletTiming(double start_t, size_t avgnum){
	return wavelet.getTiming(rate, start_t, avgnum);
}


//----------	storing		------------
template< class TySrcy_, class TySrcx_>
inline void NormalData::setAsIs(TySrcy_& y_src, TySrcx_& x_src, size_t size_src){
	//---	memory allocation	---
	MemoryControl(size_src);
	for(size_t i=0; i<size_src; ++i) data[i] = y_src[i];
	rate = 	(1.0*size_src-1.0)/(x_src[size_src-1] - x_src[0]);
	offset = x_src[0];
}
template< class TySrcy_>
inline void NormalData::setAsIs(TySrcy_& y_src, double offset_, double rate_, size_t size_src){
	//---	memory allocation	---
	MemoryControl(size_src);
	for(size_t i=0; i<size_src; ++i) data[i] = y_src[i];
	rate = 	rate_;
	offset = offset_;
}
//----------	resizing	------------
template< class TySrcY_, class TySrcX_, class Tydest_ >
inline void NormalData::set(TySrcY_& y_src, TySrcX_& x_src, Tydest_& x_dest, size_t size_src, size_t size_dest){
	//---	memory allocation	---
	MemoryControl(size_dest);

	//---	linear interpolation	---
	size_t j_floor=0, j_ceil=1;
	for(size_t i=0; i<size_dest; ++i){
		for(size_t j=j_ceil; j<size_src; ++j){
			if(x_dest[i] <= x_src[j]){ 
				j_floor=j-1; j_ceil=j;
				break;
			}
		}
		if(j_floor==0)	data[i] = 0.0;
		else if(j_floor==size_src-1) data[i] = 0.0;
		else{
			double ratio_floor =  (x_src[j_ceil] -x_dest[i])/(x_src[j_ceil]-x_src[j_floor]);
			double ratio_ceil  = -(x_src[j_floor]-x_dest[i])/(x_src[j_ceil]-x_src[j_floor]);
			data[i] = ratio_floor*y_src[j_floor] + ratio_ceil*y_src[j_ceil];
		}
	}
	rate = (1.0*size_dest-1.0)/(x_dest[size_dest-1] - x_dest[0]);
	offset = 0;
}

template< class TySrcY_, class TySrcX_, class Tydest_ >
inline void NormalData::set_averaged(TySrcY_& y_src, TySrcX_& x_src, Tydest_& x_dest, size_t size_src, size_t size_dest){
	//---	memory allocation	---
	MemoryControl(size_dest);

	size_t j_tmp=0;
	for(size_t i=0; i<size_dest-1; ++i){
		double sum = 0.0;
		size_t num_sum = 0;
		for(size_t j=j_tmp; j<size_src; ++j){
			if(x_src[j] >= x_dest[i+1]) { if(j>1)j_tmp = j-1; break;	}
			sum += y_src[j];
			++num_sum;
		}
		//---	averaging	----	
		if(num_sum != 0) data[i] = sum/num_sum;
		else if(x_src[0] > x_dest[i]) data[i] = 0.0;
		//---	linear interpolation	---
		else if(j_tmp<x_src.size()){
			double ratio_floor =  (x_src[j_tmp+1] -x_dest[i])/(x_src[j_tmp+1]-x_src[j_tmp]);
			double ratio_ceil  = -(x_src[j_tmp]-x_dest[i])/(x_src[j_tmp+1]-x_src[j_tmp]);
			data[i] = ratio_floor*y_src[j_tmp] + ratio_ceil*y_src[j_tmp+1];
		}
	}
	rate = (1.0*size_dest-1.0)/(x_dest[size_dest-1] - x_dest[0]);
	offset = x_dest[0];
}

//----------Tip functions	---------------
inline void NormalData::VectorAdd(std::vector<double>& dest, const std::vector<double>& src){ VectorAdd(dest, src, 1.0);}
inline void NormalData::VectorAdd(std::vector<double>& dest, const std::vector<double>& src, double scale){
	if(dest.size() != src.size()){
		dest.resize(src.size());
		for(size_t i=0; i<dest.size(); ++i) dest[i]=0.0;
	}
	for(size_t i=0; i<dest.size(); ++i) dest[i] += src[i]*scale;
}
