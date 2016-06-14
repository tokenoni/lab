//---------------------------------------------------------------//
//																 //
//	base functions for sharing static variables and functions	 //
//																 //
//---------------------------------------------------------------////---	constructor and destructor	---
inline FFT2base::FFT2base(void){	Allocated = false; size = 0;}
inline FFT2base::~FFT2base(void){	
	Free();
}
inline FFT2base::FFT2base(const FFT2base& obj){
	Allocated = false; size = 0;
	if(obj.IsExist()){
		MemoryControl(obj.Size());
		for(size_t i=0; i<obj.Size(); ++i) f[i] = obj.f[i];
	}
	else{	Free();	}
}
inline FFT2base& FFT2base::operator = (const FFT2base& obj){
	if(obj.IsExist()){
		MemoryControl(obj.Size());
		for(size_t i=0; i<obj.Size(); ++i) f[i] = obj.f[i];
	}
	else{	Free();	}
	return *this;
}
//---	Functions for memory control	---
inline void FFT2base::MemoryControl(size_t size_tmp){
	if(Allocated && size!= size_tmp) Free();
	if(!Allocated) Allocate(size_tmp);
}
inline void FFT2base::Allocate(size_t size_tmp){
	size = size_tmp;
	if(!Allocated) f = new double[size];
	Allocated = true;
}
inline void FFT2base::Free(){
	if(Allocated) delete[] f;
	Allocated = false;
	size = 0;
}
inline void FFT2base::set_zero(){
	for(size_t i=0; i<size; ++i) f[i]=0.0;
}
inline void FFT2base::set(const FFT2base& src){
	for(size_t i=0; i<size; ++i) f[i]=src.f[i];
}
inline void FFT2base::div(double val){
	for(size_t i=0; i<size; ++i) f[i] /= val;
}

//inline void FFT2base::div(const FFT2& src){}
inline void FFT2base::divReal(const FFT2base& src){
	f[0] /= src.f[0];
	for(size_t i=1; i<size/2; ++i){
		f[     i] /= src.f[i];
		f[size-i] /= src.f[i];
	}
	f[size/2] /= src.f[size/2];
}
inline void FFT2base::divAmp(const FFT2base& src){
	f[0] /= src.f[0];
	double amp;
	for(size_t i=1; i<size/2; ++i){
		amp =  sqrt(src.f[i]*src.f[i] + src.f[size-i]*src.f[size-i]);
		f[     i] /= amp;
		f[size-i] /= amp;
	}
	f[size/2] /= src.f[size/2];
}
inline void FFT2base::add(const FFT2base& src){
	for(size_t i=0; i<size; ++i) f[i] += src.f[i];
}
inline void FFT2base::addReal(const FFT2base& src){
	for(size_t i=0; i<=size/2; ++i) f[i] += src.f[i];
}
inline void FFT2base::mul(const double val){
	for(size_t i=0; i<size; ++i) f[i] *= val;
}

inline void FFT2base::mulAmp(const FFT2base& src){	//	src0.f(conjugate()) * src1.f
	f[0] *= src.f[0];
	for(size_t i=1; i<size/2; ++i)
		f[     i] *= sqrt(src.f[i]*src.f[i] + src.f[size-i]*src.f[size-i]);
	f[size/2] *= src.f[size/2];
}

inline void FFT2base::addmul(const FFT2base& src0, const FFT2base& src1){	//	src0.f(conjugate()) * src1.f
//	f[0] += src0.f[0]*src1.f[0];
	for(size_t i=1; i<size/2; ++i){
		f[     i] += src0.f[     i]*src1.f[i] + src0.f[size-i]*src1.f[size-i];
		f[size-i] +=-src0.f[size-i]*src1.f[i] + src0.f[     i]*src1.f[size-i];
	}
//	f[size/2] += src0.f[size/2]*src1.f[size/2];
}
inline void FFT2base::addmuldagger(const FFT2base& src0, const FFT2base& src1){	//	src0.f(conjugate()) * src1.f
//	f[0] += src0.f[0]*src1.f[0];
	for(size_t i=1; i<size/2; ++i){
		f[     i] += src0.f[     i]*src1.f[i] - src0.f[size-i]*src1.f[size-i];
		f[size-i] +=-src0.f[size-i]*src1.f[i] + src0.f[     i]*src1.f[size-i];
	}
//	f[size/2] += src0.f[size/2]*src1.f[size/2];
}
inline void FFT2base::amp(const FFT2base& src){
	f[0] = src.f[0];
	for(size_t i=1; i<size/2; ++i)
		f[     i] = sqrt(src.f[i]*src.f[i] + src.f[size-i]*src.f[size-i]);
	f[size/2] = src.f[size/2];
}
inline void FFT2base::addamp(const FFT2base& src){
	f[0] += src.f[0];
	for(size_t i=1; i<size/2; ++i)
		f[     i] += sqrt(src.f[i]*src.f[i] + src.f[size-i]*src.f[size-i]);
	f[size/2] += src.f[size/2];
}
inline std::vector<double> FFT2base::power(){
	std::vector<double> result(size/2+1);
	result[0] = f[0]*f[0];
	for(size_t i=1; i<size/2; ++i)
		result[i] = f[i]*f[i] + f[size-i]*f[size-i];
	result[size/2] = f[size/2]*f[size/2];
	return result;
}
inline std::vector<double> FFT2base::amplitude(){
	std::vector<double> result(size/2+1);
	result[0] = f[0];
	for(size_t i=1; i<size/2; ++i)
		result[i] = sqrt(f[i]*f[i] + f[size-i]*f[size-i]);
	result[size/2] = f[size/2];
	return result;
}

inline std::vector<double> FFT2base::phase(){
	std::vector<double> result(size/2+1);
	result[0] = 0.0;
	for(size_t i=1; i<size/2; ++i){
		result[i] = atan(f[size-i]/f[i]);
		if(f[i] < 0.0 && f[size-i] > 0.0) result[i] += M_PI;
		if(f[i] < 0.0 && f[size-i] < 0.0) result[i] -= M_PI;
	}
	result[size/2] = 0.0;
	return result;
}

//-------------------------------------------//
//											 //
//	core class for Fast Fourier Transform	 //
//											 //
//-------------------------------------------//
//---	Aplication function		---
//---	Power Spectrum	---
template< class Tyv_ >
inline void FFT2::Power(Tyv_& power){
	double T = rate/size;
	power[0] = f[0]*f[0];
	for(size_t i=1; i<Size()/2; ++i)
		power[i] = f[i]*f[i]/T + f[Size()-i]*f[Size()-i]/T;
	power[size/2] = f[Size()/2]*f[Size()/2]/T;
}

inline std::vector<double> FFT2::Power(){
	std::vector<double> amp(Size()/2+1);
	Power(amp);
	return amp;
}
inline std::vector<double> FFT2::Phase(){
}
//---	Cross Power Spectrum	---
template< class Tyv_ >
inline void FFT2::CrossPower(Tyv_& crosspower, const FFT2& src){
	double T = rate/size;
	crosspower[0] = f[0]*src.f[0]/T;
	double real, imag;
	for(size_t i=1; i<size/2; ++i){
		real = f[i]*src.f[i] + f[size-i]*src.f[size-i];
		imag = f[i]*src.f[size-i] - f[size-i]*src.f[i];
		crosspower[i] = sqrt(real*real + imag*imag)/T;
	}
	crosspower[size/2] = f[size/2]*src.f[size/2]/T;
}
//---	Coherency	---
template< class Tyv_ >
inline void FFT2::Coherency(Tyv_& coherency, FFT2& src){
	double T = rate/size;
	coherency[0] = 1.0;
	double real, amp;
	for(size_t i=1; i<size/2; ++i){
		real = f[i]*src.f[i] + f[size-i]*src.f[size-i];
		amp = sqrt((f[i]*f[i]+f[size-i]*f[size-i])*(src.f[i]*src.f[i]+src.f[size-i]*src.f[size-i]));
		coherency[i] = real/amp;
	}
	coherency[size/2] = 1.0;
}
//---	CrossPhase	---
template< class Tyv_ >
inline void FFT2::CrossPhase(Tyv_& crossphase, const FFT2& src){
	double T = rate/size;
	crossphase[0] = 0.0;
	double real, imag;
	for(size_t i=1; i<size/2; ++i){
		real =  f[i]*src.f[i] + f[size-i]*src.f[size-i];
		imag =  f[i]*src.f[size-i] - f[size-i]*src.f[i];
		crossphase[i] = atan(imag/real);
	}
	crossphase[size/2] = 0.0;
}


//---	core functions	---
template < class Tyv_ >
inline void FFT2::get(std::vector< Tyv_ >& stdv_raw_data, double rate_, FFTwindow window){
	get(stdv_raw_data, stdv_raw_data.size(), rate_, window);
}

template < class Tyv_ >
void FFT2::get(Tyv_& raw_data, size_t size_, double rate_, FFTwindow window){
	set(raw_data, size_, rate_, window);
}
template < class Tyv_ >
inline void FFT2::set(Tyv_& raw_data, size_t size_, double rate_, FFTwindow window){
	set(raw_data, size_, 0, rate_, window);
}
template < class Tyv_ >
inline void FFT2::set(Tyv_& raw_data, size_t size_, size_t offset, double rate_, FFTwindow window){
	rate = rate_;
	size_t size_tmp = getPower2NumFloor(size_);
	MemoryControl(size_tmp);
	//---	copy data to the pointer	---
	for(size_t i=0; i<size; ++i)	f[i] = 1.0*raw_data[i+offset];
	//---	applying the window functions	---
	if(window == Hann)		HannWindow();
	else if(window == Hamming)	HammingWindow();
	else if(window == Blackman)	BlackmanWindow();
//	else if(window == Rectangular)	RectangularWindow();
	//	fft
	gsl_fft_real_radix2_transform(f, 1, size);
	for(size_t i=0; i<size; ++i) f[i] *= rate;
}



//----------------------------------------------------------
//
//			FFT2 correlation class
//
//----------------------------------------------------------
/*inline void set

template < class Tyv0_, class Tyv1_ >
inline void FFT2Cor::get(std::vector< Tyv0_ >& src0, std::vector< Tyv1_>& src1, double rate_, size_t avgnum, FFT2::FFTwindow window){
	rate =rate_;
	size_t size_all = src0.size();
	size_t size_fft = FFT2::getPower2NumFloor(size_all);
	while(size_all - avgnum < size_fft)
		size_fft /=2;
	size_t offset = (size_t)floor(1.0*(size_all - size_fft)/avgnum);

	Allocate(size_fft);
	std::vector< Tyv0_ > src0_sub(size_fft);
	std::vector< Tyv1_ > src1_sub(size_fft);

	set_zero();
	for(size_t i=0; i<avgnum; ++i){
		for(size_t j=0; j<size_fft; ++j){
			src0_sub[j] = src0[offset*i+j];
			src1_sub[j] = src1[offset*i+j];
		}
		//	core fft 
		fft0.get(src0_sub, rate, window);
		fft1.get(src1_sub, rate, window);
		//	summation for average
		addCorrelation(fft0, fft1);
	}
	div(1.0*avgnum);
//	div(fft0);
//	div(fft1);
}

inline void FFT2Cor::add(FFT2& src0, FFT2& src1){
	double T = rate/size;
	f[0] += src0.f[0]*src1.f[0]/T;
	for(size_t i=1; i<size/2; ++i){
		f[0]      += src0.f[i]*src1.f[i]/T      + src0.f[size-i]*src1.f[size-i]/T;
		f[size-i] += src0.f[i]*src1.f[size-i]/T - src0.f[size-i]*src1.f[i]/T;
	}
	f[size/2] += src0.f[0]*src1.f[0]/T;
}

inline void FFT2Cor::addCorrelation(FFT2& src0, FFT2& src1){
	f[0] += 1.0;
	for(size_t i=1; i<size/2; ++i){
		double amp = sqrt((src0.f[i]*src0.f[i]+src0.f[size-1]*src0.f[size-1])*
						  (src1.f[i]*src1.f[i]+src1.f[size-1]*src1.f[size-1]));
		f[0]      += (src0.f[i]*src1.f[i]      + src0.f[size-i]*src1.f[size-i])/amp;
		f[size-i] += (src0.f[i]*src1.f[size-i] - src0.f[size-i]*src1.f[i])/amp;
	}
	f[size/2] += 1.0;
}

*/

//----------------------------------------------------------------------//
//																		//
//	high level class for short time FFT, and its time development		//
//																		//
//----------------------------------------------------------------------//
inline FFT2development::FFT2development(void){
	avgnum = 1;
	fft_size=0;
	fft_num =0;
}
inline void FFT2development::Free(){
	CrossPower.Free();
	Amplitude0.Free();
	Amplitude1.Free();
	CrossPhase.Free();
	fftv.clear();
}

inline FFT2development::~FFT2development(void){
	CrossPower.Free();
	Amplitude0.Free();
	Amplitude1.Free();
	CrossPhase.Free();
}

template < class Tyv0_ >
inline void FFT2development::set
(Tyv0_& src0, size_t data_size, size_t fft_size_, size_t fft_num_, double rate_, FFT2::FFTwindow window)
{
	rate = rate_;
	// number of data point for fft2 
	fft_size = FFT2::getPower2Num(fft_size_);
	//	number of fft2 analysis
	fft_num = fft_num_;
	if(fft_num > data_size-fft_size) fft_num = data_size-fft_size;
	fftv.resize(fft_num);
	//	each step for short time fft2 
	increment = (size_t)floor(1.0*(data_size-fft_size)/fft_num);
	for(size_t i=0; i<fft_num; ++i)	fftv[i].set(src0, fft_size, increment*i, rate, window);

	//	1st fft timing(point)
	start_p = fft_size/2;
}
inline std::vector<std::vector<double>> FFT2development::getPower(size_t avgnum){
	size_t timing_num = (size_t)floor(1.0*fft_num/avgnum);
	std::vector<std::vector<double>> result(timing_num);

	Amplitude0.MemoryControl(fft_size);
	for(size_t t=0; t<timing_num; ++t){
		Amplitude0.set_zero();
		for(size_t t1=0; t1< avgnum; ++t1)
			Amplitude0.addamp(fftv[t*avgnum+t1]);
		Amplitude0.mul(1.0/avgnum);
		result[t] = Amplitude0.power();
	}
	return result;
}

inline std::vector<std::vector<double>> FFT2development::getCrossPower
(FFT2development& src, size_t avgnum)
{
	CrossPower.MemoryControl(fft_size);
	Amplitude0.MemoryControl(fft_size);
	Amplitude1.MemoryControl(fft_size);

	size_t timing_num = (size_t)floor(1.0*fft_num/avgnum);
	std::vector<std::vector<double>> result(timing_num);
	if(fft_num != src.FFTnum() || fft_size != src.FFTsize()) return result;

	for(size_t t=0; t<timing_num; ++t){
		CrossPower.set_zero();
		Amplitude1.set_zero();
		for(size_t t1=0; t1< avgnum; ++t1){
			CrossPower.addmul(fftv[t*avgnum+t1], src.fftv[t*avgnum+t1]);
			Amplitude0.set_zero();
			Amplitude0.amp(fftv[t*avgnum+t1]);
//			Amplitude0.mulAmp(    fftv[t*avgnum+t1]);
			Amplitude0.mulAmp(src.fftv[t*avgnum+t1]);
			Amplitude1.addReal(Amplitude0);
		}
		CrossPower.divReal(Amplitude1);
		result[t]=CrossPower.power();
	}
	return result;
}

inline std::vector<std::vector<double>> FFT2development::getCrossPhase(FFT2development& src, size_t avgnum_){
	avgnum = avgnum_;
	CrossPhase.MemoryControl(fft_size);

	size_t timing_num = (size_t)floor(1.0*fft_num/avgnum);
	std::vector<std::vector<double>> result(timing_num);
	if(fft_num != src.FFTnum() || fft_size != src.FFTsize()) return result;

	for(size_t t=0; t<timing_num; ++t){
		CrossPhase.set_zero();
		for(size_t t1=0; t1< avgnum; ++t1){
			CrossPhase.addmuldagger(fftv[t*avgnum+t1], src.fftv[t*avgnum+t1]);
		}
		result[t]=CrossPhase.phase();
	}
	return result;
}

inline std::vector<double> FFT2development::getFrequency(){
	std::vector<double> frequency(fft_size/2+1);
	for(size_t i=0; i<fft_size/2+1; ++i)
		frequency[i] = rate*i/fft_size;
	return frequency;
}

inline std::vector<double> FFT2development::getFFTtiming(double offset){
	std::vector<double> timing((size_t)floor(1.0*fft_num/avgnum));
	for(size_t i=0; i<timing.size(); ++i)
		timing[i] = offset+(1.0*start_p+1.0*increment*avgnum/2.0*i)/rate;
	return timing;
}
