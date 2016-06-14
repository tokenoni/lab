//----------------------------------------------------------
//
//			FFT2 Core class
//
//----------------------------------------------------------
//---	Aplication function		---
//---	Power Spectrum	---
template< class Tyv_ >
inline void FFT2::Power(Tyv_& power){
	double T = rate/num;
	power[0] = f[0]*f[0];
	for(size_t i=1; i<num/2; ++i)
		power[i] = f[i]*f[i]/T + f[num-i]*f[num-i]/T;
	power[num/2] = f[num/2]*f[num/2]/T;
}

inline std::vector<double> FFT2::Power(){
	std::vector<double> amp(num/2+1);
	Power(amp);
	return amp;
}
inline std::vector<double> FFT2::Phase(){
}
//---	returnning the complex value ---
inline double FFT2::real(size_t i){
	if(i<=num/2)	return f[i]*num/rate;
	else			return f[num-i]*num/rate;
};
inline double FFT2::imag(size_t i){
	if(i==0 || i==num/2) return 0.0;
	if(i<num/2)	return f[num-i]*num/rate;
	else			return -f[i]*num/rate;
};


//---	Cross Power Spectrum	---
template< class Tyv_ >
inline void FFT2::CrossPower(Tyv_& crosspower, const FFT2& src){
	double T = rate/num;
	crosspower[0] = f[0]*src.f[0]/T;
	double real, imag;
	for(size_t i=1; i<num/2; ++i){
		real = f[i]*src.f[i] + f[num-i]*src.f[num-i];
		imag = f[i]*src.f[num-i] - f[num-i]*src.f[i];
		crosspower[i] = sqrt(real*real + imag*imag)/T;
	}
	crosspower[num/2] = f[num/2]*src.f[num/2]/T;
}
//---	Coherency	---
template< class Tyv_ >
inline void FFT2::Coherency(Tyv_& coherency, FFT2& src){
	double T = rate/num;
	coherency[0] = 1.0;
	double real, amp;
	for(size_t i=1; i<num/2; ++i){
		real = f[i]*src.f[i] + f[num-i]*src.f[num-i];
		amp = sqrt((f[i]*f[i]+f[num-i]*f[num-i])*(src.f[i]*src.f[i]+src.f[num-i]*src.f[num-i]));
		coherency[i] = real/amp;
	}
	coherency[num/2] = 1.0;
}
//---	CrossPhase	---
template< class Tyv_ >
inline void FFT2::CrossPhase(Tyv_& crossphase, const FFT2& src){
	double T = rate/num;
	crossphase[0] = 0.0;
	double real, imag;
	for(size_t i=1; i<num/2; ++i){
		real =  f[i]*src.f[i] + f[num-i]*src.f[num-i];
		imag =  f[i]*src.f[num-i] - f[num-i]*src.f[i];
		crossphase[i] = atan(imag/real);
	}
	crossphase[num/2] = 0.0;
}


inline void FFT2::set_zero(){
	for(size_t i=0; i<num; ++i) f[i]=0.0;
}
inline void FFT2::div(double val){
	for(size_t i=0; i<num; ++i) f[i] /= val;
}

inline void FFT2::div(const FFT2& src){
	double T = rate/num;
	f[0] /= sqrt(src.f[0]*src.f[0]/T);
	for(size_t i=1; i<num/2; ++i)
		f[i] /= sqrt(src.f[i]*src.f[i]/T + src.f[num-i]*src.f[num-i]/T);
	f[num/2] /= sqrt(src.f[num/2]*src.f[num/2]/T);
}
//---	core functions	---
template < class Tyv_ >
inline void FFT2::set(std::vector< Tyv_ >& stdv_raw_data, double rate_, FFTwindow window){
	rate = rate_;
	size_t num_tmp = getPower2NumFloor(stdv_raw_data.size());
	Allocate(num_tmp);
	//---	copy data to the pointer	---
	for(size_t i=0; i<num; ++i)	f[i] = 1.0*stdv_raw_data[i];
	//---	applying the window functions	---
	if(window == Hann)		HannWindow();
	else if(window == Hamming)	HammingWindow();
	else if(window == Blackman)	BlackmanWindow();
//	else if(window == Rectangular)	RectangularWindow();
	//	fft
	gsl_fft_real_radix2_transform(f, 1, num);
	for(size_t i=0; i<num; ++i) f[i] *= rate;
}

template < class Tyv_ >
inline void FFT2::set( Tyv_ & raw_data, size_t num_, double rate_, FFTwindow window){
	rate = rate_;
	size_t num_tmp = getPower2NumFloor(num_);
	Allocate(num_tmp);
	//---	copy data to the pointer	---
	for(size_t i=0; i<num; ++i)	f[i] = 1.0*raw_data[i];
	//---	applying the window functions	---
	if(window == Hann)		HannWindow();
	else if(window == Hamming)	HammingWindow();
	else if(window == Blackman)	BlackmanWindow();
//	else if(window == Rectangular)	RectangularWindow();
	//	fft
	gsl_fft_real_radix2_transform(f, 1, num);
	for(size_t i=0; i<num; ++i) f[i] *= rate;
}
//---	check the memory allocation	---
inline void FFT2::Allocate(size_t& num_tmp){
	if(!Allocated){
		f = new double[num_tmp];
		Allocated = true;
	}else if(num != num_tmp){
		delete[] f;
		f = new double[num_tmp];
	}
	num = num_tmp;
}
//---	constructor and destructor	---
inline FFT2::FFT2(void){	Allocated = false;}
inline FFT2::~FFT2(void){	if(Allocated) delete[] f;}


//----------------------------------------------------------
//
//			class for Fourier transform for complex values
//
//----------------------------------------------------------
template < class Tyv_ >
inline void FFT2Complex::set(std::vector< Tyv_ >& stdv_real, std::vector< Tyv_ >& stdv_imag, double rate, FFTwindow window){
	set(stdv_real, rate, window);
	fft2imag.set(stdv_imag, rate, window);
}

template < class Tyv_ >
inline void FFT2Complex::set( Tyv_ & raw_real, Tyv_ & raw_imag, size_t num_, double rate, FFTwindow window){
	FFT2::set(raw_real, num_, rate, window);
	fft2imag.set(raw_imag, num_, rate, window);
}

inline double FFT2Complex::real(size_t i){	return FFT2::real(i)-fft2imag.imag(i);};
inline double FFT2Complex::imag(size_t i){	return FFT2::imag(i)+fft2imag.real(i);};

//----------------------------------------------------------
//
//			FFT2 correlation class
//
//----------------------------------------------------------


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
	for(size_t i=0; i<size(); ++i){
		f[i]		= src0.real(i)*src1.real(i) + src0.imag(i)*src1.imag(i);
		f[size()-i] = src0.real(i)*src1.imag(i) - src0.imag(i)*src1.real(i);
	}
}

inline void FFT2Cor::addCorrelation(FFT2& src0, FFT2& src1){
	f[0] += 1.0;
	for(size_t i=1; i<size()/2; ++i){
		double amp = sqrt((src0.f[i]*src0.f[i]+src0.f[size()-1]*src0.f[size()-1])*
						  (src1.f[i]*src1.f[i]+src1.f[size()-1]*src1.f[size()-1]));
		f[0]      += (src0.f[i]*src1.f[i]      + src0.f[size()-i]*src1.f[size()-i])/amp;
		f[size()-i] += (src0.f[i]*src1.f[size()-i] - src0.f[size()-i]*src1.f[i])/amp;
	}
	f[size()/2] += 1.0;
}

