//-----------------------------------------------------------------//
//           		constructor and destructors					   //
//-----------------------------------------------------------------//
inline InterpolatedData::InterpolatedData(void){	allocated = false;	num = 0; }

inline InterpolatedData::~InterpolatedData(void){	FreeMemory();}

inline InterpolatedData::InterpolatedData(const InterpolatedData& src){
	if(!src.allocated) allocated = false;
	else set(src.x, src.y, src.num, src.Type);
}

inline InterpolatedData InterpolatedData::operator=(const InterpolatedData& src){
	InterpolatedData dest(src);
	return dest;
}

//-----------------------------------------------------------------//
//           		memory controlling functions				   //
//-----------------------------------------------------------------//
inline void InterpolatedData::MemoryControl(size_t size_, InterpType Type_){
	if(allocated && size_ != num)  FreeMemory();
	if(allocated && Type_ != Type) FreeMemory();
	num = size_;
	Type = Type_;
	if(!allocated) Allocate(num, Type);
};

inline void InterpolatedData::Allocate(size_t size_, InterpType Type_){
	num = size_;
	if(!allocated){
		x = new double [size_];
		y = new double [size_];
		acc = gsl_interp_accel_alloc();

		if(Type_ == Linear)	
			interp = gsl_interp_alloc(gsl_interp_linear, num);
		else if(Type_ == CubicSpline) 
			interp = gsl_interp_alloc(gsl_interp_cspline, num);
	}
	allocated = true;
}

inline void InterpolatedData::FreeMemory(){
	if(allocated){
		delete [] x;	delete [] y;
		gsl_interp_free(interp);
		gsl_interp_accel_free(acc);
	}
	allocated = false;
}

//-----------------------------------------------------------------//
//   		initializing functions for inteprolation			   //
//-----------------------------------------------------------------//
template < class Tyx_, class Tyy_ >
inline void InterpolatedData::set(std::vector< Tyx_ >& xsrc, std::vector< Tyy_ >& ysrc, InterpType type){
	set(xsrc, ysrc, xsrc.size(), type);
}

template < class Tyx_, class Tyy_ >
inline void InterpolatedData::set(Tyx_& xsrc, Tyy_& ysrc, size_t size_, InterpType type){
	MemoryControl(size_, type);
	for(size_t i=0; i<num; ++i){
		x[i] = xsrc[i];	y[i] = ysrc[i];
	}
	gsl_interp_accel_reset(acc);
	gsl_interp_init(interp, x, y, num);
}


//-----------------------------------------------------------------//
//		   		obtain the interpolated data					   //
//-----------------------------------------------------------------//
inline double InterpolatedData::get(double x0){
	return gsl_interp_eval(interp, x, y, x0, acc);
};
inline double InterpolatedData::getderiv(double x0){
	return gsl_interp_eval_deriv(interp, x, y, x0, acc);
}
inline double InterpolatedData::getderiv2(double x0){
	return gsl_interp_eval_deriv2(interp, x, y, x0, acc);
}

inline double InterpolatedData::getinteg(double xfrom, double xto){
	return gsl_interp_eval_integ(interp, x, y, xfrom, xto, acc);
}