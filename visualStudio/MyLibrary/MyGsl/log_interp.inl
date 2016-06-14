namespace mygsl{
//-------------------------------------------------//
//												   //
//---		1 dimensional interpolation			---//
//							in log(x) scale		   //
//												   //
//-------------------------------------------------//
	template < class Tyv_ >
	inline void interp_logx::set(const Tyv_& x, const Tyv_& y, InterpMethod method_){
		setSize(x.size(), method_);
		for(size_t i=0; i<x.size(); ++i)	set(i, x[i], y[i]);
		setInterp();
	}

	inline void interp_logx::set(const size_t i, const double xi, const double yi){
		interp::set(i, log(xi), yi);	
	}
	inline void interp_logx::set_x(const size_t i, const double xi){
		interp::set_x(i, log(xi));
	}

	inline std::vector<double> interp_logx::get_original_x(){
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = exp(x_array[i]);
		return tmp;
	}
	inline double interp_logx::get(double x)      {
		return interp::get(log(x));
	};
	inline double interp_logx::getderiv(double x) {
		return interp::getderiv(log(x))/x;
	};
	inline double interp_logx::getderiv2(double x){
		double x_ = log(x);
		return (-interp::getderiv(x_) + interp::getderiv2(x_))/(x*x);
	}
	inline double interp_logx::getinteg(double a, double b){
		switch(method){
		case Linear:
			if(     a < b)	 return  getinteg_linear(log(a), log(b));
			else if(b > a)   return -getinteg_linear(log(b), log(a));
			else return 0.0;
			break;
		case Cspline:
			break;
			if(     a < b)	 return  getinteg_spline(log(a), log(b));
			else if(b > a)   return -getinteg_spline(log(b), log(a));
			else return 0.0;
		}
	}

	inline double interp_logx::getinteg_linear(double a, double b){
		double val = 0.0;
		double exp_a = exp(a), exp_b = exp(b);
		val -= exp_a*interp::get(a);
		
		size_t i, imin = gsl_interp_accel_find(acc, x_array, num, a);
		double exp_x1 = exp(x_array[imin]);
		val += (exp_x1-exp_a) * interp::getderiv(0.5*x_array[imin] + 0.5*a);

		double exp_x0;
		for(i=imin; x_array[i+1] < b && i+1 < num ; ++i){
			exp_x0 = exp_x1;
			exp_x1 = exp(x_array[i+1]);
			val += (exp_x1 - exp_x0) * interp::getderiv(0.5*x_array[i]+0.5*x_array[i+1]);
		}

		val += (exp_b - exp_x1) * interp::getderiv(0.5*x_array[i] + 0.5*b);
		val += exp_b * interp::get(b);
		return val;
	}

	inline double interp_logx::getinteg_spline(double a, double b){
		double val = 0.0;
		double exp_a = exp(a), exp_b = exp(b);
		val -= exp_a*(interp::get(a) + interp::getderiv(a));
		
		size_t i, imin = gsl_interp_accel_find(acc, x_array, num, a);
		double exp_x1 = exp(x_array[imin]);
		val += (exp_x1-exp_a) * interp::getderiv2(0.5*x_array[imin] + 0.5*a);

		double exp_x0;
		for(i=imin; x_array[i+1] < b || i+1 < num ; ++i){
			exp_x0 = exp_x1;
			exp_x1 = exp(x_array[i+1]);
			val += (exp_x1 - exp_x0) * interp::getderiv2(0.5*x_array[i]+0.5*x_array[i+1]);
		}

		val += (exp_b - exp_x1) * interp::getderiv2(0.5*x_array[i] + 0.5*b);
		val += exp_b * (interp::get(b) + interp::getderiv(b));
		return val;
	}

//-------------------------------------------------//
//							in log(y) scale		   //
//-------------------------------------------------//
	template < class Tyv_ >
	inline void interp_logy::set(const Tyv_& x, const Tyv_& y, InterpMethod method_){
		setSize(x.size(), method_);
		for(size_t i=0; i<x.size(); ++i)	set(i, x[i], log(y[i]));
		setInterp();
	}

	inline void interp_logy::set(const size_t i, const double xi, const double yi){
		interp::set(i, xi, log(yi));	
	}
	inline void interp_logy::set_y(const size_t i, const double yi){
		interp::set_y(i, log(yi));
	}

	inline std::vector<double> interp_logy::get_original_y(){
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = exp(data[i]);
		return tmp;
	}

	inline double interp_logy::get(double x){
		return exp(interp::get(x));
	}

	inline double interp_logy::getderiv(double x){
		return get(x)*interp::getderiv(x);
	}

	inline double interp_logy::getderiv2(double x){
		double dy_dx = interp::getderiv(x);
		return get(x)*(dy_dx*dy_dx + interp::getderiv2(x));
	}

	inline double interp_logy::getinteg(double a, double b){
		switch(method){
		case Linear:
			if(     a < b)	 return  getinteg_linear(log(a), log(b));
			else if(b > a)   return -getinteg_linear(log(b), log(a));
			else return 0.0;
			break;
		case Cspline:
			//---	now under construction	----
			break;
			if(     a < b)	 return  getinteg_spline(log(a), log(b));
			else if(b > a)   return -getinteg_spline(log(b), log(a));
			else return 0.0;
		}
	}

	inline double interp_logy::getinteg_linear(double a_, double b_){
		size_t i = gsl_interp_accel_find(acc, x_array, size(), a_);
		double value = 0.0;
		double c[4];
		while(i<size()-1 && x_array[i+1] < b_){
			getCoef(i,c);
			value += (exp(c[0]*x_array[i+1]) - exp(c[0]*a_))*exp(c[1])/c[0];
			++i;
			a_ = x_array[i+1];
		}
		getCoef(i,c);
		value += (exp(c[0]*b_) - exp(c[0]*a_))*exp(c[1])/c[0];
	}
//-------------------------------------------------//
//					in log(x)-log(y) scale		   //
//-------------------------------------------------//
	template < class Tyv_ >
	inline void interp_logxy::set(const Tyv_& x, const Tyv_& y, InterpMethod method_){
		setSize(x.size(), method_);
		for(size_t i=0; i<x.size(); ++i)	set(i, log(x[i]), log(y[i]));
		setInterp();
	}

	inline void interp_logxy::set(const size_t i, const double xi, const double yi){
		interp::set(i, log(xi), log(yi));	
	}
	inline void interp_logxy::set_x(const size_t i, const double xi){
		interp::set_x(i, log(xi));
	}
	inline void interp_logxy::set_y(const size_t i, const double yi){
		interp::set_y(i, log(yi));
	}
	inline std::vector<double> interp_logxy::get_original_x(){
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = exp(x_array[i]);
		return tmp;
	}
	inline std::vector<double> interp_logxy::get_original_y(){
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = exp(data[i]);
		return tmp;
	}

	inline double interp_logxy::get(double x){
		return exp(interp::get(log(x)));
	}

/*		they are under construction	
inline double interp_logy::getderiv(double x){
		return get(x)*interp::getderiv(x);
	}

	inline double interp_logy::getderiv2(double x){
		double dy_dx = interp::getderiv(x);
		return get(x)*(dy_dx*dy_dx + interp::getderiv2(x));
	}

	inline double interp_logy::getinteg(double a, double b){
		switch(method){
		case Linear:
			if(     a < b)	 return  getinteg_linear(log(a), log(b));
			else if(b > a)   return -getinteg_linear(log(b), log(a));
			else return 0.0;
			break;
		case Cspline:
			//---	now under construction	----
			break;
			if(     a < b)	 return  getinteg_spline(log(a), log(b));
			else if(b > a)   return -getinteg_spline(log(b), log(a));
			else return 0.0;
		}
	}

	inline double interp_logy::getinteg_linear(double a_, double b_){
		size_t i = gsl_interp_accel_find(acc, x_array, size(), a_);
		double value = 0.0;
		double c[4];
		while(i<size()-1 && x_array[i+1] < b_){
			getCoef(i,c);
			value += (exp(c[0]*x_array[i+1]) - exp(c[0]*a_))*exp(c[1])/c[0];
			++i;
			a_ = x_array[i+1];
		}
		getCoef(i,c);
		value += (exp(c[0]*b_) - exp(c[0]*a_))*exp(c[1])/c[0];
	}
	*/
//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//							in log(x) scale		   //
//												   //
//-------------------------------------------------//
	inline void interp2d_logx::set(const std::vector<double>& x_, const std::vector<double>& y_,
							  const std::vector<std::vector<double>>& data_, InterpMethod method_){
		std::vector<double> logx(x_.size());
		for(size_t i=0; i<x_.size(); ++i)
			logx[i] = log(x_[i]);
		interp2d::set(logx, y_, data_, method_);
	}
	inline double interp2d_logx::get(const double x, const double y){
		return interp2d::get(log(x), y);
	}
	inline double interp2d_logx::getderiv_x(const double x, const double y){
		return interp2d::getderiv_x(log(x), y)/x;
	}
	inline double interp2d_logx::getderiv_xy(const double x, const double y){
		return interp2d::getderiv_xy(log(x), y)/x;
	}
	inline interp_logx interp2d_logx::getAlong_x_log(const double y, interp::InterpMethod method_){
		interp_logx interp1d(size1(), method_);
		for(size_t i=0; i<size1(); ++i)
			interp1d.set(i, xv[i], get(xv[i], y));
		interp1d.setInterp();
		return interp1d;
	}

//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//							in log(y) scale		   //
//												   //
//-------------------------------------------------//
	inline void interp2d_logy::set(const std::vector<double>& x_, const std::vector<double>& y_,
							  const std::vector<std::vector<double>>& data_, InterpMethod method_){
		std::vector<double> logy(y_.size());
		for(size_t i=0; i<y_.size(); ++i)
			logy[i] = log(y_[i]);
		interp2d::set(x_, logy, data_, method_);
	}
	inline double interp2d_logy::get(const double x, const double y){
		return interp2d::get(x, log(y));
	}
	inline double interp2d_logy::getderiv_y(const double x, const double y){
		return interp2d::getderiv_y(x, log(y))/y;
	}
	inline double interp2d_logy::getderiv_xy(const double x, const double y){
		return interp2d::getderiv_xy(x, log(y))/y;
	}
	inline interp_logx interp2d_logy::getAlong_y_log(const double x, interp::InterpMethod method_){
		interp_logx interp1d(size2(), method_);
		for(size_t i=0; i<size2(); ++i)
			interp1d.set(i, yv[i], get(x, yv[i]));
		interp1d.setInterp();
		return interp1d;
	}
//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//          		in log(x)-log(y) scale		   //
//												   //
//-------------------------------------------------//
	inline void interp2d_logxy::set(const std::vector<double>& x_, const std::vector<double>& y_,
							  const std::vector<std::vector<double>>& data_, InterpMethod method_){
		std::vector<double> logx(x_.size());
		for(size_t i=0; i<x_.size(); ++i)	logx[i] = log(x_[i]);
		std::vector<double> logy(y_.size());
		for(size_t i=0; i<y_.size(); ++i)	logy[i] = log(y_[i]);
		interp2d::set(logx, logy, data_, method_);
	}
	inline double interp2d_logxy::get(const double x, const double y){
		return interp2d::get(log(x), log(y));
	}
	inline double interp2d_logxy::getderiv_x(const double x, const double y){
		return interp2d::getderiv_y(log(x), log(y))/x;
	}
	inline double interp2d_logxy::getderiv_y(const double x, const double y){
		return interp2d::getderiv_y(log(x), log(y))/y;
	}
	inline double interp2d_logxy::getderiv_xy(const double x, const double y){
		return interp2d::getderiv_xy(log(x), log(y))/(x*y);
	}
	inline interp_logx interp2d_logxy::getAlong_x_log(const double y, interp::InterpMethod method_){
		interp_logx interp1d(size1(), method_);
		for(size_t i=0; i<size1(); ++i)
			interp1d.set(i, xv[i], get(xv[i], y));
		interp1d.setInterp();
		return interp1d;
	}
	inline interp_logx interp2d_logxy::getAlong_y_log(const double x, interp::InterpMethod method_){
		interp_logx interp1d(size2(), method_);
		for(size_t i=0; i<size2(); ++i)
			interp1d.set(i, yv[i], get(x, yv[i]));
		interp1d.setInterp();
		return interp1d;
	}
};