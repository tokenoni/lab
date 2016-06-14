//------------------	definitions	------------------

namespace mygsl{
//-------------------------------------------------//
//												   //
//---		1 dimensional interpolation			---//
//												   //
//-------------------------------------------------//

	inline interp::interp(void){
		allocated = false;
	}
	inline interp::interp(const size_t size_, InterpMethod method_){
		method = method_;
		allocated = false;
		Allocate(size_);
	}

	inline interp::interp(const interp&obj){
		allocated = false;
		FreeMemory();
		if(obj.isallocated()){
			method = obj.method;
			Allocate(obj.size());
			for(size_t i=0; i<obj.size(); ++i){
				x_array[i] = obj.x_array[i];	data[i]= obj.data[i];
			}
			if(obj.issetInterped())	setInterp();
		}
	}
	inline interp& interp::operator=(const interp& obj){
		FreeMemory();
		method = obj.method;
		if(obj.isallocated()){
			Allocate(obj.size());
			for(size_t i=0; i<obj.size(); ++i){
				x_array[i] = obj.x_array[i];	data[i]= obj.data[i];
			}
			if(obj.issetInterped())	setInterp();
		}
		return *this;
	}
	inline interp::~interp(void){FreeMemory();	}

	template < class Tyv_ >
	inline void interp::set(const Tyv_& x, const Tyv_& y, InterpMethod method_){
		method = method_;
		FreeMemory();
		Allocate(x.size());

/*		size_t *p = new size_t[num];
		for(size_t i=0; i<num; ++i)	x_array[i] = x[i];
		gsl_sort_index(p, x_array, 1, num);
		for(size_t i=0; i<num; ++i){
			x_array[i] = x[p[i]];	data[i] = y[p[i]];
		}
		gsl_interp_accel_reset(acc);
		setInterp();
		delete [] p;
*/
		for(size_t i=0; i<num; ++i){
			x_array[i] = x[i];	data[i] = y[i];
		}
		sort();
		setInterp();
	}

	inline double interp::get(double x){
/*		if(method == CsplineFlatEdge){
			if(x<x_array[0]) return data[0];
			else if (x>x_array[size()-1]) return data[size()-1];
		}
*/		return gsl_interp_eval(interp_obj, x_array, data, x, acc);
	};
	inline double interp::get(double x, gsl_interp_accel* acc_obj)const{
		return gsl_interp_eval(interp_obj, x_array, data, x, acc_obj);
	};
	inline double interp::getderiv(double x){
		return gsl_interp_eval_deriv(interp_obj, x_array, data, x, acc);
	}
	inline double interp::getderiv(double x, gsl_interp_accel* acc_obj)const{
		return gsl_interp_eval_deriv(interp_obj, x_array, data, x, acc_obj);
	};
	inline double interp::getderiv2(double x){
		return gsl_interp_eval_deriv2(interp_obj, x_array, data, x, acc);
	}
	inline double interp::getderiv2(double x, gsl_interp_accel* acc_obj)const{
		return gsl_interp_eval_deriv2(interp_obj, x_array, data, x, acc_obj);
	};
	inline double interp::getinteg(double a, double b){
		if(a<b)	return gsl_interp_eval_integ(interp_obj, x_array, data, a, b, acc);
		if(b<a)	return -gsl_interp_eval_integ(interp_obj, x_array, data, b, a, acc);
		return 0.0;
	}
	inline double interp::getinteg(double a, double b, gsl_interp_accel* acc_obj)const{
		if(a<b)	return gsl_interp_eval_integ(interp_obj, x_array, data, a, b, acc_obj);
		if(b<a)	return -gsl_interp_eval_integ(interp_obj, x_array, data, b, a, acc_obj);
		return 0.0;
	}
	inline void interp::clear(){
		FreeMemory();
	}
	inline void interp::FreeMemory(){
		if(allocated){
			delete [] x_array;	delete [] data;	
			gsl_interp_free(interp_obj);
			gsl_interp_accel_free(acc);
		}
		allocated = false;
		setInterped = false;
	}

	inline void interp::setSize(const size_t size_, InterpMethod method_){
		if(isallocated() && (size() != size_ || method_ != method)){
			method = method_;
			FreeMemory();
			Allocate(size_);
		}
		if(!isallocated()){
			method = method_;
			Allocate(size_);
		}
	}
	inline void interp::set(const size_t i, const double xi, const double yi){
		x_array[i] = xi;
		data[i] = yi;
	}
	inline void interp::set_x(const size_t i, const double xi){
		x_array[i] = xi;
	}
	inline void interp::set_y(const size_t i, const double yi){
		data[i] = yi;
	}
	inline void interp::setInterp(bool isSorted){
		if(!isSorted) sort();

		gsl_interp_accel_reset(acc);
/*		if(method == CsplineFlatEdge)
			CsplineClampedInit(interp_obj, x_array, data, num,0.0, 0.0);
		else
*/		gsl_interp_init(interp_obj, x_array, data, num);
		
		setInterped = true;
	}
	inline void interp::sort(){
		double *xtmp = new double [num];	
		double *ytmp = new double [num];
		size_t *p = new size_t[num];
		
		for(size_t i=0; i<num; ++i){
			xtmp[i] = x_array[i];
			ytmp[i] = data[i];
		}
		gsl_sort_index(p, x_array, 1, num);
		for(size_t i=0; i<num; ++i){
			x_array[i] = xtmp[p[i]];	data[i] = ytmp[p[i]];
		}
		gsl_interp_accel_reset(acc);

		//---	checking number of same x values	---
		size_t count = 0;
		for(size_t i=0; i<num; ++i){
			xtmp[i] = x_array[i];
			ytmp[i] = data[i];
			if(i<num-1 && xtmp[i] == xtmp[i+1]) ++count; 
		}
		//---	averaging same-x values	---
		if(count !=0){
			FreeMemory();
			Allocate(num-count);
		
			size_t* number_of_summation = new size_t [num];
			for(size_t i=0; i<num; ++i){
				number_of_summation[i]=0;
				data[i]=0.0;
			}
			size_t index =0;
			for(size_t i=0; i<num; ++i){
				data[i] += ytmp[index];
				x_array[i] = xtmp[index];
				++number_of_summation[i];
				if(index<num+count-1 && xtmp[index] == xtmp[index+1]){
					--i;
				}
				++index;
			}
			for(size_t i=0; i<num; ++i)data[i] /= number_of_summation[i];
			delete[] number_of_summation;
		}
		delete [] p;
		delete [] xtmp;
		delete [] ytmp;
	}
	inline void interp::set_y_zero(){
		for(size_t i=0; i<size();++i) data[i]=0.0;
	}
	inline void interp::Allocate(size_t size_){
		num = size_;
		if(!allocated){
			x_array = new double [size_];
			data = new double [size_];
			acc = gsl_interp_accel_alloc();
			switch(method){
			case Cspline:
				interp_obj = gsl_interp_alloc(gsl_interp_cspline, num);
				break;
/*			case CsplineFlatEdge:
				interp_obj = gsl_interp_alloc(gsl_interp_cspline, num);
				break;
*/			case Linear:
				interp_obj = gsl_interp_alloc(gsl_interp_linear, num);
				break;
			}
		}
		allocated = true;
		setInterped = false;
	}

	inline std::vector<double> interp::get_original_x()const{
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = x_array[i];
		return tmp;
	}
	inline std::vector<double> interp::get_original_y()const{
		std::vector<double> tmp(size());
		for(size_t i=0; i<size(); ++i)
			tmp[i] = data[i];
		return tmp;
	}
/*
	//	merge the two x_arrayition vectors
	inline void interp::merge(const double* obj, const size_t size_obj){
		//	if the current x vector is wide
		if(obj[0] >= x_array[0] && obj[size_obj-1] <= x_array[size()-1]){
			//	if the current x vector is large
			if(size_obj <= size())	return;
			//	if the current x vector is small
			else{
				double* x_array_tmp  = new double [size_obj];
				double* data_tmp = new double [size_obj];
				memcpy(x_array_tmp, obj, size_obj);
				for(size_t i=0; i<size(); ++i)
					data_tmp[i] = get(obj[i]);

			}
		}
	}
*/

	inline void interp::add(interp& obj){		//	summation
		for(size_t i=0; i<size(); ++i) data[i] += obj.get(x_array[i]);
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::sub(interp& obj){		//	subtraction
		for(size_t i=0; i<size(); ++i) data[i] -= obj.get(x_array[i]);
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::mul(interp& obj){		//	multiplication
		for(size_t i=0; i<size(); ++i) data[i] *= obj.get(x_array[i]);
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::div(interp& obj){		//	devide
		for(size_t i=0; i<size(); ++i) data[i] /= obj.get(x_array[i]);
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::add(const double obj){		//	
		for(size_t i=0; i<size(); ++i) data[i] += obj;
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::sub(const double obj){		//	
		for(size_t i=0; i<size(); ++i) data[i] -= obj;
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::mul(const double obj){		//	
		for(size_t i=0; i<size(); ++i) data[i] *= obj;
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::div(const double obj){		//
		for(size_t i=0; i<size(); ++i) data[i] /= obj;
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}
	inline void interp::inverse(){		//
		for(size_t i=0; i<size(); ++i) data[i] = 1.0/data[i];
		gsl_interp_accel_reset(acc);
		gsl_interp_init(interp_obj, x_array, data, num);
	}

	//---	get the coefficient in [x(i) ~ x(i+1)]	---
	inline void interp::getCoef(const size_t index, double c[4]){
		switch(method){
		case Linear:
			return getCoefLinear(index, c);
			break;
		case Cspline:
			return getCoefCspline(index, c);
			break;
/*		case CsplineFlatEdge:
			return getCoefCspline(index, c);
*/			break;
		}return;
	}
	
	inline void interp::getCoefCspline(const size_t index, double c[4]){
	 /* evaluate */
		double dx = x_array[index+1] - x_array[index];
		if (dx > 0.0){
			const double dy = data[index+1] - data[index];
			
			cspline_state_t *state = (cspline_state_t *)(interp_obj->state);

			const double c_i =   state->c[index];
			const double c_ip1 = state->c[index + 1];
			c[0] = data[index];
			c[1] = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
			c[2] = c_i;
			c[3] = (c_ip1 - c_i) / (3.0 * dx);
			return;
		}else{
			c[0]=0;	c[1]=0;	c[2]=0;	c[3]=0;
			return;
		}
	}
	inline void interp::getCoefLinear(const size_t index, double c[4]){
		c[0] = data[index];
		c[1] = (data[index+1]-data[index])-((x_array[index+1]-x_array[index]));
		c[2]=0.0;	c[3]=0.0;
	}
	//---	initialize for flat edge splines	---
	inline int interp::CsplineClampedInit(gsl_interp * interp, const double xa[], const double ya[],
		size_t size, const double dydx0, const double dydxn)
	{
		//---	gsl_interp_init	---
		if (size != interp->size)	GSL_ERROR ("data must match size of interpolation object", GSL_EINVAL);
				
		for (size_t i = 1; i < size; i++) {
		    if (!(xa[i-1] < xa[i])) 
		        GSL_ERROR ("x values must be monotonically increasing", GSL_EINVAL);
		}

		interp->xmin = xa[0];
		interp->xmax = xa[size - 1];
		
		//---	modified from gsl_cspline_init	---
		cspline_state_t *state = (cspline_state_t *) (interp->state);
		
		size_t i;
		size_t num_points = size;
		size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
		size_t sys_size = max_index - 1;    /* linear system is (sys_size+2) x (sys_size+2) */
				
		double *offdiag2 = new double[max_index];
		for (i = 0; i < sys_size; i++)
		{
			const double h_i       = xa[i + 1] - xa[i];
			const double h_ip1     = xa[i + 2] - xa[i + 1];
			const double ydiff_i   = ya[i + 1] - ya[i];
			const double ydiff_ip1 = ya[i + 2] - ya[i + 1];
			const double g_i       = (h_i   != 0.0) ? 1.0 / h_i   : 0.0;
			const double g_ip1     = (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
			state->offdiag[i+1] = h_ip1;
			state->diag[i+1] = 2.0 * (h_ip1 + h_i);
			state->g[i+1] = 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i);
			offdiag2[i+1] = h_ip1;
		}
		//---	for index 0	---
		const double h_0  = xa[1] - xa[0];
		const double g_0  = (h_0 != 0.0) ? 1.0 / h_0 : 0.0;
		state->offdiag[0] = h_0;
		offdiag2[0]		  = h_0;
		state->diag[0]    = 2.0*h_0;
		state->g[0]       = 3.0 * ((ya[1]-ya[0]) * g_0 -  3.0 * dydx0);
		//---	for index n_1	---
		const double h_n  = xa[max_index  ] - xa[max_index-1];
		const double h_n_1= xa[max_index-1] - xa[max_index-2];
		const double g_n  = (h_n   != 0.0) ? 1.0 / h_n   : 0.0;
		const double g_n_1= (h_n_1 != 0.0) ? 1.0 / h_n_1 : 0.0;
		state->offdiag[max_index-2] = h_n;
		offdiag2[max_index-2]   = 3.0*h_n - 4.0*h_n_1;
		state->diag[max_index-1]= 9.0*h_n - 2.0*h_n_1;
		state->g[max_index-1]   = 9.0*g_n  *(ya[max_index  ] - ya[max_index-1])
								 -6.0*g_n_1*(ya[max_index-1] - ya[max_index-2])
								 -3.0*dydxn;
			
		if (sys_size == 1)
		{
			state->c[1] = state->g[1] / state->diag[1];
			delete [] offdiag2;
			return GSL_SUCCESS;
		}
		else
		{
			gsl_vector_view g_vec = gsl_vector_view_array(         state->g,       max_index);
			gsl_vector_view diag_vec = gsl_vector_view_array(      state->diag,    max_index);
			gsl_vector_view offdiag_vec = gsl_vector_view_array(   state->offdiag, max_index-1);
			gsl_vector_view offdiag2_vec= gsl_vector_view_array(   offdiag2,       max_index-1);
			gsl_vector_view solution_vec = gsl_vector_view_array ((state->c),      max_index);
			
			int status = gsl_linalg_solve_tridiag(&diag_vec.vector, 
			                                      &offdiag_vec.vector, 
			                                      &offdiag2_vec.vector, 
			                                      &g_vec.vector, 
			                                      &solution_vec.vector);
			
			state->c[max_index] = -2.0*state->c[max_index-1] 
			 +  3.0*g_n*(g_n  *( ya[  max_index] - ya[max_index-1])
					   -(g_n_1*((ya[max_index-1] - ya[max_index-2])) 
				              + (-h_n_1/3.0+h_n)*state->c[max_index-1] + (-2.0*h_n_1/3.0+h_n)*state->c[max_index-2]));


			delete [] offdiag2;
			return status;
		}
	}

	//----	convolution	----
	inline interp interp::convolute(interp& obj){
		interp convoluted(size(), method);
		for(size_t i=0; i<size(); ++i)
			convoluted.set(i, x_array[i], IntegConv(obj, x_array[i]));
		convoluted.setInterp();
		return convoluted;
	}
	


	//--- integration used in a convolution function	---
	// integ[ f(x-tau)*g(tau) ]dtau
	inline double interp::IntegConv(interp& obj, const double x){
		double val = 0.0;
		double tau = obj.x_array[0];
		int i=size()-2;
		size_t j=0;

		if(x-x_array[size()-1] > obj.x_array[0]){
			tau = x-x_array[size()-1];
			j = gsl_interp_accel_find(obj.acc, obj.x_array, obj.size(), tau);
		}else{
			i = (int)gsl_interp_accel_find(    acc,     x_array,     size(), x-tau);
			val -= IntegConv1Interval( obj, j, x, i, x-obj.x_array[j]);
		}
		while(i>0){
			val += IntegConv1Interval( obj, j, x, i, x_array[i+1]);
			if(j==obj.size()-1) return val; 
			--i;
		}
		return val;
	}
	
	//	integ[ (a0+a1*x+a2*x^2+a3*x^3)*(b0+b1*x+b2*x^2+b3*x^3) ](0,dx)
	inline double interp::IntegCoef(const double a[4], const double b[4], const double dx){
		double val=0.0;
		double pow_dx[7]; pow_dx[0]=dx;	// pow_dx[i] = dx^(i+1)
		for(size_t i=1; i<=6; ++i) pow_dx[i] = pow_dx[i-1]*dx;
		for(size_t i=0; i<=3; ++i){
			for(size_t j=0; j<=3; ++j){
				val += a[i]*b[j]*pow_dx[i+j]/(i+j+1); 
			}
		}
		return val;
	}

	inline double interp::IntegConv1Interval
		(interp& obj, size_t& j, const double x, const size_t i, const double xmax)
	{
		double val = 0.0;
		double coef_g[4];	//	coefficients of g(x-tau) at x[i]
		double coef_f[4], coef_f_new[4];
		double x_diff = x - x_array[i];
		getCoef(i, coef_f);
		
		//---	integration in [ x(i), y(j+1) ]---
		changeCoef(coef_f, coef_f_new, x-x_array[i]-obj.x_array[j]);

		obj.getCoef(j, coef_g);
		val =  -IntegCoef(coef_f_new, coef_g, x-x_array[i] -obj.x_array[j]);

		//---	integration in [ y(j), y(j+1) ]---
		while(obj.x_array[j+1]<x-xmax){
			changeCoef(coef_f, coef_f_new, x_diff-obj.x_array[j]);
			val += IntegCoef(coef_f_new, coef_g, obj.x_array[j+1]-obj.x_array[j]);
			++j;

			if(j+1==obj.size())	return val;	//	end of the range of obj.x
		}
		//---	integration in [ y(j), x(i+1) ]---
		val +=  IntegCoef(coef_f, coef_g, x-xmax-obj.x_array[j]);
		return val;
	}

	inline void interp::changeCoef(const double c[4], double coef[4], const double s){
		coef[3] =        -c[3];
		coef[2] = 3.0*s*  c[3] +       c[2];
		coef[1] =-3.0*s*s*c[3] - 2.0*s*c[2] -   c[1];
		coef[0] =   s*s*s*c[3] +   s*s*c[2] + s*c[1] + c[0];
	}


//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//												   //
//-------------------------------------------------//
	//---	constructors	---
	inline interp2d::interp2d(void){allocated = false;}
	inline interp2d::interp2d(const interp2d&obj){
		allocated = false;
		FreeMemory();
		if(obj.isallocated()){
			Allocate(obj.size1(), obj.size2());
			
			for(size_t j=0; j<obj.size2();++j) yv[j] = obj.yv[j];
			for(size_t i=0; i<obj.size1(); ++i){
				xv[i] = obj.xv[i];
				for(size_t j=0; j<obj.size2();++j){
					data[i][j] = obj.data[i][j];
				}
			}
			gsl_interp_accel_reset(acc_x);
			gsl_interp_accel_reset(acc_y);
		}
	}

	inline interp2d& interp2d::operator=(const interp2d& obj){
		FreeMemory();
		if(obj.isallocated()){
			Allocate(obj.size1(), obj.size2());
			
			for(size_t j=0; j<obj.size2();++j) yv[j] = obj.yv[j];
			for(size_t i=0; i<obj.size1(); ++i){
				xv[i] = obj.xv[i];
				for(size_t j=0; j<obj.size2();++j){
					data[i][j] = obj.data[i][j];
				}
			}
			gsl_interp_accel_reset(acc_x);
			gsl_interp_accel_reset(acc_y);
		}
		return *this;
	}

	inline interp2d::~interp2d(void){FreeMemory();	}

	inline void interp2d::clear(){
		FreeMemory();
	}
	inline void interp2d::FreeMemory(){
		if(allocated){
			delete [] xv;	delete [] yv;	
			gsl_interp_accel_free(acc_x);
			gsl_interp_accel_free(acc_y);
			for(size_t i=0; i<size1(); ++i)
				delete[] data[i];
			delete[] data;
		}
		allocated = false;
	}
	inline void interp2d::Allocate(size_t size1_, size_t size2_){
		num1 = size1_;	num2 = size2_;
		if(!allocated){
			xv = new double [size1_];
			yv = new double [size2_];
			data = (double **)(new double* [size1_]);
			for(size_t i=0; i<size1_; ++i)	data[i] = new double [size2_];

			acc_x = gsl_interp_accel_alloc();
			acc_y = gsl_interp_accel_alloc();
		}
		allocated = true;
	}


	inline void interp2d::set(const std::vector<double>& x_, const std::vector<double>& y_,
							  const std::vector<std::vector<double>>& data_, InterpMethod method_){
		method = method_;
		if(x_.size() < 4 || y_.size() < 4) method = Bilinear;	//	bicubic needs more than 4 points
		FreeMemory();
		Allocate(x_.size(), y_.size());

	//---	sorting		---
		size_t *px = new size_t[num1];
		size_t *py = new size_t[num2];
		for(size_t i=0; i<num1; ++i)	xv[i] = x_[i];
		gsl_sort_index(px, xv, 1, num1);
		for(size_t i=0; i<num2; ++i)	yv[i] = y_[i];
		gsl_sort_index(py, yv, 1, num2);

	//---	storing the data as a pointer	---
		for(size_t j=0; j<num2; ++j)	yv[j] = y_[py[j]];
		for(size_t i=0; i<num1; ++i){
			xv[i] = x_[px[i]];
			for(size_t j=0; j<num2; ++j)
				data[i][j] = data_[px[i]][py[j]];
		}

		gsl_interp_accel_reset(acc_x);
		gsl_interp_accel_reset(acc_y);
		delete [] px;
		delete [] py;

	}
		
	inline double interp2d::get(const double x, const double y){
		switch(method){
		case Bilinear:
			return getBilinear(x, y);
			break;
		case Bicubic:
			return getBicubic(x,y);
			break;
		};
		return 0.0;
	}
/*	inline void interp2d::getGradient(const double x, const double y, double& df_dx, double& df_dy){
		switch(method){
		case Bilinear:
			return getBilinear_gradient(x, y);
			break;
		case Bicubic:
			return getBicubic_gradient(x,y);
			break;
		return;
		};
	}
	*/
	inline double interp2d::getderiv_x(const double x, const double y){
		switch(method){
		case Bilinear:
			return getBilinear_dx(x, y);
			break;
		case Bicubic:
			return getBicubic_dx(x,y);
			break;
		return 0.0;
		};
		return 0.0;
	};

	inline double interp2d::getderiv_y(const double x, const double y){
		switch(method){
		case Bilinear:
			return getBilinear_dy(x, y);
			break;
		case Bicubic:
			return getBicubic_dy(x,y);
			break;
			return 0.0;
		};
		return 0.0;
	};

	inline double interp2d::getderiv_xy(const double x, const double y){
		switch(method){
		case Bilinear:
			return getBilinear_dxy(x, y);
			break;
		case Bicubic:
			return getBicubic_dxy(x,y);
			break;
			return 0.0;//	under construction	
		};
	};

	//---	integration along line:	line pararell to x axis	---
	inline double interp2d::getIntegAlong_x(  const double y, const double xfrom, const double xto){
		switch(method){
		case Bilinear: 
			if(     xfrom<xto) return  getBilinearIntegAlong_x(y, xfrom, xto);
			else if(xfrom>xto) return -getBilinearIntegAlong_x(y, xto, xfrom);
			else return 0.0;
			break;
		case Bicubic:
			return 0.0;	//	under construction	
		}
		return 0.0;
	}
	//---	integration along line:	line pararell to y axis	---
	inline double interp2d::getIntegAlong_y(  const double x, const double yfrom, const double yto){
		switch(method){
		case Bilinear: 
			if(     yfrom<yto)	return  getBilinearIntegAlong_y(x, yfrom, yto);
			else if(yfrom>yto)	return -getBilinearIntegAlong_y(x, yto, yfrom);
			else return 0.0;
			break;
		case Bicubic:
			return 0.0;	//	under construction	
		}
		return 0.0;
	}
	//---	integration along line:	arbitrary line	---
	inline double interp2d::getIntegAlongLine(const double xfrom, const double xto, const double yfrom, const double yto){
	switch(method){
		case Bilinear: 
			if(     xfrom == xto) return fabs(getIntegAlong_y(xfrom, yfrom, yto));
			else if(yfrom == yto) return fabs(getIntegAlong_x(yfrom, xfrom, xto));
			else if(xfrom<xto && yfrom<yto)	return  getBilinearIntegAlongLine(xfrom, xto, yfrom, yto);
			else if(xfrom>xto && yfrom<yto)	return  getBilinearIntegAlongLine(xto, xfrom, yfrom, yto);
			else if(xfrom<xto && yfrom>yto)	return  getBilinearIntegAlongLine(xfrom, xto, yto, yfrom);
			else if(xfrom>xto && yfrom>yto)	return  getBilinearIntegAlongLine(xto, xfrom, yto, yfrom);
			else return 0.0;
			break;
		case Bicubic:
			return 0.0;	//	under construction	
		}
		return 0.0;
	}

	inline interp interp2d::getAlong_x(const double y, interp::InterpMethod method){
		interp interp1d(num1, method);
		for(size_t i=0; i<num1; ++i)
			interp1d.set(i, xv[i], get(xv[i], y));
		interp1d.setInterp();
		return interp1d;
	}

	inline interp interp2d::getAlong_y(const double x, interp::InterpMethod method){
		interp interp1d(num2, method);
		for(size_t i=0; i<num2; ++i)
			interp1d.set(i, yv[i], get(x, yv[i]));
		interp1d.setInterp();
		return interp1d;
	}

	inline interp interp2d::getAlongLine
		(const double x0, const double y0, const double x1, const double y1,
		 const size_t num, interp::InterpMethod method)
	{
		interp interp1d(num, method);
		double rmax = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
		for(size_t i=0; i<num; ++i){
			double ri = rmax/(num-1)*i;
			double xi = x0 + (x1-x0)/(num-1)*i;
			double yi = y0 + (y1-y0)/(num-1)*i;
			interp1d.set(i, ri, get(xi, yi));
		}
		interp1d.setInterp();
		return interp1d;
	}

//---	algorithms for bilinear interpolation	---

	inline double interp2d::getBilinear(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor == num1-1) ix_floor --;
		if(iy_floor == num2-1) iy_floor --;

		double divcoef = (xv[ix_floor+1]-xv[ix_floor])*(yv[iy_floor+1]-yv[iy_floor]);
		return (data[ix_floor  ][iy_floor  ] * (xv[ix_floor+1]-x)*(yv[iy_floor+1]-y)
			  - data[ix_floor+1][iy_floor  ] * (xv[ix_floor  ]-x)*(yv[iy_floor+1]-y)
			  - data[ix_floor  ][iy_floor+1] * (xv[ix_floor+1]-x)*(yv[iy_floor  ]-y)
			  + data[ix_floor+1][iy_floor+1] * (xv[ix_floor  ]-x)*(yv[iy_floor  ]-y)) / divcoef;
	}

	inline double interp2d::getBilinear_dx(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor == num1-1) ix_floor --;
		if(iy_floor == num2-1) iy_floor --;

		double divcoef = (xv[ix_floor+1]-xv[ix_floor])*(yv[iy_floor+1]-yv[iy_floor]);
		return (data[ix_floor  ][iy_floor  ] * (-1.0)*(yv[iy_floor+1]-y)
			  - data[ix_floor+1][iy_floor  ] * (-1.0)*(yv[iy_floor+1]-y)
			  - data[ix_floor  ][iy_floor+1] * (-1.0)*(yv[iy_floor  ]-y)
			  + data[ix_floor+1][iy_floor+1] * (-1.0)*(yv[iy_floor  ]-y)) / divcoef;
	}

	inline double interp2d::getBilinear_dy(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor == num1-1) --ix_floor ;
		if(iy_floor == num2-1) --iy_floor ;

		double divcoef = (xv[ix_floor+1]-xv[ix_floor])*(yv[iy_floor+1]-yv[iy_floor]);
		return (data[ix_floor  ][iy_floor  ] * (xv[ix_floor+1]-x)*(-1.0)
			  - data[ix_floor+1][iy_floor  ] * (xv[ix_floor  ]-x)*(-1.0)
			  - data[ix_floor  ][iy_floor+1] * (xv[ix_floor+1]-x)*(-1.0)
			  + data[ix_floor+1][iy_floor+1] * (xv[ix_floor  ]-x)*(-1.0)) / divcoef;
	}

	inline double interp2d::getBilinear_dxy(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor == num1-1) --ix_floor ;
		if(iy_floor == num2-1) --iy_floor ;

		double divcoef = (xv[ix_floor+1]-xv[ix_floor])*(yv[iy_floor+1]-yv[iy_floor]);
		return (data[ix_floor  ][iy_floor  ] * (-1.0)*(-1.0)
			  - data[ix_floor+1][iy_floor  ] * (-1.0)*(-1.0)
			  - data[ix_floor  ][iy_floor+1] * (-1.0)*(-1.0)
			  + data[ix_floor+1][iy_floor+1] * (-1.0)*(-1.0)) / divcoef;
	}
/*	inline void interp2d::getBilinear_gradient(const double x, const double y, double& df_dx, double& df_dy){

	}
*/
	inline double interp2d::getBilinearIntegAlong_x(  const double y, const double xfrom, const double xto){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, xfrom);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		double integ1=0.0, integ2=0.0;

		double x1 = xfrom, x2;
		double A, B;
		while(xv[ix_floor] < xto){
			if(ix_floor+1>= num1)         x2 = xto;		//	if xto is out of range, the data is linearly edtrapolated
			else if(xv[ix_floor+1] > xto) x2 = xto;
			else x2 = xv[ix_floor+1];

			A = data[ix_floor  ][iy_floor  ] + data[ix_floor+1][iy_floor  ];
			B = data[ix_floor  ][iy_floor+1] + data[ix_floor+1][iy_floor+1];

			integ1 += A*( - (xv[ix_floor+1]-x2)*(xv[ix_floor+1]-x2) + (xv[ix_floor+1]-x1)*(xv[ix_floor+1]-x1))/(xv[ix_floor+1]-xv[ix_floor]);
			integ2 += B*(   (x2-xv[ix_floor  ])*(x2-xv[ix_floor  ]) - (x1-xv[ix_floor  ])*(x1-xv[ix_floor  ]))/(xv[ix_floor+1]-xv[ix_floor]);

			++ix_floor;
			if(ix_floor>=num1) break;					//	if xto is out of range, the integration is already finished
			x1 = xv[ix_floor];
		}
		return 0.5*(integ1*(yv[iy_floor+1]-y)+integ2*(y-yv[iy_floor]))/(yv[iy_floor+1]-yv[iy_floor]);
	}

	inline double interp2d::getBilinearIntegAlong_y(  const double x, const double yfrom, const double yto){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, yfrom);
		double integ1=0.0, integ2=0.0;

		double y1 = yfrom, y2;
		double A, B;
		while(yv[iy_floor] < yto){
			if(iy_floor+1>= num2)         y2 = yto;		//	if yto is out of range, the data is linearly edtrapolated
			else if(yv[iy_floor+1] > yto) y2 = yto;
			else y2 = yv[iy_floor+1];

			A = data[ix_floor  ][iy_floor  ] + data[ix_floor  ][iy_floor+1];
			B = data[ix_floor+1][iy_floor  ] + data[ix_floor+1][iy_floor+1];
			
			integ1 += A*( - (yv[iy_floor+1]-y2)*(yv[iy_floor+1]-y2) + (yv[iy_floor+1]-y1)*(yv[iy_floor+1]-y1))/(yv[iy_floor+1]-yv[iy_floor]);
			integ2 += B*(   (y2-yv[iy_floor  ])*(y2-yv[iy_floor  ]) - (y1-yv[iy_floor  ])*(y1-yv[iy_floor  ]))/(yv[iy_floor+1]-yv[iy_floor]);

			++iy_floor;
			if(iy_floor>=num2) break;					//	if yto is out of range, the integration is already finished...
			y1 = yv[iy_floor];
		}
		return 0.5*(integ1*(xv[ix_floor+1]-x)+integ2*(x-xv[ix_floor]))/(xv[ix_floor+1]-xv[ix_floor]);
	}

	inline double interp2d::getBilinearIntegAlongLine(const double xfrom, const double xto, const double yfrom, const double yto){
		int ix_width = (int)gsl_interp_accel_find(acc_x, xv, num1, xto) - (int)gsl_interp_accel_find(acc_x, xv, num1, xfrom);
		int iy_width = (int)gsl_interp_accel_find(acc_y, yv, num2, yto) - (int)gsl_interp_accel_find(acc_y, yv, num2, yfrom);
		if(ix_width >= iy_width)			//	if the slope is small, 
			 return getBilinearIntegAlongLine_x_integ(xfrom, xto, yfrom, yto);		//	integrate by means of x
		else return getBilinearIntegAlongLine_y_integ(xfrom, xto, yfrom, yto);		//	integrate by means of y
	}

	inline double interp2d::getBilinearIntegAlongLine_x_integ(const double xfrom, const double xto, const double yfrom, const double yto){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, xfrom);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, yfrom);

		double alpha, beta, gamma;
		double delta_x = xto-xfrom,	delta_y = yto-yfrom;
		double slope = delta_y/delta_x;
		double y1, y2 = yv[iy_floor];
		double x1 = xfrom, x2 = xfrom;
		double x1c, x2c = 1.0/slope*(yv[iy_floor]-yfrom)+xfrom;		//	cross point of the line and y=y2 
		double dxc;
		double integ = 0.0;
		//---	for 1 division of y	---
		while(x2 < xto){
			
			y1 = y2;
			y2 = yv[iy_floor+1];					//	if still in the range
			x1c = x2c;								//	cross point of the line and y=y1 
			x2c = 1.0/slope*(y2-yfrom)+xfrom;		//	cross point of the line and y=y2 
			dxc = x2c - x1c;

			//---	for 1 division of x	---
			while(x2 < xto){
				if(ix_floor+2>= num1){//	if xto is out of range, the data is linearly edtrapolated
					x2 = xto;		
					--ix_floor;
				}
				else x2 = xv[ix_floor+1];
				//				A														B													C														D
				alpha = data[ix_floor  ][iy_floor  ]                      - data[ix_floor+1][iy_floor  ]                      - data[ix_floor  ][iy_floor+1]                      + data[ix_floor+1][iy_floor+1];
				beta  =-data[ix_floor  ][iy_floor  ]*(xv[ix_floor+1]+x2c) + data[ix_floor+1][iy_floor  ]*(xv[ix_floor  ]+x2c) + data[ix_floor  ][iy_floor+1]*(xv[ix_floor+1]+x1c) - data[ix_floor+1][iy_floor+1]*( xv[ix_floor  ]+x1c);
				gamma = data[ix_floor  ][iy_floor  ]*(xv[ix_floor+1]*x2c) - data[ix_floor+1][iy_floor  ]*(xv[ix_floor  ]*x2c) - data[ix_floor  ][iy_floor+1]*(xv[ix_floor+1]*x1c) + data[ix_floor+1][iy_floor+1]*( xv[ix_floor  ]*x1c);

				//	if still in the region
				if(x2 >= xto){
					x2 = xto;
					integ += (alpha*(x2*x2*x2 - x1*x1*x1)/3.0 + beta*(x2*x2 - x1*x1)/2.0 + gamma*(x2-x1))
								/(dxc*(xv[ix_floor+1]-xv[ix_floor]));
					break;
				}
				else if(x2 >= x2c && iy_floor+1 < num2){
					x2 = x2c;
					integ += (alpha*(x2*x2*x2 - x1*x1*x1)/3.0 + beta*(x2*x2 - x1*x1)/2.0 + gamma*(x2-x1))
								/(dxc*(xv[ix_floor+1]-xv[ix_floor]));
					x1 = x2;
					break;
				}
				
				//	if end of the region
				else{
					integ += (alpha*(x2*x2*x2 - x1*x1*x1)/3.0 + beta*(x2*x2 - x1*x1)/2.0 + gamma*(x2-x1))
								/(dxc*(xv[ix_floor+1]-xv[ix_floor]));
					++ix_floor;
					x1 = x2;
				}
			}
			++iy_floor;
		}
		return integ*sqrt(delta_y*delta_y + delta_x*delta_x)/delta_x;
	}
	inline double interp2d::getBilinearIntegAlongLine_y_integ(const double xfrom, const double xto, const double yfrom, const double yto){
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, yfrom);
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, xfrom);

		double alpha, beta, gamma;
		double delta_y = yto-yfrom,	delta_x = xto-xfrom;
		double slope = delta_x/delta_y;
		double x1, x2 = xv[ix_floor];
		double y1 = yfrom, y2 = yfrom;
		double y1c, y2c = 1.0/slope*(xv[ix_floor]-xfrom)+yfrom;		//	cross point of the line and x=x2 
		double dyc;
		double integ = 0.0;
		//---	for 1 division of x	---
		while(y2 < yto){
			
			x1 = x2;
			x2 = xv[ix_floor+1];					//	if still in the range
			y1c = y2c;								//	cross point of the line and x=x1 
			y2c = 1.0/slope*(x2-xfrom)+yfrom;		//	cross point of the line and x=x2 
			dyc = y2c - y1c;

			//---	for 1 division of y	---
			while(y2 < yto){
				if(iy_floor+2>= num2){//	if yto is out of range, the data is linearlx edtrapolated
					y2 = yto;		
					--iy_floor;
				}
				else y2 = yv[iy_floor+1];
				//				A														B													C														D
				alpha = data[ix_floor  ][iy_floor  ]                      - data[ix_floor  ][iy_floor+1]                      - data[ix_floor+1][iy_floor  ]                      + data[ix_floor+1][iy_floor+1];
				beta  =-data[ix_floor  ][iy_floor  ]*(yv[iy_floor+1]+y2c) + data[ix_floor  ][iy_floor+1]*(yv[iy_floor  ]+y2c) + data[ix_floor+1][iy_floor  ]*(yv[iy_floor+1]+y1c) - data[ix_floor+1][iy_floor+1]*( yv[iy_floor  ]+y1c);
				gamma = data[ix_floor  ][iy_floor  ]*(yv[iy_floor+1]*y2c) - data[ix_floor  ][iy_floor+1]*(yv[iy_floor  ]*y2c) - data[ix_floor+1][iy_floor  ]*(yv[iy_floor+1]*y1c) + data[ix_floor+1][iy_floor+1]*( yv[iy_floor  ]*y1c);

				//	if still in the region
				if(y2 >= yto){
					y2 = yto;
					integ += (alpha*(y2*y2*y2 - y1*y1*y1)/3.0 + beta*(y2*y2 - y1*y1)/2.0 + gamma*(y2-y1))
								/(dyc*(yv[iy_floor+1]-yv[iy_floor]));
					break;
				}
				else if(y2 >= y2c && ix_floor+1 < num1){
					y2 = y2c;
					integ += (alpha*(y2*y2*y2 - y1*y1*y1)/3.0 + beta*(y2*y2 - y1*y1)/2.0 + gamma*(y2-y1))
								/(dyc*(yv[iy_floor+1]-yv[iy_floor]));
					y1 = y2;
					break;
				}
				
				//	if end of the region
				else{
					integ += (alpha*(y2*y2*y2 - y1*y1*y1)/3.0 + beta*(y2*y2 - y1*y1)/2.0 + gamma*(y2-y1))
								/(dyc*(yv[iy_floor+1]-yv[iy_floor]));
					++iy_floor;
					y1 = y2;
				}
			}
			++ix_floor;
		}
		return integ*sqrt(delta_y*delta_y + delta_x*delta_x)/delta_y;
	}

//---	algorithms for bicubic interpolation	---

	inline void interp2d::getBicubicCoef(size_t ix_floor, size_t iy_floor){
		f[ 0] = data[ix_floor  ][iy_floor  ];		//f00   
		f[ 1] = data[ix_floor+1][iy_floor  ];		//f10   
		f[ 2] = data[ix_floor  ][iy_floor+1];		//f01   
		f[ 3] = data[ix_floor+1][iy_floor+1];		//f11   

		//---	derivatives are calculated with a central difference	--- 
		f[ 4] = (data[ix_floor+1][iy_floor  ]-data[ix_floor-1][iy_floor  ])/(2.0*(xv[ix_floor+1]-xv[ix_floor-1]));		//	dfdx00 
		f[ 5] = (data[ix_floor+2][iy_floor  ]-data[ix_floor  ][iy_floor  ])/(2.0*(xv[ix_floor+2]-xv[ix_floor  ]));		//	dfdx10 
		f[ 6] = (data[ix_floor+1][iy_floor+1]-data[ix_floor-1][iy_floor+1])/(2.0*(xv[ix_floor+1]-xv[ix_floor-1]));		//	dfdx01 
		f[ 7] = (data[ix_floor+2][iy_floor+1]-data[ix_floor  ][iy_floor+1])/(2.0*(xv[ix_floor+2]-xv[ix_floor  ]));		//	dfdx11 
											 
		f[ 8] = (data[ix_floor  ][iy_floor+1]-data[ix_floor  ][iy_floor-1])/(2.0*(yv[iy_floor+1]-yv[iy_floor-1]));	//	dfdy00
		f[ 9] = (data[ix_floor+1][iy_floor+1]-data[ix_floor+1][iy_floor-1])/(2.0*(yv[iy_floor+1]-yv[iy_floor-1]));	//	dfdy10
		f[10] = (data[ix_floor  ][iy_floor+2]-data[ix_floor  ][iy_floor  ])/(2.0*(yv[iy_floor+2]-yv[iy_floor  ]));	//	dfdy01
		f[11] = (data[ix_floor+1][iy_floor+2]-data[ix_floor+1][iy_floor  ])/(2.0*(yv[iy_floor+2]-yv[iy_floor  ]));	//	dfdy11
											 
		f[12] = (data[ix_floor+1][iy_floor+1]-data[ix_floor-1][iy_floor-1])/(2.0*(xv[ix_floor+1]-xv[ix_floor-1])*(yv[iy_floor+1]-yv[iy_floor-1]));	//	ddfdxdy00
		f[13] = (data[ix_floor+2][iy_floor+1]-data[ix_floor  ][iy_floor-1])/(2.0*(xv[ix_floor+2]-xv[ix_floor  ])*(yv[iy_floor+1]-yv[iy_floor-1]));	//	ddfdxdy10
		f[14] = (data[ix_floor+1][iy_floor+2]-data[ix_floor-1][iy_floor  ])/(2.0*(xv[ix_floor+1]-xv[ix_floor-1])*(yv[iy_floor+2]-yv[iy_floor  ]));	//	ddfdxdy01
		f[15] = (data[ix_floor+2][iy_floor+2]-data[ix_floor  ][iy_floor  ])/(2.0*(xv[ix_floor+2]-xv[ix_floor  ])*(yv[iy_floor+2]-yv[iy_floor  ]));	//	ddfdxdy11

	//---	matrix for bicubic interpolation	---
		a[0][0] = f[ 0];
		a[1][0] = f[ 4];
		a[2][0] = -3*f[ 0]+3*f[ 1]-2*f[ 4]-f[ 5];
		a[3][0] =  2*f[ 0]-2*f[ 1]+  f[ 4]+f[ 5];
		a[0][1] = f[ 8];
		a[1][1] = f[12];
		a[2][1] = -3*f[ 8]+3*f[ 9]-2*f[12]-f[13];
		a[3][1] =  2*f[ 8]-2*f[ 9]+  f[12]+f[13];
		a[0][2] = -3*f[ 0]+3*f[ 2]-2*f[ 8]-f[10];
		a[1][2] = -3*f[ 4]+3*f[ 6]-2*f[12]-f[14];
		a[2][2] =  9*f[ 0]-9*f[ 1]-9*f[ 2]+9*f[ 3]+6*f[ 4]+3*f[ 5]-6*f[ 6]-3*f[ 7]+6*f[ 8]-6*f[ 9]+3*f[10]-3*f[11]+4*f[12]+2*f[13]+2*f[14]+1*f[15];
		a[3][2] = -6*f[ 0]+6*f[ 1]+6*f[ 2]-6*f[ 3]-3*f[ 4]-3*f[ 5]+3*f[ 6]-3*f[ 7]-4*f[ 8]+4*f[ 9]-2*f[10]+2*f[11]-2*f[12]-2*f[13]-1*f[14]-1*f[15];
		a[0][3] =  2*f[ 0]-2*f[ 2]+1*f[ 8]+1*f[10];
		a[1][3] =  2*f[ 4]-2*f[ 6]+1*f[12]+1*f[14];
		a[2][3] = -6*f[ 0]+6*f[ 1]+6*f[ 2]-6*f[ 3]-4*f[ 4]-2*f[ 5]+4*f[ 6]+2*f[ 7]-3*f[ 8]+3*f[ 9]-3*f[10]+3*f[11]-2*f[12]-1*f[13]-2*f[14]-1*f[15];
		a[3][3] =  4*f[ 0]-4*f[ 1]-4*f[ 2]+4*f[ 3]-2*f[ 4]+2*f[ 5]-2*f[ 6]-2*f[ 7]+2*f[ 8]-2*f[ 9]+2*f[10]-2*f[11]+1*f[12]+1*f[13]+1*f[14]+1*f[15];
	}

	inline double interp2d::getBicubic(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor <= 1    ) ix_floor = 1;
		if(ix_floor >=num1-2) ix_floor = num1-3;
		if(iy_floor <= 1    ) iy_floor = 1;
		if(iy_floor >=num2-2) iy_floor = num2-3;
/*		if((ix_floor < 1 || ix_floor >num1-3) || (iy_floor < 1 || iy_floor >num2-3))
			return getBilinear(x, y);
*/		getBicubicCoef(ix_floor, iy_floor);

		double xtmp = (x - xv[ix_floor])/(xv[ix_floor+1]-xv[ix_floor]);
		double ytmp = (y - yv[iy_floor])/(yv[iy_floor+1]-yv[iy_floor]);

		double val=0.0;
		for(size_t i=0; i<4; ++i)	
			for(size_t j=0; j<4; ++j)
				val += a[i][j] * pow(xtmp, 1.0*i) * pow(ytmp, 1.0*j);
		return val;
	}

	inline double interp2d::getBicubic_dx(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor <= 1    ) ix_floor = 2;
		if(ix_floor >=num1-2) ix_floor = num1-3;
		if(iy_floor <= 1    ) iy_floor = 2;
		if(iy_floor >=num2-2) iy_floor = num2-3;
		getBicubicCoef(ix_floor, iy_floor);

		double xtmp = (x - xv[ix_floor])/(xv[ix_floor+1]-xv[ix_floor]);
		double ytmp = (y - yv[iy_floor])/(yv[iy_floor+1]-yv[iy_floor]);

		double val=0.0;
		for(size_t i=1; i<4; ++i)	
			for(size_t j=0; j<4; ++j)
				val += a[i][j] * i * pow(xtmp, 1.0*(i-1))/(xv[ix_floor+1]-xv[ix_floor]) * pow(ytmp, 1.0*j);
		return val;
	}

	inline double interp2d::getBicubic_dy(const double x, const double y){
		size_t ix_floor = gsl_interp_accel_find(acc_x, xv, num1, x);
		size_t iy_floor = gsl_interp_accel_find(acc_y, yv, num2, y);
		
		if(ix_floor <= 1    ) ix_floor = 2;
		if(ix_floor >=num1-2) ix_floor = num1-3;
		if(iy_floor <= 1    ) iy_floor = 2;
		if(iy_floor >=num2-2) iy_floor = num2-3;
		getBicubicCoef(ix_floor, iy_floor);

		double xtmp = (x - xv[ix_floor])/(xv[ix_floor+1]-xv[ix_floor]);
		double ytmp = (y - yv[iy_floor])/(yv[iy_floor+1]-yv[iy_floor]);

		double val=0.0;
		for(size_t i=0; i<4; ++i)	
			for(size_t j=1; j<4; ++j)
				val += a[i][j]  * pow(xtmp, 1.0*i) * j * pow(ytmp, 1.0*(j-1)/(yv[iy_floor+1]-yv[iy_floor]));
		return val;
	}
	
	

};