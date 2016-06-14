template < class TyE_ >
inline void EdgeSpline::set(const TyE_& y_src, size_t size_){
	set(y_src, size_, Free);
}
template < class TyE_ >
inline void EdgeSpline::set(const TyE_& y_src, size_t size_, EndPoint endpoint_){
	MemoryControl(size_);	
	for(size_t i=0; i<size; ++i){
		x[i] = 1.0*i;	y[i] = y_src[i];		dx[i]= 1.0;
		if(i != 0) x_edge[i-1] = 0.5*(x[i-1]+x[i]);
	}
	endpoint = endpoint_;
	setScale();
	sort();
	solve();
}
template < class TyE_ >
inline void EdgeSpline::set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, double dx_, size_t size_){
	set(x_src, y_src, edge_src, dx_, Free);
}
template < class TyE_ >
inline void EdgeSpline::set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, double dx_, size_t size_, EndPoint endpoint_){
	MemoryControl(size_);	
	for(size_t i=0; i<size; ++i){
		x[i] = x_src[i];	y[i] = y_src[i];		dx[i]= dx_;
		if(i != size-1) x_edge[i] = edge_src[i];
	}
	setScale();
	sort();
	solve();
}
template < class TyE_ >
void EdgeSpline::set(const TyE_& x_src, const TyE_& y_src, double dx_, size_t size_){
	set(x_src, y_src, dx_, size_, Free);
}
template < class TyE_ >
void EdgeSpline::set(const TyE_& x_src, const TyE_& y_src, double dx_, size_t size_, EndPoint endpoint_){
	endpoint = endpoint_;
	MemoryControl(size_);	
	for(size_t i=0; i<size; ++i){
		x[i] = x_src[i];	y[i] = y_src[i];		dx[i]= dx_;
		if(i != 0) x_edge[i-1] = 0.5*(x[i-1]+x[i]);
	}
	setScale();
	sort();
	solve();
}

void EdgeSpline::setScale(){ scale = Linear;}
	
double EdgeSpline::get(double x_interp){
	//---	search the region which x_interp belong to	---
	size_t index;
	if(x_interp <= x_edge[0]) index =0;
	else if(x_interp > x_edge[size-2]) index = size-1;
	else{
		for(size_t i=0; i<size-1; ++i){
			if(x_interp <= x_edge[i]){	index = i; break;}
		}
	}
	//---	return the value	---
	if(scale == Linear || scale == LinearArea){
		//	debug	
		double a_tmp = a[index];
		double b_tmp = b[index];
		double c_tmp = c[index];
		return a[index]+b[index]*x_interp+c[index]*x_interp*x_interp;
	}
	else if(scale == Log || scale == LogArea){
		double val = exp(a[index]+b[index]*x_interp+c[index]*x_interp*x_interp);
		if(!(val>0.0)) return 0.0;
		else return val;
	}
	else return 0.0;
}

inline void EdgeSpline::solve(){
	gsl_matrix_set_zero(X);
	gsl_vector_set_zero(coef);
	gsl_vector_set_zero(v);
	//---	condition (I)	---
	//---	interpolated value at the x[i] should be same to the y[i]	---
	if(scale == Linear){
		for(size_t i=0; i<size; ++i){
			gsl_matrix_set(X, i, 3*i,		  1.0);
			gsl_matrix_set(X, i, 3*i+1,		 x[i]);
			gsl_matrix_set(X, i, 3*i+2,	x[i]*x[i]);
			gsl_vector_set(v, i,		y[i]);
		}
	}else if(scale == Log || scale == LogArea){
		for(size_t i=0; i<size; ++i){
			gsl_matrix_set(X, i, 3*i,		  1.0);
			gsl_matrix_set(X, i, 3*i+1,		 x[i]);
			gsl_matrix_set(X, i, 3*i+2,	x[i]*x[i]);
			gsl_vector_set(v, i,		log(y[i]));
		}
	}else if(scale == LinearArea){
		for(size_t i=0; i<size; ++i){
			gsl_matrix_set(X, i, 3*i,		dx[i]);
			gsl_matrix_set(X, i, 3*i+1, x[i]*dx[i]);
			gsl_matrix_set(X, i, 3*i+2,	x[i]*x[i]*dx[i] + dx[i]*dx[i]*dx[i]/12.0);
			gsl_vector_set(v, i,			y[i]);
		}
	}
	//---	condition (II)	---
	//---	values at the x_edge[i] are same between two functions	---
	size_t offset;
	offset = size;
	for(size_t i=0; i<size-1; ++i){
		gsl_matrix_set(X, offset+i, 3*i,	1.0);					gsl_matrix_set(X, offset+i, 3*i+3, -1.0);
		gsl_matrix_set(X, offset+i, 3*i+1,	x_edge[i]);				gsl_matrix_set(X, offset+i, 3*i+4, -x_edge[i]);
		gsl_matrix_set(X, offset+i, 3*i+2,	x_edge[i]*x_edge[i]);	gsl_matrix_set(X, offset+i, 3*i+5, -x_edge[i]*x_edge[i]);
	}
	//---	condition (III)	---
	//---	1t derivatives at the x_edge[i] are same between two functions	---
	offset = 2*size-1;
	for(size_t i=0; i<size-1; ++i){
		gsl_matrix_set(X, offset+i, 3*i+1,	1.0);			gsl_matrix_set(X, offset+i, 3*i+4, -1.0);
		gsl_matrix_set(X, offset+i, 3*i+2, 2.0*x_edge[i]);	gsl_matrix_set(X, offset+i, 3*i+5, -2.0*x_edge[i]);
	}
	//---	condition (VI)	---
	//---	2nd derivatives of the end region is the same to the next one	---
	offset = 3*size-2;
	if(endpoint == Free){
		gsl_matrix_set(X, offset  , 2,	1.0);			gsl_matrix_set(X, offset  , 3+2, -1.0);
		gsl_matrix_set(X, offset+1, 3*size-4,	1.0);	gsl_matrix_set(X, offset+1, 3*size-1, -1.0);
	}else{
		gsl_matrix_set(X, offset  , 2,	1.0);
		gsl_matrix_set(X, offset+1, 3*size-1,	1.0);	
	}
	//---	solving the system of the equation	---
	int signum=0; 
	gsl_linalg_LU_decomp(X, p, &signum);
	gsl_linalg_LU_solve( X, p, v, coef);

	//---	keep the derived value in the pointer a, b, c	---
	for(size_t i=0; i<size; ++i){
		a[i] = get_a(coef, i);
		b[i] = get_b(coef, i);
		c[i] = get_c(coef, i);
	}
}
	
inline void EdgeSpline::sort(){
	// under construction	
}

inline void EdgeSpline::MemoryControl(size_t size_){
	if(size_ < 3){	exit(8);}
	//---	memory is already allocated, but size is different	---
	if(Allocated && size!= size_){		//	once free 
		FreeMemory();
	}
	size = size_;
	if(!Allocated) Allocation();
}
inline void EdgeSpline::Allocation(){
	//---	memory is not allocated	---
	x = new double[size];		y = new double[size];
	dx = new double[size];		x_edge = new double[size-1];
	a = new double[size];		b = new double[size];		c = new double[size];
	v = gsl_vector_calloc(size*3);
	coef = gsl_vector_calloc(size*3);
	X = gsl_matrix_calloc(size*3, size*3);
	p = gsl_permutation_calloc(size*3);
	Allocated = true;
}
inline void EdgeSpline::FreeMemory(){
	delete[] x;			delete[] y;
	delete[] x_edge;	delete[] dx;
	delete[] a;			delete[] b;			delete[] c;
	gsl_vector_free(v);	gsl_vector_free(coef);
	gsl_matrix_free(X);
	gsl_permutation_free(p);
	Allocated = false;
}

inline EdgeSpline::EdgeSpline(void){
	Allocated = false;
	endpoint  = Free;
}

inline EdgeSpline::~EdgeSpline(void){	if(Allocated) FreeMemory();}

//-------------------	EdgeLogSpline Class	--------------------
inline void EdgeLogSplineArea::solve(){
	//---	initial guess	---
	EdgeSpline::solve();
	for(size_t i=0; i<size; ++i) 
		for(size_t j=0; j<3; ++j)
			gsl_vector_set(coef, 3*i+j, gsl_vector_get(coef, 3*i+j)*dx[i]);
	//-------------------------

	size_t iter = 0;
	FDF.f  = &EdgeLogSpline_Multiroot::f;
	FDF.n  = size*3;
	FDF.params = (void*)this;
	gsl_multiroot_fsolver_set(s, &FDF, coef);
	//---	rootine	---
	int status;
	do{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		std::cerr << "iteration = "<<iter << std::endl;
		if(status) break;
		status = gsl_multiroot_test_residual(s->f, 1.0e-7);	//<--- residual should be determined numerically
	}while(status==GSL_CONTINUE && iter < 100);
	std::cerr<<"status = "<< gsl_strerror(status);
	//---	keep the result into a, b, c vectors	---
	for(size_t i=0; i<size; ++i){
		a[i] = get_a(s->x, i);
		b[i] = get_b(s->x, i);
		c[i] = get_c(s->x, i);
	}
}

void EdgeLogSplineArea::f_multiroot(const gsl_vector *coef_, gsl_vector *f){
	gsl_vector_set_zero(f);
	size_t offset;
	//---	condition (I)	---
	//---	interpolated function's integral from x[i]-dx[i]/2 to x[i]+dx[i]/2 should be same to the y[i]	---
	offset = 0;
	for(size_t i=0; i<size; ++i){
		double exp_abc = exp(get_a(coef_,i)+get_b(coef_,i)*x[i]+get_c(coef_,i)*x[i]*x[i]);
		double a_ = get_a(coef_,i);
		double b_ = get_b(coef_,i);
		double c_ = get_c(coef_,i);
		double ftmp_real = -y[i]
				+ exp_abc*dx[i]
				+1.0/6.0*(2.0*c_+(b_+2.0*c_*x[i])*(b_+2.0*c_*x[i]))*exp_abc*(dx[i]*dx[i]*dx[i]/4.0);
		gsl_vector_set(f, i, ftmp_real);
	}
	offset = size;
	for(size_t i=0; i<size-1; ++i){
		gsl_vector_set(f, offset+i, 
			  get_a(coef_, i  ) + get_b(coef_, i  )*x_edge[i] + get_c(coef_, i  )*x_edge[i]*x_edge[i]
			-(get_a(coef_, i+1) + get_b(coef_, i+1)*x_edge[i] + get_c(coef_, i+1)*x_edge[i]*x_edge[i]));
	}
	//---	condition (III)	---
	//---	1t derivatives at the x_edge[i] are same between two functions	---
	offset = 2*size-1;
	for(size_t i=0; i<size-1; ++i){
		gsl_vector_set(f, offset+i, 
			  get_b(coef_, i  ) + 2.0*get_c(coef_, i  )*x_edge[i]
			-(get_b(coef_, i+1) + 2.0*get_c(coef_, i+1)*x_edge[i]));
	}
	//---	condition (VI)	---
	//---	2nd derivatives of the end region is the same to the next one	---
	offset = 3*size-2;
	if(endpoint == Free){
		gsl_vector_set(f, offset  , 	get_c(coef_, 1)      - get_c(coef_, 0)     );
		gsl_vector_set(f, offset+1, 	get_c(coef_, size-2) - get_c(coef_, size-1));
	}else{
		gsl_vector_set(f, offset  , 	get_c(coef_, 0)      );
		gsl_vector_set(f, offset+1, 	get_c(coef_, size-1) );
	}
}


inline void EdgeLogSplineArea::Allocation(){
	f = gsl_vector_calloc(size*3);
	s = gsl_multiroot_fsolver_alloc( T, size*3);
	EdgeSpline::Allocation();
}
inline void EdgeLogSplineArea::FreeMemory(){
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(f);
	EdgeSpline::FreeMemory();
}

const gsl_multiroot_fsolver_type* EdgeLogSplineArea::T = gsl_multiroot_fsolver_hybrids;
