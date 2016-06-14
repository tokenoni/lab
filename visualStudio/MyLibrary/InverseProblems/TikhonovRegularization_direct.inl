inline TikhonovRegularization_direct::TikhonovRegularization_direct
(const gsl_matrix * A_, const gsl_matrix * L_,
const gsl_vector * d_,
const gsl_vector * f_inf_,
const double ratio_stop_, const size_t max_iter_)
:
A(A_), L(L_), d(d_), f_inf(f_inf_),
ratio_stop(ratio_stop_), max_iter(max_iter_),
smat(gsl_matrix_calloc(f_inf->size, 3)),
cov(gsl_matrix_calloc(f_inf->size, f_inf->size)),
r(gsl_vector_calloc(A->size1))
{
}


inline TikhonovRegularization_direct::~TikhonovRegularization_direct()
{
	gsl_matrix_free(smat);
	gsl_matrix_free(cov);
	gsl_vector_free(r);
}

//	with specified Lambda gradient
//inline bool TikhonovRegularization_direct::run
//	(gsl_vector* finit,
//	const double lambdaGradient,
//	const double lambdaMin, const double lambdaMax, const size_t lambdaNum)
//{
//	return true;
//}

inline bool TikhonovRegularization_direct::run_maximum_Curvature
(gsl_vector* finit, double& lambda_init, 
 const double lambda_min, const double lambda_max, 
 const size_t iterMax, const double lambda_epsrel)
{
	int status;
	int iter = 0; 
/*	const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);

	gsl_function F;
	F.function = &(TikhonovRegularization_direct::getCurvature);
	F.params = this;

	double a = lambda_min;
	double b = lambda_max;
	gsl_min_fminimizer_set(s, &F, lambda_init, a, b);
	do {
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		lambda_init = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);
		status = gsl_min_test_interval(a, b, 0.0, lambda_epsrel);
		if (status == GSL_SUCCESS) printf("Converged:\n");
		printf("%5d %10.7f %10.7f\n",
			iter, lambda_init, a - b);
	} while (status == GSL_CONTINUE && iter < max_iter);
	gsl_min_fminimizer_free(s);
*/	
	const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 1);
	gsl_vector *ss = gsl_vector_calloc(1), *x = gsl_vector_calloc(1);
	gsl_vector_set(x, 0, log(lambda_init));
	gsl_vector_set(ss, 0, 10.0);
	
	gsl_multimin_function minex_func;
	minex_func.n = 1;
	minex_func.f = &(TikhonovRegularization_direct::getCurvature);
	minex_func.params = this;
	
	double size;
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, lambda_epsrel);
		if (status == GSL_SUCCESS) printf("converged to minimum at\n");
		printf("%5d %.5f %10.5e\n", iter, exp(gsl_vector_get(s->x, 0)), s->f);
	} while (status == GSL_CONTINUE && iter < 100);

	lambda_init = exp(gsl_vector_get(s->x, 0));

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	run(finit, lambda_init);
	if (iter < iterMax) return true;
	else return false;
}



inline bool TikhonovRegularization_direct::run_specified_LGradient
	(gsl_vector* finit, const double L_gradient_, double& lambda_init, const size_t iterMax, const double lambda_epsrel)
{
	L_gradient = L_gradient_;
	int status;
	int iter = 0;
	const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
	gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc(T);
	gsl_function_fdf FDF;
	
	FDF.f = &(TikhonovRegularization_direct::getGradient_f);
	FDF.df =  &(TikhonovRegularization_direct::getGradient_df);
	FDF.fdf = &(TikhonovRegularization_direct::getGradient_fdf);
	FDF.params = this;

	gsl_root_fdfsolver_set(s, &FDF, lambda_init);
	double lambda_new = lambda_init; 
	do {
		iter++;
		status = gsl_root_fdfsolver_iterate(s);
		lambda_init = lambda_new;
		lambda_new = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta (lambda_new, lambda_init, 0, lambda_epsrel);
		printf("%5d %10.7e %10.7e\n",
			iter, lambda_init, lambda_new - lambda_init);
		if (status == GSL_SUCCESS) printf("Converged:\n");
	} while (status == GSL_CONTINUE && iter < iterMax);

	lambda_init = lambda_new;
	run(finit, lambda_init);
	if(iter < iterMax) return true;
	else return false;
}

//----------------------------------------------
//	static functions for lambda estimation
//----------------------------------------------
inline double TikhonovRegularization_direct::getGradient_f(double x, void * params){
	TikhonovRegularization_direct* p = (TikhonovRegularization_direct*) params;
	gsl_vector *f = gsl_vector_alloc(p->f_inf->size);
	p->run(f, x);
	gsl_vector_free(f);
	return - p->rho/(p->xi*x*x) - p->L_gradient;
}
inline double TikhonovRegularization_direct::getGradient_df(double x, void * params){
	TikhonovRegularization_direct* p = (TikhonovRegularization_direct*) params;
	gsl_vector *f = gsl_vector_alloc(p->f_inf->size);
	p->runWithCurvature(f, x);
	double rho = p->rho;
	double eta = p->xi;
	gsl_vector_free(f);
	return 2.0*rho / (x*x*eta) + (rho/(x*x*eta*eta) + 1.0/eta)/p->dxi_dlambda;
}
inline void TikhonovRegularization_direct::getGradient_fdf(double x, void * params, double * f_, double * df_){
	TikhonovRegularization_direct* p = (TikhonovRegularization_direct*) params;
	gsl_vector *f = gsl_vector_alloc(p->f_inf->size);
	p->runWithCurvature(f, x);
	double rho = p->rho;
	double eta = p->xi;
	*f_ = - p->rho/(p->xi*x*x) - p->L_gradient;
	*df_ = 2.0*rho / (x*x*eta) + (rho/(x*x*eta*eta) + 1.0/eta)/p->dxi_dlambda;
	gsl_vector_free(f);
	return;
}
inline double TikhonovRegularization_direct::getCurvature(double x, void * params){
	TikhonovRegularization_direct* p = (TikhonovRegularization_direct*)params;
	gsl_vector *f = gsl_vector_alloc(p->f_inf->size);
	double curvature = p->runWithCurvature(f, fabs(x));
	gsl_vector_free(f);
	return curvature;
}
inline double TikhonovRegularization_direct::getCurvature(const gsl_vector *x, void * params){
	TikhonovRegularization_direct* p = (TikhonovRegularization_direct*)params;
	gsl_vector *f = gsl_vector_alloc(p->f_inf->size);
	double curvature = p->runWithCurvature(f, exp(gsl_vector_get(x, 0)));
	gsl_vector_free(f);
	return curvature;
}

inline double TikhonovRegularization_direct::getNCP()const
{
	gsl_vector* data = gsl_vector_calloc(r->size);
	gsl_fft_real_wavetable * const real = gsl_fft_real_wavetable_alloc(r->size);
	gsl_fft_real_workspace * const work = gsl_fft_real_workspace_alloc(r->size);

	gsl_vector_memcpy(data, r);
	double ncp = 0.0;
	gsl_fft_real_transform(data->data, data->stride, data->size, real, work);

	size_t i=0, n = r->size;
	size_t nhalf = ceil(n/2)+1;
	gsl_vector * c = gsl_vector_calloc(nhalf-1);
	for(i=1; i<n-i; ++i){
		double hc_real = data->data[(2*i-1) * data->stride];
		double hc_imag = data->data[(2*i  ) * data->stride];
		double amp2 = hc_real *hc_real  + hc_imag *hc_imag ;
		gsl_vector_set(c, i-1, amp2);
	}
	if(i==n-i){
		double hc_real = data->data[(n-1) * data->stride];
		double hc_imag = 0.0;
		double amp2 = hc_real *hc_real  + hc_imag *hc_imag ;
		gsl_vector_set(c, i-1, amp2);
	}

//	IGORdata::write_itx(c, "fft_v.itx", "fft_v");
//	IGORdata::write_itx(r, "res_v.itx", "res_v");

	double sum = gsl_vector_get(c, 0);
	for(size_t i=1; i<c->size; ++i){
		gsl_vector_set(c, i, gsl_vector_get(c, i-1) + gsl_vector_get(c, i));
	}
	
	gsl_vector_scale(c, 1.0/gsl_vector_get(c, c->size-1));
//	IGORdata::write_itx(c, "cvec.itx", "cvec");

	gsl_vector_free(data);
	gsl_vector_free(c);
	gsl_fft_real_workspace_free(work);
	gsl_fft_real_wavetable_free(real);
	return ncp;
}

/*
inline double TikhonovRegularization_direct::runWithCurvature(gsl_vector* fnow, const double lambda, const double dlambda){
	double lambda_tmp = lambda * (1.0 + dlambda);
	run(fnow, lambda_tmp);
	evaluatePrefRes(fnow);
	double xi_prev = pref*pref;

	run(fnow, lambda);
	evaluatePrefRes(fnow);
	double xi = pref*pref;
	double rho = res*res;
	double xi_prime = (xi_prev - xi) / (lambda*dlambda);

	double lambda2 = lambda*lambda;
	double lambda4 = lambda2*lambda2;
	double curvature = 2.0 * xi * rho / xi_prime * (lambda2 * xi_prime * rho + 2.0 * lambda * xi * rho + lambda4 * xi * xi_prime) / pow(lambda4 * xi * xi + rho*rho, 1.5);
	return curvature;
}
*/

inline double TikhonovRegularization_direct::runWithCurvature(gsl_vector* fnow, const double lambda){
	gsl_matrix * Pinv_CC = gsl_matrix_calloc(fnow->size, fnow->size);
	gsl_vector * Cy_Px = gsl_vector_calloc(fnow->size);
	gsl_permutation * p = gsl_permutation_calloc(fnow->size);

	//	Pinv_CC = (lambda^2) L^t + L
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, lambda*lambda, L, L, 0.0, Pinv_CC);

	//	Cy_Px = At * d
	gsl_blas_dgemv(CblasTrans, 1.0, A, d, 0.0, Cy_Px);
	//	Cy_Px = At * d + Pinv * x_inf
	gsl_blas_dgemv(CblasNoTrans, 1.0, Pinv_CC, f_inf, 1.0, Cy_Px);

	//	PinvCC = A^t * A + Pinv
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 1.0, Pinv_CC);


	//	LU decomposition of Pinv_CC
	int signum;
	gsl_linalg_LU_decomp(Pinv_CC, p, &signum);

	gsl_linalg_LU_invert(Pinv_CC, p, cov);

	//	x = cov * Cy_Px
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov, Cy_Px, 0.0, fnow);

	evaluatePrefRes(fnow);
	//--------------------------------
	//	for determining the curvature
	//--------------------------------
	gsl_vector * Ax_d = gsl_vector_alloc(d->size);
	gsl_vector_memcpy(Ax_d, d);
	//	Ax_d = A*x - d
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, -1.0, Ax_d);

	gsl_vector * At_Ax_d = gsl_vector_calloc(fnow->size);
	gsl_blas_dgemv(CblasTrans, 1.0, A, Ax_d, 0.0, At_Ax_d);

	gsl_vector * z = gsl_vector_calloc(fnow->size);
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov, At_Ax_d, 0.0, z);

	//	Pinv_CC = L^t * L
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L, L, 0.0, Pinv_CC);
	//	Cy_Px = L^t * L * z
	gsl_blas_dgemv(CblasNoTrans, 1.0, Pinv_CC, z, 0.0, Cy_Px);

	gsl_blas_ddot(fnow, Cy_Px, &dxi_dlambda);
	dxi_dlambda*= 4.0 / lambda;

	xi = pref*pref;
	rho = res*res;

	double lambda2 = lambda*lambda;
	double lambda4 = lambda2*lambda2;
	double curvature = 2.0 * xi * rho / dxi_dlambda * (lambda2 * dxi_dlambda * rho + 2.0 * lambda * xi * rho + lambda4 * xi * dxi_dlambda) / pow(lambda4 * xi * xi + rho*rho, 1.5);

	gsl_matrix_free(Pinv_CC);
	gsl_vector_free(Cy_Px);
	gsl_permutation_free(p);
	gsl_vector_free(Ax_d);
	gsl_vector_free(At_Ax_d);
	gsl_vector_free(z);

	return curvature;
}

inline bool TikhonovRegularization_direct::run(gsl_vector* fnow, const double lambda){
	gsl_matrix * Pinv_CC = gsl_matrix_calloc(fnow->size, fnow->size);
	gsl_vector * Cy_Px = gsl_vector_calloc(fnow->size);
	gsl_permutation * p = gsl_permutation_calloc(fnow->size);

	//	Pinv_CC = (lambda^2) L^t + L
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, lambda*lambda, L, L, 0.0, Pinv_CC);

	//	Cy_Px = At * d
	gsl_blas_dgemv(CblasTrans, 1.0, A, d, 0.0, Cy_Px);
	//	Cy_Px = At * d + Pinv * x_inf
	gsl_blas_dgemv(CblasNoTrans, 1.0, Pinv_CC, f_inf, 1.0, Cy_Px);

	//	PinvCC = A^t * A + Pinv
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 1.0, Pinv_CC);


	//	LU decomposition of Pinv_CC
	int signum;
	gsl_linalg_LU_decomp(Pinv_CC, p, &signum);

	gsl_linalg_LU_invert(Pinv_CC, p, cov);

	//	x = cov * Cy_Px
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov, Cy_Px, 0.0, fnow);

	gsl_matrix_free(Pinv_CC);
	gsl_vector_free(Cy_Px);
	gsl_permutation_free(p);
	evaluatePrefRes(fnow);
	return true;
}
inline void TikhonovRegularization_direct::evaluatePrefRes(const gsl_vector* fnow){
	//---	results	---
	gsl_vector* tmp = gsl_vector_calloc(fnow->size);
	gsl_vector* tmp2 = gsl_vector_calloc(L->size1);
	//	tmp = f - finf
	gsl_vector_memcpy(tmp, fnow);
	gsl_vector_sub(tmp, f_inf);
	//	s1 = L * (f-finf)
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, tmp, 0.0, tmp2);
	pref = gsl_blas_dnrm2(tmp2);

	//	tmp = A*f - d
	gsl_vector_memcpy(r, d);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, -1.0, r);
	res = gsl_blas_dnrm2(r);

	xi = pref*pref;
	rho = res*res;

	gsl_vector_free(tmp);
	gsl_vector_free(tmp2);

	return;
}

/*
inline bool TikhonovRegularization_direct::run(gsl_vector* fnow, const double lambda){
	double ratio_pararrel = run_first(fnow, lambda);

	size_t iter = 0;
	while (ratio_pararrel > ratio_stop){
		iter++;
		ratio_pararrel = run_next(fnow, lambda);
		if (iter > max_iter) return false;
#ifdef DEBUG
		gsl_vector* tmp2 = gsl_vector_calloc(d->size);
		IGORdata::write_itx(smat, "smat.itx", "smat");
		IGORdata::write_itx(fnow, "fnow.itx", "fnow");
		IGORdata::write_itx(d, "dd.itx", "dd");
		IGORdata::write_itx(A, "AA.itx", "AA");
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, -1.0, tmp2);
		res = gsl_blas_dnrm2(tmp2);
		IGORdata::write_itx(tmp2, "tmp2.itx", "tmp2");
		gsl_vector_free(tmp2);
#endif
	}
	evaluatePrefRes(fnow);
}


inline double TikhonovRegularization_direct::run_first(gsl_vector* fnow, const double lambda){
	gsl_vector*c = gsl_vector_calloc(2);
	gsl_vector* tmp = gsl_vector_calloc(fnow->size);
	gsl_vector* tmp2 = gsl_vector_calloc(d->size);
	gsl_vector* tmp3 = gsl_vector_calloc(L->size1);

	gsl_vector_view s1_view = gsl_matrix_column(smat, 0);
	gsl_vector_view s2_view = gsl_matrix_column(smat, 1);
	gsl_vector_view s3_view = gsl_matrix_column(smat, 2);

	//	tmp = f - finf
	gsl_vector_memcpy(tmp, fnow);
	gsl_vector_sub(tmp, f_inf);
	//	tmp3 = L * (f-finf)
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, tmp, 0.0, tmp3);
	pref = gsl_blas_dnrm2(tmp3);
	//	s1 = Lt * L * (f-finf)
	gsl_blas_dgemv(CblasTrans, 1.0, L, tmp3, 0.0, &s1_view.vector);

	//	tmp = A*f - d
	gsl_vector_memcpy(tmp2, d);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, -1.0, tmp2);
	res = gsl_blas_dnrm2(tmp2);
	//	s2 = At * (Af - d)
	gsl_blas_dgemv(CblasTrans, 1.0, A, tmp2, 0.0, &s2_view.vector);

	gsl_matrix_view s12_view = gsl_matrix_submatrix(smat, 0, 0, fnow->size, 2);
	getCoefVector(fnow, lambda, &s12_view.matrix, c);

	//	s3 = s * c 
	gsl_blas_dgemv(CblasNoTrans, 1.0, &s12_view.matrix, c, 0.0, &s3_view.vector);

	gsl_vector_add(fnow, &s3_view.vector);

	gsl_vector_free(c);
	gsl_vector_free(tmp);
	gsl_vector_free(tmp2);
	gsl_vector_free(tmp3);


//#ifdef DEBUG
//	IGORdata::write_itx(smat, "smat.itx", "smat");
//#endif

	double dot = 0.0;
	gsl_blas_ddot(&s1_view.vector, &s2_view.vector, &dot);
	dot = fabs(dot);
	double norm2 = gsl_blas_dnrm2(&s1_view.vector)*gsl_blas_dnrm2(&s2_view.vector);

//#ifdef DEBUG
//	std::cerr << "dot = " << dot << "  norm" << sqrt(gsl_blas_dnrm2(&s1_view.vector)*gsl_blas_dnrm2(&s2_view.vector)) << std::endl;
//#endif

	return (norm2 - dot) / norm2;
}

inline double TikhonovRegularization_direct::run_next(gsl_vector* fnow, const double lambda){
	gsl_vector*c = gsl_vector_calloc(3);
	gsl_vector* tmp = gsl_vector_calloc(fnow->size);
	gsl_vector* tmp2 = gsl_vector_calloc(d->size);
	gsl_vector* tmp3 = gsl_vector_calloc(L->size1);

	gsl_vector_view s1_view = gsl_matrix_column(smat, 0);
	gsl_vector_view s2_view = gsl_matrix_column(smat, 1);
	gsl_vector_view s3_view = gsl_matrix_column(smat, 2);

	//	tmp = f - finf
	gsl_vector_memcpy(tmp, fnow);
	gsl_vector_sub(tmp, f_inf);
	//	tmp3 = L * (f-finf)
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, tmp, 0.0, tmp3);
	pref = gsl_blas_dnrm2(tmp3);
	//	s1 = Lt * L * (f-finf)
	gsl_blas_dgemv(CblasTrans, 1.0, L, tmp3, 0.0, &s1_view.vector);

	//	tmp = A*f - d
	gsl_vector_memcpy(tmp2, d);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, -1.0, tmp2);
	res = gsl_blas_dnrm2(tmp2);
	//	s2 = At * (Af - d)
	gsl_blas_dgemv(CblasTrans, 1.0, A, tmp2, 0.0, &s2_view.vector);

	getCoefVector(fnow, lambda, smat, c);

	//	s3 = s * c 
	gsl_blas_dgemv(CblasNoTrans, 1.0, smat, c, 0.0, tmp);
	gsl_vector_memcpy(&s3_view.vector, tmp);
	gsl_vector_add(fnow, tmp);

	gsl_vector_free(c);
	gsl_vector_free(tmp);
	gsl_vector_free(tmp2);
	gsl_vector_free(tmp3);

	double dot = 0.0;
	gsl_blas_ddot(&s1_view.vector, &s2_view.vector, &dot);
	dot = fabs(dot);
	double norm2 = gsl_blas_dnrm2(&s1_view.vector)*gsl_blas_dnrm2(&s2_view.vector);
	return (norm2 - dot) / norm2;
}

inline bool TikhonovRegularization_direct::getCoefVector(const gsl_vector * fnow, const double lambda, const gsl_matrix * s, gsl_vector * c)
{

	bool success_flag = true;
//----	derive H~_omega	---
	//	Ls = L * s
	gsl_matrix * Ls = gsl_matrix_calloc(L->size1, s->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, s, 0.0, Ls);

	//	Htilde_w = (Ls)t * Ls
	gsl_matrix * Htilde_w = gsl_matrix_calloc(s->size2, s->size2);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Ls, Ls, 0.0, Htilde_w);

//----	derive G~_omega	---
	//	fnow_minus_finf  = fnow - finf
	gsl_vector * fnow_minus_finf = gsl_vector_calloc(fnow->size);
	gsl_vector_memcpy(fnow_minus_finf, fnow);
	gsl_vector_sub(fnow_minus_finf, f_inf);

	//	L_fnow_minus_finf = L * (fnow - finf)
	gsl_vector * L_fnow_minus_finf = gsl_vector_calloc(L->size1);
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, fnow_minus_finf, 0.0, L_fnow_minus_finf);

	//	Gtilde_w = -(Ls)t * L(fnow - finf)
	gsl_vector * Gtilde_w = gsl_vector_calloc(s->size2);
	gsl_blas_dgemv(CblasTrans, -1.0, Ls, L_fnow_minus_finf, 0.0, Gtilde_w);

//----	derive H~_C---
	//	As = A * s1
	gsl_matrix * As = gsl_matrix_calloc(A->size1, s->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, s, 0.0, As);

	//	Htilde_c = (As)t * As
	gsl_matrix * Htilde_c = gsl_matrix_calloc(s->size2, s->size2);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, As, As, 0.0, Htilde_c);

//----	derive G~_C---
	//	d_minus_Afnow = d - A * fnow
	gsl_vector * d_minus_Afnow = gsl_vector_calloc(A->size1);
	gsl_vector_memcpy(d_minus_Afnow, d);
	gsl_blas_dgemv(CblasNoTrans, -1.0, A, fnow, 1.0, d_minus_Afnow);

	//	Gtilde_c = (As)t * (d - A * fnow)
	gsl_vector * Gtilde_c = gsl_vector_calloc(s->size2);
	gsl_blas_dgemv(CblasTrans, 1.0, As, d_minus_Afnow, 0.0, Gtilde_c);


#ifdef DEBUG
	IGORdata::write_itx(Htilde_w, "Htilde_w.itx", "Htilde_w");
	IGORdata::write_itx(s, "ss.itx", "ss");
#endif
//---	derive c	---
	//	Htilde_w = lambda^2 * Htilde_w + Htilde_c
	gsl_matrix_scale(Htilde_w, lambda*lambda);
	gsl_matrix_add(Htilde_w, Htilde_c);

	//	Gtilde_w = lambda^2 * Gtilde_w + Gtilde_c
	gsl_vector_scale(Gtilde_w, lambda*lambda);
	gsl_vector_add(Gtilde_w, Gtilde_c);
	
//	gsl_permutation*p = gsl_permutation_alloc(s->size2);
//	int signum;
//	gsl_linalg_LU_decomp(Htilde_w, p, &signum);
//	gsl_linalg_LU_solve(Htilde_w, p, Gtilde_w, c);
//	gsl_permutation_free(p);

	//	gsl_error handler
	gsl_error_handler_t * gsl_error_handler = gsl_set_error_handler_off();
	int error = gsl_linalg_cholesky_decomp(Htilde_w);
	gsl_set_error_handler(gsl_error_handler);
	if (error == GSL_EDOM)
		success_flag = false;
	else
	{	//	c = Htilde_w ^-1 Gtilde_w
		gsl_linalg_cholesky_solve(Htilde_w, Gtilde_w, c);
	}
//	gsl_linalg_cholesky_decomp(Htilde_w);
//	gsl_linalg_cholesky_solve(Htilde_w, Gtilde_w, c);

	gsl_matrix_free(Ls);
	gsl_matrix_free(Htilde_w);
	gsl_vector_free(fnow_minus_finf);
	gsl_vector_free(L_fnow_minus_finf);
	gsl_vector_free(Gtilde_w);

	gsl_matrix_free(As);
	gsl_matrix_free(Htilde_c);
	gsl_vector_free(d_minus_Afnow);
	gsl_vector_free(Gtilde_c);

	
	return success_flag;
}

*/
//----------old codes----------------
/*
bool TikhonovRegularization_direct::run1(gsl_vector* fnow, const double lambda){
const size_t snum = 4;
gsl_matrix * s = gsl_matrix_calloc(fnow->size, snum);

gsl_vector_view s1_view = gsl_matrix_column(s, 0);
gsl_vector_view s2_view = gsl_matrix_column(s, 2);
gsl_vector_view s3_view = gsl_matrix_column(s, 3);
gsl_vector_view s4_view = gsl_matrix_column(s, 4);

//	gsl_matrix_view s1_mview = gsl_matrix_submatrix(s, 0, 0, fnow->size, 1);
get1stEstimate(fnow, lambda, &s1_view.vector);


get1stEstimate(fnow, lambda, &s1_view.vector);
double c1;
get1stCoef(fnow, lambda, &s1_view.vector, &c1);

gsl_vector_scale(&s1_view.vector, c1);

//	cs1 = fnow + c1 * s1
gsl_vector_add(fnow, &s1_view.vector);
get1stEstimate(fnow, lambda, &s2_view.vector);

gsl_vector * c = gsl_vector_calloc(2);
bool success_flag = getCoefVector(fnow, lambda, s, c);

//	fnow = s * c + fnow
if(success_flag)
gsl_blas_dgemv(CblasNoTrans, 1.0, s, c, 1.0, fnow);

gsl_matrix_free(s);
gsl_vector_free(c);
return success_flag;
}

bool TikhonovRegularization_direct::get1stEstimate(const gsl_vector * fnow, const double lambda, gsl_vector * s1){
//	fnow_minus_finf  = fnow - finf
gsl_vector * fnow_minus_finf = gsl_vector_calloc(fnow->size);
gsl_vector_memcpy(fnow_minus_finf, fnow);
gsl_vector_sub(fnow_minus_finf, f_inf);

//	L_fnow_minus_finf = L * (fnow - finf)
gsl_vector * L_fnow_minus_finf = gsl_vector_calloc(L->size1);
gsl_blas_dgemv(CblasNoTrans, 1.0, L, fnow_minus_finf, 0.0, L_fnow_minus_finf);

//	d_minus_Af = d - A * f
gsl_vector * d_minus_Af = gsl_vector_calloc(A->size1);
gsl_vector_memcpy(d_minus_Af, d);
gsl_blas_dgemv(CblasNoTrans, -1.0, A, fnow, 1.0, d_minus_Af);

//	At_d_minus_Af = At (d - Af)
gsl_vector * At_d_minus_Af = gsl_vector_calloc(fnow->size);
gsl_blas_dgemv(CblasTrans, 1.0, A, d_minus_Af, 1.0, At_d_minus_Af);

//	s1 = -lambda^2 * Lt * L * (f - finf) + At (d - Af)
gsl_vector_memcpy(s1, At_d_minus_Af);
gsl_blas_dgemv(CblasTrans, -lambda*lambda, L, L_fnow_minus_finf, 1.0, s1);


#ifdef DEBUG
IGORdata::write_itx(s1, "s1tmp.itx", "s1tmp");
#endif

gsl_vector_free(fnow_minus_finf);
gsl_vector_free(L_fnow_minus_finf);
gsl_vector_free(d_minus_Af);
gsl_vector_free(At_d_minus_Af);

return true;
}

bool TikhonovRegularization_direct::get1stCoef(const gsl_vector * fnow, const double lambda, const gsl_vector * s1, double* c1){
//----	derive H~_omega	---
//	Ls = L * s1
gsl_vector * Ls = gsl_vector_calloc(L->size1);
gsl_blas_dgemv(CblasNoTrans, 1.0, L, s1, 0.0, Ls);

//	Htilde_w = (Ls)t * Ls
double Htilde_w;
gsl_blas_ddot(Ls, Ls, &Htilde_w);

//----	derive G~_omega	---
//	fnow_minus_finf  = fnow - finf
gsl_vector * fnow_minus_finf = gsl_vector_calloc(fnow->size);
gsl_vector_memcpy(fnow_minus_finf, fnow);
gsl_vector_sub(fnow_minus_finf, f_inf);

//	L_fnow_minus_finf = L * (fnow - finf)
gsl_vector * L_fnow_minus_finf = gsl_vector_calloc(L->size1);
gsl_blas_dgemv(CblasNoTrans, 1.0, L, fnow_minus_finf, 0.0, L_fnow_minus_finf);

//	Gtilde_w = -(Ls)t * L(fnow - finf)
double Gtilde_w;
gsl_blas_ddot(Ls, L_fnow_minus_finf, &Gtilde_w);
Gtilde_w *= -1.0;

//----	derive H~_C---
//	As = A * s1
gsl_vector * As = gsl_vector_calloc(A->size1);
gsl_blas_dgemv(CblasNoTrans, 1.0, A, s1, 0.0, As);

//	Htilde_c = (As)t * As
double Htilde_c;
gsl_blas_ddot(As, As, &Htilde_c);

//----	derive G~_C---
//	d_minus_Afnow = d - A * fnow
gsl_vector * d_minus_Afnow = gsl_vector_calloc(A->size1);
gsl_vector_memcpy(d_minus_Afnow, d);
gsl_blas_dgemv(CblasNoTrans, -1.0, A, fnow, 1.0, d_minus_Afnow);

//	Htilde_c = (As)t * As
double Gtilde_c;
gsl_blas_ddot(As, d_minus_Afnow, &Gtilde_c);

*c1 = (lambda*lambda * Gtilde_w + Gtilde_c) / (lambda*lambda * Htilde_w + Htilde_c);

gsl_vector_free(Ls);
gsl_vector_free(fnow_minus_finf);
gsl_vector_free(L_fnow_minus_finf);
gsl_vector_free(As);
gsl_vector_free(d_minus_Afnow);

return true;
}
*/

