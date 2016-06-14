inline RobustRegularization::RobustRegularization
(const gsl_matrix * A_, const gsl_matrix * L_,
const gsl_vector * d_,
const gsl_vector * f_inf_,
const double nu_,
const double ratio_stop_, const size_t max_iter_)
:
A(A_), L(L_), d(d_), f_inf(f_inf_),
nu(nu_),
ratio_stop(ratio_stop_), max_iter(max_iter_),
size_i(A->size1), size_j(f_inf->size),
cov(gsl_matrix_calloc(f_inf->size, f_inf->size)),
jac(gsl_vector_calloc(f_inf->size)),
hes(gsl_matrix_calloc(f_inf->size, f_inf->size)),
LtL(gsl_matrix_calloc(f_inf->size, f_inf->size))
{
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L, L, 0.0, LtL);
}


inline RobustRegularization::~RobustRegularization()
{
	gsl_matrix_free(cov);
	gsl_vector_free(jac);
	gsl_matrix_free(hes);
	gsl_matrix_free(LtL);
}

//	with specified Lambda
inline bool RobustRegularization::run(gsl_vector* finit, const double lambda, const size_t iterMax, const double epsrel)
{
	//	seeking f value which gives zero Jacobian
	for (size_t t = 0; t < max_iter; ++t){
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			return true;
	}
	return false;
}

inline double RobustRegularization::update_f(gsl_vector* fnow, const double lambda){
	gsl_vector* df = gsl_vector_calloc(fnow->size);
	set_Jac_and_Hes(fnow, lambda);

#ifdef DEBUG
	gsl_vector* Af = gsl_vector_calloc(d->size);
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
	IGORdata::write_itx(fnow, "fnow.itx", "fnow");
	IGORdata::write_itx(LtL, "LtL.itx", "LtL");
#endif

	gsl_vector_scale(jac, -1.0);
	
	//	solve 
//	gsl_linalg_cholesky_decomp(hes);
//	gsl_linalg_cholesky_solve(hes, jac, df);
	
	int signum;
	gsl_permutation* p = gsl_permutation_alloc(hes->size1);
	gsl_linalg_LU_decomp(hes, p, &signum);
	gsl_linalg_LU_solve(hes, p, jac, df);
	gsl_permutation_free(p);
	gsl_vector_add(fnow, df);

#ifdef DEBUG
	IGORdata::write_itx(df, "df.itx", "df");
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, 0.0, Af);
	IGORdata::write_itx(Af, "Af.itx", "Af");
	gsl_vector_free(Af);
#endif

//	for (size_t j = 0; j < size_j; ++j){
//		double fj = gsl_vector_get(fnow, j);
//		double dfj = gsl_vector_get(df, j);
//		if (fj > 0.1){
//			if (dfj / fj < 0.0)
//				gsl_vector_set(fnow, j, fj * (exp(5.0*dfj / fj) - 1.0)*0.1+fj);
//			else
//				gsl_vector_set(fnow, j, fj + fj * 0.1*log(1.0 + 5.0*dfj / fj));
//		}
//		else 
//			gsl_vector_set(fnow, j, fj + 0.5 * dfj);
//	}
	double norm = gsl_blas_dnrm2(df);
	gsl_vector_free(df);
	return norm / gsl_blas_dnrm2(fnow);
}

inline double RobustRegularization::get_Sigma(const double x)const{
	return 3000.0;
}
inline double RobustRegularization::get_dSigma_dx(const double x)const{
	return 0.0;
}
inline double RobustRegularization::get_d2Sigma_dx2(const double x)const{
	return 0.0;
}

//	calculate jacobian and hessian
inline bool RobustRegularization::set_Jac_and_Hes(const gsl_vector* fnow, const double lambda){
	gsl_vector_set_zero(jac);
	gsl_matrix_set_zero(hes);

	for (size_t i = 0; i < size_i; ++i){
		//	calculate value at i
		double t = 0.0;
		for (size_t j = 0; j < size_j; ++j)
			t += gsl_matrix_get(A, i, j) * gsl_vector_get(fnow, j);
		
		//	calculate eta at i
		double eta = get_Sigma(t);

		double u = gsl_vector_get(d, i);
		double coef = 1.0 / (1.0 + eta * (u - t)*(u - t));

		for (size_t k = 0; k < size_j; ++k){
			double deta_dfk = get_dSigma_dx(t) * gsl_matrix_get(A, i, k);
			double dt_dfk = gsl_matrix_get(A, i, k);
			double parenthesis_k = (2.0*eta*(t - u)*dt_dfk + (t - u)*(t - u)*deta_dfk);
			double jac1 =
				- 0.5 * 1.0 / eta * deta_dfk +
				(0.5*nu + 0.5) * coef * parenthesis_k;
			
			gsl_vector_set(jac, k, gsl_vector_get(jac, k) + jac1);

			for (size_t j = k; j < size_j; ++j){
				double deta_dfj = get_dSigma_dx(t) * gsl_matrix_get(A, i, j);
				double d2eta_dfj_dfk = get_d2Sigma_dx2(t) * gsl_matrix_get(A, i, k) * gsl_matrix_get(A, i, j);
				double dt_dfj = gsl_matrix_get(A, i, j);
				double parenthesis_j = (2.0*eta*(t - u)*dt_dfj + (t - u)*(t - u)*deta_dfj);

				double term1 =  0.5 / (eta*eta)*deta_dfj*deta_dfk;
				double term2 = -0.5 / eta * d2eta_dfj_dfk;
				
				double hes1 = term1 + term2 + 
					(0.5*nu + 0.5) * coef * 
						(
						- coef * parenthesis_j * parenthesis_k
						- 2.0*(t-u)*(dt_dfk*deta_dfj + dt_dfj*deta_dfk) 
						+ 2.0*eta*dt_dfj*dt_dfk
						+ (t - u)*(t-u)*d2eta_dfj_dfk
						);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + hes1);
			}
		}
	}

#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif


	//	regularization term
	for (size_t k = 0; k < size_j; ++k){
		//	for jacobian		
		double ft_LtL = 0.0;
		for (size_t j = 0; j < size_j; ++j)
			ft_LtL += gsl_vector_get(fnow, j) * gsl_matrix_get(LtL, j, k);
		
		gsl_vector_set(jac, k, gsl_vector_get(jac, k) + lambda*ft_LtL);

		//	for hessian
		for (size_t j = k; j < size_j; ++j){
			gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + lambda * gsl_matrix_get(LtL, j, k));
			gsl_matrix_set(hes, k, j, gsl_matrix_get(hes, j, k));
		}
		
	}
#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif
	return true;
}
