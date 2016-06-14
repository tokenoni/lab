
//------------------------------------------//
//											//
//				Robust Fit					//
//											//
//------------------------------------------//

inline RobustFit::RobustFit(const size_t size_f_, const size_t size_d_)
:
size_d(size_d_), size_f(size_f_),
A(gsl_matrix_calloc(size_d_, size_f_)),
d(gsl_vector_calloc(size_d_)),
f(gsl_vector_calloc(size_f_)),
jac(gsl_vector_calloc(size_f_)),
hes(gsl_matrix_calloc(size_f_, size_f_)),
cov(gsl_matrix_calloc(size_f_, size_f_)),
LtL(gsl_matrix_calloc(size_f_, size_f_))
{
}


inline RobustFit::~RobustFit()
{
	gsl_matrix_free(A);
	gsl_vector_free(d);
	gsl_vector_free(f);
	gsl_vector_free(jac);
	gsl_matrix_free(hes);
	gsl_matrix_free(cov);
	gsl_matrix_free(LtL);
}

inline void RobustFit::setL(const gsl_matrix* L){
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, L, L, 0.0, LtL);
}

void inline RobustFit::setData(const gsl_vector* data){
	gsl_vector_memcpy(d, data);
}

//	with specified Lambda
inline bool RobustFit::run
(gsl_vector* finit, const double lambda,
const size_t iterMax, const double epsrel)
{
	//	seeking f value which gives zero Jacobian
	for (size_t t = 0; t < iterMax; ++t){
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			return true;
	}
	return false;
}
inline void RobustFit::calculateCov(const gsl_vector* fnow){
	int signum;
	gsl_permutation* p = gsl_permutation_alloc(hes->size1);
	gsl_linalg_LU_decomp(hes, p, &signum);
	gsl_linalg_LU_invert(hes, p, cov);
	gsl_permutation_free(p);
}

inline double RobustFit::update_f(gsl_vector* fnow, const double lambda){
	gsl_vector* df = gsl_vector_calloc(fnow->size);
	set_Jac_and_Hes(fnow, lambda);

#ifdef DEBUG
	gsl_vector* Af = gsl_vector_calloc(d->size);
	IGORdata::write_itx(d, "d.itx", "d");
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
	gsl_vector_add(fnow, df);
	gsl_permutation_free(p);

#ifdef DEBUG
	IGORdata::write_itx(df, "df.itx", "df");
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, fnow, 0.0, Af);
	IGORdata::write_itx(Af, "Af.itx", "Af");
	IGORdata::write_itx(cov, "cov.itx", "cov");
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

inline std::vector<double> RobustFit::getResult(const gsl_vector * frslt)const{
	gsl_vector* Af = gsl_vector_calloc(d->size);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, frslt, 0.0, Af);

	std::vector<double> rslt(size_d);
	for (size_t i = 0; i < size_d; ++i)
		rslt[i] = gsl_vector_get(Af, i);
	gsl_vector_free(Af);
	return rslt;
}

inline std::vector<double> RobustFit::getRes(const gsl_vector * frslt)const{
	gsl_vector* Af = gsl_vector_alloc(d->size);
	gsl_vector_memcpy(Af, d);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, frslt, -1.0, Af);

	std::vector<double> rslt(size_d);
	for (size_t i = 0; i < size_d; ++i)
		rslt[i] = gsl_vector_get(Af, i);
	gsl_vector_free(Af);
	return rslt;
}

//	regularization term
inline bool RobustFit::add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda){
	for (size_t k = 0; k < size_f; ++k){
		//	for jacobian		
		double ft_LtL = 0.0;
		for (size_t j = 0; j < size_f; ++j)
			ft_LtL += gsl_vector_get(fnow, j) * gsl_matrix_get(LtL, j, k);

		double lambda_weight = lambda * get_lambda_weight(k);
		gsl_vector_set(jac, k, gsl_vector_get(jac, k) + lambda_weight*ft_LtL);

		//	for hessian
		for (size_t j = k; j < size_f; ++j){
			gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + lambda_weight * gsl_matrix_get(LtL, j, k));
			gsl_matrix_set(hes, k, j, gsl_matrix_get(hes, j, k));
		}
	}
#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif
	return true;
}


//-----------------		inheriting classes		----------------------

//------------------------------------------------------------------//
//																	//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is fixed						//
//------------------------------------------------------------------//
inline RobustFit_GaussianPDF::RobustFit_GaussianPDF(const size_t size_f_, const size_t size_d_)
:
RobustFit(size_f_, size_d_),
alpha(gsl_vector_calloc(size_d_))
{}

inline RobustFit_GaussianPDF::~RobustFit_GaussianPDF()
{
	gsl_vector_free(alpha);
};

inline void RobustFit_GaussianPDF::setSigma(const gsl_vector* sigma){
	for (size_t i = 0; i < size_d; ++i){
		double tmp = 1.0 / gsl_vector_get(sigma, i);
		gsl_vector_set(alpha, i, tmp*tmp);
	}
}
inline void RobustFit_GaussianPDF::setSigma(const std::vector<double>& sigma){
	for (size_t i = 0; i < size_d; ++i){
		double tmp = 1.0 / sigma[i];
		gsl_vector_set(alpha, i, tmp*tmp);
	}
}
inline double RobustFit_GaussianPDF::get_alpha(const double y, const size_t i)const{
	return gsl_vector_get(alpha, i);
}

//	calculate jacobian and hessian
inline bool RobustFit_GaussianPDF::set_Jac_and_Hes(const gsl_vector* fnow, const double lambda){
	gsl_vector_set_zero(jac);
	gsl_matrix_set_zero(hes);

	for (size_t i = 0; i < size_d; ++i){
		//	calculate value at i
		double t = 0.0;
		for (size_t j = 0; j < size_f; ++j)
			t += gsl_matrix_get(A, i, j) * gsl_vector_get(fnow, j);

		//	data
		double y = gsl_vector_get(d, i);

		//	calculate alpha at i
		double alpha = get_alpha(t, i);
		for (size_t j = 0; j < size_f; ++j){
			//	Jacobian
			double jac_ij = 0.0;
			jac_ij += 2.0 * alpha * (y - t) * gsl_matrix_get(A, i, j);
			gsl_vector_set(jac, j, gsl_vector_get(jac, j) + jac_ij);

			//	Hessian
			for (size_t k = j; k < size_f; ++k){
				double hes_ijk = 0.0;
				hes_ijk -= 2.0 * alpha *gsl_matrix_get(A, i, k)*gsl_matrix_get(A, i, j);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + hes_ijk);
			}
		}
	}

	gsl_vector_scale(jac, -1.0);
	gsl_matrix_scale(hes, -1.0);
	//	regularization term
	for (size_t k = 0; k < size_f; ++k){
		for (size_t j = k; j < size_f; ++j)
			gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, k, j));
	}
	add_Jac_and_Hes_reg(fnow, lambda);
	return true;
}

//------------------------------------------------------------------//
//																	//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is function of y				//
//------------------------------------------------------------------//

inline double RobustFit_GaussianPDF_flexible::get_alpha(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	return 1.0 / (2.0*sigma*sigma);
}
inline double RobustFit_GaussianPDF_flexible::get_dalpha_dy(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	return -1.0 / (sigma*sigma*sigma)*get_dSigma_dy(y, i);
}
inline double RobustFit_GaussianPDF_flexible::get_d2alpha_dy2(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	double sigma2 = sigma*sigma;
	double dsigma_dy = get_dSigma_dy(y, i);
	return 3.0 / (sigma2*sigma2)*dsigma_dy*dsigma_dy - 1.0 / (sigma2*sigma) * get_d2Sigma_dy2(y, i);
}
inline double RobustFit_GaussianPDF_flexible::get_dalpha_dy_over_alpha(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	return -2.0 / sigma * get_dSigma_dy(y, i);
}
inline double RobustFit_GaussianPDF_flexible::get_d2alpha_dy2_over_alpha2(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	double sigma2 = sigma*sigma;
	double dsigma_dy = get_dSigma_dy(y, i);
	return 12.0 * dsigma_dy*dsigma_dy - 2.0 / (sigma) * get_d2Sigma_dy2(y, i);
}

//	calculate jacobian and hessian
inline bool RobustFit_GaussianPDF_flexible::set_Jac_and_Hes(const gsl_vector* fnow, const double lambda){
	gsl_vector_set_zero(jac);
	gsl_matrix_set_zero(hes);

	for (size_t i = 0; i < size_d; ++i){
		//	calculate value at i
		double t = 0.0;
		for (size_t j = 0; j < size_f; ++j)
			t += gsl_matrix_get(A, i, j) * gsl_vector_get(fnow, j);

		//	data
		double y = gsl_vector_get(d, i);

		//	calculate alpha at i
		double alpha = get_alpha(t, i);
		for (size_t j = 0; j < size_f; ++j){
			double dalpha_dfj = get_dalpha_dy(t, i) * gsl_matrix_get(A, i, j);
			double dalpha_dfj_over_alpha = get_dalpha_dy_over_alpha(t, i) * gsl_matrix_get(A, i, j);

			//	Jacobian
			double jac_ij = 0.0;
			jac_ij += 0.5 * dalpha_dfj_over_alpha;
			jac_ij -= (y - t)*(y - t) * dalpha_dfj;
			jac_ij += 2.0 * alpha * (y - t) * gsl_matrix_get(A, i, j);
			gsl_vector_set(jac, j, gsl_vector_get(jac, j) + jac_ij);


			//	Hessian
			for (size_t k = j; k < size_f; ++k){
				double dalpha_dfk = get_dalpha_dy(t, i) * gsl_matrix_get(A, i, k);
				double d2alpha_dfj_dfk = get_d2alpha_dy2(t, i) * gsl_matrix_get(A, i, j) *gsl_matrix_get(A, i, k);
				double dalpha_dfk_over_alpha = get_dalpha_dy_over_alpha(t, i) * gsl_matrix_get(A, i, k);
				double d2alpha_dfj_dfk_over_alpha = get_d2alpha_dy2_over_alpha2(t, i) * gsl_matrix_get(A, i, j) *gsl_matrix_get(A, i, k);

				double hes_ijk = 0.0;
				hes_ijk -= 0.5 * dalpha_dfj_over_alpha*dalpha_dfk_over_alpha;
				hes_ijk += 0.5 * d2alpha_dfj_dfk_over_alpha;
				hes_ijk -= (y - t)*(y - t)*d2alpha_dfj_dfk;
				hes_ijk += 2.0 * (y - t)*gsl_matrix_get(A, i, j)*dalpha_dfk;
				hes_ijk += 2.0 * (y - t)*gsl_matrix_get(A, i, k)*dalpha_dfj;
				hes_ijk -= 2.0 * alpha *gsl_matrix_get(A, i, k)*gsl_matrix_get(A, i, j);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + hes_ijk);
			}
		}
	}

	gsl_vector_scale(jac, -1.0);
	gsl_matrix_scale(hes, -1.0);
	//	regularization term
	for (size_t k = 0; k < size_f; ++k){
		for (size_t j = k; j < size_f; ++j)
			gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, k, j));
	}

	add_Jac_and_Hes_reg(fnow, lambda);
	return true;
}



inline bool RobustFit_students_tPDF::set_Jac_and_Hes(const gsl_vector* fnow, const double lambda){
	gsl_vector_set_zero(jac);
	gsl_matrix_set_zero(hes);

	for (size_t i = 0; i < size_d; ++i){
		//	calculate value at i
		double t = 0.0;
		for (size_t j = 0; j < size_f; ++j)
			t += gsl_matrix_get(A, i, j) * gsl_vector_get(fnow, j);
		
		//	calculate eta at i
		double eta = get_eta(t, i);

		double y = gsl_vector_get(d, i);
		double coef = 1.0 / (1.0 + eta * (y - t)*(y - t));

		for (size_t k = 0; k < size_f; ++k){
			double deta_dfk = get_deta_dy(t, i) * gsl_matrix_get(A, i, k);
			double dt_dfk = gsl_matrix_get(A, i, k);
			double parenthesis_k = (-2.0*eta*(y - t)*dt_dfk + (y - t)*(y - t)*deta_dfk);
			double jac1 =
				+  0.5 / eta * deta_dfk
				- (0.5*nu + 0.5) * coef * parenthesis_k;
			
			gsl_vector_set(jac, k, gsl_vector_get(jac, k) + jac1);

			for (size_t j = k; j < size_f; ++j){
				double deta_dfj = get_deta_dy(t, i) * gsl_matrix_get(A, i, j);
				double d2eta_dfj_dfk = get_d2eta_dy2(t, i) * gsl_matrix_get(A, i, k) * gsl_matrix_get(A, i, j);
				double dt_dfj = gsl_matrix_get(A, i, j);
				double parenthesis_j = (-2.0*eta*(y - t)*dt_dfj + (y - t)*(y - t)*deta_dfj);

				double term1 = -0.5 / (eta*eta)*deta_dfj*deta_dfk;
				double term2 =  0.5 / eta * d2eta_dfj_dfk;
				
				double hes1 = term1 + term2 + 
					(0.5*nu + 0.5) * coef * 
						(
						+ coef * parenthesis_j * parenthesis_k
						+ 2.0*(y-t)*(dt_dfk*deta_dfj + dt_dfj*deta_dfk) 
						- 2.0*eta*dt_dfj*dt_dfk
						- (y-t)*(y-t)*d2eta_dfj_dfk
						);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + hes1);
			}
		}
	}

#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif

	gsl_vector_scale(jac, -1.0);
	gsl_matrix_scale(hes, -1.0);
	//	regularization term
	for (size_t k = 0; k < size_f; ++k){
		for (size_t j = k; j < size_f; ++j)
			gsl_matrix_set(hes, k, j, gsl_matrix_get(hes, j, k));
	}

	add_Jac_and_Hes_reg(fnow, lambda);
#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif
	return true;
}

inline double RobustFit_students_tPDF::get_eta(const double y, const size_t i)const{
	return (nu - 2.0)*get_Sigma(y, i);
}
inline double RobustFit_students_tPDF::get_deta_dy(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	return 2.0*(nu - 2.0)*get_dSigma_dy(y, i);
}
inline double RobustFit_students_tPDF::get_d2eta_dy2(const double y, const size_t i)const{
	double sigma = get_Sigma(y, i);
	double dsigma_dy = get_dSigma_dy(y, i);
	return 2.0*(nu - 2.0)*(dsigma_dy*dsigma_dy + sigma * get_d2Sigma_dy2(y, i));
}

inline bool RobustFit_OffsetGaussianPDF::set_Jac_and_Hes(const gsl_vector* fnow, const double lambda){
	gsl_vector_set_zero(jac);
	gsl_matrix_set_zero(hes);

	std::vector<double> Bj(size_f, 0.0), ej(size_f, 0.0);

	for (size_t i = 0; i < size_d; ++i){
		//	calculate value at i

		double t = 0.0;
		for (size_t j = 0; j < size_f; ++j)
			t += gsl_matrix_get(A, i, j) * gsl_vector_get(fnow, j);

		//	calculate eta at i
		double alpha = get_alpha(t, i);
		double B = sqrt(alpha / M_PI);
		double y_minus_t = gsl_vector_get(d, i) - t;
		double E = exp(-alpha*y_minus_t*y_minus_t);
		double p = get_p(t, i);
		double epsilon_i = (1.0 - p)*B*E + p*deltaYinv;

		for (size_t k = 0; k < size_f; ++k){
			double dalpha_dfk = get_dalpha_dy(t, i) * gsl_matrix_get(A, i, k);
			Bj[k] = 1.0 / (2.0*sqrt(alpha * M_PI)) * dalpha_dfk;
			ej[k] = -dalpha_dfk * y_minus_t * y_minus_t
				+ 2.0*alpha * y_minus_t * gsl_matrix_get(A, i, k);
		}

		for (size_t k = 0; k < size_f; ++k){
			double Ek = E*ej[k];
			double dt_dfk = gsl_matrix_get(A, i, k);
			double dalpha_dfk = get_dalpha_dy(t, i) * dt_dfk;
			double dp_dfk = get_dp_dy(t, i) * dt_dfk;
			double jac1 = (1.0 - p) * (Bj[k] * E + B*Ek) - dp_dfk * (B * E - deltaYinv);

			gsl_vector_set(jac, k, gsl_vector_get(jac, k) + jac1 / epsilon_i);

			for (size_t j = k; j < size_f; ++j){
				double dt_dfj = gsl_matrix_get(A, i, j);
				double d2alpha_dfk_dfj = get_d2alpha_dy2(t, i) * dt_dfj*dt_dfk;
				double dalpha_dfj = get_dalpha_dy(t, i) * dt_dfj;
				double Ej = E*ej[j];
				double Bjk = 1.0 / (2.0*sqrt(alpha * M_PI)) * d2alpha_dfk_dfj - 1.0 / (4.0*alpha*sqrt(alpha*M_PI)) * dalpha_dfj*dalpha_dfk;
				double Ejk = Ek * ej[j]
					+ E*(-d2alpha_dfk_dfj * y_minus_t * y_minus_t
					+ 2.0 * y_minus_t * (dalpha_dfj * dt_dfk + dalpha_dfk * dt_dfj)
					- 2.0 * alpha * dt_dfj*dt_dfk);

				double dp_dfj = get_dp_dy(t, i) * dt_dfj;
				double d2p_dfj_dfk = get_d2p_dy2(t, i) * dt_dfk * dt_dfj;

				double term1 = (1.0 - p)*(Bjk * E + Bj[j] * Ek + Bj[k] * Ej + B*Ejk);
				double term2 = -dp_dfj * (Bj[j] * E + B * Ej)
					- dp_dfk * (Bj[k] * E + B * Ek)
					- d2p_dfj_dfk * B * E
					+ d2p_dfj_dfk * deltaYinv;
				double term3 = ((1.0 - p) * (Bj[k] * E + B*Ek) - dp_dfk * (B * E - deltaYinv))
					* ((1.0 - p) * (Bj[j] * E + B*Ej) - dp_dfj * (B * E - deltaYinv));

				double hes1 = (term1 + term2) / epsilon_i - term3 / (epsilon_i*epsilon_i);

				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + hes1);
			}
		}
	}

#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif

	gsl_vector_scale(jac, -1.0);
	gsl_matrix_scale(hes, -1.0);
	//	regularization term
	for (size_t k = 0; k < size_f; ++k){
		for (size_t j = k; j < size_f; ++j)
			gsl_matrix_set(hes, k, j, gsl_matrix_get(hes, j, k));
	}

	add_Jac_and_Hes_reg(fnow, lambda);
	return true;
}

//------------------------------------------------------------------//
//		L1 regularization															//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is fixed						//
//------------------------------------------------------------------//
inline RobustFit_GaussianPDF_LASSO::RobustFit_GaussianPDF_LASSO
(const size_t size_f_, const size_t size_d_, const size_t size_L)
:
RobustFit_GaussianPDF(size_f_, size_d_),
reg_weight(gsl_vector_calloc(size_L)),
L(gsl_matrix_calloc(size_L, size_f)), Lf(gsl_vector_calloc(size_L))
{
	gsl_vector_set_all(reg_weight, 1.0);
}
inline RobustFit_GaussianPDF_LASSO::~RobustFit_GaussianPDF_LASSO()
{
	gsl_vector_free(reg_weight);
	gsl_vector_free(Lf);
	gsl_matrix_free(L);
}

//	with specified Lambda
inline bool RobustFit_GaussianPDF_LASSO::run
(gsl_vector* finit, const double lambda,
const size_t iterMax, const double epsrel)
{

	size_t t0;
	update_f(finit, lambda);
	for (t0 = 0; t0 < iterMax; ++t0){
		//	seeking f value which gives zero Jacobian
		update_Reg_Weight(finit);
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			break;
	}
	return false;
}
//	regularization term
inline bool RobustFit_GaussianPDF_LASSO::add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda){
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, fnow, 0.0, Lf);

	for (size_t i = 0; i < Lf->size; ++i){
		double lambda_weight = lambda * gsl_vector_get(reg_weight, i);
		double x = gsl_vector_get(Lf, i);
		//	for jacobian		
		for (size_t j = 0; j < size_f; ++j){
			double Jj_tmp = gsl_matrix_get(L, i, j) * x;
			gsl_vector_set(jac, j, gsl_vector_get(jac, j) + lambda_weight * Jj_tmp);
		}

		//	for hessian
		for (size_t j = 0; j < size_f; ++j){
			for (size_t k = 0; k < size_f; ++k){
				double Hjk_tmp = gsl_matrix_get(L, i, j) * gsl_matrix_get(L, i, k);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + lambda_weight * Hjk_tmp);
			}
		}
	}
	return true;
}

inline void RobustFit_GaussianPDF_LASSO::update_Reg_Weight(const gsl_vector * finit){
	gsl_permutation *p = gsl_permutation_alloc(Lf->size);

	gsl_blas_dgemv(CblasNoTrans, 1.0, L, finit, 0.0, Lf);

	for (size_t i = 0; i < Lf->size; ++i){
		gsl_vector_set(Lf, i, fabs(gsl_vector_get(Lf, i)));
	}
	gsl_sort_vector_index(p, Lf);
	double MAD = gsl_vector_get(Lf, gsl_permutation_get(p, Lf->size / 2));

	get_lambda_weight_from_fLLf(Lf, 1.4826 * MAD);

#ifdef DEBUG
	IGORdata::write_itx(Lf, "LtLf.itx", "LtLf");
	IGORdata::write_itx(reg_weight, "reg_weight.itx", "reg_weight");
#endif
	gsl_permutation_free(p);
}
//--------------------------------------------------//
//													//
//				Robust Regression					//
//													//
//--------------------------------------------------//
inline RobustFit_Mestimator::RobustFit_Mestimator(const size_t size_f_, const size_t size_d_, const double a)
:
RobustFit_GaussianPDF(size_f_, size_d_),
weight(gsl_vector_calloc(size_d_)),
a_inv(1.0 / a)
{
	gsl_vector_set_all(weight, 1.0);
}

inline RobustFit_Mestimator::~RobustFit_Mestimator()
{
	gsl_vector_free(weight);
}
inline double RobustFit_Mestimator::get_alpha(const double y, const size_t i)const{
	return RobustFit_GaussianPDF::get_alpha(y, i) * gsl_vector_get(weight, i);
}

//	with specified Lambda
inline bool RobustFit_Mestimator::run
(gsl_vector* finit, const double lambda,
const size_t iterMax, const double epsrel)
{
	//	seeking f value which gives zero Jacobian
	update_f(finit, lambda);
	for (size_t t = 0; t < iterMax; ++t){
		updateWeight(finit);
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			break;
	}
	return false;
}

inline void RobustFit_Mestimator::updateWeight(const gsl_vector * finit){
	gsl_vector* res = gsl_vector_alloc(d->size);
	gsl_permutation *p = gsl_permutation_alloc(d->size);

	gsl_vector_memcpy(res, d);
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, finit, -1.0, res);
	
	for (size_t i = 0; i < size_d; ++i){
		if (gsl_vector_get(res, i) < 0.0)
			gsl_vector_set(res, i, fabs(gsl_vector_get(res, i)));
	}

	//	update scale
	gsl_sort_vector_index(p, res);
	double MAD = gsl_vector_get(res, gsl_permutation_get(p, size_d / 2));
	double scale_inv = 1.0 / (1.4826 * MAD);

	//	update weight
	for (size_t i = 0; i < size_d; ++i){
		gsl_vector_set(weight, i, getWeight(gsl_vector_get(res, i) * scale_inv * a_inv));
	}
	gsl_vector_free(res);
	gsl_permutation_free(p);

#ifdef DEBUG
	IGORdata::write_itx(weight, "weight.itx", "weight");
#endif
}

inline double RobustFit_Mestimator::getWeight(double u_over_a){
	double one_minus_u_over_a2 = 1.0 - u_over_a*u_over_a;
	if (one_minus_u_over_a2 > 0.0)
		return one_minus_u_over_a2 * one_minus_u_over_a2;
	else return 0.0;
}

//----------------------------------------
//		for LASSO regularization
//----------------------------------------

inline RobustFit_Mestimator_LASSO::RobustFit_Mestimator_LASSO
(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L)
:
RobustFit_Mestimator(size_f_, size_d_, a),
reg_weight(gsl_vector_calloc(size_L)),
L(gsl_matrix_calloc(size_L, size_f)), Lf(gsl_vector_calloc(size_L))
{
	gsl_vector_set_all(reg_weight, 1.0);
}

inline RobustFit_Mestimator_LASSO::~RobustFit_Mestimator_LASSO(){
	gsl_vector_free(reg_weight);
	gsl_vector_free(Lf);
	gsl_matrix_free(L);
}

//	with specified Lambda
inline bool RobustFit_Mestimator_LASSO::run
(gsl_vector* finit, const double lambda,
const size_t iterMax, const double epsrel)
{
	size_t t0;
	update_f(finit, lambda);
	for (t0 = 0; t0 < iterMax; ++t0){
		//	seeking f value which gives zero Jacobian
		for (size_t t1 = 0; t1 < iterMax; ++t1){
			updateWeight(finit);
			double norm = update_f(finit, lambda);
			if (norm < epsrel)
				break;
		}
		update_Reg_Weight(finit);
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			break;
	}
	return false;
}

inline bool RobustFit_Mestimator_LASSO::
runWithoutMestimate(gsl_vector* finit, const double lambda, const size_t iterMax, const double epsrel)
{
	size_t t0;
	update_f(finit, lambda);
	for (t0 = 0; t0 < iterMax; ++t0){
		//	seeking f value which gives zero Jacobian
		update_Reg_Weight(finit);
		double norm = update_f(finit, lambda);
		if (norm < epsrel)
			break;
	}
	return false;
}

//	regularization term
inline bool RobustFit_Mestimator_LASSO::add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda){
	gsl_blas_dgemv(CblasNoTrans, 1.0, L, fnow, 0.0, Lf);

	for (size_t i = 0; i < Lf->size; ++i){
		double lambda_weight = lambda * gsl_vector_get(reg_weight, i);
		double x = gsl_vector_get(Lf, i);
		//	for jacobian		
		for (size_t j = 0; j < size_f; ++j){
			double Jj_tmp = gsl_matrix_get(L, i, j) * x;
			gsl_vector_set(jac, j, gsl_vector_get(jac, j) + lambda_weight * Jj_tmp);
		}

		//	for hessian
		for (size_t j = 0; j < size_f; ++j){
			for (size_t k = 0; k < size_f; ++k){
				double Hjk_tmp = gsl_matrix_get(L, i, j) * gsl_matrix_get(L, i, k);
				gsl_matrix_set(hes, j, k, gsl_matrix_get(hes, j, k) + lambda_weight * Hjk_tmp);
			}
		}
	}
#ifdef DEBUG
	IGORdata::write_itx(jac, "jac.itx", "jac");
	IGORdata::write_itx(hes, "hes.itx", "hes");
#endif
	return true;
}

inline void RobustFit_Mestimator_LASSO::update_Reg_Weight(const gsl_vector * finit){
	gsl_permutation *p = gsl_permutation_alloc(Lf->size);

	gsl_blas_dgemv(CblasNoTrans, 1.0, L, finit, 0.0, Lf);
	
	for (size_t i = 0; i < Lf->size; ++i){
		gsl_vector_set(Lf, i, fabs(gsl_vector_get(Lf, i)));
	}
	gsl_sort_vector_index(p, Lf);
	double MAD = gsl_vector_get(Lf, gsl_permutation_get(p, Lf->size / 2));

	get_lambda_weight_from_fLLf(Lf, 1.4826 * MAD);

#ifdef DEBUG
	IGORdata::write_itx(Lf, "LtLf.itx", "LtLf");
	IGORdata::write_itx(reg_weight, "reg_weight.itx", "reg_weight");
#endif
	gsl_permutation_free(p);
}

inline void RobustFit_Mestimator_LASSO::get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD)
{
	for (size_t i = 0; i < Lf->size; ++i){
		gsl_vector_set(reg_weight, i, 
			MAD / sqrt(gsl_vector_get(fLLf, i)*gsl_vector_get(fLLf, i) + MAD*MAD)
			);
	}
}

//---------------------------------------------//
//			for L0.5 regularization			   //
//---------------------------------------------//
inline RobustFit_Mestimator_Lhalf::RobustFit_Mestimator_Lhalf
(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L)
:
RobustFit_Mestimator_LASSO(size_f_, size_d_, a, size_L){}

inline void RobustFit_Mestimator_Lhalf::get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD)
{
	for (size_t i = 0; i < size_f; ++i){
		gsl_vector_set(reg_weight, i, MAD / (pow(gsl_vector_get(fLLf, i), 0.75) + MAD));
	}
}

//---------------------------------------------//
//		for Uniform Penalty regularization	   //
//---------------------------------------------//
inline RobustFit_Mestimator_UniformPenalty::RobustFit_Mestimator_UniformPenalty
(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L)
:
RobustFit_Mestimator_LASSO(size_f_, size_d_, a, size_L){}

inline void RobustFit_Mestimator_UniformPenalty::get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD)
{
	for (size_t i = 0; i < size_f; ++i){
		gsl_vector_set(reg_weight, i, MAD / (gsl_vector_get(fLLf, i) + MAD));
	}
}
