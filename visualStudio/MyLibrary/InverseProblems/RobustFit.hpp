#pragma once
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_math.h>
#include <vector>

//#define DEBUG
#ifdef DEBUG
#include <IGOR/IGORitx.hpp>
#else
#define GSL_RANGE_CHECK_OFF
#endif

//--------------------------------------------------//
//													//
//		base class for the Robust fit				//
//													//
//--------------------------------------------------//
class RobustFit{
public:
	RobustFit(const size_t size_f, const size_t size_d);
	virtual ~RobustFit();

	//	for initializing
	void setL(const gsl_matrix* L);
	void setA(const gsl_matrix* A_){ gsl_matrix_memcpy(A, A_); }
	void setData(const gsl_vector* data);

	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);

	void calculateCov(const gsl_vector* fnow);
protected:
	//	member functions
	virtual double update_f(gsl_vector* fnow, const double lambda);
	std::vector<double> getRes(const gsl_vector * frslt)const;
	std::vector<double> getResult(const gsl_vector * frslt)const;

	//	member variables 
	gsl_matrix * const A;		//	the equation of system
	gsl_vector * const d;		//	the answer of the matrix
	gsl_vector * const f;
	gsl_vector * const jac;	//	Jacobian vector
	gsl_matrix * const hes;	//	Hessian matrix of the solution
	gsl_matrix * const cov;	//	covariant matrix of the solution
	gsl_matrix * const LtL;	//	L^t * L
	const size_t size_d;	//	number of data point
	const size_t size_f;	//	number of variables

	//	virtual fuinctions
	//	function for the relation between Af value and sigma
	virtual double get_Sigma(const double y, const size_t i)const = 0;
	virtual double get_dSigma_dy(const double y, const size_t i)const = 0;
	virtual double get_d2Sigma_dy2(const double y, const size_t i)const = 0;
	//	Jacobian and Hessian
	virtual bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda) = 0;
	//	regularization term
	virtual bool add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda);

	//	for regularization term
	virtual double get_lambda_weight(const size_t k)const{ return 1.0; }

private:
	RobustFit(void);
	RobustFit(const RobustFit& obj);
	
public:
	//---	get functions ---
	const size_t get_size_d()const{ return size_d; }
	const size_t get_size_f()const{ return size_f; }
	const gsl_matrix* get_matrix_A()const{ return A; }
	const gsl_vector* get_vector_f()const{ return f; }
	const gsl_vector* get_vector_d()const{ return d; }
};

//------------------------------------------------------------------//
//																	//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is fixed						//
//------------------------------------------------------------------//
//	The noise amplitude of the signal is assumed to be given.
//	Initially, the noise amplitude is set to be 1.
//	It should be modified through "setSigma" or "setAlpha" functions
class RobustFit_GaussianPDF : public RobustFit{
public:
	RobustFit_GaussianPDF(const size_t size_f_, const size_t size_d_);
	virtual ~RobustFit_GaussianPDF();

	void setSigma(const gsl_vector* sigma);
	void setSigma(const std::vector<double>& sigma);
	void setAlpha(const gsl_vector* alpha_){gsl_vector_memcpy(alpha, alpha_);}
protected:
	//	function 
	virtual bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda);
	//	a
	virtual double get_alpha(const double y, const size_t i)const;
	gsl_vector * const alpha;

private:	//	unused functions
	double get_Sigma(const double y, const size_t i)const{	return 0.0;	};
	double get_dSigma_dy(const double y, const size_t i)const{ return 0.0; };
	double get_d2Sigma_dy2(const double y, const size_t i)const{ return 0.0; };

public:
	//---	get functions ---
	const gsl_vector* get_vector_alpha()const{ return alpha; }
};

//------------------------------------------------------------------//
//																	//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is function of y				//
//------------------------------------------------------------------//
//	The noise amplitude of the signal is given as a function of y.
//	function getAlpha and functions calculating its derivatives should be inherited
class RobustFit_GaussianPDF_flexible : public RobustFit{
public:
	RobustFit_GaussianPDF_flexible(const size_t size_f_, const size_t size_d_)
		: RobustFit(size_f_, size_d_){};
	virtual ~RobustFit_GaussianPDF_flexible(){};
protected:
	//	function 
	virtual bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda);
	//	a
	virtual double get_alpha(const double y, const size_t i)const;
	//	da/dy
	virtual double get_dalpha_dy(const double y, const size_t i)const;
	//	d^2a/dy^2
	virtual double get_d2alpha_dy2(const double y, const size_t i)const;

	//	da_dy / a
	virtual double get_dalpha_dy_over_alpha(const double y, const size_t i)const;
	//	d^2a_dy^2 / a^2
	virtual double get_d2alpha_dy2_over_alpha2(const double y, const size_t i)const;

};
//	when the signal distribution is given by Student's t-distribution
class RobustFit_students_tPDF : public RobustFit{
public:
	RobustFit_students_tPDF(const size_t size_f_, const size_t size_d_, const double nu_) :
		nu(nu_), RobustFit(size_f_, size_d_){};
	virtual ~RobustFit_students_tPDF(){};
protected:

	//	function 
	bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda);
	double get_eta(const double y, const size_t i)const;
	double get_deta_dy(const double y, const size_t i)const;
	double get_d2eta_dy2(const double y, const size_t i)const;
	const double nu;
};

//	when the signal is composed with the randomely distributed outliers;
class RobustFit_OffsetGaussianPDF : public RobustFit_GaussianPDF_flexible{
public:
	RobustFit_OffsetGaussianPDF
		(const size_t size_f_, const size_t size_d_, const double range_)
		: RobustFit_GaussianPDF_flexible(size_f_, size_d_), deltaYinv(1.0 / range_){};
	virtual ~RobustFit_OffsetGaussianPDF(){};

protected:
	//	function 
	virtual bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda);

	//	outlier probability at y index i and its derivatives against y
	virtual double get_p(const double y, const size_t i)const = 0;
	virtual double get_dp_dy(const double y, const size_t i)const = 0;
	virtual double get_d2p_dy2(const double y, const size_t i)const = 0;

	const double deltaYinv;
};

//------------------------------------------------------------------//
//		L1 regularization															//
//		The signal distribution is assumed as a Gauss function		//
//							with sigma is fixed						//
//------------------------------------------------------------------//
//	The noise amplitude of the signal is assumed to be given.
//	Initially, the noise amplitude is set to be 1.
//	It should be modified through "setSigma" or "setAlpha" functions
//	(L1 regularization)
class RobustFit_GaussianPDF_LASSO : public RobustFit_GaussianPDF{
public:
	RobustFit_GaussianPDF_LASSO
		(const size_t size_f_, const size_t size_d_, const size_t size_L);
	virtual ~RobustFit_GaussianPDF_LASSO();
	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);

protected:
	//	regularization term
	virtual bool add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda);
	//	reweighted least square
	virtual void update_Reg_Weight(const gsl_vector * finit);
	//	for regularization term
	virtual void get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD);
	virtual double get_lambda_weight(const size_t k)const{ return gsl_vector_get(reg_weight, k); }

	//	reguralization weight for LASSO regularization
	gsl_vector * const reg_weight;
	gsl_vector * const Lf;
	//	regularization matrix
	gsl_matrix * const L;
public:
	//---	get functions ---
	const gsl_vector * get_vector_reg_weight()const{ return reg_weight; }
	const gsl_matrix * get_matrix_L()const{ return L; }
	const size_t get_size_L()const{ return L->size1; }
};


//	M-estimators
//	(L2 regularization)
class RobustFit_Mestimator : public RobustFit_GaussianPDF{
public:
	RobustFit_Mestimator(const size_t size_f_, const size_t size_d_, const double a);
	virtual ~RobustFit_Mestimator();

	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);

protected:
	//	reweighted least square
	virtual void updateWeight(const gsl_vector * finit);
	//	overriding function 
	double get_alpha(const double y, const size_t i)const;

	gsl_vector * const weight;
	const double a_inv;

	static double getWeight(double u_over_a);
public:
	//---	get functions ---
	const double get_a_inv()const{ return a_inv; }
	const gsl_vector * get_vector_weight()const{ return weight; }
};

//	M-estimators (L1 regularization)
class RobustFit_Mestimator_LASSO : public RobustFit_Mestimator{
public:
	RobustFit_Mestimator_LASSO
		(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L);
	virtual ~RobustFit_Mestimator_LASSO();
	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);
	//	with specified Lambda
	bool runWithoutMestimate(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);

protected:
	//	regularization term
	virtual bool add_Jac_and_Hes_reg(const gsl_vector* fnow, const double lambda);

	//	reweighted least square
	virtual void update_Reg_Weight(const gsl_vector * finit);
	//	for regularization term
	virtual void get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD);
	virtual double get_lambda_weight(const size_t k)const{ return gsl_vector_get(reg_weight, k); }

	//	reguralization weight for LASSO regularization
	gsl_vector * const reg_weight;
	gsl_vector * const Lf;
	//	regularization matrix
	gsl_matrix * const L;
public:
	//---	get functions ---
	const gsl_vector * get_vector_reg_weight()const{ return reg_weight; }
	const gsl_matrix * get_matrix_L()const{ return L; }
	const size_t get_size_L()const{ return L->size1; }
};



//	(L0.5 regularization)
class RobustFit_Mestimator_Lhalf : public RobustFit_Mestimator_LASSO{
public:
	RobustFit_Mestimator_Lhalf
		(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L);
	virtual ~RobustFit_Mestimator_Lhalf(){};

protected:
	//	for regularization term
	virtual void get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD);
};


//	(L0 regularization)
class RobustFit_Mestimator_UniformPenalty : public RobustFit_Mestimator_LASSO{
public:
	RobustFit_Mestimator_UniformPenalty
		(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L);
	virtual ~RobustFit_Mestimator_UniformPenalty(){};

protected:
	//	for regularization term
	virtual void get_lambda_weight_from_fLLf(const gsl_vector* fLLf, const double MAD);
};



//	S-estimators (under construction)
class RobustFit_Sestimator : public RobustFit_Mestimator{
public:
	RobustFit_Sestimator
		(const size_t size_f_, const size_t size_d_, const double a, const size_t size_L);
	virtual ~RobustFit_Sestimator();

protected:
	//	reweighted least square
	virtual void updateWeight(const gsl_vector * finit);
	double a_inv;
};


#include "RobustFit.inl"
