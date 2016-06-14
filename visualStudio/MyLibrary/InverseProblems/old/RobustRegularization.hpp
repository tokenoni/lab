#pragma once
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_roots.h>

#define DEBUG
#ifdef DEBUG
#include <IGOR/IGORitx.hpp>
#endif

class RobustRegularization
{
public:
	RobustRegularization
		(const gsl_matrix * A_, const gsl_matrix * L_,
		 const gsl_vector * d_,
		 const gsl_vector * f_inf_,
		 const double nu_ = 1.0,
		 const double ratio_stop_ = 1.0e-6, const size_t max_iter_ = 100
		);
	~RobustRegularization();

	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda, const size_t iterMax, const double lambda_epsrel);

protected:
	//	member functions
	bool set_Jac_and_Hes(const gsl_vector* fnow, const double lambda);
	double update_f(gsl_vector* fnow, const double lambda);

	//	member variables 
	const gsl_matrix * const A;		//	the equation of system
	const gsl_vector * const d;		//	the answer of the matrix
	const gsl_matrix * const L;		//	regularization matrix
	const gsl_vector * const f_inf;	//	ideal (expected) solution
	gsl_matrix* const cov;	//	covariant matrix of the solution
	gsl_vector* const jac;	//	Jacobian vector
	gsl_matrix* const hes;	//	Hessian matrix of the solution
	gsl_matrix* const LtL;	//	L^t * L

	const double nu;	//	parameter for the Student's t-distribution
	const double ratio_stop;
	const size_t max_iter;
	const size_t size_i, size_j;


	//	function for the relation between Af value and sigma
	double get_Sigma(const double x)const;
	double get_dSigma_dx(const double x)const;
	double get_d2Sigma_dx2(const double x)const;
};

#include "RobustRegularization.inl"
