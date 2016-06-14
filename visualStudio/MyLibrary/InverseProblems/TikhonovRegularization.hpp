#pragma once
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_roots.h>
//#include <IGOR/IGORitx.hpp>

//#define DEBUG
#ifdef DEBUG
#include <IGOR/IGORitx.hpp>
#endif

class TikhonovRegularization
{
public:
	TikhonovRegularization
		(const gsl_matrix * A_, const gsl_matrix * L_,
		const gsl_vector * d_,
		const gsl_vector * f_inf_,
		const double ratio_stop_ = 1.0e-6, const size_t max_iter_ = 100
		);
	~TikhonovRegularization();

	//	iterative method. This method is useful for large f->size case
	//	with specified Lambda gradient
	bool run
		(gsl_vector* finit,
		const double lambdaGradient,
		const double lambdaMin, const double lambdaMax, const double dlambdaRatio);

	//	with specified Lambda
	bool run(gsl_vector* finit, const double lambda);


	//	Direct method. This method is useful for small f->size case
	//	returns the curvature 
	double runDirectWithCurvature(gsl_vector* finit, const double lambda);
	double runDirectWithCurvature(gsl_vector* finit, const double lambda, const double dlambda);
	//	with specified Lambda
	bool runDirect(gsl_vector* finit, const double lambda);
	//	with specified L-gradient 
	bool runDirect(gsl_vector* finit, const double L_gradient, double& lambda_init, const size_t iterMax = 10, const double lambda_epsrel = 0.1);
	//	get covariance
	const gsl_matrix* getCov()const{ return cov; }

	void evaluatePrefRes(const gsl_vector* fnow);
	double getPref()const{ return pref; }
	double getRes()const{ return res; }
	double getNCP()const;

protected:
	double run_first(gsl_vector* finit, const double lambda);
	double run_next(gsl_vector* finit, const double lambda);
	bool getCoefVector(const gsl_vector * fnow, const double lambda, const gsl_matrix * s, gsl_vector * c);

	const gsl_matrix * const A;		//	the equation of system
	const gsl_vector * const d;		//	the answer of the matrix
	const gsl_matrix * const L;		//	regularization matrix
	const gsl_vector * const f_inf;	//	ideal (expected) solution
	gsl_matrix* const cov;	//	covariant matrix of the solution
	gsl_matrix* const smat;
	gsl_vector* const r;	//	residual vector
	double pref, res;
	double L_gradient;		//	specified L-gradient for the corner detection

	const double ratio_stop;
	const size_t max_iter;

public:
	static double getGradient_f(double x, void * params);
	static double getGradient_df(double x, void * params);
	static void getGradient_fdf(double x, void * params, double * f, double * df);
};

class TikhonovRegularizationWIthNCP: public TikhonovRegularization{
public:
	TikhonovRegularizationWIthNCP
		(const gsl_matrix * A_, const gsl_matrix * L_,
		const gsl_vector * d_,
		const gsl_vector * f_inf_,
		const double ratio_stop_ = 1.0e-6, const size_t max_iter_ = 100
		);
	~TikhonovRegularizationWIthNCP();

	double getNCP()const;
private:
	gsl_fft_real_wavetable * const real;
	gsl_fft_halfcomplex_wavetable * const hc;
	gsl_fft_real_workspace * const work;
};

#include "TikhonovRegularization.inl"
