#ifndef __MYGSL_GAUSS_FIT_HPP__
#define __MYGSL_GAUSS_FIT_HPP__

#include <gsl\gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include "constants.hpp"

namespace mygsl{
	//-------------------------
	//	fit by a Gauss function with
	//	y0 + A/(sqrt(2pi) * sigma) * exp( - (x-x0)^2/(2*width^2))
	//
	//--------------------------
	class GaussFit{
	public:
		GaussFit(const size_t dataSize);
		~GaussFit();
		const size_t size;
		int fit(const std::vector<double>& x_, const std::vector<double>& data_, const std::vector<double>& sigma, const std::vector<double>& initial_guess, double epsabs = 1.0e-4, double epsrel = 1.0e-4);
		int fit(const std::vector<double>& x_, const std::vector<double>& data_, const std::vector<double>& initial_guess, double epsabs = 1.0e-4, double epsrel = 1.0e-4);

		//	if fixedFlag == 1, the parameter is fixed
		void setFixedParameters(const std::vector<int>& fixedFlag_);

		//	returns the result
		std::vector<double> getFitResult();
		double getY0()const{ return gsl_vector_get(params, 0); }
		double getA()const{ return gsl_vector_get(params, 1); }
		double getX0()const{ return gsl_vector_get(params, 2); }
		double getWidth()const{ return gsl_vector_get(params, 3); }
		double getChi()const{ return chi; }
	protected:
		int fitCore(double epsabs, double epsrel);
		gsl_matrix* cov;
		gsl_multifit_fdfsolver* s;
		const gsl_multifit_fdfsolver_type* T;
		double chi;
		gsl_multifit_function_fdf func;
	public:	
		gsl_vector* params;	//	size: 4
		std::vector<double> x;
		std::vector<double> data;
		std::vector<double> sigma_inv;
		std::vector<int> fixedFlag;
		std::vector<double> result, result_err;
		bool isSigmaProvided;
	public:
		static int f(const gsl_vector*v, void* params, gsl_vector* f);
		static int df(const gsl_vector*v, void* params, gsl_matrix* J);
		static int fdf(const gsl_vector*v, void* params, gsl_vector* f, gsl_matrix* J);
	};

	//-------------------------
	//	fit by a Airy function with
	//	y0 + A/(pi * sigma) * ( sin((x-x0)/width)/((x-x0)/width) )^2
	//
	//--------------------------
	class AiryFit : public GaussFit{
	public:
		AiryFit(const size_t dataSize);
		std::vector<double> getFitResult();
		static int f(const gsl_vector*v, void* params, gsl_vector* f);
		static int df(const gsl_vector*v, void* params, gsl_matrix* J);
		static int fdf(const gsl_vector*v, void* params, gsl_vector* f, gsl_matrix* J);
	};

	class LorentzFit : public GaussFit{
	public:
		static int f(const gsl_vector*v, void* params, gsl_vector* f);
		static int df(const gsl_vector*v, void* params, gsl_matrix* J);
		static int fdf(const gsl_vector*v, void* params, gsl_vector* f, gsl_matrix* J);
	};

};

#include "GaussFit.inl"

#endif