#ifndef __MULTIMIN_HPP__
#define __MULTIMIN_HPP__

#include <string>
#include <vector>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

namespace mygsl{
	class Multimin
	{
	public:
		enum{
			Simplex2,
			Simplex
		};

		enum IsPrint{
			Print,
			NoPrint
		}isprint;

		std::vector<double> params_result;
		std::vector<double> params_err;
		std::vector<std::vector<double>> ScaledCovMatrix;

	//---	Constructor	---
		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err, IsPrint isprint_=Print);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, double step_ratio = 1.0e-4, double relative_err= 1.0e-4, IsPrint isprint_=Print);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err, std::vector<std::string> params_name_, IsPrint isprint_=Print);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, double step_ratio, double relative_err, std::vector<std::string> params_name_, IsPrint isprint_=Print);

		~Multimin(void){};

	private:
		std::vector<std::string> params_name;
		std::vector<double> sigma;
	
		double MAXerr, MAXderr;
		size_t MAXiter;

		template< class TyF_ >
		int process( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err);


		//---	tip functions	---
		void print_state(size_t iter, gsl_multimin_fminimizer *s, double size);
	
		static void gsl2std_vector(const gsl_vector *gsl_v, std::vector<double>& std_v){
			if(gsl_v->size != std_v.size()) std_v.resize(gsl_v->size);
			for(size_t i=0; i<gsl_v->size; i++)	std_v[i] = gsl_vector_get(gsl_v, i);}
		static std::vector<double> gsl2std_vector(const gsl_vector *gsl_v){
			std::vector<double> std_v(gsl_v->size);
			for(size_t i=0; i<gsl_v->size; i++)	std_v[i] = gsl_vector_get(gsl_v, i);
			return std_v;}

		static void std2gsl_vector(const std::vector<double>& std_v, gsl_vector *gsl_v){
			for(size_t i=0; i<gsl_v->size; i++)	gsl_vector_set(gsl_v, i, std_v[i]);}
		static gsl_vector*	std2gsl_vector(const std::vector<double>& std_v){
			gsl_vector* gsl_v = gsl_vector_calloc(std_v.size());
			for(size_t i=0; i<std_v.size(); i++) gsl_vector_set(gsl_v, i, std_v[i]);
			return gsl_v;}

	};

	template < class TyF_ >
	class FunctionAndParameters{
	public:
		FunctionAndParameters(TyF_ * functor_, const size_t num){
			params.resize(num); func = functor_;}

		static double call_back(const gsl_vector* pgsl, void* func_and_params){
			FunctionAndParameters< TyF_ > *fap= reinterpret_cast< FunctionAndParameters< TyF_ >  * >(func_and_params);
			for(size_t i=0; i<pgsl->size; ++i) fap->params[i] = gsl_vector_get(pgsl, i);
			return fap->func->operator()(fap->params);
		}
	private:
		TyF_ * func;
		std::vector<double> params;
	};
};
#include "Multimin.inl"

#endif