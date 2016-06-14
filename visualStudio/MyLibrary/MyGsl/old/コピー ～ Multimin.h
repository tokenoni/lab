#ifndef __MULTIMIN_H__
#define __MULTIMIN_H__

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

		std::vector<double> params_result;
		std::vector<double> params_err;
		std::vector<std::vector<double>> ScaledCovMatrix;

	//---	Constructor	---
		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, double step_ratio = 1.0e-4, double relative_err= 1.0e-4);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err, std::vector<std::string> params_name_);

		template<class TyF_ >
		Multimin
			( TyF_& function, std::vector<double>& params, double step_ratio, double relative_err, std::vector<std::string> params_name_);

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

	//---	class for adopting a functor to C function	---
	template < class TyF_ >
	class Cfunction{
	public:
		typedef double (*C_FUNCTION)(const gsl_vector*, void*);	//	address of the functor
		Cfunction(TyF_& function_, size_t num_params){	fAddress = &function_;	std_v.resize(num_params);}
		C_FUNCTION adopt(){ return &glue;}
	private:
		static TyF_* fAddress;
		static std::vector<double> std_v;
		static double glue(const gsl_vector* gsl_v, void* dummy){ 
	//		std::vector<double> std_v(gsl_v->size);
			for(size_t i=0; i<gsl_v->size; i++)	std_v[i] = gsl_vector_get(gsl_v, i);
			return (*fAddress)(std_v);
		}
	};

template< class TyF_ > TyF_* Cfunction< TyF_ >::fAddress=0;
template< class TyF_ > std::vector<double> Cfunction< TyF_ >::std_v(0);
};
#include "Multimin.inl"

#endif