#pragma once
#ifndef __MULTIFIT_HPP__
#define __MULTIFIT_HPP__
#include <vector>
#include <string>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <iostream>

namespace mygsl{
	class Multifit
	{
	public:
		enum Fixed_NotFixed{
			Fixed,
			NotFixed
		};

		std::vector<double> params_result;
		std::vector<double> params_err;
		std::vector<std::vector<double>> ScaledCovMatrix;

		//---	several types of the constructors are available	---
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, size_t NumData);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<double>& sigma_);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, size_t NumData, std::vector<Fixed_NotFixed> FixedParams_);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData, std::vector<Fixed_NotFixed> FixedParams_);
		template< class TyF_ >
		Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_);


		~Multifit(void){};

	private:
	//---	data set for gsl_proceedure	---
		size_t NumData;
		std::vector<double> sigma;
		std::vector<Fixed_NotFixed> FixedParams;
		double MAXerr, MAXderr;
		size_t MAXiter;

		template< class TyF_ >
		void Process( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData, std::vector<Fixed_NotFixed> FixedParams_);
		template< class TyF_ >
		void Process( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_);
	//---	core functions	---
		template< class TyF_ >
		static int calc_function_f  (const gsl_vector * params, void *pFunc, gsl_vector * f);
		template< class TyF_ >
		static int calc_function_df (const gsl_vector * params, void *pFunc, gsl_matrix * J);
		template< class TyF_ >
		static int calc_function_fdf(const gsl_vector * params, void *pFunc, gsl_vector * f, gsl_matrix * J);
		template< class TyF_ >
		void CoreMultifit(TyF_& function, std::vector<double>& params);
	//---	data for output	----
		double chi, chisq;
		std::vector<std::vector<double>> params_cov;
		std::vector<std::string> params_name;


	//---	tip functions	---
		void print_state(size_t iter, gsl_multifit_fdfsolver * s);

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

	template <class TyF_ > 
	class FuncParamData{
	public:	
		TyF_* function;
		std::vector<double> param, data, sigma;
		std::vector<double> param_dif, data_plus, data_minus;
		std::vector<Multifit::Fixed_NotFixed> FixedParams;
		FuncParamData(TyF_ &function_, std::vector<double> &param_, std::vector<double> &sigma_, std::vector<Multifit::Fixed_NotFixed> FixedParams_){
			function = &function_; param=param_; sigma = sigma_; FixedParams = FixedParams_;
			param_dif.resize(param.size());
			data.resize(sigma.size());
			data_plus.resize(sigma.size());	data_minus.resize(sigma.size());
		};
	};
};
#include "Multifit.inl"
#endif