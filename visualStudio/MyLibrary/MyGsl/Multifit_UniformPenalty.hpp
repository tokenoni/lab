#ifndef __MULTIFIT_UNIFORM_PENALTY_HPP__
#define __MULTIFIT_UNIFORM_PENALTY_HPP__
#include <gsl/gsl_multimin.h>
#include <vector>
namespace mygsl{
	class Multifit_UniformPenalty{
	public:
		bool initialize();
	private:
		//	evaluation value of the 1st and 2nd term
		double E1st, E2nd; 
		Multifit_UniformPenalty_params params;
	
	public:
	};

	class Multifit_UniformPenalty_params{
	public:
		Multifit_UniformPenalty_params(void){};
		~Multifit_UniformPenalty_params(void){};

		bool initialize(const std::vector<double>& data_, const std::vector<double>& sigma_, 
						const std::vector<std::vector<std::vector<double>>>& phi_ijk_,
			            const std::vector<std::vector<double>>& bik_, const std::vector<double>& delta2_, const std::vector<double>& lambda_);
		void free();
		
		size_t num_j, num_i, num_k;

		double* data;
		double* sigma;
		double*** phi_ijk;
		double** bik;
		double* delta2;
		double* lambda;

		bool isLogarithmic, isPositiveConstrained;
	};

};
#include "Multifit_UniformPenalty.inl"
#endif
