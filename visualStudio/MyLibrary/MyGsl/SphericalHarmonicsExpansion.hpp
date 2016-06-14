#ifndef __MYGSL_SPHERICAL_HARMONICS_HPP__
#define __MYGSL_SPHERICAL_HARMONICS_HPP__

#include <gsl\gsl_sf.h>
#include <gsl\gsl_sf_legendre.h>
#include <mygsl\complex.hpp>
#include <MyGsl\constants.hpp>

namespace mygsl{
	//	class for spherical harmonics expansion: 
	//	Currently, the maximum rank is l = 4
	class SphericalHarmonicsExpansion{
	public:
		SphericalHarmonicsExpansion(const size_t rank_);
		~SphericalHarmonicsExpansion(void);
		const size_t rank_max;
		const size_t rank;
		const size_t size;

		void run(const double theta, const double phi);
		void run(const double x, const double y, const double z);
		//	x, y, z shoul be a unit vector
		const complex* getValues()const{return sh_values;}
	private:
		complex* sh_values;
		
		// l=0
		const double coef_l0;
		// l=1
		const double coef_l1_m0;
		const double coef_l1_m1;
		// l=2
		const double coef_l2_m0;
		const double coef_l2_m1;
		const double coef_l2_m2;
		// l=3
		const double coef_l3_m0;
		const double coef_l3_m1;
		const double coef_l3_m2;
		const double coef_l3_m3;
		// l=4
		const double coef_l4_m0;
		const double coef_l4_m1;
		const double coef_l4_m2;
		const double coef_l4_m3;
		const double coef_l4_m4;

	//	static functions
	public:
		static size_t rankToSize(const size_t rank){return (rank+1)*(rank+1);};
	};
};
#include "SphericalHarmonicsExpansion.inl"
#endif