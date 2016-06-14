#ifndef __MYGSL_GAUSSIAN_BASIS_HPP__
#define __MYGSL_GAUSSIAN_BASIS_HPP__

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <cfloat>

namespace mygsl{
	class Gaussian;	//	just a definition 
	class Overlap;
//-----------------------------------------------------------------------//
//																		 //
//				class for One Gaussian Functions						 //
//																		 //
//-----------------------------------------------------------------------//
	
	//---	gaussian basis:	x^l y^m z^n exp(-alpha(x^2+y^2+z^2)	---
	class OneGaussian{
	public:
		//	constructors	
		OneGaussian(void){ alpha=0.0; coef=0.0;}
		~OneGaussian(void){};
		OneGaussian(const OneGaussian& ob);
		OneGaussian& operator = (const OneGaussian& obj);

		//	overriden operators
		double operator * (const OneGaussian& obj)const;
		double operator * (const Gaussian& obj)const;
		bool operator == (const OneGaussian& obj)const;
		bool operator != (const OneGaussian& obj)const;

		//	initializing funcitons
		void set(double exponent, size_t xi_, size_t yi_, size_t zi_=0);	
		void set(double exponent, size_t xi_=0);
		
		double get(double x, double y, double z)const;	//	3 dimensional
		double get(double x, double y)const;				//	2 dimensional
		double get(double x)const;						//	1 dimensional
		double getArea()const{return area;};
		double getMoment_x()const{return moment[0];};
		double getMoment_y()const{return moment[1];};
		double getMoment_z()const{return moment[2];};
		size_t l()const{	return xi;}
		size_t m()const{	return yi;}
		size_t n()const{	return zi;}
		double a()const{return alpha;}
		double c()const{return coef;}

	private:
		double alpha;			//	exponent of tha basis
		double coef;			//	coefficient
		size_t xi, yi, zi;		//	power index for the x, y, z
		double area;
		double moment[3];
		void calculateCoef();
	
	public:
		static double OneIntegral(size_t i, double alpha_);				//	returns the value of integ[x^i exp(-a x^2)]dx
	//	returns the value of <phi1| x^i y^j z^k |phi2>
		static double OneIntegral(  const OneGaussian& phi1, const OneGaussian& phi2,
									const size_t pow_x, const size_t pow_y, size_t pow_z);
	};

//-----------------------------------------------------------------------//
//																		 //
//				class for ContractedGaussian Functions					 //
//																		 //
//-----------------------------------------------------------------------//
	class Gaussian{
	public:
		//---	setting functions for Non-contracted Gaussian basis	---
		void set(double exponent, size_t xi_, size_t yi_, size_t zi_=0);	
		void set(double exponent, size_t xi_=0);
		//---	setting functions for Contracted Gaussian basis	---
		void set(const std::vector<OneGaussian>& gaussian_set_, const std::vector<double>& ratio_);

		double operator * (const OneGaussian& obj)const;
		double operator * (const Gaussian& obj)const;

		double get(double x, double y, double z)const;	//	3 dimensional
		double get(double x, double y)const;			//	2 dimensional
		double get(double x)const;						//	1 dimensional
		double getArea()const{return area;};
		double getMoment_x()const{return moment[0];};
		double getMoment_y()const{return moment[1];};
		double getMoment_z()const{return moment[2];};
		
		const OneGaussian& operator[] (size_t i)const{return gaussian_set[i];}
		double getCoef(size_t i)const{return ratio[i];}
		size_t size()const {return ratio.size();}

	private:
		std::vector<OneGaussian> gaussian_set;
		std::vector<double>   ratio;
		double area;
		double moment[3];
		void caluclateCoef();
		//---	protected functions	---		
		size_t l();size_t m(); size_t n();	double a(); double c();

	public:
	//	returns the value of <phi1| x^i y^j z^k |phi2>
		static double OneIntegral(  const Gaussian& phi1, const Gaussian& phi2,
									const size_t pow_x, const size_t pow_y, size_t pow_z);		//	returns the value of <phi1|phi2>
	};

//-----------------------------------------------------------------------//
//																		 //
//					class for Overlap Matricies							 //
//																		 //
//-----------------------------------------------------------------------//
	//---	overlap matrix for gaussian basis set	---
	class Overlap{
	public:
		//---	constructors	---
		Overlap(void);
		Overlap(const Overlap& obj);
		Overlap& operator= (const Overlap& obj);
		void Copy(const Overlap& obj);
		~Overlap(void);

		//---	setting functions	---
		void set(std::vector<Gaussian>& Phi_);
		void resize(size_t size_);

		//---	core functions for calculating the overlap	---
		void calculate();

		//---	functions for accessing to the contained data	---
		size_t size() const {return Phi.size();}
		const std::vector<Gaussian>& getPhi()const {return Phi;}
		const std::vector<Gaussian>& getBasis()const {return Phi;}
		const Gaussian& operator [] (size_t i)const{return Phi[i];}

		gsl_matrix* operator()() const{return delta;}			//	return <phi|phi>
		gsl_matrix* getdelta() const{return delta;}				//	return <phi|phi>
		gsl_matrix* x()const{return delta_x;}					//	return <phi|x|phi>
		gsl_matrix* y()const{return delta_y;}					//	return <phi|x|phi>
		gsl_matrix* z()const{return delta_z;}					//	return <phi|x|phi>
		gsl_matrix* inv()const{return delta_inv;}				//	return <phi|phi>^-1
		gsl_matrix* invsqrt()const{return delta_invsqrt;}		//	return <phi|phi>^-0.5

		Overlap getTransformedBasisSet(const gsl_matrix* transform_matrix)const;
		void BasisTransformation(const gsl_matrix* transform_matrix);

	private:
		std::vector<Gaussian> Phi;		//	|phi>
		gsl_matrix* delta;							//	<phi|phi>
		gsl_matrix* delta_x;						//	<phi|x|phi>
		gsl_matrix* delta_y;						//	<phi|y|phi>
		gsl_matrix* delta_z;						//	<phi|z|phi>
		gsl_matrix* delta_inv;						//	inverse matrix of the delta
		gsl_matrix* delta_invsqrt;					//	inverse sqrt matrix of the delta
		
		//	for memory control
		bool allocated;
		void MemoryControl(size_t size_);
		void Allocate(size_t size_);
		void Free();
	};
};

#include "gaussian_basis.inl"

#endif
