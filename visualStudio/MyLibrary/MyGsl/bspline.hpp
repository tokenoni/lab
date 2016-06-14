#ifndef __MYGSL_BSPLINE_HPP__
#define __MYGSL_BSPLINE_HPP__

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <MyGsl/Multimin.hpp>
#include "interp.hpp"

//for debug
#include <IGOR\IGORitx.hpp>
namespace mygsl{
//-----------------------------------------------------//
//													   //
//		class for basis functions interpolation		   //
//													   //
//-----------------------------------------------------//
	class bspline{
	public:
		enum KnotsGettingMethod{
			Uniform,		//	uniformally distributed knots
			Greville,		//	determined by grebille method
			Automatic		//	automatically determined by data
		};

/*
//-------------------------------------
//			under construction 
//-------------------------------------
*/
		enum HoldingMethod{
			NoHoldings,
			Flat_at_FirstEdge,
			Zero_at_FirstEdge,
			ZeroAndFlat_at_FirstEdge,
			Flat_at_LastEdge,
			Zero_at_LastEdge,
			ZeroAndFlat_at_LastEdge,
			Flat_at_BothEdge,
			Zero_at_BothEdge,
			ZeroAndFlat_at_BothEdge,
			ZeroAndFlat_at_FirstEdge_Zero_at_LastEdge,
			ZeroAndFlat_at_FirstEdge_Flat_at_LastEdge,
			ZeroAndFlat_at_LastEdge_Zero_at_FirstEdge,
			ZeroAndFlat_at_LastEdge_Flat_at_FirstEdge
		};


		bspline(void);
		bspline(const size_t size_, const size_t nbreak_, const size_t k_ = 4);
		bspline(const bspline& obj);
		bspline& operator = (const bspline& obj);

		double get(const double x)const;
		double operator() (const double x)const{return get(x);}
		double getderiv(const double x)const;
		//--- setting function	---
		template < class Tyv_ >
		void set(const Tyv_& x, const Tyv_& y, const size_t nbreak_, const size_t k_ = 4, KnotsGettingMethod knots_method_ = Automatic);
		
		//---	core calculation function for bspline interpolation
		void setInterp(KnotsGettingMethod knots_method_ = Automatic, HoldingMethod holdingMethod_ = NoHoldings);
		void setInterp(std::vector<double> breakpoints, HoldingMethod holdingMethod_ = NoHoldings);
		double setBspline(HoldingMethod holdingMethod_ = NoHoldings);

		//---	for directly setting the contained data	---
		void setSize(const size_t size_, const size_t nbreak_, const size_t k_ = 4)	{	MemoryControl(size_, nbreak_, k_);}
		void set(const size_t i, const double xi, const double yi, const double sigma_=1.0);
		void set_x(const size_t i, const double xi)									{	pos[i] = xi;}
		void set_y(const size_t i, const double yi)									{	data[i] = yi;}
		void set_y_zero();

		//---	tip functions	---
		void clear()					{	FreeMemory();	FreeMemoryWorkspace();}
		size_t size()const				{	return num;}
		size_t order()const				{	return k;}
		bool isAllocated()const			{	return allocated;}
		bool isAllocatedWorkspace()const{	return allocated_workspace;}
		double getMin_x()const{if(allocated)return pos[0];else return 0.0;};
		double getMax_x()const{if(allocated)return pos[size()-1]; else return 0.0;};
		std::vector<double> getOriginal_x()const;
		std::vector<double> getOriginal_y()const;

	private:

		gsl_bspline_workspace *bw;
		gsl_bspline_deriv_workspace *dbw;
		gsl_vector *Bk;		//	spline coefficient (contracted)
		gsl_matrix *dBk;
		gsl_vector *c;		//	coefficient	for Bspline. It is derived by linear multifitting
		size_t nbreak;		//	nubmer of break point 
		size_t k;			//	order of the bspline functions. default is k=4 (3rd order)
		
		bool allocated;
		bool allocated_workspace;
		double *pos;
		double *data;
		double *w;		//	weight (defined as 1/sigma)
		size_t num;

		KnotsGettingMethod knots_method;
		HoldingMethod holdingMethod;

		void sort();
		void getKnotsUniform();
		void getKnotsGreville();
		double setBspline_NoHoldings();
		double setBspline_Flat_at_FirstEdge();
		double setBspline_Zero_at_FirstEdge();
		double setBspline_ZeroAndFlat_at_FirstEdge();
		double setBspline_Flat_at_LastEdge();
		double setBspline_Zero_at_LastEdge();
		double setBspline_ZeroAndFlat_at_LastEdge();
		double setBspline_Flat_at_BothEdge();
		double setBspline_Zero_at_BothEdge();
		double setBspline_ZeroAndFlat_at_BothEdge();
		double setBspline_ZeroAndFlat_at_FirstEdge_Zero_at_LastEdge();
		double setBspline_ZeroAndFlat_at_FirstEdge_Flat_at_LastEdge();
		double setBspline_ZeroAndFlat_at_LastEdge_Zero_at_FirstEdge();
		double setBspline_ZeroAndFlat_at_LastEdge_Flat_at_FirstEdge();

		double setBsplineAutomatic();
		size_t coef_size()const{return nbreak - 2 + k;}

		//---	memory controling functions		---
		void MemoryControl(const size_t size_, const size_t nbreak_, const size_t k_ = 4);
		void Allocate(const size_t size_);
		void FreeMemory();
		void AllocateWorkspace(const size_t nbreak_, const size_t k_ = 4);
		void FreeMemoryWorkspace();
		void copyData(const bspline& obj);
		void copyWorkspace(const bspline& obj);
	};


//---	class only for used in the automatically knots dentemining in Bspline class	---
	class BsplineAutomaticKnotsDetermine{
	public:
		BsplineAutomaticKnotsDetermine
			(const size_t data_size_, const size_t coef_size_, const size_t break_size_, 
			const double *x, const double *y, const double *w_, 
			gsl_vector *Bk, gsl_vector *c_, gsl_bspline_workspace *bw);
		~BsplineAutomaticKnotsDetermine(void);
		double operator () (std::vector<double>& coef);

	private:
		const size_t data_size, coef_size, nbreak;
		const double *x, *y, *w;
		gsl_vector *Bk, *c;
		gsl_bspline_workspace *bw;
		BsplineAutomaticKnotsDetermine(void);
		gsl_matrix *X;
		gsl_matrix *cov;
		gsl_multifit_linear_workspace *mw;
		gsl_vector *breakpts;
		gsl_permutation *p;
	};
};

#include "bspline.inl"
#endif