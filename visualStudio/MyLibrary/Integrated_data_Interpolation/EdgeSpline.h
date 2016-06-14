#pragma once
#ifndef __EDGESPLINE_H__
#define __EDGESPLINE_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <iostream>

class EdgeSpline
{
public:
	EdgeSpline(void);
	~EdgeSpline(void);

	//---	constraint method at x ---
	enum Scale{
		Linear,
		Log,
		LinearArea,
		LogArea
	};
	Scale scale;
	enum EndPoint{
		Free,	//	default is Free (2nd derivetive is same to the next point) 
		Natural
	};
	EndPoint endpoint;
//---	various type of set function	---
	template < class TyE_ >
	void set(const TyE_& y_src, size_t size_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, double dx_, size_t size_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, size_t size_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, double dx_, size_t size_);
	template < class TyE_ >
	void set(const TyE_& y_src, size_t size_, EndPoint endpoint_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, double dx_, size_t size_, EndPoint endpoint_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, size_t size_, EndPoint endpoint_);
	template < class TyE_ >
	void set(const TyE_& x_src, const TyE_& y_src, const TyE_& edge_src, double dx_, size_t size_, EndPoint endpoint_);

//---	get the interpolated value at x_interp	---
	double get(double x_interp);

//---	valiables contains the coefficients of the interpolated function	---	
	double* x;	double* y;
	double* dx;			//	length to calculate the area
	double* x_edge;		//	point the function would be continus
	size_t size;
	double* a;	double* b;	double* c;	//	variables for spline interpolation
	void sort();
//---	variables for solving the systems of the equation
	gsl_matrix* X;	
	gsl_vector* v;
	gsl_vector* coef;
	gsl_permutation* p;
	virtual void solve();
	virtual void setScale();
	inline double get_a(const gsl_vector* coef_, size_t i){ return coef_->data[  3*i  * coef_->stride];} 
	inline double get_b(const gsl_vector* coef_, size_t i){ return coef_->data[(3*i+1)* coef_->stride];} 
	inline double get_c(const gsl_vector* coef_, size_t i){ return coef_->data[(3*i+2)* coef_->stride];} 

//---	variables for solving the nonlinear system of the equation	---
//---	variables and function for memory control	---
private:
	bool Allocated;
//	size_t size;
	void MemoryControl(size_t size_);
public:
	virtual void Allocation();
	virtual void FreeMemory();
};

class EdgeSplineArea: public EdgeSpline{
public:
	void setScale(){ scale = LinearArea;}
};

class EdgeLogSpline: public EdgeSpline{
public:
	void setScale(){ scale = Log;}
};

class EdgeLogSplineArea: public EdgeSpline{
private:
	gsl_vector* f;
	static const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_function FDF;
	gsl_multiroot_fsolver* s;
	void setScale(){ scale = LogArea;}
	void solve();

public:	
	void f_multiroot( const gsl_vector *coef_, gsl_vector* f);

//---	function for memory control	---
	void Allocation();
	void FreeMemory();
private:
};

class EdgeLogSpline_Multiroot{
public:
	static int f(const gsl_vector* coef, void* p, gsl_vector* f){
		((EdgeLogSplineArea *)p)->f_multiroot(coef, f);
		return GSL_SUCCESS;
	}
private:
};

#include "EdgeSpline.inl"

#endif