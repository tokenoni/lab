#ifndef __IntegratedInterpolation_H__
#define __IntegratedInterpolation_H__
#include <vector>
#include <gsl/gsl_spline.h>
#include <math.h>

class IntegSpline{
public:
	IntegSpline(void);
	~IntegSpline(void);

	enum Scale{
		Log,
		Linear
	};
	Scale scale;

	//---	set the values to be interpolated ---
	template< class TyI_ >
	void set(TyI_ &y_src, double dx_, size_t size_, Scale scale_ );
	template< class TyI_ >
	void set(const TyI_ &y_src, const TyI_ & x_src, size_t size_, Scale scale_ );

	//---	get the value with arbitrary x_interp	---
	double get(const double x_interp);
	//---	get the value with arbitrary x_interp	---
	double getIntegrated(const double x_interp);
private:
	void refreshSplineCoefficient(const double x_interp);

private:
	double* x;
	double* y;
	double* dx;
	bool Allocated;
	//---	workspace for spline interpolation	---
	gsl_interp_accel *acc;
	gsl_spline  *spline;
	size_t size;
	//---	keep the index, and coefficient for integration	---
	size_t index;
	double a, b, c;
	//---	function for the memory treating	---
	void Allocation(size_t size_);

	void sort();
};

#include "IntegratedInterpolation.inl"
#endif
