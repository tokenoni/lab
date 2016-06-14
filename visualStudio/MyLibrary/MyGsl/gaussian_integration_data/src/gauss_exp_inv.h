#pragma once
#include <vector>
#include <MyGsl\interp.hpp>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

static const size_t vnum = 1000;
static const double threashold = 0.0001;
static const double vmin_factor= 0.000001;

class gauss_exp_inv
{
public:
	gauss_exp_inv(void){};
	~gauss_exp_inv(void){};
	double get0(double x);
	double get1_2(double x);
	double get1(double x);
	double get2(double x);
	double get3(double x);
	double get_1(double x);

	std::vector<double> v;
	std::vector<double> f;
	double vmax, vmin;
	double x;
	mygsl::interp interp;

	double get0_gsl(double x);
	double get1_2_gsl(double x);
	double get1_gsl(double x);
	double get2_gsl(double x);
	double get3_gsl(double x);
	double getInteg(double x, double func(double, void*));
	static double func0(double x, void *params);
	static double func1(double x, void *params);
	static double func2(double x, void *params);
	static double func3(double x, void *params);
};
