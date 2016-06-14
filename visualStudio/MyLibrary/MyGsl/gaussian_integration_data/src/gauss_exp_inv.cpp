#include "gauss_exp_inv.h"


double gauss_exp_inv::get0(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = exp(-v[i]*v[i] - x/v[i]);
	}
	interp.set(v, f, mygsl::interp::Cspline);
	return interp.getinteg(0.0, vmax) + 0.5*gsl_sf_erfc(vmax)*sqrt(M_PI);
}
double gauss_exp_inv::get1_2(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = exp(-x/v[i])*0.5*(0.5*sqrt(M_PI)*gsl_sf_erfc(v[i]) + v[i]*exp(-v[i]*v[i]));
	}
	interp.set(v, f, mygsl::interp::Cspline);
	double val = interp.getinteg(0.0, vmax);
	val += -sqrt(M_PI)*0.25*vmax*gsl_sf_erfc(vmax)
		+ 0.5*(0.5*sqrt(M_PI)+1.0)*(0.5*exp(-vmax*vmax)+sqrt(M_PI)*0.25*gsl_sf_erfc(vmax));
	return val;
}

double gauss_exp_inv::get1(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = v[i]*exp(-v[i]*v[i] - x/v[i]);
	}
	interp.set(v, f, mygsl::interp::Cspline);
	return interp.getinteg(0.0, vmax) + 0.5*exp(-vmax*vmax) + 0.25*gsl_sf_erfc(vmax)*sqrt(M_PI);
}

double gauss_exp_inv::get2(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = v[i]*v[i]*exp(-v[i]*v[i] - x/v[i]);
	}
	interp.set(v, f, mygsl::interp::Cspline);
	return interp.getinteg(0.0, vmax) + 0.5*vmax*exp(-vmax*vmax) + 0.25*exp(-vmax*vmax) + 0.125*gsl_sf_erfc(vmax)*sqrt(M_PI);
}
double gauss_exp_inv::get3(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = v[i]*v[i]*v[i]*exp(-v[i]*v[i] - x/v[i]);
	}
	interp.set(v, f, mygsl::interp::Cspline);
	return interp.getinteg(0.0, vmax) + 0.5*vmax*vmax*exp(-vmax*vmax) + 0.25*vmax*exp(-vmax*vmax) + 0.125*exp(-vmax*vmax) + 0.0625* gsl_sf_erfc(vmax)*sqrt(M_PI);
}

double gauss_exp_inv::get_1(double x_){
	x = x_;
	vmax = -x/log(1.0-threashold);
	vmin = vmax*vmin_factor;
	v.resize(vnum);
	f.resize(vnum);

	v[0] = 0.0;
	f[0] = 0.0;
	for(size_t i=1; i<vnum; ++i){
		v[i] = vmin*exp(log(vmax/vmin)/(vnum-2)*(i-1));
		f[i] = 1.0/v[i]*exp(-v[i]*v[i] - x/v[i]);
	}
	interp.set(v, f, mygsl::interp::Cspline);
	return interp.getinteg(0.0, vmax);
}


double gauss_exp_inv::get0_gsl(double x_){
	return getInteg(x, gauss_exp_inv::func0);
}
double gauss_exp_inv::get1_gsl(double x_){
	return getInteg(x, gauss_exp_inv::func1);
}
double gauss_exp_inv::get2_gsl(double x_){
	return getInteg(x, gauss_exp_inv::func2);
}
double gauss_exp_inv::get3_gsl(double x_){
	return getInteg(x, gauss_exp_inv::func3);
}


double gauss_exp_inv::getInteg(double x, double func(double, void*)){
	double result, error;
	gsl_function F;
	F.function = func;
	F.params = &x;
	gsl_integration_workspace * w =gsl_integration_workspace_alloc(2000);
	gsl_integration_qagiu(&F, 0.0, 0,  1e-4,  2000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

double gauss_exp_inv::func0(double v, void *params){
	if(v==0.0) return 0.0;
	double x = *(double *)params;
	return exp(-v*v - x/v);
}
double gauss_exp_inv::func1(double v, void *params){
	if(v==0.0) return 0.0;
	double x = *(double *)params;
	return v*exp(-v*v - x/v);
}
double gauss_exp_inv::func2(double v, void *params){
	if(v==0.0) return 0.0;
	double x = *(double *)params;
	return v*v*exp(-v*v - x/v);
}
double gauss_exp_inv::func3(double v, void *params){
	if(v==0.0) return 0.0;
	double x = *(double *)params;
	return v*v*v*exp(-v*v - x/v);
}
