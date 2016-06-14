#include "gauss_exp_inv.h"
#include "gauss_exp_inv_complex.h"
#include <IGOR\IGORitx.hpp>
void calc_real();
void calc_complex();

int main(){
	calc_real();
//	calc_complex();
}

void calc_complex(){
	gauss_exp_inv_complex test;
	test.get0(mygsl::complex(0.1, 0.0));
	IGORdata::write_itx(test.interp.get_original_x(), "vel_00.itx", "v00");
	IGORdata::write_itx(test.interp.get_original_y(), "func_00.itx", "f00");

	test.get0(mygsl::complex(0.1, 0.1));
	IGORdata::write_itx(test.interp.get_original_x(),  "vel_01.itx", "v01");
	IGORdata::write_itx(test.interp.get_original_y(), "func_01.itx", "f01");

	test.get0(mygsl::complex(0.1, 0.2));
	IGORdata::write_itx(test.interp.get_original_x(),  "vel_02.itx", "v02");
	IGORdata::write_itx(test.interp.get_original_y(), "func_02.itx", "f02");

	//----	calculation	---
	const size_t num = 200;
	const size_t num1 = 201;
	const double xmax = 100.0;
	const double xmin = 0.001;
	std::vector<double> x(num);
	std::vector<double> y(num1);
	std::vector<std::vector<double>> F_r(num, num1);
	std::vector<std::vector<double>> F_i(num, num1);

	IGORdata::write_itx(x, "xx.itx", "xx");
	//--------	for F0	---------
	x[0] = 0.0;
	y[0] = 0.0;
	
	F_r[0][0] = 0.5*sqrt(M_PI);
	F_i[0][0] =  0.0;

	for(size_t i=1; i<num; ++i)
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
	for(size_t i=1; i<num1; ++i)
		y[i] = xmin*exp(log(xmax/xmin)/(num1-2)*(i-1));
	
	gauss_exp_inv_complex calc_v;
	for(size_t i=1; i<num; ++i){
		for(size_t j=1; j<num1; ++j){
			mygsl::complex tmp = calc_v.get0(mygsl::complex(x[i], y[j]));
			F_r[i][j] = tmp.real();
			F_i[i][j] = tmp.imag();
		}
	}

	IGORdata::write_itx(  F_r, "exp_inv_F0_z_real.itx", "F0_z_real");
	IGORdata::write_itx(  F_i, "exp_inv_F0_z_imag.itx", "F0_z_imag");
	IGORdata::write_edgeVector( x, "z_real.itx", "z_real");
	IGORdata::write_edgeVector( y, "z_imag.itx", "z_imag");

}

void calc_real(){
	gauss_exp_inv test;
	test.get1_2(0.001);
	IGORdata::write_itx(test.v, test.f, "func_1.itx", "v1", "f1");

	test.get1_2(10.0);
	IGORdata::write_itx(test.v, test.f, "func_10.itx", "v10", "f10");

	test.get1_2(100.0);
	IGORdata::write_itx(test.v, test.f, "func_100.itx", "v100", "f100");


	//----	calculation	---
	const size_t num = 1000;
	const double xmax = 100.0;
	const double xmin = 0.001;
	std::vector<double> x(num);
	std::vector<double> F(num);

	//--------	for F0	---------
	x[0] = 0.0;
	F[0] = 0.5*sqrt(M_PI);

	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get0(x[i]);
//		F_gsl[i] = calc_v.get0_gsl(x[i]);
	}

	IGORdata::write_itx(x, "xx.itx", "xx");
	IGORdata::write_itx( F, "exp_inv_F0_x.itx", "F0_x");

	//--------	for F1_2	---------
	x[0] = 0.0;
	F[0] = 0.5;

	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get1_2(x[i]);
//		F_gsl[i] = calc_v.get0_gsl(x[i]);
	}

	IGORdata::write_itx(x, "xx.itx", "xx");
	IGORdata::write_itx( F, "exp_inv_F1_2_x.itx", "F1_2_x");

	//--------	for F1	---------
	x[0] = 0.0;
	F[0] = 0.5;
//	F_gsl[0] = 0.5;
	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get1(x[i]);
//		F[i] = calc_v.get1(x[i]);
	}

	IGORdata::write_itx( F, "exp_inv_F1_x.itx", "F1_x");

	//--------	for F2	---------
	x[0] = 0.0;
	F[0] = sqrt(M_PI)/4.0;
	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get2(x[i]);
	}

	IGORdata::write_itx( F, "exp_inv_F2_x.itx", "F2_x");

	//--------	for F3	---------
	x[0] = 0.0;
	F[0] = 0.5;
	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get3(x[i]);
	}

	IGORdata::write_itx( F, "exp_inv_F3_x.itx", "F3_x");

	//--------	for F_1	---------
	x[0] = 0.0;
	F[0] = 0.5;
	for(size_t i=1; i<num; ++i){
		x[i] = xmin*exp(log(xmax/xmin)/(num-2)*(i-1));
		gauss_exp_inv calc_v;
		F[i] = calc_v.get_1(x[i]);
	}

	IGORdata::write_itx( F, "exp_inv_F_1_x.itx", "F_1_x");
}