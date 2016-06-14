#include "gamma_erf.hpp"
#include <IGOR/IGORitx.hpp>
#define NUMN 100
#include <complex>
#include <iostream>
#include <ccgsl/matrix.hpp>
#include <time.h>
#include <iostream>
#include "fmath.hpp"

int main(){
	size_t num_calc = 100000000;
	size_t store_ratio = 100000;
	size_t count;
	int calc_t;

	double* rv = new double [num_calc];
	double* rslt_original = new double [num_calc/store_ratio];
	double* rslt_new= new double [num_calc/store_ratio];

	double ln2 = log(2.0);
	for(size_t i=0; i<num_calc; ++i){
		rv[i] = -10.0+20.0/num_calc*i;
	}
	

	count =0;
	calc_t = clock();
	for(size_t i=0; i<num_calc; ++i){
		double rslt = fmath::expd(rv[i]);
		if(i%store_ratio == 0){
			rslt_new[count] = rslt;
			++ count;
		}
	}
	std::cout << "time elapsed " <<  clock() - calc_t << "ms"<< std::endl;

	count =0;
	calc_t = clock();
	for(size_t i=0; i<num_calc; ++i){
		double rslt = exp(rv[i]);
		if(i%store_ratio == 0){
			rslt_original[count] = rslt;
			++ count;
		}
	}
	std::cout << "time elapsed " <<  clock() - calc_t << "ms"<< std::endl;
	getchar();

	IGORdata::write_itx(rslt_original, num_calc/store_ratio, "check/exp_original.itx", "exp_original");
	IGORdata::write_itx(rslt_new, num_calc/store_ratio, "check/exp_new.itx", "exp_new");
	delete [] rv;
	delete [] rslt_original;
	delete [] rslt_new;


/*	gsl::vector z(NUMN);
	gsl::vector z_img(NUMN+1);
	gsl::vector zval(NUMN);
	
	for(size_t i=0; i<NUMN; i++){
		z[i] = -5.0+5.0/(0.5*NUMN-1.0)*i;
		z_img[i] = -5.0+5.0/(0.5*NUMN-1.0)*(i-0.5);
	}
	z_img[NUMN] = -5.0+5.0/(0.5*NUMN-1.0)*(NUMN-0.5);
	IGORdata::write_itx(z,     "check/zz.itx",    "zz");
	IGORdata::write_itx(z_img, "check/zz_img.itx",    "zz_img");
	

	for(size_t i=0; i<NUMN; i++) zval[i] = MyFunc::gamma(z[i]);
	IGORdata::write_itx(z, zval, "check/gamma_test.itx",    "zz", "gamma_test");

	for(size_t i=0; i<NUMN; i++) zval[i] = MyFunc::ln_gamma(z[i]);
	IGORdata::write_itx(z, zval, "check/ln_gamma_test.itx",    "zz", "ln_gamma_test");

	for(size_t i=0; i<NUMN; i++) zval[i] = MyFunc::erf(z[i]);
	IGORdata::write_itx(z, zval, "check/erf_test.itx",    "zz", "erf_test");

	for(size_t i=0; i<NUMN; i++) zval[i] = MyFunc::erfc(z[i]);
	IGORdata::write_itx(z, zval, "check/erfc_test.itx",    "zz", "erfc_test");

	gsl::matrix zr(NUMN, NUMN);
	gsl::matrix zi(NUMN, NUMN);

	for(size_t i=0; i<NUMN; i++) {
		for(size_t j=0; j<NUMN; j++){
			std::complex<double> arg(z[i], z[j]);
			std::complex<double> val = MyFunc::gamma(arg);
			zr[i][j] = val.real();
			zi[i][j] = val.imag();
		}
	}
	IGORdata::write_itx(zr, "check/gammac_r_test.itx", "gammac_r_test");
	IGORdata::write_itx(zi, "check/gammac_i_test.itx", "gammac_i_test");

	for(size_t i=0; i<NUMN; i++) {
		for(size_t j=0; j<NUMN; j++){
			std::complex<double> arg(z[i], z[j]);
			std::complex<double> val = MyFunc::erf(arg);
			zr[i][j] = val.real();
			zi[i][j] = val.imag();
		}
	}
	IGORdata::write_itx(zr, "check/erf_r_test.itx", "erf_r_test");
	IGORdata::write_itx(zi, "check/erf_i_test.itx", "erf_i_test");

	for(size_t i=0; i<NUMN; i++) {
		for(size_t j=0; j<NUMN; j++){
			std::complex<double> arg(z[i], z[j]);
			std::complex<double> val = MyFunc::erfc(arg);
			zr[i][j] = val.real();
			zi[i][j] = val.imag();
		}
	}
	IGORdata::write_itx(zr, "check/erfc_r_test.itx", "erfc_r_test");
	IGORdata::write_itx(zi, "check/erfc_i_test.itx", "erfc_i_test");

	
	double complex zc;
	for(size_t i=0; i<NUMN; i++){
		for(size_t j=0; j<NUMN; j++){
			u[i][j] = zval.real();
			v[i][j] = zval.imag();
		}
	}
	IGORdata::write_itx(z, "cerf_check/zz.itx", "zz");
	IGORdata::write_itx(z_img, "cerf_check/zimg.itx", "zimg");
	IGORdata::write_itx(u, "cerf_check/uu.itx", "uu");
	IGORdata::write_itx(v, "cerf_check/vv.itx", "vv");
*/}
