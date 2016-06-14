#include "cerf.hpp"
#include <IGOR/IGORitx.hpp>

#include <ccgsl/matrix.hpp>
#define NUMN 100


int main(){
	gsl::vector z(NUMN);
	gsl::vector z_img(NUMN+1);
	gsl::matrix u(NUMN,NUMN);
	gsl::matrix v(NUMN,NUMN);
	
	for(size_t i=0; i<NUMN; i++){
		z[i] = -2.5+2.5/(0.5*NUMN-1.0)*i;
//		z[i] = 2.5/(NUMN-1.0)*i;
		z_img[i] = -2.5+2.5/(0.5*NUMN-1.0)*(i-0.5);
//		z_img[i] = 2.5/(NUMN-1.0)*(i-0.5);
	}
	z_img[NUMN] = -2.5+2.5/(0.5*NUMN-1.0)*(NUMN-0.5);
//	z_img[NUMN] = 2.5/(NUMN-1.0)*(NUMN-0.5);

//	for(size_t i=40; i<NUMN; i++)
//		gsl_plus::cerf(z[i], 0.0, u[i][0], v[i][0]);
		
	for(size_t i=0; i<NUMN; i++){
		for(size_t j=0; j<NUMN; j++){
		gsl_plus::cerf(z[i], z[j], u[j][i], v[j][i]);
		}
	}
	IGORdata::write_itx(z, "cerf_check/zz.itx", "zz");
	IGORdata::write_itx(z_img, "cerf_check/zimg.itx", "zimg");
	IGORdata::write_itx(u, "cerf_check/uu.itx", "uu");
	IGORdata::write_itx(v, "cerf_check/vv.itx", "vv");
}
