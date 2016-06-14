#ifndef _CERF_hPP_
#define _CERF_hPP_

#include <math.h>
#include <complex>
#include <gsl/gsl_sf_hyperg.h>

namespace gsl_plus{
	#define PI 3.14159265358979323846264338327950288 

	int round_off(double x){ return int(x+0.5);}


	bool cerf(double xi, double yi, double& u, double& v){
/*		const double RMAXREAL = 1.0/floor(0.1);
		const double RMAXEXP  = log(RMAXREAL) - log(2.0);
		const double RMAXGONI = 0.0;
*/
		double factor = 2.0/sqrt(PI);	//	factor for the (EQuATIOn (7.1.5), P.297)
		
		double xabs = fabs(xi);
		double yabs = fabs(yi);
		double x = xabs/6.3;
		double y = yabs/4.4;
		double qrho = xabs * xabs + yabs * yabs;
		
		if(!(qrho>=0.0)) return false;
		double xabsq = xabs * xabs;
		double xquad = xabsq - yabs*yabs;
		double yquad = 2.0*xabs*yabs;

		double u2, v2;

	//--	if qrho < 0.085264														--//
	//--	the FADDEEVA-FUNCTION is evaluated using a power series					--//
	//--				(ABRAMOWITZ/STEGun, EQuATIOn (7.1.5), P.297)				--//
	//--	n is the minimum number of terms needed to obtain the required accuracy	--//
		if(qrho < 0.085264){	
			qrho = (1-0.85*y) * sqrt(qrho);
			size_t n = (size_t)round_off(6.0+72.0*qrho);
			size_t J = 2*n+1;
			double xsum = 1.0/J;
			double ysum = 0.0;
			double xaux = 0.0;
			for(size_t i=n; i>=1; i--){
				J = J-2;
				xaux = (xsum*xquad - ysum*yquad)/i;
				ysum = (xsum*yquad + ysum*xquad)/i;
				xsum = xaux + 1.0/J;
			}
			double u1 = -factor*(xsum*yabs + ysum*xabs)+1.0;
			double v1 =  factor*(xsum*xabs - ysum*yabs);
			double daux = exp(-xquad);

			u2 =  daux * cos(yquad);
			v2 = -daux * sin(yquad);

			u = u1*u2 - v1*v2;
			v = u1*v2 + v1*u2;
		}
	//--	if qrho > 1.0
	//--	the W(Z) is evaluated using the LAPLACE CONTINUED Fraction	--//
	//--	nu is the minimum number of terms needed to obtain the required accuracy	--//

	//--	if 0.085264 < qrho < 1.0
	//--	then W(z) is evaluated by a TRunLATED TAYLOR EXPANSIhON, 
	//--	where the laplace continued fraction is used to calculate the derivatives of W(z)					--//
	//--	kapn is the minimum number of terms in the TAYLOR expansion needed to obtain the required accuracy	--//
	//--	nu is the minimum number of terms of the continued fraction needed to calculate the derivatives		--//
		else{
			double h = 0.0, h2 = 0.0;
			size_t kapn = 0;
			size_t nu=0;
			if(qrho > 1.0) {
				qrho = sqrt(qrho);
				nu =   (size_t)round_off(3.0+(1442.0/(26.0*qrho+77.0)));
			}else{
				qrho = (1.0-y)*sqrt(1.0-qrho);
				h = 1.88*qrho;
				h2 = 2.0*h;
				kapn = (size_t)round_off( 7.0+34.0*qrho);
				nu =   (size_t)round_off(16.0+26.0*qrho);
			}
			
			double qlambda=0.0;
			if(h>0.0) qlambda = pow(h2, 1.0*kapn);

			double rx = 0.0;
			double ry = 0.0;
			double sx = 0.0;
			double sy = 0.0;

			for(int n = nu; n>=0; n--){
				double np1 = 1.0*n+1.0;
				double tx = yabs + h + np1*rx;
				double ty = xabs - np1*ry;
				double c  = 0.5/(tx*tx + ty*ty);
				rx = c * tx;
				ry = c * ty;
				if(h>0.0 && n <= (int)kapn){
					tx = qlambda + sx;
					sx = rx*tx - ry*sy;
					sy = ry*tx + rx*sy;
					qlambda = qlambda/h2;
				}
			}
//---	correction by fujii	for the case of (qrho > 1.0)---
/*			if(h==0.0){
				rx = 0.0;
				ry = 0.0;
				for(int n=nu; n>=2; n--){
					double n_double = 1.0*n;
					double sign = 1.0;
					if(n%2==1) sign =--1.0;
					rx = n_double-0.5 - sign*0.5*n_double*xquad;
					ry = n_double-0.5 - sign*0.5*n_double*yquad;
*/
			if(h==0.0){
				u = factor * rx;
				v = factor * ry;
			}else{
				u = factor * sx;
				v = factor * sy;
			}
			if(yabs == 0.0)	u = exp(-xabs*xabs);
		}

	//--	evaluation of w(z) in other quadrant	--
		if(yi < 0.0){
			if(qrho < 0.085264){
				u2 = 2.0 * u2;
				v2 = 2.0 * v2;
			}else{
				xquad = -xquad;
				double w1 = 2.0*exp(xquad);
				u2 =  w1 * cos(yquad);
				v2 = -w1 * sin(yquad);
			}
			u = u2 - u;
			v = v2 - v;
			if(xi > 0.0)        v = -v;
			else if ( xi < 0.0) v = -v;
		}

		return true;
	}

};

#endif