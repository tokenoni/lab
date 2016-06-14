#ifndef __MYGSL_VOIGT_HPP__
#define __MYGSL_VOIGT_HPP__

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_erf.h>
#include <cmath>

namespace mygsl{
	static const double M_PI = 3.14159265358979323846264338327950288;

	static double voigt_calc(const double x,const double y);

	static double gauss_integ(const double x, const double x0, const double wg, const double dx){	//	integ[x-dx/2, x+dx/2]{gauss(x, x0, wg)}
		double we = wg / (2.0*sqrt(log(2.0)));
		double val = (0.5*gsl_sf_erf((x-x0+0.5*dx)/we) - 0.5*gsl_sf_erf((x-x0-0.5*dx)/we))/dx;
		if(val==val) return val;
		else 0.0;
	}

	static double gauss(const double x, const double x0, const double wg)
	{
		double we = wg / (2.0*sqrt(2.0*log(2.0)));
		double ag = 1.0/(sqrt(2.0*M_PI)*we);
		return ag * exp(-(x-x0)*(x-x0)/(2.0*we*we));
	}

	static double lorentz(const double x, const double x0, const double wl)
	{		
		double gamma = 0.5*wl;	// HWHM of Lorentz function
		return 1.0/M_PI*gamma/((x-x0)*(x-x0) + gamma*gamma);
	}

	static double voigt(const double x, const double x0, const double wge, const double wl)
	/*****		voigt function		*******
	*	x0 : center of spectrum
	*	wge : 1/e width of Gaussian width	(σ in usual meaning)
	*	wl : FWHM of Lorentzian width		(γ 
	*/
	{
		double w2=1.0/(2.0*sqrt(2.0)*wge);
		double w4=w2*wl;
	
		double xx=((x-x0)*w2);
		double yy=w4;

		return voigt_calc(xx,yy)*w2/sqrt(M_PI);
	}

	static double voigt_calc(const double x,const double y)
	{
		#define ADDc(a,b) (gsl_complex_add(a,b))
		#define SUBc(a,b) (gsl_complex_sub(a,b))
		#define MULc(a,b) (gsl_complex_mul(a,b))
		#define DIVc(a,b) (gsl_complex_div(a,b))
		#define ADDr(a,b) (gsl_complex_add_real(a,b))
		#define SUBr(a,b) (gsl_complex_sub_real(a,b))
		#define MULr(a,b) (gsl_complex_mul_real(a,b))
		#define DIVr(a,b) (gsl_complex_div_real(a,b))
		#define REAL(a)   (gsl_complex_rect(a,0.0))
		gsl_complex t = gsl_complex_rect(y, -x);
		gsl_complex u;
		gsl_complex w;
		double s = fabs(x) + y;
	
	//	Region I		W= T*0.5641896/(0.5+T*T)
		if(s >= 15.0)
		{	w = DIVc( MULr(t, 0.5641896), ADDr( MULc(t,t) , 0.5));}
	
	//	Region II		U= T*T
	// 					 W= T*(U*0.5641896+1.410474)/(U*(3+U)+0.75)
		else if( s>=5.5)
		{	u = MULc(t,t);
			w = DIVc(MULc(t, ADDr( MULr(u,0.5641896), 1.410474)), ADDr( MULc(u, ADDr(u, 3.0)), 0.75));
		}
	
	//	Region III		W  = (T*(T*(T*(T* 0.5642236+3.778987+11.96482)+20.20933)+16.4955)	
	//	                W /= (T*(T*(T*(T*(T+6.699398)+21.69274)+39.27121)+38.82363)+16.4955)
		else if( y >= (0.195*fabs(x)-0.176))
		{	w =          ADDr( MULc(t, ADDr( MULc(t, ADDr( MULc(t, ADDr( MULr(t,0.5642236) ,3.778987)),11.96482)), 20.20933)), 16.4955);
			w = DIVc( w, ADDr( MULc(t, ADDr( MULc(t, ADDr( MULc(t, ADDr( MULc(t,ADDr(t, 6.699398)),21.69274)),39.27121)),38.82363)), 16.4955));
		}
	//	Region IV		U= T*T
	//	                W= T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313-U*(35.76683-U*(1.320522-U*0.56419))))))
	//	                W /= ( 32066.6-U*(24322.84-U*(9022.228-U*(2186.181-U*(364.2191-U*(61.57037-U*(1.841439-U)))))))
	//					W= cmplx(exp(real(U))*cos(imag(U)),0)-W
		else
		{	u = MULc(t,t);
			w = MULc(t, SUBc(REAL(36183.31),MULc(u,SUBc(REAL(3321.9905),MULc(u,SUBc(REAL(1540.787),MULc(u,SUBc(REAL(219.0313),
											MULc(u,SUBc(REAL( 35.76683),MULc(u,SUBc(REAL(1.320522),MULr(u,0.56419)))))))))))));
			w = DIVc(w, SUBc(REAL( 32066.6),MULc(u,SUBc(REAL( 24322.84),MULc(u,SUBc(REAL(9022.228),MULc(u,SUBc(REAL(2186.181),
											MULc(u,SUBc(REAL( 364.2191),MULc(u,SUBc(REAL(61.57037),MULc(u,SUBc(REAL(1.841439),u))))))))))))));
			w = SUBc( REAL(exp(GSL_REAL(u))*cos(GSL_IMAG(u))), w);
		}
		return GSL_REAL(w);
	}

	/* IGOR Tech note でのVoigt関数の表記を拝借．
	 *	以下，そのソース．
	 *
	 *		Function Voigt(X,Y)
	 *			variable X,Y
	 *
	 *			variable/C W,U,T= cmplx(Y,-X)
	 *			variable S =abs(X)+Y
	 *
	 *			if( S >= 15 )								//        Region I
	 *				W= T*0.5641896/(0.5+T*T)
	 *			else
	 *				if( S >= 5.5 ) 							//        Region II
	 *					U= T*T
	 *					W= T*(1.410474+U*0.5641896)/(0.75+U*(3+U))
	 *				else
	 *					if( Y >= (0.195*ABS(X)-0.176) ) 	//        Region III
	 *						W= (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*0.5642236))))
	 *						W /= (16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))))
	 *					else									//        Region IV
	 *						U= T*T
	 *						W= T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313-U*(35.76683-U*(1.320522-U*0.56419))))))
	 *						W /= (32066.6-U*(24322.84-U*(9022.228-U*(2186.181-U*(364.2191-U*(61.57037-U*(1.841439-U)))))))
	 *						W= cmplx(exp(real(U))*cos(imag(U)),0)-W
	 *					endif
	 *				endif
	 *			endif
	 *			return real(W)
	 *		end
	 *
 
	 *		Accompanying files:
	 *			New Voigt.txt'	--	text file containing the Voigt function
	 *		_____________________________________________________________________________________________
	 *
	 *		The Voigt profile is used in spectroscopy when a given line shape is neither a Gaussian nor a Lorentzian but rather is a convolution of the two.  The accompanying file 鮮ew Voigt TEXT・contains an Igor user defined function that calculates this function.  Its relative accuracy is better than 0.0001 and most of the time is much better.
	 *
	 *		To use, position the insert point in the procedure window of an experiment and use 選nsert Text...・from the file menu to read in the function. To use the function in a user defined curve fit, you will need to wrap it in another function like so:
	 *
	 *		Function MyVoigt(w,x)
	 *			Wave w
	 *			Variable x
	 *
	 *			return w[0]+w[1]*voigt(w[2]*(x-w[3]),w[4])
	 *		End
	 *
	 *		Parameter w[0] sets the DC offset, w[1] sets the amplitude, w[2]  affects the width, w[3] sets the location of the peak and w[4] adjusts the shape.
	 *
	 *		The algorithm used is patterned after the source given by J. Humlicek. The algorithm is discussed by F. Schreier.
	 *
	 *		After the fit, you can use the returned coefficients to calculate the area (a) along with the half width at half max for the Gaussian (wg), Lorentzian (wl) and the Voigt (wv). Assuming the coefficient wave is named coef:
	 *		 
	 *			a= coef[1]*sqrt(pi)/coef[2]
	 *			wg= sqrt(ln(2))/coef[2]
	 *			wl= coef[4]/coef[2] 
	 *			wv= wl/2 + sqrt( wl^2/4 + wg^2)
	 *
	 *		You can calculate the amplitude of the peak via:
	 *			amp= coef[1]*exp(coef[4]^2)*erfc(coef[4])
	 *		This was derived by integrating Eq. (2) or Eq. (18) in Armstrong with x=0.
	 *
	 *		The calculation of the area is based on the following
	 *			voigt(x,y)=(y/pi) Integral[ dt*exp(-t^2)/[y^2+(x-t)^2]],
	 *		where the integration is from -infinity to infinity.
	 *
	 *		The integral of the voigt function given by Armstrong is:
	 *			Integal[voigt(x,y) dx] =sqrt(pi),
	 *		where the integration is again from -infinity to infinity.
	 *
	 *		Therefore, the area (neglecting the DC offset) is given by:
	 *			Integal[w[1]*voigt(w[2]*x',y) dx'] =sqrt(pi)*w[1]/w[2]
	 *
	 *		The calculation of wv is an approximation given by Whiting.
	 *
	 *		*/

};
#endif