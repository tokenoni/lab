#ifndef _MYGSL_COMPLEX_LONG_DOUBLE_HPP_
#define _MYGSL_COMPLEX_LONG_DOUBLE_HPP_

#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef __GSL_COMPLEX_USE__
namespace mygsl{
	class complex_long_double_long_double : public gsl_complex_long_double_long_double{
	public:
		//---	constructor	---
		complex_long_double(){};
		complex_long_double(double real_part, long double imag_part)	{set(real_part, imag_part);	}
		complex_long_double(const gsl_complex_long_double& obj){GSL_SET_complex_long_double(this, GSL_REAL(obj), GSL_IMAG(obj));}
	
		//---	copy constructor	---
		complex_long_double(const complex_long_double& obj)	  {	GSL_SET_complex_long_double(this, GSL_REAL(obj), GSL_IMAG(obj));}
		gsl_complex_long_double& operator = (const gsl_complex_long_double& a){	GSL_SET_complex_long_double(this, GSL_REAL(a), GSL_IMAG(a));	return dynamic_cast<gsl_complex_long_double & >(*this);}
		gsl_complex_long_double& operator = (const long double a){	GSL_SET_complex_long_double(this, a, 0.0);	return dynamic_cast<gsl_complex_long_double & >(*this);}

		//---	overriden operators for mathematical expressions	---
		complex_long_double operator + (const gsl_complex_long_double& a){	return gsl_complex_long_double_add(     dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator + (const long double a){		return gsl_complex_long_double_add_real(dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator - (const gsl_complex_long_double& a){	return gsl_complex_long_double_sub(     dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator - (const long double a){		return gsl_complex_long_double_sub_real(dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator * (const gsl_complex_long_double& a){	return gsl_complex_long_double_mul(     dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator * (const long double a){		return gsl_complex_long_double_mul_real(dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator / (const gsl_complex_long_double& a){	return gsl_complex_long_double_div(     dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		complex_long_double operator / (const long double a){		return gsl_complex_long_double_div_real(dynamic_cast<gsl_complex_long_double & >(*this), a);		}
		
		complex_long_double& operator +=(const gsl_complex_long_double& a){	(*this)=(*this)+a;	return (*this);}
		complex_long_double& operator +=(const long double a){		(*this)=(*this)+a;	return (*this);}
		complex_long_double& operator -=(const gsl_complex_long_double& a){	(*this)=(*this)-a;	return (*this);}
		complex_long_double& operator -=(const long double a){		(*this)=(*this)-a;	return (*this);}
		complex_long_double& operator *=(const gsl_complex_long_double& a){	(*this)=(*this)*a;	return (*this);}
		complex_long_double& operator *=(const long double a){		(*this)=(*this)*a;	return (*this);}
		complex_long_double& operator /=(const gsl_complex_long_double& a){	(*this)=(*this)/a;	return (*this);}
		complex_long_double& operator /=(const long double a){		(*this)=(*this)/a;	return (*this);}

		bool operator == (const gsl_complex_long_double& a)	{if(real()==GSL_REAL(a) && imag()==GSL_IMAG(a)) return true; else return false;}
		complex_long_double operator - () { return (*this)*(-1.0);}

		complex_long_double conj(){ mygsl::complex_long_double val(real(), -imag()); return val;}
		complex_long_double inv() { mygsl::complex_long_double val = (mygsl::complex_long_double)gsl_complex_long_double_inverse( dynamic_cast<gsl_complex_long_double & >(*this)); return val;}
		double abs(){ return gsl_complex_long_double_abs(dynamic_cast<gsl_complex_long_double & >(*this));}

		//---	tip functions ---
		void set(double real_part, long double imag_part){	GSL_SET_complex_long_double(this, real_part, imag_part);}
		void setreal(double real_part)	{ GSL_SET_REAL(this, real_part);}
		void setimag(double imag_part)	{ GSL_SET_IMAG(this, imag_part);}
		double real()const{return GSL_complex_long_double_P_REAL(this);}
		double imag()const{return GSL_complex_long_double_P_IMAG(this);}

	};

	//---	overriden mathematical functions	---
	static complex_long_double complex_long_double_poler(const long double amp, const long double phase){	return complex_long_double(amp*cos(phase), amp*sin(phase));}
};

//---	overriden mathematical functions	---
static mygsl::complex_long_double exp(const gsl_complex_long_double& a){return (mygsl::complex_long_double)gsl_complex_long_double_exp(a);}
static mygsl::complex_long_double sqrt(const gsl_complex_long_double& a){return (mygsl::complex_long_double)gsl_complex_long_double_sqrt(a);}
static mygsl::complex_long_double inv(const gsl_complex_long_double& a) {return (mygsl::complex_long_double)gsl_complex_long_double_inverse(a);}

#else
namespace mygsl{
	class complex_long_double{
	public:
		//---	constructor	---
		complex_long_double(){};
		complex_long_double(long double real_part_, long double imag_part_)	{real_part=real_part_, imag_part=imag_part_;	}
		complex_long_double(const gsl_complex_long_double& obj){real_part = GSL_REAL(obj);	imag_part = GSL_IMAG(obj);}
	
		//---	copy constructor	---
		complex_long_double(const complex_long_double& obj)	  {real_part=obj.real_part, imag_part=obj.imag_part;	}
		complex_long_double& operator = (const complex_long_double& obj){real_part=obj.real_part, imag_part=obj.imag_part;	return (*this);}
		complex_long_double& operator = (const long double a){real_part=a; imag_part=0.0;	return (*this);}

		//---	elementary functions	---
		complex_long_double conj()const { return complex_long_double(real(), -imag()); }
		double abs2()  const { return real_part*real_part + imag_part*imag_part;}
		double abs()  const { return sqrt(abs2());}
		double arg() const  { return atan2(real_part/imag_part);}
		complex_long_double inv() const 
			{	double abs_inv = 1.0/abs2();
				return complex_long_double(real_part*abs_inv, -imag_part*abs_inv);}

		//---	overriden operators for mathematical expressions	---
		complex_long_double operator + (const complex_long_double& a)const {	return complex_long_double(real_part+a.real(),	imag_part+a.imag());}
		complex_long_double operator + (const long double a)  const {	return complex_long_double(real_part+a,			imag_part);			}
		complex_long_double operator - (const complex_long_double& a)const {	return complex_long_double(real_part-a.real(),	imag_part-a.imag());}
		complex_long_double operator - (const long double a)  const {	return complex_long_double(real_part-a,			imag_part);			}
	
		complex_long_double operator * (const long double a)  const {	return complex_long_double(real_part*a,         imag_part*a);}
		complex_long_double operator / (const long double a)  const {	return complex_long_double(real_part/a,         imag_part/a);}
		complex_long_double operator * (const complex_long_double& a)const {	return complex_long_double(real_part*a.real() - imag_part*a.imag(), real_part+a.imag() + imag_part*a.real());}
		complex_long_double operator / (const complex_long_double& a)const {	return (*this)*a.inv();}
		
		complex_long_double operator - () { return (*this)*(-1.0);}

		complex_long_double& operator +=(const complex_long_double& a){	real_part += a.real();	imag_part += a.imag();	return (*this);}
		complex_long_double& operator +=(const long double a){	real_part += a;									return (*this);}
		complex_long_double& operator -=(const complex_long_double& a){	real_part -= a.real();	imag_part -= a.imag();	return (*this);}
		complex_long_double& operator -=(const long double a){	real_part -= a;									return (*this);}
		complex_long_double& operator *=(const complex_long_double& a){	(*this)=(*this)*a;	return (*this);}
		complex_long_double& operator *=(const long double a){	real_part *= a;			imag_part *=	a;			return (*this);}
		complex_long_double& operator /=(const complex_long_double& a){	(*this)=(*this)/a;	return (*this);}
		complex_long_double& operator /=(const long double a){	real_part /= a;			imag_part /=    a;			return (*this);}

		bool operator == (const complex_long_double& a)	{if(real()==a.real() && imag()==a.imag()) return true; else return false;}


		//---	tip functions ---
		void set(long double real_part_, long double imag_part_){real_part=real_part_;	imag_part=imag_part_;}
		void setreal(long double real_part_)	{ real_part=real_part_;}
		void setimag(long double imag_part_)	{ imag_part=imag_part_;}
		long double real()const{return real_part;}
		long double imag()const{return imag_part;}
	private:
		long double real_part, imag_part;
	};
};
//---	overriden mathematical functions	---
//static mygsl::complex_long_double exp( const mygsl::complex_long_double& a){return (mygsl::complex_long_double)gsl_complex_long_double_exp(a.gsl());}
//static mygsl::complex_long_double sqrt(const mygsl::complex_long_double& a){return (mygsl::complex_long_double)gsl_complex_long_double_sqrt(a.gsl());}
//static mygsl::complex_long_double inv( const mygsl::complex_long_double& a){return (mygsl::complex_long_double)gsl_complex_long_double_inverse(a.gsl());}

#endif

static mygsl::complex_long_double exp_euler(const double& a){return mygsl::complex_long_double(cos(a), sin(a));}

#endif