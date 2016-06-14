#ifndef _MYGSL_COMPLEX_HPP_
#define _MYGSL_COMPLEX_HPP_

#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#ifdef __GSL_COMPLEX_USE__
namespace mygsl{
	class complex : public gsl_complex{
	public:
		//---	constructor	---
		complex(){};
		complex(double real_part, double imag_part)	{set(real_part, imag_part);	}
		complex(const gsl_complex& obj){GSL_SET_COMPLEX(this, GSL_REAL(obj), GSL_IMAG(obj));}
	
		//---	copy constructor	---
		complex(const complex& obj)	  {	GSL_SET_COMPLEX(this, GSL_REAL(obj), GSL_IMAG(obj));}
		gsl_complex& operator = (const gsl_complex& a){	GSL_SET_COMPLEX(this, GSL_REAL(a), GSL_IMAG(a));	return dynamic_cast<gsl_complex & >(*this);}
		gsl_complex& operator = (const double a){	GSL_SET_COMPLEX(this, a, 0.0);	return dynamic_cast<gsl_complex & >(*this);}

		//---	overriden operators for mathematical expressions	---
		complex operator + (const gsl_complex& a){	return gsl_complex_add(     dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator + (const double a){		return gsl_complex_add_real(dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator - (const gsl_complex& a){	return gsl_complex_sub(     dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator - (const double a){		return gsl_complex_sub_real(dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator * (const gsl_complex& a){	return gsl_complex_mul(     dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator * (const double a){		return gsl_complex_mul_real(dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator / (const gsl_complex& a){	return gsl_complex_div(     dynamic_cast<gsl_complex & >(*this), a);		}
		complex operator / (const double a){		return gsl_complex_div_real(dynamic_cast<gsl_complex & >(*this), a);		}
		
		complex& operator +=(const gsl_complex& a){	(*this)=(*this)+a;	return (*this);}
		complex& operator +=(const double a){		(*this)=(*this)+a;	return (*this);}
		complex& operator -=(const gsl_complex& a){	(*this)=(*this)-a;	return (*this);}
		complex& operator -=(const double a){		(*this)=(*this)-a;	return (*this);}
		complex& operator *=(const gsl_complex& a){	(*this)=(*this)*a;	return (*this);}
		complex& operator *=(const double a){		(*this)=(*this)*a;	return (*this);}
		complex& operator /=(const gsl_complex& a){	(*this)=(*this)/a;	return (*this);}
		complex& operator /=(const double a){		(*this)=(*this)/a;	return (*this);}

		bool operator == (const gsl_complex& a)	{if(real()==GSL_REAL(a) && imag()==GSL_IMAG(a)) return true; else return false;}
		complex operator - () { return (*this)*(-1.0);}

		complex conj(){ mygsl::complex val(real(), -imag()); return val;}
		complex inv() { mygsl::complex val = (mygsl::complex)gsl_complex_inverse( dynamic_cast<gsl_complex & >(*this)); return val;}
		double abs(){ return gsl_complex_abs(dynamic_cast<gsl_complex & >(*this));}

		//---	tip functions ---
		void set(double real_part, double imag_part){	GSL_SET_COMPLEX(this, real_part, imag_part);}
		void setreal(double real_part)	{ GSL_SET_REAL(this, real_part);}
		void setimag(double imag_part)	{ GSL_SET_IMAG(this, imag_part);}
		double real()const{return GSL_COMPLEX_P_REAL(this);}
		double imag()const{return GSL_COMPLEX_P_IMAG(this);}

	};

	//---	overriden mathematical functions	---
	static complex complex_poler(const double amp, const double phase){	return complex(amp*cos(phase), amp*sin(phase));}
};

//---	overriden mathematical functions	---
static mygsl::complex exp(const gsl_complex& a){return (mygsl::complex)gsl_complex_exp(a);}
static mygsl::complex sqrt(const gsl_complex& a){return (mygsl::complex)gsl_complex_sqrt(a);}
static mygsl::complex inv(const gsl_complex& a) {return (mygsl::complex)gsl_complex_inverse(a);}

#else
namespace mygsl{
	class complex{
	public:
		//---	constructor	---
		complex(){};
		complex(double real_part_, double imag_part_)	{real_part=real_part_, imag_part=imag_part_;	}
		complex(const gsl_complex& obj){real_part = GSL_REAL(obj);	imag_part = GSL_IMAG(obj);}
	
		//---	copy constructor	---
		complex(const complex& obj)	  {real_part=obj.real_part, imag_part=obj.imag_part;	}
		complex& operator = (const complex& obj){real_part=obj.real_part, imag_part=obj.imag_part;	return (*this);}
		complex& operator = (const double a){real_part=a; imag_part=0.0;	return (*this);}

		//---	elementary functions	---
		complex conj()const { return complex(real(), -imag()); }
		double abs2()  const { return real_part*real_part + imag_part*imag_part;}
		double abs()  const { return sqrt(abs2());}
		complex inv() const 
			{	double abs_inv = 1.0/abs2();
				return complex(real_part*abs_inv, -imag_part*abs_inv);}

		//---	overriden operators for mathematical expressions	---
		complex operator + (const complex& a)const {	return complex(real_part+a.real(),	imag_part+a.imag());}
		complex operator + (const double a)  const {	return complex(real_part+a,			imag_part);			}
		complex operator - (const complex& a)const {	return complex(real_part-a.real(),	imag_part-a.imag());}
		complex operator - (const double a)  const {	return complex(real_part-a,			imag_part);			}
	
		complex operator * (const double a)  const {	return complex(real_part*a,         imag_part*a);}
		complex operator / (const double a)  const {	return complex(real_part/a,         imag_part/a);}
		complex operator * (const complex& a)const {	return complex(real_part*a.real() - imag_part*a.imag(), real_part+a.imag() + imag_part*a.real());}
		complex operator / (const complex& a)const {	return (*this)*a.inv();}
		
		complex operator - () { return (*this)*(-1.0);}

		complex& operator +=(const complex& a){	real_part += a.real();	imag_part += a.imag();	return (*this);}
		complex& operator +=(const double a){	real_part += a;									return (*this);}
		complex& operator -=(const complex& a){	real_part -= a.real();	imag_part -= a.imag();	return (*this);}
		complex& operator -=(const double a){	real_part -= a;									return (*this);}
		complex& operator *=(const complex& a){	(*this)=(*this)*a;	return (*this);}
		complex& operator *=(const double a){	real_part *= a;			imag_part *=	a;			return (*this);}
		complex& operator /=(const complex& a){	(*this)=(*this)/a;	return (*this);}
		complex& operator /=(const double a){	real_part /= a;			imag_part /=    a;			return (*this);}

		bool operator == (const complex& a)	{if(real()==a.real() && imag()==a.imag()) return true; else return false;}


		//---	tip functions ---
		void set(double real_part_, double imag_part_){real_part=real_part_;	imag_part=imag_part_;}
		void setreal(double real_part_)	{ real_part=real_part_;}
		void setimag(double imag_part_)	{ imag_part=imag_part_;}
		double real()const{return real_part;}
		double imag()const{return imag_part;}
		gsl_complex gsl()const{return gsl_complex_rect(real_part, imag_part);}
	private:
		double real_part, imag_part;
	};
};
//---	overriden mathematical functions	---
static mygsl::complex exp( const mygsl::complex& a){return (mygsl::complex)gsl_complex_exp(a.gsl());}
static mygsl::complex sqrt(const mygsl::complex& a){return (mygsl::complex)gsl_complex_sqrt(a.gsl());}
static mygsl::complex inv( const mygsl::complex& a){return (mygsl::complex)gsl_complex_inverse(a.gsl());}

#endif

static mygsl::complex exp_euler(const double& a){return mygsl::complex(cos(a), sin(a));}

#endif