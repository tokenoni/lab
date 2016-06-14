#ifndef __MYGSL_INTERP_HPP__
#define __MYGSL_INTERP_HPP__

#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

namespace mygsl{
//-------------------------------------------------//
//												   //
//---		1 dimensional interpolation			---//
//												   //
//-------------------------------------------------//
	class interp{
	public:
		enum InterpMethod{
			Cspline,			//	cubic spline widh natural edges	(second derivative at the edge points is zero)
			Linear
//			CsplineFlatEdge,	//	cubic spline with flat edges	(first derivative at the edge points is zero)
//				now the "Cspline flat edge" version is not available. 
//				it should be separated into another source file 
//				like "interp_edges_controlled.hpp" for fastening the speed
		};
		//---	constructors	
		interp(void);
		interp(const size_t size_, InterpMethod method_ = Cspline);
		interp(const interp& obj);
		interp& operator = (const interp& obj);
		~interp(void);
		
		//--- setting function	---
		template < class Tyv_ >
		void set(const Tyv_& x, const Tyv_& y, InterpMethod method_ = Cspline);
		
		//---	for directly setting the contained data	---
		void setSize(const size_t size, InterpMethod method_ = Cspline);
		void setInterp(bool isSorted = false);
		void sort();
		void set(const size_t i, const double xi, const double yi);
		void set_x(const size_t i, const double xi);
		void set_y(const size_t i, const double yi);
		void set_y_zero();

		//---	operations	---
		double& operator[] (size_t i){return data[i];}
		double& get_original_x(size_t i){return x_array[i];}
		double const* get_original_x_p()const {return x_array;}
		double& get_original_y(size_t i){return data[i];}
		std::vector<double> get_original_x()const;
		std::vector<double> get_original_y()const;
		void clear();

		//---	high level functions	---
		double get(double x);					//	get value at x
		double get(double x, gsl_interp_accel *acc)const;
		double getderiv(double x);				//	get 1st derivative
		double getderiv(double x, gsl_interp_accel *acc)const;
		double getderiv2(double x);				//	get 2nd derivative
		double getderiv2(double x, gsl_interp_accel *acc)const;
		double getinteg(double a, double b);	//	get integration from [a, b]
		double getinteg(double a, double b, gsl_interp_accel *acc)const;
		
		//---	more high level functions	---
		//		the "x_array" vector of the original one (or that in first term) is used 
		//		users should check which x_array vector is appropriate for the interpolation
		void add(interp& obj);		//	summation
		void sub(interp& obj);		//	subtraction
		void mul(interp& obj);		//	multiplication
		void div(interp& obj);		//	devide
		void add(const double obj);		//	
		void sub(const double obj);		//	
		void mul(const double obj);		//	
		void div(const double obj);		//
		void inverse();					//	the contained value is inversed (data = 1.0/data)

		interp& operator += (interp& obj) {      add(obj); return *this;};
		interp& operator -= (interp& obj) {      sub(obj); return *this;};
		interp& operator *= (interp& obj) {      mul(obj); return *this;};
		interp& operator /= (interp& obj) {      div(obj); return *this;};
		interp& operator += (const double obj){  add(obj); return *this;};
		interp& operator -= (const double obj){  sub(obj); return *this;};
		interp& operator *= (const double obj){  mul(obj); return *this;};
		interp& operator /= (const double obj){  div(obj); return *this;};

		//		!!  caution  !!
		//		these functions are slow because many data are duplicated (and memory allocation)
		interp add_copy(interp&      obj){interp copy = *this; copy.add(obj); return copy;};		//	summation
		interp sub_copy(interp&      obj){interp copy = *this; copy.sub(obj); return copy;};		//	subtraction
		interp mul_copy(interp&      obj){interp copy = *this; copy.mul(obj); return copy;};		//	multiplication
		interp div_copy(interp&      obj){interp copy = *this; copy.div(obj); return copy;};		//	devide
		interp add_copy(const double obj){interp copy = *this; copy.add(obj); return copy;};		//	
		interp sub_copy(const double obj){interp copy = *this; copy.sub(obj); return copy;};		//	
		interp mul_copy(const double obj){interp copy = *this; copy.mul(obj); return copy;};		//	
		interp div_copy(const double obj){interp copy = *this; copy.div(obj); return copy;};		//
		
		interp operator + (interp& obj){     return add_copy(obj);};
		interp operator - (interp& obj){     return sub_copy(obj);};
		interp operator * (interp& obj){     return mul_copy(obj);};
		interp operator / (interp& obj){     return div_copy(obj);};
		interp operator + (const double obj){return add_copy(obj);};
		interp operator - (const double obj){return sub_copy(obj);};
		interp operator * (const double obj){return mul_copy(obj);};
		interp operator / (const double obj){return div_copy(obj);};

		interp& operator - (int){ mul(-1.0); return *this;};		//	

		interp getderiv( InterpMethod method_ = Cspline);
		interp getderiv2(InterpMethod method_ = Cspline);
		interp getinteg( const double xfrom,  InterpMethod method_ = Cspline);

		//--- convolution with a interp object	---
		// integ[ f(x-tau)*g(tau) ]dtau
		// x coordinate is same to this* object	---
		interp convolute(interp& obj);
		//---	tip functions	---
		size_t size()const{return num;}
		bool isallocated()const{return allocated;}
		bool issetInterped()const{return setInterped;}
		double get_xmin()const{return x_array[0];}
		double get_xmax()const{return x_array[size()-1];}

		//---	output the acceralation method	---
		gsl_interp_accel* getAcc(){return acc;};
		void setAcc(const gsl_interp_accel *acc_obj)
			{	acc->cache = acc_obj->cache; acc->miss_count = acc_obj->miss_count; acc->hit_count = acc_obj->hit_count;}

	protected:
		InterpMethod method;
		gsl_interp* interp_obj;
		gsl_interp_accel *acc;
		bool allocated, setInterped;
		double *x_array;
		double *data;
		size_t num;
		void FreeMemory();
		void Allocate(size_t size_);
		void merge(const double *obj, const size_t size_obj);		//	merge the two x_arrayition vectors
		
		//---	polynominal coefficients in [ x(i),x(i+1) )	---
		//		the value is 
		//		value = c[0] + c[1]*dx + c[2]*dx*dx + c[3]*dx*dx*dx
		//			(	dx = x - x[i]	)
		void getCoef(const size_t index, double c[4]);
		void getCoefLinear(const size_t index, double c[4]);
		void getCoefCspline(const size_t index, double c[4]);


		//--- integration used in a convolution function	---
		//	integ[ f(x-tau)*g(tau) ]dtau
		double IntegConv(interp& obj, const double x);
		//	integ(xi,xi+1)[
		double IntegConv1Interval(interp& obj, size_t& j, const double x, const size_t i, const double dx);
		//	integ[ (a0+a1*x+a2*x^2+a3*x^3)*(b0+b1*x+b2*x^2+b3*x^3) ](0,dx)
		static double IntegCoef(const double a[4], const double b[4], const double dx);
		//	change the polynominal coefficient as f(x) -> f(s-x)
		static void changeCoef(const double coef[4], double coef_new[4], const double s);
		//---	initialize for flat edge splines	---
		static int CsplineClampedInit(gsl_interp * interp, const double xa[], const double ya[], size_t size, const double dydx0, const double dydxn);

	};

//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//												   //
//-------------------------------------------------//
	class interp2d{
	public:
		enum InterpMethod{
			Bilinear,
			Bicubic
		};
		//	constructors	
		interp2d(void);
		interp2d(const interp2d& obj);
		interp2d& operator = (const interp2d& obj);
		~interp2d(void);

		//---	setting functions	---
		void set(const std::vector<double>& x, const std::vector<double>& y,
				 const std::vector<std::vector<double>>& data_, InterpMethod method_ = Bilinear);
		
		//---	operations	---
		double& get_original_x(size_t i){return xv[i];}
		double& get_original_y(size_t i){return yv[i];}

		//---	high level functions	---
		double get(double x, double y);
		double operator () (const double x, const double y);
		double getderiv_x(  const double x, const double y);	
		double getderiv_y(  const double x, const double y);	
		double getderiv_xy( const double x, const double y);
//		void getGradient(const double x, const double y, double& df_dx, double& df_dy);
		double getIntegAlong_x(  const double y, const double xfrom, const double xto);
		double getIntegAlong_y(  const double x, const double yfrom, const double yto);
		double getIntegAlongLine(const double xfrom, const double xto, const double yfrom, const double yto);

		//---	get interp-type data at yv=y	---
		interp getAlong_x(const double y, interp::InterpMethod method_ = interp::Cspline);
		//---	get interp-type data at xv=x	---
		interp getAlong_y(const double x, interp::InterpMethod method_ = interp::Cspline);
		//---	get interp-type data along the line determined as [x0,y0], [x1,y1]. the x coordinate of the interp is r, which is defined as the length from [x0,y0]	---
		interp getAlongLine(const double x0, const double y0, const double x1, const double y1, const size_t num, interp::InterpMethod method_ = interp::Cspline);

		void clear();
		
		//---	tip functions	---
		size_t size1()const{return num1;}
		size_t size2()const{return num2;}
		bool isallocated()const{return allocated;}

	protected:
		InterpMethod method;
		gsl_interp_accel *acc_x, *acc_y;
		bool allocated;
		double *xv, *yv;
		double **data;
		size_t num1, num2;
		void FreeMemory();
		void Allocate(size_t size1_, size_t size2_);
		double getBilinear(    const double x, const double y);
		double getBilinear_dx( const double x, const double y);
		double getBilinear_dy( const double x, const double y);
		double getBilinear_dxy(const double x, const double y);
//		void getBilinear_gradient(const double x, const double y, double );
		double getBilinearIntegAlong_x(  const double y, const double xfrom, const double xto);
		double getBilinearIntegAlong_y(  const double x, const double yfrom, const double yto);
		double getBilinearIntegAlongLine(const double xfrom, const double xto, const double yfrom, const double yto);
		double getBilinearIntegAlongLine_x_integ(const double xfrom, const double xto, const double yfrom, const double yto);
		double getBilinearIntegAlongLine_y_integ(const double xfrom, const double xto, const double yfrom, const double yto);

		double getBicubic(    const double x, const double y);
		double getBicubic_dx( const double x, const double y);
		double getBicubic_dy( const double x, const double y);
		double getBicubic_dxy(const double x, const double y);
		void getBicubicCoef(size_t ix_floor, size_t iy_floor);
		//	constants for bicubic interpolation	---
		double f[16];
		double a[4][4];

	};

//---	for gsl cspline interpolation	---
	typedef struct
	{
		double * c;
		double * g;
		double * diag;
		double * offdiag;
	} cspline_state_t;


};

#include "interp.inl"
#endif
