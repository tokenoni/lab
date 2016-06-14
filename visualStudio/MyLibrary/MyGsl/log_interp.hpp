#ifndef __MYGSL_LOG_INTERP_HPP__
#define __MYGSL_LOG_INTERP_HPP__

#include "interp.hpp"

namespace mygsl{
//-------------------------------------------------//
//												   //
//---		1 dimensional interpolation			---//
//							in log(x) scale		   //
//												   //
//-------------------------------------------------//
	class interp_logx : public interp{
	public:
		interp_logx(const size_t size_, InterpMethod method_ = Cspline): interp(size_, method_){};
		interp_logx(void){};
		~interp_logx(void){};
		//--- setting function	---
		template < class Tyv_ >
		void set(const Tyv_& x, const Tyv_& y, InterpMethod method_ = Cspline);	
		//---	for directly setting the contained data	---
		void set(const size_t i, const double xi, const double yi);
		void set_x(const size_t i, const double xi);

		//---	operations	---
		std::vector<double> get_original_x();

		//---	high level functions	---
		double get(double x)      ; 
		double getderiv(double x) ;
		double getderiv2(double x);
		double getinteg(double a, double b);

	private:
		const double* get_original_x_p()const {return x_array;}
		double& get_original_x(size_t i);

		double getinteg_linear(double a_, double b_);
		double getinteg_spline(double a_, double b_);
	};

//-------------------------------------------------//
//							in log(y) scale		   //
//-------------------------------------------------//
	class interp_logy : public interp{
	public:
		interp_logy(void){};
		interp_logy(const size_t size_, InterpMethod method_ = Cspline): interp(size_, method_){};
		//--- setting function	---
		template < class Tyv_ >
		void set(const Tyv_& x, const Tyv_& y, InterpMethod method_ = Cspline);	
		//---	for directly setting the contained data	---
		void set(const size_t i, const double xi, const double yi);
		void set_y(const size_t i, const double yi);

		//---	operations	---
		std::vector<double> get_original_y();

		//---	high level functions	---
		double get(double x)      ; 
		double getderiv(double x) ;
		double getderiv2(double x);
		double getinteg(double a, double b);
	private:
		double getinteg_linear(double a_, double b_);
		double getinteg_spline(double a_, double b_);
	};


//-------------------------------------------------//
//					in log(x)-log(y) scale		   //
//-------------------------------------------------//
	class interp_logxy : public interp{
	public:
		interp_logxy(const size_t size_, InterpMethod method_ = Cspline): interp(size_, method_){};
		interp_logxy(void):interp(){};
		//--- setting function	---
		template < class Tyv_ >
		void set(const Tyv_& x, const Tyv_& y, InterpMethod method_ = Cspline);	
		//---	for directly setting the contained data	---
		void set(const size_t i, const double xi, const double yi);
		void set_x(const size_t i, const double yi);
		void set_y(const size_t i, const double yi);

		//---	operations	---
		std::vector<double> get_original_x();
		std::vector<double> get_original_y();

		//---	high level functions	---
		double get(double x)      ; 
/*		double getderiv(double x) ;
		double getderiv2(double x);
		double getinteg(double a, double b);
*/	private:
		double getinteg_linear(double a_, double b_);
		double getinteg_spline(double a_, double b_);
	};


//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//							in log(x) scale		   //
//												   //
//-------------------------------------------------//
	class interp2d_logx : public interp2d{
	public:
		void set(const std::vector<double>& x, const std::vector<double>& y,
				 const std::vector<std::vector<double>>& data_, InterpMethod method_ = Bilinear);
		//---	high level functions	---
		double get(double x, double y);
		double operator () (const double x, const double y){return get(x, y);}
		double getderiv_x(  const double x, const double y);	
		double getderiv_xy( const double x, const double y);

		//---	get interp-type data at yv=y	---
		interp_logx getAlong_x_log(const double y, interp::InterpMethod method_ = interp::Cspline);
	private:
		double& get_original_x(size_t i){return xv[i];}
		interp getAlongLine(const double x0, const double y0, const double x1, const double y1, const size_t num, interp::InterpMethod method_ = interp::Cspline);
		interp getAlong_x(const double y, interp::InterpMethod method_ = interp::Cspline);
	};

//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//							in log(x) scale		   //
//												   //
//-------------------------------------------------//
	class interp2d_logy : public interp2d{
	public:
		void set(const std::vector<double>& x, const std::vector<double>& y,
				 const std::vector<std::vector<double>>& data_, InterpMethod method_ = Bilinear);
		//---	high level functions	---
		double get(double x, double y);
		double operator () (const double x, const double y){return get(x, y);}
		double getderiv_y(  const double x, const double y);	
		double getderiv_xy( const double x, const double y);

		//---	get interp-type data at xv=x	---
		interp_logx getAlong_y_log(const double x, interp::InterpMethod method_ = interp::Cspline);
	private:
		double& get_original_y(size_t i){return yv[i];}
		//---	get interp-type data along the line determined as [x0,y0], [x1,y1]. the x coordinate of the interp is r, which is defined as the length from [x0,y0]	---
		interp getAlongLine(const double x0, const double y0, const double x1, const double y1, const size_t num, interp::InterpMethod method_ = interp::Cspline);
		interp getAlong_y(const double x, interp::InterpMethod method_ = interp::Cspline);
	};

//-------------------------------------------------//
//												   //
//---		2 dimensional interpolation			---//
//          		in log(x)-log(y) scale		   //
//												   //
//-------------------------------------------------//
	class interp2d_logxy : public interp2d{
	public:
		void set(const std::vector<double>& x, const std::vector<double>& y,
				 const std::vector<std::vector<double>>& data_, InterpMethod method_ = Bilinear);
		//---	high level functions	---
		double get(double x, double y);
		double operator () (const double x, const double y){return get(x, y);}
		double getderiv_x(  const double x, const double y);	
		double getderiv_y(  const double x, const double y);	
		double getderiv_xy( const double x, const double y);

		//---	get interp-type data at yv=y	---
		interp_logx getAlong_x_log(const double y, interp::InterpMethod method_ = interp::Cspline);
		//---	get interp-type data at xv=x	---
		interp_logx getAlong_y_log(const double x, interp::InterpMethod method_ = interp::Cspline);
	private:
		double& get_original_x(size_t i){return xv[i];}
		double& get_original_y(size_t i){return yv[i];}
		interp getAlongLine(const double x0, const double y0, const double x1, const double y1, const size_t num, interp::InterpMethod method_ = interp::Cspline);

		interp getAlong_x(const double y, interp::InterpMethod method_ = interp::Cspline);
		interp getAlong_y(const double x, interp::InterpMethod method_ = interp::Cspline);
	};

};
#include "log_interp.inl"
#endif