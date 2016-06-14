#ifndef __MYGSL_INTERP_COMPLEX_HPP__
#define __MYGSL_INTERP_COMPLEX_HPP__

#include <MyGsl\complex.hpp>
#include <MyGsl\interp.hpp>

namespace mygsl{
	class interpComplex{
	public:
		void set(std::vector<double>& x, std::vector<double>& yreal, std::vector<double>& yimag, mygsl::interp::InterpMethod method=mygsl::interp::Cspline);
		void set(std::vector<double>& x, std::vector<mygsl::complex>& y, mygsl::interp::InterpMethod method=mygsl::interp::Cspline);
		void setSize(const size_t size, mygsl::interp::InterpMethod method=mygsl::interp::Cspline);
		void setInterp();
		void set(const size_t i, const double x, const mygsl::complex z);
		void set_x(const size_t i, const double x);
		void set_y(const size_t i, const mygsl::complex z);
		void set_y_real(const size_t i, const double y_real);
		void set_y_imag(const size_t i, const double y_imag);

		mygsl::complex get(double x);
		mygsl::complex getderiv(double x);
		mygsl::complex getderiv2(double x);
		mygsl::complex getinteg(double a, double b);
		std::vector<double> get_original_x(){ return real.get_original_x();}
		std::vector<double> get_original_y_real(){return real.get_original_y();}
		std::vector<double> get_original_y_imag(){return imag.get_original_y();}
		std::vector<mygsl::complex> get_original_y();
	private:
		interp real, imag;
	};
};

//------------------	definitions	------------------
namespace mygsl{
	inline void interpComplex::set(std::vector<double>& x, std::vector<double>& yreal, std::vector<double>& yimag, mygsl::interp::InterpMethod method){
		real.set(x, yreal, method);
		imag.set(x, yimag, method);
	}

	inline void interpComplex::set(std::vector<double>& x, std::vector<mygsl::complex>& y, mygsl::interp::InterpMethod method){
		setSize(x.size(), method);
		for(size_t i=0; i<x.size(); ++i) set(i, x[i], y[i]);
		setInterp();
	}

	inline void interpComplex::setInterp(){
		real.setInterp();
		imag.setInterp();
	}
	inline void interpComplex::setSize(const size_t size, mygsl::interp::InterpMethod method){
		real.setSize(size, method);
		imag.setSize(size, method);
	}
	inline void interpComplex::set(const size_t i, const double x, const mygsl::complex z){
		real.set(i, x, z.real());
		imag.set(i, x, z.imag());
	}
	inline void interpComplex::set_x(const size_t i, const double x){
		real.set_x(i, x);
		imag.set_x(i, x);
	}
	inline void interpComplex::set_y(const size_t i, const mygsl::complex z){
		real.set_y(i, z.real());
		imag.set_y(i, z.imag());
	}
	inline void interpComplex::set_y_real(const size_t i, const double y_real){
		real.set_y(i, y_real);
	}
	inline void interpComplex::set_y_imag(const size_t i, const double y_imag){
		imag.set_y(i, y_imag);
	}
	inline mygsl::complex interpComplex::get(double x){
		return mygsl::complex(real.get(x), imag.get(x));
	}
	inline mygsl::complex interpComplex::getderiv(double x){
		return mygsl::complex(real.getderiv(x), imag.getderiv(x));
	}
	inline mygsl::complex interpComplex::getderiv2(double x){
		return mygsl::complex(real.getderiv2(x), imag.getderiv2(x));
	}
	inline mygsl::complex interpComplex::getinteg(double a, double b){
		return mygsl::complex(real.getinteg(a, b), imag.getinteg(a, b));
	}
	inline std::vector<mygsl::complex> interpComplex::get_original_y(){
		std::vector<mygsl::complex> zv(real.size());
		for(size_t i=0; i<real.size(); ++i){
			zv[i] = mygsl::complex(real.get_original_y(i), imag.get_original_y(i));
		}
		return zv;
	};

};
#endif