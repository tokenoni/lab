#ifndef __MYGSL_GAUSSIAN_EXP_INV_HPP__
#define __MYGSL_GAUSSIAN_EXP_INV_HPP__

#include <IGOR\IGORitx.hpp>
#include <vector>
#include <string>
#include <MyGsl\interp.hpp>

namespace mygsl{

//	const std::string gaussian_exp_inv_dir = "D:\\MyDocument\\My Dropbox\\visual_studio\\MyLibrary\\MyGsl\\gaussian_integration_data\\";
	const std::string gaussian_exp_inv_dir = "C:\\Users\\KeisukeFujii\\Dropbox\\visual_studio\\MyLibrary\\MyGsl\\gaussian_integration_data\\";
		
	//	function returns the value of
	//					integ[0,inf]{ exp[-v*v -x/v]dv }
	class gaussian_integ_exp_inv_F0{
	public:
		gaussian_integ_exp_inv_F0(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F0_x.itx");
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			return data.get(x);
		};
		double get(const double alpha, const double r){
			return 1.0/sqrt(alpha)*data.get(r*sqrt(alpha));
		}
	private:
		mygsl::interp data;
	};
	
	//	function returns the value of
		//					integ[0,inf]{ exp[-x/v]dv (xprojection vMaxwellian) }
	class gaussian_integ_exp_inv_F1_2{
	public:
		gaussian_integ_exp_inv_F1_2(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F1_2_x.itx");
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			if(fabs(x)>data.get_xmax()) return 0.0;
			return data.get(fabs(x));
		};
		double get(const double alpha, const double r){
			return 1.0/sqrt(alpha)*(operator()(r*sqrt(alpha)));
		}
	private:
		mygsl::interp data;
	};

	//	function returns the value of
	//					integ[0,inf]{ v*exp[-v*v -x/v]dv }
	class gaussian_integ_exp_inv_F1{
	public:
		gaussian_integ_exp_inv_F1(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F1_x.itx");
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			return data.get(x);
		};
		double get(const double alpha, const double r){
			return 1.0/alpha*data.get(r*sqrt(alpha));



		}
	private:
		mygsl::interp data;
	};

	//	function returns the value of
	//					integ[0,inf]{ v*v*exp[-v*v -x/v]dv }
	class gaussian_integ_exp_inv_F2{
	public:
		gaussian_integ_exp_inv_F2(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F2_x.itx");
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			return data.get(x);
		};
		double get(const double alpha, const double r){
			return 2.0/(alpha*sqrt(alpha))*data.get(r*sqrt(alpha));
		}
	private:
		mygsl::interp data;
	};

	//	function returns the value of
	//					integ[0,inf]{ v*v*v*exp[-v*v -x/v]dv }
	class gaussian_integ_exp_inv_F3{
	public:
		gaussian_integ_exp_inv_F3(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F3_x.itx");
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			return data.get(x);
		};
		double get(const double alpha, const double r){
			return 2.0/(alpha*alpha)*data.get(r*sqrt(alpha));
		}
	private:
		mygsl::interp data;
	};

	//	function returns the value of
	//					integ[0,inf]{ 1.0/v * exp[-v*v -x/v]dv }
	class gaussian_integ_exp_inv_F_1{
	public:
		gaussian_integ_exp_inv_F_1(void){
			std::vector<double> x = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "xx.itx");
			std::vector<double> f = IGORdata::itx_1stdv(gaussian_exp_inv_dir + "exp_inv_F0_x.itx");
			f[0] = 0.0;	//	x=0 is dummy data. this function diverges at x=0
			xmin = x[1];
			fmin = f[1];
			data.set(x, f, mygsl::interp::Cspline);
		};
		double operator() (const double x){
			return get(x);}
		double get(const double x){
			if(x<xmin) return log(x)/log(xmin)*fmin;
			else return data.get(x);
		};
		double get(const double alpha, const double r){
			return get(r*sqrt(alpha));
		}
	private:
		mygsl::interp data;
		double xmin,fmin;
	};

};

#endif