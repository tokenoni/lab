#ifndef __LAPLACEINVERSION_H__
#define __LAPLACEINVERSION_H__
#include <vector>
#include <gsl_plus/jacobi.h>

#include <IGOR/IGORitx.hpp>

class LaplaceInversion{
public:
	template < class TyL_ >	//	s = (b + n)*d	( n=1, 2, 3, ..., order )
	void set( TyL_ & Fs, double b_, double d_, size_t order_);
	
	double get(double t);
	double operator()(double t){ return get(t);}
	static double getb(double start, double end, size_t order);
	static double getd(double start, double end, size_t order);
private:
	std::vector<double> C;
	double b, d;
	size_t order;

};

inline double LaplaceInversion::getb(double start, double end, size_t order_){
	return (start*order_ - end)/(end-start);
}
inline double LaplaceInversion::getd(double start, double end, size_t order){
	return (end- start)/(order-1);
}


template < class TyL_ >
inline void LaplaceInversion::set( TyL_ & Fs, double b_, double d_, size_t order_){
	b = b_;		d = d_;		order = order_;
	C.resize(order);
//---	derivation of C	---
	C[0] = Fs[0]*d*(b+1.0);			//<--- i=0
	double term=1.0;
	for(size_t i=1; i<order; ++i){	
		double term_sum=0.0;
		for(size_t j=0; j<i; ++j){
			term = 1.0;
			for(size_t k=0; k<j; ++k)
				term *= (1.0*i-k)/(b+1.0+i+k);
			term /= (b+1.0+i+j);
			term_sum += C[j]*term;
		}
		term = 1.0;
		for(size_t k=0; k<i; ++k)
			term *= (1.0*i-k)/(b+1.0+i+k);
		term /= (b+1.0+2.0*i);
		C[(size_t)i] = (d*Fs[(size_t)i]-term_sum)/term;
	}
}

inline double LaplaceInversion::get(double t){
	double value = 0.0;
	double x = 2.0*exp(-d*t) -1.0;
	for(size_t n=0; n<order; ++n)
		value += C[n]*jac_jacobi(x, n, 0.0, b);
	return value;
}
#endif