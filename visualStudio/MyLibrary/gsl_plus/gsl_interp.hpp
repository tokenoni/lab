#ifndef _GSL_INTERP_HPP_
#define _GSL_INTERP_HPP_

#include <gsl/gsl_interp.h>
#include <vector>
#include <gsl/gsl_heapsort.h>

class InterpolatedData{
public:
//---	type of interpolation	---
	enum InterpType{
		Linear,
		CubicSpline,
		LinearInLogScale
	};

//---	constructor, destructor ---
	InterpolatedData(void);
	~InterpolatedData(void);
	InterpolatedData(const InterpolatedData& src);
	InterpolatedData operator=(const InterpolatedData& src);
//	void Copy(const InterpolatedData& src);

//---	initializing function	---
	template < class Tyx_, class Tyy_ >
	void set(std::vector< Tyx_ >& xsrc, std::vector< Tyy_ >& ysrc, InterpType type = Linear);
	template < class Tyx_, class Tyy_ >
	void set(Tyx_& xsrc, Tyy_& ysrc, size_t size_, InterpType type = Linear);

//---	obtain the interpolated data	---
	double get(double x0);
	double getderiv(double x0);
	double getderiv2(double x0);
	double getinteg(double xfrom, double xto);

//---	obtain the original data	---
	double operator[] (size_t i){return y[i];}

//---	tip functions	---
	size_t size()const{return num;}
	bool IsAllocated()const{return allocated;}
	double* getXpointer()const{return x;};
	double* getYpointer()const{return y;};
	InterpType getType()const{return Type;};

private:
	double* x;	double* y;
	size_t num;
	bool allocated;
	gsl_interp* interp;
	gsl_interp_accel* acc;
	InterpType Type;
	void MemoryControl(size_t size_, InterpType Type_);
	void FreeMemory();
	void Allocate(size_t size_, InterpType Type_);
	
	void Sort();	//	<---	under construction

};

#include "gsl_interp.inl"

#endif
