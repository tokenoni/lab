#ifndef __RAYTRACING_WAVEOPTICS_HPP__
#define __RAYTRACING_WAVEOPTICS_HPP__

#include <vector> 
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <MyGsl\complex.hpp>
#include <MyGsl\constants.hpp>

class PointOnSurface{
public:
	PointOnSurface(void):
		pi( 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679 ),
		twopi(long double(2.0) * pi)
	{};
	~PointOnSurface(void){};
	
	//	core function	lambda in [mm]
	bool Go1Dim(PointOnSurface const* obj, const size_t num, const double lambda);
	bool Go(    PointOnSurface const* Obj, const size_t num, const double lambda);
	const long double pi, twopi;
public:
	//	position of the point
	long double x[3];
	//	resultant amplitude and phase
	mygsl::complex A;
};

class WaveOpticsBase{
public:
	//---	constructors	---
	WaveOpticsBase(const size_t num1Dim_);
	
	std::string name;		//	name of the optics
//---	center position, direction of the optics	---
	void setLocation( double  X_, double  Y_, double  Z_,		//	set center position of the Optics
					  double dX_, double dY_, double dZ_);
	void setLocation( double  X_, double  Y_, double  Z_,
					  double dX_, double dY_, double dZ_,		//	indicate the z direction of the optics
					  double dXz_, double dYz_, double dZz_);	//	only the component perpendicular to the dX vector is considered
	void setLocationCore( double  X_, double  Y_, double  Z_,
					  double dX_, double dY_, double dZ_,		//	indicate the z direction of the optics
					  double dXz_, double dYz_, double dZz_);	//	only the component perpendicular to the dX vector is considered
	
	template<class OpticsType_ >
	void copyLocation(const OpticsType_ &dest);

private:
	//	number of points on the surface
	const size_t num1Dim, num;

	//	Locations of Optics
	long double  X[3];		//	Center of the Optics
	long double dX[3];		//	Direction	"
	long double dXz[3];		
//---	physical property	---
	double diameter, size_y, size_z;	// diameter: diameter of the optics		(available if "shape" is "Circle")

	gsl_matrix *ProjectionMatrix;	//	  projection matrix to the relative coodinate
	gsl_matrix *reProjectionMatrix;	//	reprojection matrix to the absolute coodinate

public:
	enum Shape{							//	size_y, size_z: length of each side (available if "shape" is "Recutangular")
		Circle,
		Rectangular};
	enum Shape shape;
	void getProjectionMatrix(gsl_matrix *Minv, gsl_matrix *M);
	void RelativeCoordinate(double* x);		void RelativeCorrdinate(double *xrel, const double* x){ copy(xrel, x); RelativeCoordinate(xrel);}
	void AbsoluteCoordinate(double* x);		void AbsoluteCorrdinate(double *xrel, const double* x){ copy(xrel, x); AbsoluteCoordinate(xrel);}
	void RelativeDirection(double* x);		void RelativeDirection(double *xrel, const double* x){ copy(xrel, x); RelativeDirection(xrel);}
	void AbsoluteDirection(double* x);		void AbsoluteDirection(double *xrel, const double* x){ copy(xrel, x); AbsoluteDirection(xrel);}

//---	glass name	---
	enum Glass{
		BK7,	BK7Nikon,	SQ,	SF11,	SF10,	BaF10,	BaFN10, SF5, SK11, SSK8, BAL35, SF2,
		BaK1,   BaK4,	SF8
	};

	//---	unallowed constructors	---
private:
	WaveOpticsBase(void);


protected:
//---	Tip functions	---
	inline double dot(const double* xi, const double* xj){
		double inner_product=0.0; for(size_t i=0; i<3; i++) inner_product+=xi[i]*xj[i]; return inner_product;};
	inline void normalize(double* xi){
		double xi_abs = sqrt(dot(xi, xi)); for(size_t i=0; i<3; i++) xi[i] /= xi_abs;}
	inline void copy(double* xi, const double* xj){for(size_t i=0; i<3; i++) xi[i] = xj[i];};
	inline void copy(double* xi, const std::vector<double>& xj){for(size_t i=0; i<3; i++) xi[i] = xj[i];};
	inline void copy(std::vector<double>& xi, const double* xj){for(size_t i=0; i<3; i++) xi[i] = xj[i];};
	
	std::vector<double> get_stdvNAN()
		{ std::vector<double> stdvNAN(3, sqrt(-1.0));return stdvNAN;}
	std::vector<double> stdvNAN;
	double pNAN[3];
	std::vector<double> get_stdv(double* x){std::vector<double> tmp(3); tmp[0]=x[0]; tmp[1]=x[1]; tmp[2]=x[2]; return tmp;}

};


#include "WaveOptics.inl"

#endif