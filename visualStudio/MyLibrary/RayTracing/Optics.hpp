#ifndef __RAYTRACING_OPTICS_HPP__
#define __RAYTRACING_OPTICS_HPP__

#include <vector> 
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#include "Rays.hpp"


//---	very basic class for every optics	---
//	this class provide several basic funstions 
//	every optics class inherits this class. (directory of indirectory) 
class Optics
{
public:
//---	constructors	---
	Optics(void){OpticsExplicitConstructor();}
	Optics(const Optics& dest){OpticsExplicitConstructor(); copyLocation(dest);}
	~Optics(void);
	void OpticsExplicitConstructor(void);
	template<class OpticsType_ >
	OpticsType_& operator = (const OpticsType_& dest){	copyLocation(dest); return *this;}
	
	template<class OpticsType_ >
	void copyLocation(const OpticsType_ &dest);

public:
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
	
	double  X[3];		//	Center of the Optics
	double dX[3];		//	Direction	"
	double dXz[3];		
//---	physical property	---
	double diameter, size_y, size_z;	// diameter: diameter of the optics		(available if "shape" is "Circle")

	double x_tmp[3];
	double ray_pos[3];	
	Ray ray_tmp;

	std::vector<double> stdv_ray_pos;
private:
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

public:
	void getRayTrajectory(RayTrajectory& ray_traj){ray_traj.set(ray_pos);}
	std::vector<double> getRayTrajectory(){for(size_t i=0; i<3; i++) stdv_ray_pos[i] = ray_pos[i]; return stdv_ray_pos;}

//---	glass name	---
	enum Glass{
		BK7,	BK7Nikon,	SQ,	SF11,	SF10,	BaF10,	BaFN10, SF5, SK11, SSK8, BAL35, SF2,
		BaK1,   BaK4,	SF8
	};

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

//------------------	class for Creating Rays		-----------------------------------
class CreateRay :public Optics{
public:
	CreateRay(void){RaySet_Relative.resize(0);};
	CreateRay& operator =(const CreateRay& dest){
		copyLocation(dest);
		RaySet_Relative = dest.RaySet_Relative;
		RaySet_Absolute = dest.RaySet_Absolute;	return *this;
	}
	CreateRay(const CreateRay& dest) :Optics(dest){ *this = dest;};

	std::vector<Ray>& Go();
	std::vector<Ray> RaySet_Relative, RaySet_Absolute;

	void PointSource(size_t NumRay, double SpreadAngle_in_radian){SpreadingSource(NumRay,SpreadAngle_in_radian,0.0,0.0);}
	void PointGaussianSource(size_t NumRay, double FWHMAngle_in_radian, double MaxAngle_in_radian)
		{SpreadingGaussianSource(NumRay,FWHMAngle_in_radian, MaxAngle_in_radian, 0.0, 0.0);}
	void PointPrallelSource(){PointSource(0, 0.0);};
	void ParallelSource(size_t NumRay, double diameter){SpreadingSource(NumRay,0.0,diameter);}
	void ParallelSource(size_t NumRay, double size_y, double size_z){SpreadingSource(NumRay,0.0,size_y, size_z);}
	void SpreadingSource(size_t NumRay, double SpreadAngle_in_radian, double diameter);
	void SpreadingGaussianSource(const size_t NumRay, const double FWHMAngle_in_radian, const double MaxAngle_in_radian, const double diameter);
	void SpreadingSource(const size_t NumRay, const double FWHMAngle_in_radian, const double size_y, double size_z);
	void SpreadingGaussianSource(const size_t NumRay, const double FWHMAngle_in_radian, const double MaxAngle_in_radian, double size_y, double size_z);
	void PointingSource(const size_t NumRay, const double SpreadAngle_in_radian, const double FocalLength)
		{FinitePointingSource(NumRay, SpreadAngle_in_radian, FocalLength, 0.0);}
	void FinitePointingSource(const size_t NumRay, const double SpreadAngle_in_radian, const double FocalLength, const double diameter);

	void getRayTrajectory(RayTrajectory& ray_traj, size_t index){ ray_traj.set(RaySet_Absolute[index].x);}
	size_t size()const{ return RaySet_Relative.size();}

};

//------------------	class for Surface		-----------------------------------
class Surface : public Optics{
public:
	//---	constructors	---
	Surface(const Surface& dest):Optics(dest){};
	Surface& operator = (const Surface& dest);
	Surface(void):Optics(){};

	bool CheckOutside(double* x);
	bool CheckOutsideRelative(double* x);
	virtual bool CrossPoint(Ray& ray)=0;	//	this virtual function must be overriden. returning true if the ray is outside.
	bool GoForward(Ray& ray);				//	ray_pos must be update in the overriding function

	virtual double getRelativeSurface(const double y, const double z)=0;
	virtual void getRelativeNormal(double* normal, double* abs_pos)=0;
			//	normal must be normalized
	void Reflect(Ray& ray);
	void Refract(Ray& ray, double n_in, double n_out);

	double R;					//	radius of the sphere.
public:
	//---	functions for imaging	---
	std::vector<std::vector<double>> getShape();
private:
	double nv[3], perp_nv[3];
};

class FlatSurface: public Surface{
public:
	FlatSurface(void):Surface(){R=0.0;};
	FlatSurface(const FlatSurface& dest):Surface(dest){};

	bool CrossPoint(Ray& ray);	//	if the cross point is in the opsite side, return true;
	double getRelativeSurface(const double y, const double z){return 0.0;}
	void getRelativeNormal(double* normal, double* abs_pos)
	{	normal[0]=1.0; normal[1]=0.0; normal[2]=0.0;}
};

class SphericalSurface: public Surface{	// positive direction : center position -> spherical center
public:
	SphericalSurface(void):Surface(){};
	SphericalSurface(const SphericalSurface& dest) :Surface(dest){ *this = dest;};
	void operator =(const SphericalSurface& dest){copyLocation(dest); R = dest.R;}

	double Xc_x[3];		//	vector defined as Xc - x
	bool CrossPoint(Ray& ray);
	double getRelativeSurface(const double y, const double z);
	void getRelativeNormal(double* normal, double* abs_pos);
};

class CylindricalSurface: public Surface{
public:
	CylindricalSurface(void): Surface(){};
	CylindricalSurface(const CylindricalSurface& dest): Surface(dest){*this = dest;};
	void operator = (const CylindricalSurface& dest){ copyLocation(dest); R = dest.R;};

	double Xc_x[3];				//	vector defined as Xc - x
	bool CrossPoint(Ray& ray);
	double getRelativeSurface(const double y, const double z);
	void getRelativeNormal(double* normal, double* abs_pos);
};

class AsphericalSurface: public Surface{
public:
	AsphericalSurface(void):Surface(), err(1.0e-8), lmax(1.0e3){};
	AsphericalSurface(const AsphericalSurface& obj);
	AsphericalSurface& operator =(const AsphericalSurface& dest);

	void setParameter(double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7);
	bool CrossPoint(Ray& ray);
	double getRelativeSurface(const double y, const double z);
	void getRelativeNormal(double* normal, double* abs_pos);
private:
	double c;		//curvature at the vertex (1/radius of curvature)    
	double k;		// K=conic constant  
	double coef[8];	//	aspheric coefficient
	const double err;	//	err for iterative calculation
	const double lmax;	//	max value of travelling length;
};

class ParabolicSurface: public Surface{	// positive direction : center position -> spherical center
public:
	ParabolicSurface(void):Surface(){};
	ParabolicSurface(const ParabolicSurface& dest) :Surface(dest){ *this = dest;};
	void operator =(const ParabolicSurface& dest);

	double Xfrel[3];	//	focal point 

	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_){setLocation(X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);}
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_);
	bool CrossPoint(Ray& ray);
	double getRelativeSurface(const double y, const double z);
	void getRelativeNormal(double* normal, double* abs_pos);
private:
	double pcoef, Xf, Yf, Zf;	//	point on the surface holds a relation of X = p*((Y-Yf)^2+(Z-Zf)^2-Yf^2-Zf^2)
};

//------------------	class for optical component		-----------------------------------
class PlaneMirror : public FlatSurface{
public:
	PlaneMirror(void):FlatSurface(){name="PlaneMirror";}
	void setParameter(double diameter_){diameter=diameter_; shape=Circle;};
	void setParameter(double size_y_, double size_z_){size_y=size_y_; size_z=size_z_; shape=Rectangular; };
	void Go(Ray& ray){Reflect(ray);};
};

class ConcaveMirror : public SphericalSurface{
public:
	ConcaveMirror(void):SphericalSurface(){name="SphericalConcaveMirror";}
	ConcaveMirror(double f, double diameter_){setParameter(f, diameter_);}
	ConcaveMirror(double f, double size_y_, double size_z_){setParameter(f, size_y_, size_z_);}
	void setParameter(double f, double diameter_)	{R=2.0*f; diameter = diameter_; shape =Circle;};
	void setParameter(double f, double size_y_, double size_z_)	{R=2.0*f; size_y = size_y_; size_z = size_z_; shape = Rectangular;};
	void Go(Ray& ray){Reflect(ray);};
};

class ParabolicMirror : public ParabolicSurface{	// focal point is in x-y (relative coordinate) plane
public:
	ParabolicMirror(void):ParabolicSurface(){name="ParabolicMirror";}
	void setParameter(double f, double focal_angle, double diameter_)	
		{Xfrel[0]=f*cos(focal_angle); Xfrel[1]=f*sin(focal_angle); Xfrel[2]=0.0; diameter=diameter_; shape=Circle;};
	void setParameter(double f, double focal_angle, double size_y_, double size_z_)	
		{Xfrel[0]=f*cos(focal_angle); Xfrel[1]=f*sin(focal_angle); Xfrel[2]=0.0; size_y = size_y_*cos(focal_angle*180.0/M_PI); size_z = size_z_; shape = Rectangular;};
	void Go(Ray& ray){Reflect(ray);};
};
class PlaneGrating : public FlatSurface{
public:
	PlaneGrating(void):FlatSurface(){name="PlaneGrating";}
	PlaneGrating& operator =(const PlaneGrating& dest){N=dest.N; Neff=dest.Neff;return *this;}
	PlaneGrating(const PlaneGrating& dest):FlatSurface(dest){ *this = dest;};

	double N;	//	number of grooves per unit length in [/mm]
	double Neff;
	void setParameter(double N_in_mm, double diameter_)	{N=N_in_mm; Neff=N_in_mm; diameter=diameter_; shape=Circle; shape=Circle;};
	void setParameter(double N_in_mm, double size_y_, double size_z_)	{N=N_in_mm; Neff=N_in_mm; size_y = size_y_; size_z = size_z_;	shape=Rectangular;};
	void setOrder(int order){Neff = N*order;};
	void Diffraction(Ray &ray);
	void Diffraction(Ray &ray, int order){ setOrder(order);Diffraction(ray);};
};

class Slit : public FlatSurface{
public:
	Slit(void):FlatSurface(){name="Slit";}
	void setParameter(double size_y_, double size_z_, enum Shape shape_){size_y=size_y_; size_z=size_z_; shape=shape_;};
	void GO(Ray& ray){	GoForward(ray);};
};

class Aparture: public FlatSurface{
public:
	Aparture(void):FlatSurface(){name="Aparture";}
	void setParameter(const double diameter_){diameter = diameter_;shape = Circle;}
	void setParameter(const double size_y_, const double size_z_){size_y = size_y_; size_z = size_z_; shape = Rectangular;}
	void Go(Ray& ray){	GoForward(ray);};
};
class Screen: public FlatSurface{
public:
	enum Direction{
		Ydirection,
		Zdirection
	};
	Screen(void):FlatSurface(){name="Screen";}
	void setParameter(double size_y_, double size_z_){size_y = size_y_; size_z = size_z_; shape = Rectangular; };
	void Go(Ray &ray);
	double GetSigma(std::vector<Ray>& RaySet, Direction direction);
	double GetFlatness(std::vector<Ray>& RaySet, Direction direction);
	double GetSigmaFromX(std::vector<Ray>& RaySet, double Yrelative, double Zrelative);
	double GetSigmaFromXOneDirection(std::vector<Ray>& RaySet, double Xrelative, Direction direction);
	double GetNumRays(const std::vector<Ray>& RaySet);
	double GetNumRaysOutsideRegion(const std::vector<Ray>& RaySet, Direction direction, double size);
	std::vector<std::vector<double>> getPosition(std::vector<Ray>& RaySet);
	void addPosition(std::vector<Ray>& RaySet, std::vector<std::vector<double>>& Image);
};

//--------------------------------	basic class for glass optics	---------------------------//
//																							   //
//	this class has two surface class. surface0 and surface1, respectively.					   //
//	every glass optics class inherits this class. (directory of indirectory) 				   //
//																							   //
//---------------------------------------------------------------------------------------------//
template < class Surface0, class Surface1 >
class GlassOptics: public Optics{
public:
//---		constructors	---
	double nair;
	GlassOptics():Optics(){		nair = 1.000292;	}
	GlassOptics(const GlassOptics& dest):Optics(dest){*this=dest;}
	GlassOptics< Surface0, Surface1 >& operator = (const GlassOptics< Surface0, Surface1 >& dest){
		surface0 = dest.surface0;		surface1 = dest.surface1;
		glass = dest.glass;				setGlass();
		return *this;
	}
	//---	glass property	---
	Glass glass;
private:
	double a0, a1, a2, a3, a4, a5, a6, a7, a8;	//	coefficient for the refractive index of the glass; with the wavelength wl in um
	double B1, B2, B3, C1, C2, C3;	//	coefficient for the refractive index of the glass; with the wavelength wl in um
	void setGlass();
	double getn(double wavelength_in_nm);
	double getnair(double wavelength_in_nm);
	virtual void setSurfaces(){};	

protected:
	//---	defining the path ---
	enum Path{
		Path_Go,							//	ray is refracted on surface0 and again refracted on surface1
		Path_Refract0_Reflect1_Refract0,	//	ray is refracted on surface0 and reflected on surface1 and again refracted on surface0
		Path_Reflect0,						//	ray is reflected on surface0
		Path_Reverse,						//	the oposite direction 
		Path_Refract1_Reflect0_Refract1,
		Path_Reflect1
	};
	Path path;
	//-------------------------

public:
	Surface0 surface0;		//	Glass Optics class has two surface classes (surface0 and surface1)
	Surface1 surface1;		//	and additional two classes (surface0_rev, surface1)
											//	they are almost same but the direction is oposit (to consider the reverse direction)

public:
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_){setLocation(X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);}
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_);	

//----	functions for ray tracing	---
	void Go(Ray& ray);
	void Go_Refract_Refract(Ray& ray);
	void Reverse(Ray& ray);
	void GoRev_Refract_Refract(Ray& ray);
	void Go_Refract_Reflect_Refract(Ray& ray);
	void GoRev_Refract_Reflect_Refract(Ray& ray);
	void Go_Reflect(Ray& ray);
	void GoRev_Reflect(Ray& ray);

//----	functions for drawing the shape	---
	void getRayTrajectory(RayTrajectory& ray_traj);

	std::vector<std::vector<double>> getShape(){
		std::vector<std::vector<double>> surf0 = surface0.getShape();	surf0.push_back(get_stdvNAN());
		std::vector<std::vector<double>> surf1 = surface1.getShape();
		for(size_t i=0; i<surf1.size(); i++) surf0.push_back(surf1[i]);	
		return surf0;
	};
};


//---------------------	Substantial GlassOptics inheriting "GlassOptics"class	---
class SphericalSingletLens : public GlassOptics<SphericalSurface, SphericalSurface>
{
public:
	SphericalSingletLens(void): GlassOptics(){name="SphericalSingletLens";}
	SphericalSingletLens(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_);

	double diameter, center_thickness, radius0, radius1;
	void setParameter(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_){
		diameter=diameter_; center_thickness=center_thickness_;	radius0=radius0_;	radius1=radius1_;	glass=glass_;
	}
	void setSurfaces();
};
class ConvexLens : public SphericalSingletLens{};

class CylindricalSingletLens: public GlassOptics< CylindricalSurface, CylindricalSurface >{
public:
	CylindricalSingletLens(void): GlassOptics(){name="CylindricalSingletLens";}
	CylindricalSingletLens(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_);
	CylindricalSingletLens(double size_y_, double size_z_, double center_thickness_, double radius0_, double radius1_, Glass glass_);

	double diameter, center_thickness, radius0, radius1;
	void setParameter(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_);
	void setParameter(double size_y_, double size_z_, double center_thickness_, double radius0_, double radius1_, Glass glass_);
	void setSurfaces();
};

class AsphericalSingletLens: public GlassOptics<AsphericalSurface, FlatSurface>{
public:
	AsphericalSingletLens(void):GlassOptics(){name="AsphericalSingletLens";}
	AsphericalSingletLens(double diameter_, double center_thickness_, double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7, Glass glass_); 
	double diameter, center_thickness;
	void setParmeter(double diameter_, double center_thickness_, double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7, Glass glass_); 
	void setSurfaces();
};

class Prism: public GlassOptics<FlatSurface, FlatSurface>
{
public:
	Prism(void):GlassOptics(){name="Prism";}
	double angle, size_y, size_z;
	void setParameter(double angle_, double size_y_, double size_z_, Glass glass_){
		angle=angle_;	size_y=size_y_; size_z=size_z_;	glass=glass_;
	}
	void setSurfaces();
};



class SphericalDoubletLens : public Optics
{
public:
	//---	constructors	---
	SphericalDoubletLens(void):Optics(){name="SphericalDoubletLens";}
	SphericalDoubletLens& operator=(const SphericalDoubletLens& dest);
	SphericalDoubletLens(const SphericalDoubletLens& dest): Optics(dest){*this=dest;}
	SphericalDoubletLens(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, 
		 Glass glass0_, Glass glass1_);
	SphericalDoubletLens(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, double radius3_,double separation,
		 Glass glass0_, Glass glass1_);

/*	SphericalDoubletLens OppositeDirectionOptics(){
		SphericalDoubletLens tmp;
		tmp.lens0 = lens1.OppositeDirectionOptics();	tmp.lens1 = lens0.OppositeDirectionOptics();	
		tmp.xcenter0 = xcenter1;	tmp.xcenter1 = xcenter0;
		return tmp;
	}
*/
	enum Path{
		Path_Go,							//	ray is refracted on surface0 and again refracted on surface1
		Path_Reverse
	};

private:
	SphericalSingletLens lens0;
	SphericalSingletLens lens1;
	Path path;
	double xcenter0, xcenter1;	//	relative location of both lens 
	
public:
	void setParameter
		(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, 
		 Glass glass0_, Glass glass1_);
	void setParameter
		(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, double radius3_,double separation,
		 Glass glass0_, Glass glass1_);

	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_){setLocation(X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);}
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_);

	//---	function for ray tracing	---
	void Go(Ray& ray);
	void Go_Refract_Refract_Refract(Ray& ray);
	void Reverse(Ray& ray);
	void GoRev_Refract_Refract_Refract(Ray& ray);



	//---	function for imaging	---
	void getRayTrajectory(RayTrajectory& ray_traj){
		if(path==Path_Go)		{	lens0.getRayTrajectory(ray_traj);	lens1.getRayTrajectory(ray_traj);}
		if(path==Path_Reverse)	{	lens1.getRayTrajectory(ray_traj);	lens0.getRayTrajectory(ray_traj);}
	}

	std::vector<std::vector<double>> getShape(){
		std::vector<std::vector<double>> imag0 = lens0.getShape();	imag0.push_back(get_stdvNAN());
		std::vector<std::vector<double>> imag1 = lens1.getShape();
		for(size_t i=0; i<imag1.size(); i++) imag0.push_back(imag1[i]);	return imag0;
	};

};


#include "Optics.inl"
#endif