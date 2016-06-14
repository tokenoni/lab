#ifndef __OPTICS_HPP__
#define __OPTICS_HPP__

#include <vector> 
#include <string>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

//---	class for describing one ray	---
//	this class contains the coordinate, direction, wavelength, and out flag
class Ray{
public:
	double x[3];	//	position of the ray at the each step
	double dx[3];	//	direction	"
	double wavelength;		//	wavelength in [nm]
	bool Out;				//	true if ray is outside the optics
};

//---	class for storing the ray trajectories	---
//	the trajectories are stored as a matrix with 3 column, x, y, z coordinate.
//	NAN is used as a separator.
class RayTrajectory{
public:
	std::vector<std::vector<double>> x;
	std::vector<double> stdv_x;
	std::vector<double> wavelength;
private:
	double wavelength_tmp;
	std::vector<double> stdvNAN;
public:
	//	constructor	
	RayTrajectory(void){
		x.resize(0);
		stdv_x.resize(3);
		stdvNAN.resize(3);
		for(size_t i=0; i<3;i++) stdvNAN[i]=sqrt(-1.0);
	}
	//	adds one trajectory
	void set(double *x_add){
		for(size_t i=0; i<3; i++) stdv_x[i] = x_add[i];
		x.push_back(stdv_x);
		wavelength.push_back(wavelength_tmp);
	}
	void addSeparator(){x.push_back(stdvNAN); wavelength.push_back(wavelength_tmp);}
	void setWavelength(double wl){wavelength_tmp = wl;}
	void clear(void){ x.clear(); wavelength.clear();}
};

//---	very basic class for every optics	---
//	this class provide several basic funstions 
//	every optics class inherits this class. (directory of indirectory) 
class Optics
{
public:
	static const double PI;
	std::string name;		//	name of the optics

	Optics(void){OpticsExplicitConstructor();}
	void OpticsExplicitConstructor(void);
	~Optics(void);

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

	double x_tmp[3];
	double ray_pos[3];	
	std::vector<double> stdv_ray_pos;
	gsl_matrix *ProjectionMatrix;	//	  projection matrix to the relative coodinate
	gsl_matrix *reProjectionMatrix;	//	reprojection matrix to the absolute coodinate

//---	physical property	---
	double diameter, size_y, size_z;	// diameter: diameter of the optics		(available if "shape" is "Circle")
	enum Shape{							//	size_y, size_z: length of each side (available if "shape" is "Recutangular")
		Circle,
		Rectangular};
	enum Shape shape;

	void getProjectionMatrix(gsl_matrix *Minv, gsl_matrix *M);
	void RelativeCoordinate(double* x);		void RelativeCorrdinate(double *xrel, const double* x){ copy(xrel, x); RelativeCoordinate(xrel);}
	void AbsoluteCoordinate(double* x);		void AbsoluteCorrdinate(double *xrel, const double* x){ copy(xrel, x); AbsoluteCoordinate(xrel);}
	void RelativeDirection(double* x);		void RelativeDirection(double *xrel, const double* x){ copy(xrel, x); RelativeDirection(xrel);}
	void AbsoluteDirection(double* x);		void AbsoluteDirection(double *xrel, const double* x){ copy(xrel, x); AbsoluteDirection(xrel);}

	void getRayTrajectory(RayTrajectory& ray_traj){ray_traj.set(ray_pos);}
	std::vector<double> getRayTrajectory(){for(size_t i=0; i<3; i++) stdv_ray_pos[i] = ray_pos[i]; return stdv_ray_pos;}

//---	glass name	---
	enum Glass{
		BK7,	BK7Nikon,	SQ,	SF11,	SF10,	BaF10,	BaFN10
	};

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

//---	overloading functions for copy	---
	template<class OpticsType_ >
	void copyLocation(const OpticsType_ &dest){
		setLocation(dest.X[0], dest.X[1], dest.X[2], dest.dX[0], dest.dX[1], dest.dX[2], dest.dXz[0], dest.dXz[1], dest.dXz[2] );
		shape = dest.shape;
		diameter=dest.diameter; size_y=dest.size_y;	size_z=dest.size_z;
		name=name;
	};
	template<class OpticsType_ >
	OpticsType_& operator = (const OpticsType_& dest){	copyLocation(dest); return *this;}
	Optics(const Optics& dest){OpticsExplicitConstructor(); copyLocation(dest);}
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
	void ParallelSource(size_t NumRay, double size_y, double size_z);
	void SpreadingSource(size_t NumRay, double SpreadAngle_in_radian, double diameter);
	void SpreadingGaussianSource(size_t NumRay, double FWHMAngle_in_radian, double MaxAngle_in_radian, double diameter);
	void SpreadingSource(size_t NumRay, double FWHMAngle_in_radian, double size_y, double size_z);
	void SpreadingGaussianSource(size_t NumRay, double FWHMAngle_in_radian, double MaxAngle_in_radian, double size_y, double size_z);
	
	void getRayTrajectory(RayTrajectory& ray_traj, size_t index){ ray_traj.set(RaySet_Absolute[index].x);}

};

//------------------	class for Surface		-----------------------------------
class Surface : public Optics{
private:
	double nv[3], perp_nv[3];
public:
	Surface(const Surface& dest):Optics(dest){};
	Surface(void):Optics(){};

	bool CheckOutside(double* x);
	virtual bool CrossPoint(Ray& ray)=0;	//	this virtual function must be overriden. returning true if the ray is outside.
	bool GoForward(Ray& ray);				//	ray_pos must be update in the overriding function

	virtual double getRelativeSurface(const double y, const double z)=0;
	virtual void getRelativeNormal(double* normal, double* abs_pos)=0;
			//	normal must be normalized
	void Reflect(Ray& ray);
	void Refract(Ray& ray, double n_in, double n_out);

	//---	functions for imaging	---
	std::vector<std::vector<double>> getShape();

};

class FlatSurface: public Surface{
public:
	FlatSurface(void):Surface(){};
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

	double R;					//	radius of the sphere.
	double Xc_x[3];		//	vector defined as Xc - x
	bool CrossPoint(Ray& ray);
	double getRelativeSurface(const double y, const double z);
	void getRelativeNormal(double* normal, double* abs_pos);
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
		{Xfrel[0]=f*cos(focal_angle); Xfrel[1]=f*sin(focal_angle); Xfrel[2]=0.0; size_y = size_y_*cos(focal_angle*180.0/PI); size_z = size_z_; shape = Rectangular;};
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
	double GetSigmaFromX(std::vector<Ray>& RaySet, double Yrelative, double Zrelative);
	std::vector<std::vector<double>> getPosition(std::vector<Ray>& RaySet);
	void addPosition(std::vector<Ray>& RaySet, std::vector<std::vector<double>>& Image);
};

//---	basic class for glass optics	---
//	this class has two surface class. surface0 and surface1, respectively.
//	every glass optics class inherits this class. (directory of indirectory) 
template < class Surface0, class Surface1 >
class GlassOptics: public Optics{
public:
	double nair;
	GlassOptics():Optics(){		nair = 1.000292;	}
	GlassOptics(const GlassOptics& dest):Optics(dest){*this=dest;}
	GlassOptics< Surface0, Surface1 >& operator = (const GlassOptics< Surface0, Surface1 >& dest){
		surface0 = dest.surface0;		surface1 = dest.surface1;
		glass = dest.glass;				setGlass();
		return *this;
	}
/*	GlassOptics< Surface1, Surface0 > OppositeDirectionOptics(){		
		GlassOptics< Surface1, Surface0 > tmp;
		tmp.surface0 = surface1;	tmp.surface1 = surface0;	tmp.glass = glass;	tmp.setGlass();
		return tmp;
	}
*/	//---	glass property	---
	Glass glass;
private:
	double a0, a1, a2, a3, a4, a5, a6, a7, a8;	//	coefficient for the refractive index of the glass; with the wavelength wl in um
	void setGlass(){
		if(glass==SF11)		{a0=3.05304325E+00;	a1=-1.27339910E-02; a2=0.00000000E+00;	a3=3.99774262E-02;	a4=3.16619134E-03;	a5=-5.02824259E-04;	a6=1.22491876E-04;	a7=-1.25325941E-05;	a8=6.19354223E-07;}
		if(glass==SF10)		{a0=2.87916509E+00;	a1=-1.19049122E-02;	a2=0.00000000E+00;	a3=3.28054585E-02;	a4=2.70047713E-03;	a5=-4.76826023E-04;	a6=1.07927203E-04;	a7=-1.07672748E-05;	a8=5.00986227E-07;}
		if(glass==BaF10)	{a0=2.72808119E+00;	a1=-9.30210914E-03;	a2=-9.30210914E-03;	a3=2.08031569E-02;	a4=4.57311835E-04;	a5=-2.96273778E-06;	a6=1.63114030E-06;	a7=0.00000000E+00;	a8=0.00000000E+00;}
		if(glass==BaFN10)	{a0=2.72876635;		a1=-0.010198898;	a2=-9.41367796e-05;	a3=0.0207671289;	a4=0.000402540351;	a5=1.39206695e-06;	a6=3.73422611e-06;	a7=-4.01939581e-07;	a8=2.09956028e-08;}
		if(glass==BK7Nikon ){a0=2.27109726E+00;	a1=-9.47304881E-03; a2=-8.91871520E-05;	a3=1.09352525E-02;	a4=1.36527555E-04;	a5=1.68617824E-06;	a6=5.85391298E-08;	a7=0.00000000E+00;	a8=0.00000000E+007;}
		if(glass==BK7)		{a0=2.27487993;		a1=-0.0104911551;   a2=-5.23238996e-05;	a3=0.00846123789;	a4=0.000701086887;	a5=-5.12306651e-05;	a6=1.64858795e-06;	a7=0.0;				a8=0.0;}
		if(glass==SQ)		{a0=2.10224128;		a1=-0.00875050668;	a2=-0.000120370139;	a3=0.0102370232;	a4=-0.000252447731;	a5=2.68476597e-05;	a6=-5.42703333e-07;	a7=0.0;				a8=0.0;}
	}
	double getn(double wavelength_in_nm){
		double l = wavelength_in_nm*1.0e-3;	//	wavelength in um;
		double l2 = l*l;	double l4 = l2*l2;
		return sqrt(a0 + a1*l2 + a2*l4 + a3/l2 + a4/l4 + a5/(l2*l4) + a6/(l4*l4) + a7/(l2*l4*l4) + a8/(l4*l4*l4));
	};
	virtual void setSurfaces(){};	

public:
	//---	defining the pass ---
	enum Pass{
		Pass_Go,							//	ray is refracted on surface0 and again refracted on surface1
		Pass_Refract0_Reflect1_Refract0,	//	ray is refracted on surface0 and reflected on surface1 and again refracted on surface0
		Pass_Reflect0,						//	ray is reflected on surface0
		Pass_Reverse,						//	the oposite direction 
		Pass_Refract1_Reflect0_Refract1,
		Pass_Reflect1,
	};
	Pass pass;
	//-------------------------

	Surface0 surface0;		//	Glass Optics class has two surface classes (surface0 and surface1)
	Surface1 surface1;		//	and additional two classes (surface0_rev, surface1)
											//	they are almost same but the direction is oposit (to consider the reverse direction)
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_){setLocation(X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);}
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_){
		setLocationCore(X_, Y_, Z_, dX_, dY_, dZ_, dXz_, dYz_, dZz_);
		setSurfaces();
		AbsoluteCoordinate(surface0.X); AbsoluteDirection(surface0.dX);
		surface0.setLocation(surface0.X[0],surface0.X[1],surface0.X[2],surface0.dX[0],surface0.dX[1],surface0.dX[2]);
		AbsoluteCoordinate(surface1.X); AbsoluteDirection(surface1.dX);
		surface1.setLocation(surface1.X[0],surface1.X[1],surface1.X[2],surface1.dX[0],surface1.dX[1],surface1.dX[2]);
		setGlass();
	}
		
//----	functions for ray tracing	---
	void Go(Ray& ray){Go_Refract_Refract(ray);}
	void Go_Refract_Refract(Ray& ray){
		surface0.Refract(ray, nair, getn(ray.wavelength));
		surface1.Refract(ray, getn(ray.wavelength), nair); 
		pass = Pass_Go;
	}
	void Reverse(Ray& ray){GoRev_Refract_Refract(ray);}
	void GoRev_Refract_Refract(Ray& ray){
		surface1.Refract(ray, nair, getn(ray.wavelength));
		surface0.Refract(ray, getn(ray.wavelength), nair); 
		pass = Pass_Reverse;
	}
	void Go_Refract_Reflect_Refract(Ray& ray){
		surface0.Refract(ray, nair, getn(ray.wavelength));	for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i];
		surface1.Reflect(ray); 
		surface0.Refract(ray, getn(ray.wavelength), nair);
		pass = Pass_Refract0_Reflect1_Refract0;
	}
	void GoRev_Refract_Reflect_Refract(Ray& ray){
		surface1.Refract(ray, nair, getn(ray.wavelength));	for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i];
		surface0.Reflect(ray); 
		surface1.Refract(ray, getn(ray.wavelength), nair);
		pass = Pass_Refract1_Reflect0_Refract1;
	}
	void Go_Reflect(Ray& ray){
		surface0.Reflect(ray); 
		pass = Pass_Reflect0;
	}
	void GoRev_Reflect(Ray& ray){
		surface1.Reflect(ray); 
		pass = Pass_Reflect1;
	}
//----	functions for drawing the shape	---
	void getRayTrajectory(RayTrajectory& ray_traj){
		if(pass == Pass_Go							){	surface0.getRayTrajectory(ray_traj); surface1.getRayTrajectory(ray_traj);}
		if(pass == Pass_Reverse						){	surface1.getRayTrajectory(ray_traj); surface0.getRayTrajectory(ray_traj);}
		if(pass == Pass_Refract0_Reflect1_Refract0	){	ray_traj.set(ray_pos); surface1.getRayTrajectory(ray_traj); surface0.getRayTrajectory(ray_traj);}
		if(pass == Pass_Refract1_Reflect0_Refract1	){	ray_traj.set(ray_pos); surface0.getRayTrajectory(ray_traj); surface1.getRayTrajectory(ray_traj);}
		if(pass == Pass_Reflect0					){	surface0.getRayTrajectory(ray_traj);}
		if(pass == Pass_Reflect1					){	surface1.getRayTrajectory(ray_traj);}
	}
	std::vector<std::vector<double>> getShape(){
		std::vector<std::vector<double>> surf0 = surface0.getShape();	surf0.push_back(get_stdvNAN());
		std::vector<std::vector<double>> surf1 = surface1.getShape();
		for(size_t i=0; i<surf1.size(); i++) surf0.push_back(surf1[i]);	return surf0;};
	
};

class SphericalSingletLens : public GlassOptics<SphericalSurface, SphericalSurface>
{
public:
	SphericalSingletLens(void): GlassOptics(){name="SphericalSingletLens";}

	double diameter, center_thickness, radius0, radius1;
	void setParameter(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_){
		diameter=diameter_; center_thickness=center_thickness_;	radius0=radius0_;	radius1=radius1_;	glass=glass_;
	}
	void setSurfaces();
};
class ConvexLens : public SphericalSingletLens{};

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
	SphericalDoubletLens(void):Optics(){name="SphericalDoubletLens";}
	SphericalDoubletLens(const SphericalDoubletLens& dest): Optics(dest){*this=dest;}
	void operator=(const SphericalDoubletLens& dest)
	{	lens0 = dest.lens0;	lens1 = dest.lens1;	xcenter0 = dest.xcenter0;	xcenter1 = dest.xcenter1;}

/*	SphericalDoubletLens OppositeDirectionOptics(){
		SphericalDoubletLens tmp;
		tmp.lens0 = lens1.OppositeDirectionOptics();	tmp.lens1 = lens0.OppositeDirectionOptics();	
		tmp.xcenter0 = xcenter1;	tmp.xcenter1 = xcenter0;
		return tmp;
	}
*/
	SphericalSingletLens lens0;
	SphericalSingletLens lens1;

	enum Pass{
		Pass_Go,							//	ray is refracted on surface0 and again refracted on surface1
		Pass_Reverse
	};
	Pass pass;
	double xcenter0, xcenter1;	//	relative location of both lens 
	
public:
	void setParameter
		(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, 
		 Glass glass0_, Glass glass1_){
			lens0.setParameter(diameter_, center_thickness0_, radius0_, radius1_, glass0_);
			xcenter0 =  -0.5*fabs(center_thickness1_) -0.5e-6;
			lens1.setParameter(diameter_, center_thickness1_, radius1_, radius2_, glass1_);
			xcenter1 =   0.5*fabs(center_thickness0_) +0.5e-6;
	}

	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_){setLocation(X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);}
	void setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_);

	//---	function for ray tracing	---
	void Go(Ray& ray){Go_Refract_Refract_Refract(ray);}
	void Go_Refract_Refract_Refract(Ray& ray){
		lens0.Go_Refract_Refract(ray);
		lens1.Go_Refract_Refract(ray);
		pass=Pass_Go;
	}
	void Reverse(Ray& ray){GoRev_Refract_Refract_Refract(ray);}
	void GoRev_Refract_Refract_Refract(Ray& ray){
		lens1.GoRev_Refract_Refract(ray);
		lens0.GoRev_Refract_Refract(ray);
		pass=Pass_Reverse;
	}

	//---	function for imaging	---
	void getRayTrajectory(RayTrajectory& ray_traj){
		if(pass==Pass_Go)		{	lens0.getRayTrajectory(ray_traj);	lens1.getRayTrajectory(ray_traj);}
		if(pass==Pass_Reverse)	{	lens1.getRayTrajectory(ray_traj);	lens0.getRayTrajectory(ray_traj);}
	}

	std::vector<std::vector<double>> getShape(){
		std::vector<std::vector<double>> imag0 = lens0.getShape();	imag0.push_back(get_stdvNAN());
		std::vector<std::vector<double>> imag1 = lens1.getShape();
		for(size_t i=0; i<imag1.size(); i++) imag0.push_back(imag1[i]);	return imag0;};

};

#include "Optics.inl"
#endif