#ifndef __RAYTRACING_RAYS_HPP__
#define __RAYTRACING_RAYS_HPP__

//---	class for describing one ray	---
//	this class contains the coordinate, direction, wavelength, and out flag
class Ray{
public:
	Ray& operator=(const Ray& obj){
		 x[0] = obj.x[0];	x[1] = obj.x[1];	x[2] = obj.x[2];
		dx[0] = obj.dx[0];	dx[1] = obj.dx[1];	dx[2] = obj.dx[2];
		wavelength = obj.wavelength;
		Out = obj.Out;
		return *this;
	}
	Ray(void){};
	~Ray(void){};
	Ray(const Ray& obj){	*this = obj;	}

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



#endif