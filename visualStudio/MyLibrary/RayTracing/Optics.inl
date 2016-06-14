//------------------	Class Optics	---------------------

//------	constructors and related functions	-----
inline void Optics::OpticsExplicitConstructor(void){
	stdv_ray_pos.resize(3);
	stdvNAN = get_stdvNAN();
	copy(pNAN, stdvNAN);
	ProjectionMatrix   = gsl_matrix_calloc(3,3);
	reProjectionMatrix = gsl_matrix_calloc(3,3);
	X[0]=0.0;	X[1]=0.0;	X[2]=0.0;	
	dX[0]=1.0;	dX[1]=0.0;	dX[2]=0.0;	
	dXz[0]=1.0;	dXz[1]=0.0;	dXz[2]=1.0;	
	diameter=0.0; size_y=0.0; size_z=0.0;
};

inline Optics::~Optics(void){
	gsl_matrix_free(ProjectionMatrix);
	gsl_matrix_free(reProjectionMatrix);
};

inline void Optics::setLocation(double X_, double Y_, double Z_, 
						  double dX_, double dY_, double dZ_){
	setLocationCore( X_, Y_, Z_, dX_, dY_, dZ_, 0.0, 0.0, 1.0);
}
inline void Optics::setLocation( double X_, double Y_, double Z_, 
						  double dX_, double dY_, double dZ_,
						  double dXz_, double dYz_, double dZz_){
	setLocationCore( X_, Y_, Z_, dX_, dY_, dZ_, dXz_, dYz_, dZz_);
}
inline void Optics::setLocationCore( double X_, double Y_, double Z_, 
						  double dX_, double dY_, double dZ_,
						  double dXz_, double dYz_, double dZz_){
	X[0] = X_; X[1] = Y_; X[2] = Z_;
	dX[0] = dX_; dX[1] = dY_; dX[2] = dZ_;
	normalize(dX);
	
//---	get the Absolute coordinate of the "Relative axis"	---
	dXz[0] = dXz_; dXz[1] = dYz_, dXz[2] = dZz_;
	double b = dot(dX, dXz);
	for(size_t i=0; i<3; i++) dXz[i] -= b*dX[i];
	normalize(dXz);

	double dXy[3];
	dXy[0] = dXz[1]*dX[2] - dX[1]*dXz[2];
	dXy[1] = dXz[2]*dX[0] - dX[2]*dXz[0];
	dXy[2] = dXz[0]*dX[1] - dX[0]*dXz[1];
	for(size_t i=0; i<3; i++) gsl_matrix_set(reProjectionMatrix, 0, i, dX[i]);
	for(size_t i=0; i<3; i++) gsl_matrix_set(reProjectionMatrix, 1, i, dXy[i]);
	for(size_t i=0; i<3; i++) gsl_matrix_set(reProjectionMatrix, 2, i, dXz[i]);
	getProjectionMatrix(ProjectionMatrix, reProjectionMatrix);
}

inline void Optics::getProjectionMatrix(gsl_matrix *Minv, gsl_matrix *M){
	gsl_permutation *p = gsl_permutation_alloc(3);	//	calculating an inverse matrix
	gsl_matrix *LU = gsl_matrix_calloc(3,3);
	gsl_matrix_memcpy(LU, M);
	int signum;
	gsl_linalg_LU_decomp(LU, p, &signum);
	gsl_linalg_LU_invert(LU, p, Minv);
	gsl_matrix_free(LU);
	gsl_permutation_free(p);
};

//---	overloading functions for copy	---
template<class OpticsType_ >
inline void Optics::copyLocation(const OpticsType_ &dest){
	setLocation(dest.X[0], dest.X[1], dest.X[2], dest.dX[0], dest.dX[1], dest.dX[2], dest.dXz[0], dest.dXz[1], dest.dXz[2] );
	shape = dest.shape;
	diameter=dest.diameter; size_y=dest.size_y;	size_z=dest.size_z;
	name=dest.name;
};

inline void Optics::RelativeCoordinate(double* x){
	for(size_t i=0; i<3; i++){
		x_tmp[i] =0.0;
		for(size_t j=0;j<3; j++)
//			x_tmp[i] += (x[j]-X[j])*gsl_matrix_get(ProjectionMatrix, j, i);
			x_tmp[i] += (x[j]-X[j]) * ProjectionMatrix->data[j * ProjectionMatrix->tda + i];
	}
	copy(x, x_tmp);
};

inline void Optics::AbsoluteCoordinate(double* x){
	for(size_t i=0; i<3; i++){
		x_tmp[i] =X[i];
		for(size_t j=0;j<3; j++)
//			x_tmp[i] += x[j] * gsl_matrix_get(reProjectionMatrix, j, i);
			x_tmp[i] += x[j] * reProjectionMatrix->data[j * reProjectionMatrix->tda + i];
	}
	copy(x, x_tmp);
};

inline void Optics::RelativeDirection(double* x){
	for(size_t i=0; i<3; i++){
		x_tmp[i] =0.0;
		for(size_t j=0;j<3; j++)
//			x_tmp[i] += x[j] * gsl_matrix_get(ProjectionMatrix, j, i);
			x_tmp[i] += x[j] * ProjectionMatrix->data[j * ProjectionMatrix->tda + i];
	}
	copy(x, x_tmp);
};

inline void Optics::AbsoluteDirection(double* x){
	for(size_t i=0; i<3; i++){
		x_tmp[i] = 0.0;
		for(size_t j=0;j<3; j++)
//			x_tmp[i] += x[j] * gsl_matrix_get(reProjectionMatrix, j, i);
			x_tmp[i] += x[j] * reProjectionMatrix->data[j * reProjectionMatrix->tda + i];
	}
	copy(x, x_tmp);
};

//------------------	Creating Rays		-----------------------------------
inline std::vector<Ray>& CreateRay::Go(){
	RaySet_Absolute.resize(RaySet_Relative.size());

	for(size_t i=0; i<RaySet_Relative.size();i++){
		RaySet_Absolute[i] = RaySet_Relative[i];
		AbsoluteCoordinate(RaySet_Absolute[i].x);
		AbsoluteDirection(RaySet_Absolute[i].dx);
	}
	return RaySet_Absolute;
}


inline void CreateRay::SpreadingSource(const size_t NumRay, const double SpreadAngle, const double diameter){
	RaySet_Relative.resize(NumRay);
	double theta_max = 0.5*SpreadAngle;
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		do{
			RaySet_Relative[i].x[1]=diameter*(0.5-1.0*rand()/RAND_MAX);
			RaySet_Relative[i].x[2]=diameter*(0.5-1.0*rand()/RAND_MAX);
		}while(RaySet_Relative[i].x[1]*RaySet_Relative[i].x[1]+RaySet_Relative[i].x[2]*RaySet_Relative[i].x[2] > 0.25*diameter*diameter);
		double x     = 1.0*rand()/RAND_MAX;
		double theta = sqrt(x*theta_max*theta_max);
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}

inline void CreateRay::SpreadingSource(const size_t NumRay, const double SpreadAngle, const double size_y, const double size_z){
	RaySet_Relative.resize(NumRay);
	double theta_max = 0.5*SpreadAngle;
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		RaySet_Relative[i].x[1]=size_y*(0.5-1.0*rand()/RAND_MAX);
		RaySet_Relative[i].x[2]=size_z*(0.5-1.0*rand()/RAND_MAX);
		double x     = 1.0*rand()/RAND_MAX;
		double theta = sqrt(x*theta_max*theta_max);
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}

inline void CreateRay::FinitePointingSource(const size_t NumRay, const double SpreadAngle, const double FocalLength, const double diameter){
	SpreadingSource(NumRay, SpreadAngle, diameter);
	for(size_t i=0; i<NumRay; i++){
		for(size_t j=0; j<3; ++j)
			RaySet_Relative[i].x[j] -= RaySet_Relative[i].dx[j] * FocalLength;
		RaySet_Relative[i].x[0] += FocalLength;
	}
}


/*
inline void CreateRay::SpreadingGaussianSource(size_t NumRay, double FWHMAngle_in_radian, double MaxAngle_in_radian, double diameter)
{
	RaySet_Relative.resize(NumRay);
	double SigmaAngle = FWHMAngle_in_radian /(2.0*sqrt(2.0*log(2.0)));
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		do{
			RaySet_Relative[i].x[1]=diameter*(0.5-1.0*rand()/RAND_MAX);
			RaySet_Relative[i].x[2]=diameter*(0.5-1.0*rand()/RAND_MAX);
		}while(RaySet_Relative[i].x[1]*RaySet_Relative[i].x[1]+RaySet_Relative[i].x[2]*RaySet_Relative[i].x[2] > 0.25*diameter*diameter);
		double theta;
		do{
			double x = 
			theta = sqrt(-2.0*log(1.0*rand()/RAND_MAX))*sin(2.0*M_PI*rand()/RAND_MAX)*SigmaAngle;
		}while(!(fabs(theta) < fabs(0.5*MaxAngle_in_radian)));
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}

inline void CreateRay::SpreadingGaussianSource(size_t NumRay, double FWHMAngle_in_radian, double MaxAngle_in_radian, double size_y, double size_x)
{
	RaySet_Relative.resize(NumRay);
	double SigmaAngle = FWHMAngle_in_radian /(2.0*sqrt(2.0*log(2.0)));
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		RaySet_Relative[i].x[1]=diameter*(0.5-1.0*rand()/RAND_MAX);
		RaySet_Relative[i].x[2]=diameter*(0.5-1.0*rand()/RAND_MAX);
		double theta;
		do{
			theta = sqrt(-2.0*log(1.0*rand()/RAND_MAX))*sin(2.0*M_PI*rand()/RAND_MAX)*SigmaAngle;
		}while(!(fabs(theta) < fabs(0.5*MaxAngle_in_radian)));
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}
*/
//----------	surface class	---


inline bool Surface::CheckOutside(double* x){
	bool Outflag = false;
	if(shape == Circle){
		double distance = 0.0;
		for(size_t i=0; i<3; i++) distance += (x[i]-X[i])*(x[i]-X[i]);
		if(sqrt(distance) > 0.5*diameter) Outflag = true;
	}
	else if(shape == Rectangular){
		RelativeCoordinate(x);
		if(fabs(x[1])>0.5*fabs(size_y) || fabs(x[2])> 0.5*fabs(size_z))
			Outflag = true;
		else if(R != 0.0 && fabs(x[0]) > fabs(R)) Outflag = true;
		AbsoluteCoordinate(x);			
	}
	return Outflag;
}

//----------	surface class	---
inline Surface& Surface::operator=(const Surface& obj){
	Optics::copyLocation(obj); 
	R = obj.R;
	shape = obj.shape;
	return *this;
}
inline bool Surface::CheckOutsideRelative(double* x){
	bool Outflag = false;
	if(shape == Circle){
		if(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] > 0.25*diameter*diameter) 
			Outflag = true;
	}
	else if(shape == Rectangular){
		if(fabs(x[1])>0.5*fabs(size_y) || fabs(x[2])> 0.5*fabs(size_z)){
			Outflag = true;
		}else if(R != 0.0 && fabs(x[0]) > fabs(R)) Outflag = true;
	}
	return Outflag;
}


inline bool Surface::GoForward(Ray &ray){
	if(ray.Out)			 {copy(ray_pos, pNAN);	return true;}	//	already out of optics
	if(CrossPoint(ray))	 {
		copy(ray_pos, pNAN);	ray.Out = true; return true;}	//	cross point is in opsite side
	if(CheckOutside(ray.x)) {ray.Out=true;	return true;}	//	cross point is outside the optics
	return false;
}

inline void Surface::Reflect(Ray& ray){
	if( GoForward(ray) ) return;
	getRelativeNormal(nv, ray.x);
	AbsoluteDirection(nv);
	double l = -2.0*dot(ray.dx, nv);
	for(size_t i=0; i<3; i++) ray.dx[i] = ray.dx[i] + l*nv[i];
}

inline void Surface::Refract(Ray& ray, double n_in, double n_out){
	if( GoForward(ray) ) return;
	getRelativeNormal(nv, ray.x);
	AbsoluteDirection(nv);
	normalize(ray.dx);

	double coef_para = dot(ray.dx, nv);		//	component pararell to the normal

	for(size_t i=0; i<3; i++) perp_nv[i] = ray.dx[i] - coef_para*nv[i];
	double coef_perp = sqrt(dot(perp_nv, perp_nv));	//			perpendiculer to the normal

	for(size_t i=0; i<3; i++) perp_nv[i] = n_in/n_out*perp_nv[i];
	coef_perp = sqrt(dot(perp_nv, perp_nv));
//---	complete reflection condition	---
	if(fabs(coef_perp) >= 1.0)	{ray.Out= true; return;}
	if(coef_para>=0.0 ) coef_para=sqrt(1.0-coef_perp*coef_perp);
	else coef_para = -sqrt(1.0-coef_perp*coef_perp);
	for(size_t i=0; i<3; i++) ray.dx[i] = perp_nv[i]+coef_para*nv[i];
}

inline std::vector<std::vector<double>> Surface::getShape(){
	std::vector<std::vector<double>> ShapeXYZ(0);
	double Shapetmp[3];
	ShapeXYZ.push_back(get_stdv(X));	//	center position of the surface
	ShapeXYZ.push_back(stdvNAN);

	if(shape == Circle){
		double dia = diameter;
		for(size_t j=0; j<3; j++){
			for(size_t i=0; i<13; i++){
				double phi = 2.0*M_PI/12.0*i;
				Shapetmp[1] = 0.5*dia*cos(phi);	Shapetmp[2] = 0.5*dia*sin(phi);	
				Shapetmp[0] = getRelativeSurface(Shapetmp[1], Shapetmp[2]);
				AbsoluteCoordinate(Shapetmp);
				ShapeXYZ.push_back(get_stdv(Shapetmp));
			}
			dia *= 0.66;
			ShapeXYZ.push_back(stdvNAN);
		}
	}
	else if(shape == Rectangular){
		double dia_y = sqrt(2.0)*size_y;
		double dia_z = sqrt(2.0)*size_z;
		for(size_t j=0; j<3; j++){
			for(size_t i=0; i<5; i++){
				double phi = 2.0*M_PI/4.0*(0.5+i);
				Shapetmp[1] = 0.5*dia_y*cos(phi);	Shapetmp[2] = 0.5*dia_z*sin(phi);	
				Shapetmp[0] = getRelativeSurface(Shapetmp[1], Shapetmp[2]);
				AbsoluteCoordinate(Shapetmp);
				ShapeXYZ.push_back(get_stdv(Shapetmp));
			}
			dia_y *= 0.66;
			dia_z *= 0.66;
			ShapeXYZ.push_back(stdvNAN);
		}
	}
	return ShapeXYZ;
};

//------------------	FlatSurface		-----------------------------------
inline bool FlatSurface::CrossPoint(Ray& ray){
	double l = (dot(X, dX) - dot(ray.x, dX))/dot(ray.dx, dX);
	if(l<0.0) {ray.Out= true; return true;}		//	direction is oposite
	for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
	for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i];
	return false;				//	reaching to the optics
}

//------------------	SphericalSurface		-----------------------------------
inline bool SphericalSurface::CrossPoint(Ray& ray){
//---	if the surface is flat	---
	if(R==0.0){
		double l = (dot(X, dX) - dot(ray.x, dX))/dot(ray.dx, dX);
		if(l<0.0) {ray.Out= true; return true;}		//	direction is oposite
		for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
		for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i];
		return false;				//	reaching to the optics
	}
//---	surface is spherical	---
	double dx2 = dot(ray.dx, ray.dx);
	for(size_t i=0; i<3; i++) Xc_x[i] =X[i] + R*dX[i]-ray.x[i];
	double dxXc_x = dot(ray.dx,Xc_x);

	double l=-1.0;
	//---	Cross point is on the concave surface	---
	l = (dxXc_x+sqrt((dxXc_x*dxXc_x) - (dot(Xc_x, Xc_x) - R*R)*dx2))/dx2;
	for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i] + l*ray.dx[i];
	
	if(l>=0.0 && !CheckOutside(ray_pos)){	//---	Cross point is on the concave surface	---
		for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
		return false;
	}
	else{ // on the convex surface
		l = (dxXc_x-sqrt((dxXc_x*dxXc_x) - (dot(Xc_x, Xc_x) - R*R)*dx2))/dx2;
		if(l>=0.0){
			for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
			for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i];
			return false;				//	reaching the optics
		}
		else return true;		//	no solution exist!
	}	
	return false;
}

inline double SphericalSurface::getRelativeSurface(const double y, const double z){
	if(R>0.0)		return R-sqrt(R*R-y*y-z*z);
	else if(R<0.0)	return R+sqrt(R*R-y*y-z*z);
	//--- if the surface is flat	---
	else return 0.0;
}

inline void SphericalSurface::getRelativeNormal(double* normal, double* abs_pos){
	//--- if the surface is flat	---
	if(R==0.0){		normal[0]=1.0;	normal[1]=0.0;	normal[2]=0.0;	}
	//---	surface is spherical	---
	else{
		RelativeCoordinate(abs_pos);
		normal[0] = R - abs_pos[0];
		normal[1] =   - abs_pos[1];
		normal[2] =   - abs_pos[2];
		normalize(normal);
		AbsoluteCoordinate(abs_pos);
	}
}


//------------------	CylindricalSurface		-----------------------------------
inline bool CylindricalSurface::CrossPoint(Ray& ray){
//---	if the surface is flat	---
	if(R==0.0){
		double l = (dot(X, dX) - dot(ray.x, dX))/dot(ray.dx, dX);
		if(l<0.0) {ray.Out= true; return true;}		//	direction is oposite
		for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
		copy(ray_pos, ray.x);
		return false;				//	reaching to the optics
	}

//---	surface is cylindrical	---
	RelativeCoordinate(ray.x);
	RelativeDirection(ray.dx);
	double a = (ray.dx[0]*ray.dx[0]+ ray.dx[1]*ray.dx[1]);
	double b = (ray.x[0]-R)*ray.dx[0] + ray.x[1]*ray.dx[1];
	double c =  ray.x[0]*ray.x[0] + ray.x[1]*ray.x[1] - 2.0*R*ray.x[0];
	double l = (- b + sqrt(b*b - a*c))/a;
	
	for(size_t i=0; i<3; ++i) ray_pos[i] = ray.x[i] + l*ray.dx[i];
	if(l >= 0.0 && !CheckOutsideRelative(ray_pos)){
		AbsoluteCoordinate(ray_pos);
		copy(ray.x, ray_pos);
		AbsoluteDirection(ray.dx);
		return false;
	}
	l = (- b - sqrt(b*b - a*c))/a;
	for(size_t i=0; i<3; ++i) ray_pos[i] = ray.x[i] + l*ray.dx[i];
	if(l >= 0.0 && !CheckOutsideRelative(ray_pos)){
		AbsoluteCoordinate(ray_pos);
		copy(ray.x, ray_pos);
		AbsoluteDirection(ray.dx);
		return false;
	}
	else{
		AbsoluteCoordinate(ray.x);
		copy(ray_pos, ray.x);
		AbsoluteDirection(ray.dx);
		return true;
	}
}

inline double CylindricalSurface::getRelativeSurface(const double y, const double z){
	if(R>0.0)		return R-sqrt(R*R-y*y);
	else if(R<0.0)	return R+sqrt(R*R-y*y);
	//--- if the surface is flat	---
	else return 0.0;
}

inline void CylindricalSurface::getRelativeNormal(double* normal, double* abs_pos){
	//--- if the surface is flat	---
	if(R==0.0){		normal[0]=1.0;	normal[1]=0.0;	normal[2]=0.0;	}
	else{
		RelativeCoordinate(abs_pos);
		normal[0] = R - abs_pos[0];	
		normal[1] =   - abs_pos[1];	
		normal[2] =   -  0.0;
		normalize(normal);
		AbsoluteCoordinate(abs_pos);
	}
}

//------------------	AsphericalSurface		-----------------------------------
inline AsphericalSurface::AsphericalSurface(const AsphericalSurface& obj): Surface(obj), err(1.0e-8), lmax(1.0e3){
	for(size_t i=0; i<8; ++i) coef[i] = obj.coef[i];
	c=obj.c;	k=obj.k;
}
inline void AsphericalSurface::setParameter(double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7){
	c=c_;	k=k_;
	coef[0] = coef0;	coef[1]=coef1;	coef[2]=coef2;	coef[3]=coef3;
	coef[4] = coef4;	coef[5]=coef5;	coef[6]=coef6;	coef[7]=coef7;
}

inline AsphericalSurface& AsphericalSurface::operator =(const AsphericalSurface& obj){
	copyLocation(obj);
	for(size_t i=0; i<8; ++i) coef[i] = obj.coef[i];
	c=obj.c;	k=obj.k;
	shape = obj.shape;
	return *this;
}

inline double AsphericalSurface::getRelativeSurface(const double y, const double z){
	double r2 = y*y + z*z;

	return c*r2/(1.0+sqrt(1.0-(1.0+k)*c*c*r2)) 
		+ coef[0]*r2 + coef[1]*r2*r2 + coef[2]*r2*r2*r2 + coef[3]*r2*r2*r2*r2 + coef[4]*r2*r2*r2*r2*r2
		+ coef[5]*r2*r2*r2*r2*r2*r2 + coef[6]*r2*r2*r2*r2*r2*r2*r2 + coef[7]*r2*r2*r2*r2*r2*r2*r2*r2;

}
inline void AsphericalSurface::getRelativeNormal(double* normal, double* abs_pos){
	RelativeCoordinate(abs_pos);
	double r2 = abs_pos[1]*abs_pos[1] + abs_pos[2]*abs_pos[2];
	double r = sqrt(r2);
	double sqrtv = sqrt(1.0-(1.0+k)*c*c*r2);
	double dzdr =
		(2.0*c*r*(1.0+sqrtv) + c*r2*(1.0+k)*c*c*r/sqrtv )/((1.0+sqrtv)*(1.0+sqrtv))
		+ (	 2.0*coef[0] + 4.0*coef[1]*r2 + 6.0*coef[2]*r2*r2 + 8.0*coef[3]*r2*r2*r2 + 10.0*coef[4]*r2*r2*r2*r2
		   +12.0*coef[5]*r2*r2*r2*r2*r2 + 14.0*coef[6]*r2*r2*r2*r2*r2*r2 + 16.0*coef[7]*r2*r2*r2*r2*r2*r2*r2)*r;
	
	normal[0] = r;
	normal[1] = -abs_pos[1]*dzdr;
	normal[2] = -abs_pos[2]*dzdr;
	normalize(normal);
	AbsoluteCoordinate(abs_pos);
}

inline bool AsphericalSurface::CrossPoint(Ray& ray){
	RelativeCoordinate(ray.x);
	RelativeDirection(ray.dx);
	double dx = getRelativeSurface(ray.x[1], ray.x[2]) - ray.x[0], dx_next, dx_new;
	bool upper_flag = dx > 0.0;
	double l = 0.5*dx;
	if(!(upper_flag)) l= -0.5*dx;

	if((ray.dx[0] > 0.0) ^ upper_flag){
		double tmp = getRelativeSurface(ray.x[1], ray.x[2]) - ray.x[0];
		AbsoluteCoordinate(ray.x);
		AbsoluteDirection(ray.dx);
		copy(ray_pos, ray.x);
		return true;
	}

	//	the solution is within [dx, dx_next]
	do{
		l *= 2.0;
		dx_next = getRelativeSurface(ray.x[1] + l*ray.dx[1], ray.x[2] + l*ray.dx[2]) - (ray.x[0] + l*ray.dx[0]);
		if(fabs(l) > lmax){
			AbsoluteCoordinate(ray.x);
			AbsoluteDirection(ray.dx);
			copy(ray_pos, ray.x);
			return true;
		}
	}while(dx_next* dx>0.0);
	//	the solution is derived with intermediate-value theorem
	while(fabs(dx)>err){
		l *= fabs(dx/(dx - dx_next));
		ray.x[0] += l*ray.dx[0];
		ray.x[1] += l*ray.dx[1];
		ray.x[2] += l*ray.dx[2];
		dx_new = getRelativeSurface(ray.x[1], ray.x[2]) - ray.x[0];
		//	if the solution is in [dx, dx_new]
		if(dx * dx_new <0.0){
			l*=-1.0;
			dx_next = dx; dx = dx_new;
		}
		//	if the solution is in [dx_new, dx_next]
		else{
			l *= fabs(dx_next/dx);
			dx = dx_new;
		}
	}
	/*
	//--- iterative calculation	---
	while(fabs(dx)>err){
		bool tmp1 = (dx>=0.0);
		bool tmp2 = tmp1^upper_flag;
		dx = getRelativeSurface(ray.x[1] + l*ray.dx[1], ray.x[2] + l*ray.dx[2]) - (ray.x[0] + l*ray.dx[0]);

		if(fabs(dx) < fabs(l)) {
			ray.x[0] += l*ray.dx[0];
			ray.x[1] += l*ray.dx[1];
			ray.x[2] += l*ray.dx[2];
			l = dx;
			if(dx>0.0) upper_flag=true;
			else upper_flag = false;
		}
		else if(!((dx>0.0) ^ upper_flag)) l*= 2.0;
		else l*=0.5;
	}
	*/
	AbsoluteCoordinate(ray.x);
	AbsoluteDirection(ray.dx);
	copy(ray_pos, ray.x);
	return false;
}


//------------------	ParabolicSurface		-----------------------------------
inline void ParabolicSurface::setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_){
	setLocationCore(X_, Y_, Z_, dX_, dY_, dZ_, dXz_, dYz_, dZz_);
	
	Xf=Xfrel[0];	Yf=Xfrel[1];	Zf=Xfrel[2];
	pcoef = (-Xf+sqrt(Xf*Xf+Yf*Yf+Zf*Zf))/(2.0*(Yf*Yf+Zf*Zf));
}
inline bool ParabolicSurface::CrossPoint(Ray& ray){
	RelativeCoordinate(ray.x);
	RelativeDirection(ray.dx);
	double a = ray.dx[1]*ray.dx[1] + ray.dx[2]*ray.dx[2];
	double b = 2.0*(ray.x[1]*ray.dx[1] - Yf*ray.dx[1] + ray.x[2]*ray.dx[2] - Zf*ray.dx[2]) - ray.dx[0]/pcoef;
	double c = (ray.x[1]-2.0*Yf)*ray.x[1] + (ray.x[2]-2.0*Zf)*ray.x[2] - ray.x[0]/pcoef;

	//---	Cross point is on the concave surface	---
	double l = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);
	for(size_t i=0; i<3; i++) ray_pos[i] = ray.x[i] + l*ray.dx[i];
	AbsoluteCoordinate(ray_pos);
	if(l>=0.0 && !CheckOutside(ray_pos)){	//---	Cross point is on the concave surface	---
		copy(ray.x, ray_pos);
		AbsoluteDirection(ray.dx);
		return false;
	}
	else{								//---	Cross point is on the convex surface	---
		l = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
		if(l<0.0) {ray.Out= true; return true;}		//	direction is oposite
		for(size_t i=0; i<3; i++) ray.x[i] = ray.x[i] + l*ray.dx[i];
		AbsoluteCoordinate(ray.x);
		AbsoluteDirection(ray.dx);
		copy(ray_pos, ray.x);
		return false;				//	reaching the optics
	}
}

inline double ParabolicSurface::getRelativeSurface(const double y, const double z){
	return pcoef*((y-Yf)*(y-Yf) + (z-Zf)*(z-Zf) -Yf*Yf-Zf*Zf);
}
inline void ParabolicSurface::getRelativeNormal(double* normal, double* abs_pos){
	RelativeCoordinate(abs_pos);
	double dXdy = 2.0*pcoef*(abs_pos[1]-Yf);
	double dXdz = 2.0*pcoef*(abs_pos[2]-Zf);
	normal[0]=1.0;	normal[1] = -dXdy;	normal[2] = -dXdz;
	normalize(normal);
	AbsoluteCoordinate(abs_pos);
}

inline void ParabolicSurface::operator=(const ParabolicSurface& dest){
	setLocation(dest.X[0], dest.X[1], dest.X[2], dest.dX[0], dest.dX[1], dest.dX[2], dest.dXz[0], dest.dXz[1], dest.dXz[2]);
	Xfrel[0] = dest.Xfrel[0];	Xfrel[1] = dest.Xfrel[1];	Xfrel[2] = dest.Xfrel[2];	
}

//------------------		class Plane Grating	: Optics		--------------------
inline void PlaneGrating::Diffraction(Ray &ray){
	if( GoForward(ray) ) return;
	RelativeDirection(ray.dx);
	double epsilon = atan(ray.dx[2]/sqrt(ray.dx[0]*ray.dx[0]+ray.dx[1]*ray.dx[1]));
	double alpha   = atan(ray.dx[1]/ray.dx[0]);
	double beta    = asin(Neff*1.0e-6*ray.wavelength/cos(epsilon) - sin(alpha));
	if(!(-0.5*M_PI<=beta && beta<=0.5*M_PI)){
		ray.Out=true;return;}
	ray.dx[0] = cos(epsilon)*cos(beta);
	ray.dx[1] = cos(epsilon)*sin(beta);
	ray.dx[2] = sin(epsilon);
	AbsoluteDirection(ray.dx);
}

inline void Screen::Go(Ray& ray){
	GoForward(ray);
}

inline double Screen::GetSigma(std::vector<Ray>& RaySet, Direction direction){
	size_t direction_index;
	if     (direction == Ydirection)	direction_index = 1;
	else if(direction == Zdirection)	direction_index = 2;
	else return stdvNAN[0];

	double Avg=0.0;
	double Sigma=0.0;
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			num++;
			RelativeCoordinate(RaySet[i].x);
			Avg+=RaySet[i].x[direction_index];
		}
	}
	Avg /= num;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			Sigma+=(RaySet[i].x[direction_index] - Avg)*(RaySet[i].x[direction_index] - Avg);
			AbsoluteCoordinate(RaySet[i].x);
		}
	}
	return sqrt(Sigma/num);
}

inline double Screen::GetNumRays(const std::vector<Ray>& RaySet){
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++)
		if(!RaySet[i].Out) ++num;
	return num;
}

inline double Screen::GetFlatness(std::vector<Ray>& RaySet, Direction direction){
	size_t direction_index;
	if     (direction == Ydirection)	direction_index = 1;
	else if(direction == Zdirection)	direction_index = 2;
	else return stdvNAN[0];

	double Avg=0.0;
	double Sigma=0.0;
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			num++;
			RelativeCoordinate(RaySet[i].x);
			Avg+=RaySet[i].x[direction_index];
		}
	}
	Avg /= num;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			Sigma+=(RaySet[i].x[direction_index] - Avg)*(RaySet[i].x[direction_index] - Avg)
				  *(RaySet[i].x[direction_index] - Avg)*(RaySet[i].x[direction_index] - Avg);
			AbsoluteCoordinate(RaySet[i].x);
		}
	}
	return sqrt(sqrt(Sigma/num));
}

inline double Screen::GetNumRaysOutsideRegion(const std::vector<Ray>& RaySet, Direction direction, double size){
	size_t direction_index;
	if     (direction == Ydirection)	direction_index = 1;
	else if(direction == Zdirection)	direction_index = 2;
	else return stdvNAN[0];

	double Avg=0.0;
	double Sigma=0.0;
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++){
		if((RaySet[i].Out) ||
			(RaySet[i].x[direction_index] < Avg-0.5*size || Avg+0.5*size < RaySet[i].x[direction_index]))
				++num;
	}
	return 1.0*num;
}
inline double Screen::GetSigmaFromXOneDirection(std::vector<Ray>& RaySet, double Xrelative, Direction direction){
	size_t direction_index;
	if     (direction == Ydirection)	direction_index = 1;
	else if(direction == Zdirection)	direction_index = 2;
	else return stdvNAN[0];
	double Sigma=0.0;
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			RelativeCoordinate(RaySet[i].x);
			Sigma+=(RaySet[i].x[direction_index] - Xrelative)*(RaySet[i].x[direction_index] - Xrelative);
			AbsoluteCoordinate(RaySet[i].x);
			num++;
		}
	}
	return sqrt(Sigma/num);
}

inline double Screen::GetSigmaFromX(std::vector<Ray>& RaySet, double Yrelative, double Zrelative){
	double Sigma=0.0;
	size_t num=0;
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			RelativeCoordinate(RaySet[i].x);
			Sigma+=(RaySet[i].x[1] - Yrelative)*(RaySet[i].x[1] - Yrelative);
			Sigma+=(RaySet[i].x[2] - Zrelative)*(RaySet[i].x[2] - Zrelative);
			AbsoluteCoordinate(RaySet[i].x);
			num++;
		}
	}
	return sqrt(Sigma/num);
}

inline std::vector<std::vector<double>> Screen::getPosition(std::vector<Ray>& RaySet){
	std::vector<std::vector<double>> Position(0);
	for(size_t i=0; i<RaySet.size(); i++){
		if(!RaySet[i].Out){
			RelativeCoordinate(RaySet[i].x);
			copy(stdv_ray_pos, RaySet[i].x);
			Position.push_back(stdv_ray_pos);
			AbsoluteCoordinate(RaySet[i].x);
		}
	}
	if(Position.size()==0) Position.push_back(stdvNAN);
	return Position;
};

inline void Screen::addPosition(std::vector<Ray>& RaySet, std::vector<std::vector<double>>& Image){
	std::vector<std::vector<double>> adding_image = getPosition(RaySet);
	for(size_t i=0; i<adding_image.size(); i++)	Image.push_back(adding_image[i]);
}

//---------	Class for Glass Optics	------------
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::setGlass(){
/*
	Sellmeier dispersion formula
	n^2(l) –1 = B1 l^2/ (l^2– C1) + B2 l^2/ (l^2– C2) + B3 l^2 / (l^2– C3)
*/
	if(glass==SF10)		{B1=1.62153902;    B2=0.256287842;    B3=1.64447552;    C1=0.0122241457;  C2=0.0595736775;  C3=147.468793;}
	if(glass==SF11)		{B1=1.73759695;	   B2=0.313747346;	  B3=1.89878101;    C1=0.013188707;	  C2=0.0623068142;  C3=155.23629;}
	if(glass==BaF10)    {B1=1.5851495;	   B2=0.143559385;    B3=1.08521269;    C1=0.00926681282; C2=0.0424489805;  C3=105.613573;}
	if(glass==BaFN10)   {B1=0.0; B2=0.0; B3=0.0; C1=0.0;C2=0.0; C3=0.0;a0=2.72876635;		a1=-0.010198898;	a2=-9.41367796e-05;	a3=0.0207671289;	a4=0.000402540351;	a5=1.39206695e-06;	a6=3.73422611e-06;	a7=-4.01939581e-07;	a8=2.09956028e-08;}
	if(glass==BK7)      {B1=1.03961212;    B2=0.231792344;    B3=1.01046945;    C1=0.00600069867; C2=0.0200179144;  C3=103.560653;}
	if(glass==SQ)		{B1=0.6961663;     B2=0.4079426;      B3=0.8974794;     C1=0.0684043*0.0684043;    C2=0.0684043*0.0684043;     C3=  9.896161* 9.896161; }
	if(glass==SF5)      {B1=1.52481889;    B2=0.187085527;    B3=1.42729015;    C1=0.011254756;	  C2=0.0588995392;  C3=129.141675;}
	if(glass==SK11)     {B1=1.17963631;    B2=0.229817295;    B3=0.935789652;   C1=0.00680282081; C2=0.0219737205;  C3=101.513232;}
	if(glass==SSK8)     {B1=1.44857867;    B2=0.117965926;    B3=1.06937528;    C1=0.00869310149; C2=0.0421566593;  C3=111.300666;}
	if(glass==SK11)     {B1=1.17963631;    B2=0.229817295;    B3=0.935789652;   C1=0.00680282081; C2=0.0219737205;  C3=101.513232;}
	if(glass==BAL35)    {B1=9.41357273E-01;B2=5.46174895E-01; B3=1.16168917E+00;C1=1.40333996E-02;C2=9.06635683E-04;C3=1.14163758E+02;}
	if(glass==SF2)      {B1=1.47343127;    B2=0.163681849;    B3=1.36920899;    C1=0.0109019098;  C2=0.0585683687;  C3=127.404933;}
	if(glass==BaK1)     {B1=1.12365662;    B2=0.309276848;    B3=0.881511957;   C1=0.00644742752; C2=0.0222284402;  C3=107.297751;}
	if(glass==BaK4)     {B1=1.28834642;    B2=0.132817724;    B3=0.945395373;   C1=0.00779980626; C2=0.0315631177;  C3=105.965875;}
	if(glass==SF8)      {B1=1.55075812;    B2=0.209816918;    B3=1.46205491;    C1=0.0114338344;  C2=0.0582725652;  C3=133.24165;}

	
	
	
	/*	Nikkon type dispersion formula	
	if(glass==SF10)		{a0=2.87916509E+00;	a1=-1.19049122E-02;	a2=0.00000000E+00;	a3=3.28054585E-02;	a4=2.70047713E-03;	a5=-4.76826023E-04;	a6=1.07927203E-04;	a7=-1.07672748E-05;	a8=5.00986227E-07;}
	if(glass==SF11)		{a0=3.05304325E+00;	a1=-1.27339910E-02; a2=0.00000000E+00;	a3=3.99774262E-02;	a4=3.16619134E-03;	a5=-5.02824259E-04;	a6=1.22491876E-04;	a7=-1.25325941E-05;	a8=6.19354223E-07;}
	if(glass==BaF10)	{a0=2.72808119E+00;	a1=-9.30210914E-03;	a2=-9.30210914E-03;	a3=2.08031569E-02;	a4=4.57311835E-04;	a5=-2.96273778E-06;	a6=1.63114030E-06;	a7=0.00000000E+00;	a8=0.00000000E+00;}
	if(glass==BaFN10)	{a0=2.72876635;		a1=-0.010198898;	a2=-9.41367796e-05;	a3=0.0207671289;	a4=0.000402540351;	a5=1.39206695e-06;	a6=3.73422611e-06;	a7=-4.01939581e-07;	a8=2.09956028e-08;}
	if(glass==BK7Nikon ){a0=2.27109726E+00;	a1=-9.47304881E-03; a2=-8.91871520E-05;	a3=1.09352525E-02;	a4=1.36527555E-04;	a5=1.68617824E-06;	a6=5.85391298E-08;	a7=0.00000000E+00;	a8=0.00000000E+00;}
	if(glass==BK7)		{a0=2.27487993;		a1=-0.0104911551;   a2=-5.23238996e-05;	a3=0.00846123789;	a4=0.000701086887;	a5=-5.12306651e-05;	a6=1.64858795e-06;	a7=0.0;				a8=0.0;}
	if(glass==SQ)		{a0=2.10224128;		a1=-0.00875050668;	a2=-0.000120370139;	a3=0.0102370232;	a4=-0.000252447731;	a5=2.68476597e-05;	a6=-5.42703333e-07;	a7=0.0;				a8=0.0;}
	if(glass==SF5)      {a0=2.71072072E+00;	a1=-1.02160186E-02; a2=-9.06763794E-05;	a3=2.88337808E-02;	a4=5.57561753E-04;	a5=1.33564048E-04;	a6=-1.34358407E-05;	a7=1.19202152E-06;	a8=0.00000000E+00;}
	if(glass==SK11)     {a0=2.40941529E+00;	a1=-9.65092890E-05; a2=-9.06763794E-05;	a3=1.31170272E-02;	a4=1.53988355E-04;	a5=4.69136387E-06;	a6=-2.59660236E-08;	a7=0.0           ;	a8=0.00000000E+00;}
	if(glass==SSK8)     {a0=2.56658096E+00;	a1=-9.72847347E-03; a2=-9.45439785E-05;	a3=1.74935076E-02;	a4=3.71433240E-04;	a5=-4.00752907E-06;	a6=-1.64198401E-06;	a7=0.0           ;	a8=0.00000000E+00;}
	if(glass==BAL35)	{a0=2.4882622     ;	a1=-0.011313747   ; a2=-0.00010227986 ;	a3=0.013587683   ;	a4=0.00019272228 ;	a5=-1.8539064e-06 ;	a6=1.1705456e-06  ;	a7=-1.3593238e-07;	a8=6.0345089e-09 ;}
	if(glass==SF2)      {a0=2.63659143    ; a1=-0.0105396966  ; a2=-0.000114667615; a3=0.0262121018  ;  a4=0.000468911167;  a5=8.63713212e-05 ; a6=9.98379505e-07 ; a7=-7.166621e-07 ;  a8=8.07923897e-08;}
*/
}

template < class Surface0, class Surface1 >
inline double GlassOptics< Surface0, Surface1 >::getn(double wavelength_in_nm)
{
/*
	Sellmeier dispersion formula
	n^2(l) –1 = B1 l^2/ (l^2– C1) + B2 l^2/ (l^2– C2) + B3 l^2 / (l^2– C3)
*/
	double l = wavelength_in_nm*1.0e-3;	//	wavelength in um;
	double l2 = l*l;	
	double n2_1 = B1*l2 / (l2-C1) + B2*l2 / (l2 - C2) + B3*l2 / (l2 - C3);
	if(n2_1!=0.0)
		return sqrt(n2_1 +1.0);
	else{
		double l4 = l2*l2;
		return sqrt(a0 + a1*l2 + a2*l4 + a3/l2 + a4/l4 + a5/(l2*l4) + a6/(l4*l4) + a7/(l2*l4*l4) + a8/(l4*l4*l4));
	}
};

template < class Surface0, class Surface1 >
inline double GlassOptics< Surface0, Surface1 >::getnair(double wavelength_in_nm)
{
	return 1.000292;
};
	
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_){
	setLocationCore(X_, Y_, Z_, dX_, dY_, dZ_, dXz_, dYz_, dZz_);
	setSurfaces();
	AbsoluteCoordinate(surface0.X); AbsoluteDirection(surface0.dX);
	surface0.setLocation(surface0.X[0],surface0.X[1],surface0.X[2],surface0.dX[0],surface0.dX[1],surface0.dX[2], dXz_, dYz_, dZz_);
	AbsoluteCoordinate(surface1.X); AbsoluteDirection(surface1.dX);
	surface1.setLocation(surface1.X[0],surface1.X[1],surface1.X[2],surface1.dX[0],surface1.dX[1],surface1.dX[2], dXz_, dYz_, dZz_);
	setGlass();
}

template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::Go(Ray& ray){
	ray_tmp = ray;
	Go_Refract_Refract(ray_tmp);
	if(!ray_tmp.Out){
		ray = ray_tmp;
	}else{
		GoRev_Refract_Refract(ray);
	}
}

template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::Go_Refract_Refract(Ray& ray){
	surface0.Refract(ray, getnair(ray.wavelength), getn(ray.wavelength));
	surface1.Refract(ray, getn(ray.wavelength), getnair(ray.wavelength)); 
	path = Path_Go;
}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::Reverse(Ray& ray){GoRev_Refract_Refract(ray);}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::GoRev_Refract_Refract(Ray& ray){
	surface1.Refract(ray, getnair(ray.wavelength), getn(ray.wavelength));
	surface0.Refract(ray, getn(ray.wavelength), getnair(ray.wavelength)); 
	path = Path_Reverse;
}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::Go_Refract_Reflect_Refract(Ray& ray){
	surface0.Refract(ray, getnair(ray.wavelength), getn(ray.wavelength));	copy(ray_pos, ray.x);
	surface1.Reflect(ray); 
	surface0.Refract(ray, getn(ray.wavelength), getnair(ray.wavelength));
	path = Path_Refract0_Reflect1_Refract0;
}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::GoRev_Refract_Reflect_Refract(Ray& ray){
	surface1.Refract(ray, getnair(ray.wavelength), getn(ray.wavelength));	copy(ray_pos, ray.x);
	surface0.Reflect(ray); 
	surface1.Refract(ray, getn(ray.wavelength), getnair(ray.wavelength));
	path = Path_Refract1_Reflect0_Refract1;
}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::Go_Reflect(Ray& ray){
	surface0.Reflect(ray); 
	path = Path_Reflect0;
}
template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::GoRev_Reflect(Ray& ray){
	surface1.Reflect(ray); 
	path = Path_Reflect1;
}

template < class Surface0, class Surface1 >
inline void GlassOptics< Surface0, Surface1 >::getRayTrajectory(RayTrajectory& ray_traj){
	if(path == Path_Go							){	surface0.getRayTrajectory(ray_traj); surface1.getRayTrajectory(ray_traj);}
	if(path == Path_Reverse						){	surface1.getRayTrajectory(ray_traj); surface0.getRayTrajectory(ray_traj);}
	if(path == Path_Refract0_Reflect1_Refract0	){	ray_traj.set(ray_pos); surface1.getRayTrajectory(ray_traj); surface0.getRayTrajectory(ray_traj);}
	if(path == Path_Refract1_Reflect0_Refract1	){	ray_traj.set(ray_pos); surface0.getRayTrajectory(ray_traj); surface1.getRayTrajectory(ray_traj);}
	if(path == Path_Reflect0					){	surface0.getRayTrajectory(ray_traj);}
	if(path == Path_Reflect1					){	surface1.getRayTrajectory(ray_traj);}
}


inline void Prism::setSurfaces(){
	surface0.X[0]  = 0.5*size_z*sin(0.5*angle);		surface0.X[1] = 0.0;	surface0.X[2] = 0.0;
	surface0.dX[0] = cos(0.5*angle);				surface0.dX[1] = 0.0;	surface0.dX[2] = sin(0.5*angle);
	surface0.shape = Optics::Rectangular;			surface0.size_y = size_y;	surface0.size_z = size_z;

	surface1.X[0]  =-0.5*size_z*sin(0.5*angle);		surface1.X[1] = 0.0;	surface1.X[2] = 0.0;
	surface1.dX[0] = cos(0.5*angle);				surface1.dX[1] = 0.0;	surface1.dX[2] =-sin(0.5*angle);
	surface1.shape = Optics::Rectangular;			surface1.size_y = size_y;	surface1.size_z = size_z;
}
//---	SphericalSingletLens	---
inline SphericalSingletLens::SphericalSingletLens(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_){
	setParameter(diameter_, center_thickness_, radius0_, radius1_, glass_);
}

inline void SphericalSingletLens::setSurfaces(){
	surface0.X[0]=-0.5*fabs(center_thickness);	surface0.X[1]=0.0;	surface0.X[2]=0.0;
	surface0.dX[0]= 1.0;					surface0.dX[1]=0.0;	surface0.dX[2]=0.0;
	surface0.shape = Surface::Circle;		surface0.diameter=fabs(diameter);
	surface0.R = radius0;					//surface0.surface_type=surface0.Convex;

	surface1.X[0]= 0.5*fabs(center_thickness);	surface1.X[1]=0.0;	surface1.X[2]=0.0;
	surface1.dX[0]=1.0;						surface1.dX[1]=0.0;	surface1.dX[2]=0.0;
	surface1.shape = Surface::Circle;		surface1.diameter=fabs(diameter);
	surface1.R = radius1;					//surface1.surface_type=surface1.Concave;
}

//---	CylindricalSingletLens	---
inline CylindricalSingletLens::CylindricalSingletLens(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_){
	setParameter(diameter_, center_thickness_, radius0_, radius1_, glass_);
}
inline CylindricalSingletLens::CylindricalSingletLens(double size_y_, double size_z_, double center_thickness_, double radius0_, double radius1_, Glass glass_){
	setParameter(size_y_, size_z_, center_thickness_, radius0_, radius1_, glass_);
}

inline void CylindricalSingletLens::setParameter(double diameter_, double center_thickness_, double radius0_, double radius1_, Glass glass_)
{
	diameter=diameter_; center_thickness=center_thickness_;	radius0=radius0_;	radius1=radius1_;	glass=glass_;
	surface0.shape = Surface::Circle;		surface0.diameter=fabs(diameter);
	surface1.shape = Surface::Circle;		surface1.diameter=fabs(diameter);
}

inline void CylindricalSingletLens::setParameter(double size_y_, double size_z_, double center_thickness_, double radius0_, double radius1_, Glass glass_)
{
	size_y = size_y_; size_z = size_z_; center_thickness=center_thickness_;	radius0=radius0_;	radius1=radius1_;	glass=glass_;
	surface0.shape = Surface::Rectangular;		surface0.size_y=fabs(size_y);	surface0.size_z=fabs(size_z);
	surface1.shape = Surface::Rectangular;		surface1.size_y=fabs(size_y);	surface1.size_z=fabs(size_z);
}

inline void CylindricalSingletLens::setSurfaces(){
	surface0.X[0]=-0.5*fabs(center_thickness);	surface0.X[1]=0.0;	surface0.X[2]=0.0;
	surface0.dX[0]= 1.0;					surface0.dX[1]=0.0;	surface0.dX[2]=0.0;
	surface0.R = radius0;	

	surface1.X[0]= 0.5*fabs(center_thickness);	surface1.X[1]=0.0;	surface1.X[2]=0.0;
	surface1.dX[0]=1.0;						surface1.dX[1]=0.0;	surface1.dX[2]=0.0;
	surface1.R = radius1;	
}

//---	Aspherical singlet lens	---
inline AsphericalSingletLens::AsphericalSingletLens(double diameter_, double center_thickness_, double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7, Glass glass_){
	setParmeter(diameter_, center_thickness_, c_, k_, coef0, coef1, coef2, coef3, coef4, coef5, coef6, coef7, glass_);
}

inline void AsphericalSingletLens::setParmeter(double diameter_, double center_thickness_, double c_, double k_, double coef0, double coef1, double coef2, double coef3, double coef4, double coef5, double coef6, double coef7, Glass glass_){
	diameter = diameter_;	center_thickness = center_thickness_;	glass=glass_;
	surface0.setParameter(c_, k_, coef0, coef1, coef2, coef3, coef4, coef5, coef6, coef7);
	surface0.shape = Surface::Circle;		surface0.diameter = diameter;
	surface1.shape = Surface::Circle;		surface1.diameter = diameter;
}

inline void AsphericalSingletLens::setSurfaces(){
	surface0.X[0]=-0.5*fabs(center_thickness);	surface0.X[1]=0.0;	surface0.X[2]=0.0;
	surface0.dX[0]= 1.0;					surface0.dX[1]=0.0;	surface0.dX[2]=0.0;
	surface1.X[0]= 0.5*fabs(center_thickness);	surface1.X[1]=0.0;	surface1.X[2]=0.0;
	surface1.dX[0]=1.0;						surface1.dX[1]=0.0;	surface1.dX[2]=0.0;
};

//---	SphericalDoubletLens	---
inline SphericalDoubletLens::SphericalDoubletLens
	(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, 
		 Glass glass0_, Glass glass1_){
	setParameter(diameter_, center_thickness0_, center_thickness1_, radius0_, radius1_, radius2_, glass0_, glass1_);
}
inline SphericalDoubletLens::SphericalDoubletLens
	(double diameter_, double center_thickness0_, double center_thickness1_, 
		 double radius0_, double radius1_, double radius2_, double radius3_, double separation,
		 Glass glass0_, Glass glass1_){
	setParameter(diameter_, center_thickness0_, center_thickness1_, radius0_, radius1_, radius2_, radius3_, separation, glass0_, glass1_);
}

inline void SphericalDoubletLens::setLocation( double  X_, double  Y_, double  Z_, double dX_, double dY_, double dZ_,double dXz_, double dYz_, double dZz_){
	setLocationCore(X_, Y_, Z_, dX_, dY_, dZ_, dXz_, dYz_, dZz_);
	
	lens0.setLocation(xcenter0, 0.0, 0.0, 1.0, 0.0, 0.0);	//	relative location	
	AbsoluteCoordinate(lens0.X);	AbsoluteDirection(lens0.dX);
	lens0.setLocation(lens0.X[0], lens0.X[1], lens0.X[2], lens0.dX[0], lens0.dX[1], lens0.dX[2]);
															//	absolute location
	lens1.setLocation(xcenter1, 0.0, 0.0, 1.0, 0.0, 0.0);	//	relative location
	AbsoluteCoordinate(lens1.X);	AbsoluteDirection(lens1.dX);
	lens1.setLocation(lens1.X[0], lens1.X[1], lens1.X[2], lens1.dX[0], lens1.dX[1], lens1.dX[2]);
															//	absolute location
}

inline SphericalDoubletLens& SphericalDoubletLens::operator = (const SphericalDoubletLens& dest){
	lens0 = dest.lens0;	lens1 = dest.lens1;	xcenter0 = dest.xcenter0;	xcenter1 = dest.xcenter1;
	return *this;
}

inline void SphericalDoubletLens::setParameter
	(double diameter_, double center_thickness0_, double center_thickness1_, 
	 double radius0_, double radius1_, double radius2_, 
	 Glass glass0_, Glass glass1_)
{
	lens0.setParameter(diameter_, center_thickness0_, radius0_, radius1_, glass0_);
	xcenter0 =  -0.5*fabs(center_thickness1_) -0.5e-6;
	lens1.setParameter(diameter_, center_thickness1_, radius1_, radius2_, glass1_);
	xcenter1 =   0.5*fabs(center_thickness0_) +0.5e-6;
}

inline void SphericalDoubletLens::setParameter
	(double diameter_, double center_thickness0_, double center_thickness1_, 
	 double radius0_, double radius1_, double radius2_, double radius3_, double separation,
	 Glass glass0_, Glass glass1_)
{
	lens0.setParameter(diameter_, center_thickness0_, radius0_, radius1_, glass0_);
	xcenter0 =  -0.5*fabs(center_thickness1_) -0.5*fabs(separation);
	lens1.setParameter(diameter_, center_thickness1_, radius2_, radius3_, glass1_);
	xcenter1 =   0.5*fabs(center_thickness0_) +0.5*fabs(separation);
}

inline void SphericalDoubletLens::Go(Ray& ray){
	ray_tmp = ray;
	Go_Refract_Refract_Refract(ray_tmp);
	if(!ray_tmp.Out){
		ray = ray_tmp;
	}else
		GoRev_Refract_Refract_Refract(ray);
}
inline void SphericalDoubletLens::Go_Refract_Refract_Refract(Ray& ray){
	lens0.Go_Refract_Refract(ray);
	lens1.Go_Refract_Refract(ray);
	path=Path_Go;
}
inline void SphericalDoubletLens::Reverse(Ray& ray){GoRev_Refract_Refract_Refract(ray);}
inline void SphericalDoubletLens::GoRev_Refract_Refract_Refract(Ray& ray){
	lens1.GoRev_Refract_Refract(ray);
	lens0.GoRev_Refract_Refract(ray);
	path=Path_Reverse;
}

