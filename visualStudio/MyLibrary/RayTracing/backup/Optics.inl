
//------------------	Class Optics	---------------------
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

inline void CreateRay::SpreadingSource(size_t NumRay, double SpreadAngle, double diameter){
	RaySet_Relative.resize(NumRay);
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		do{
			RaySet_Relative[i].x[1]=diameter*(0.5-1.0*rand()/RAND_MAX);
			RaySet_Relative[i].x[2]=diameter*(0.5-1.0*rand()/RAND_MAX);
		}while(RaySet_Relative[i].x[1]*RaySet_Relative[i].x[1]+RaySet_Relative[i].x[2]*RaySet_Relative[i].x[2] > 0.25*diameter*diameter);
		double theta = 0.5*SpreadAngle*rand()/RAND_MAX;
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}

inline void CreateRay::SpreadingSource(size_t NumRay, double SpreadAngle, double size_y, double size_z){
	RaySet_Relative.resize(NumRay);
	for(size_t i=0; i<NumRay; i++){
		RaySet_Relative[i].x[0]=0.0;
		RaySet_Relative[i].x[1]=size_y*(0.5-1.0*rand()/RAND_MAX);
		RaySet_Relative[i].x[2]=size_z*(0.5-1.0*rand()/RAND_MAX);
		double theta = 0.5*SpreadAngle*rand()/RAND_MAX;
		double phi   = 2.0*M_PI*rand()/RAND_MAX;
		RaySet_Relative[i].dx[0]=cos(theta);
		RaySet_Relative[i].dx[1]=sin(theta)*sin(phi);
		RaySet_Relative[i].dx[2]=sin(theta)*cos(phi);
		RaySet_Relative[i].Out = false;
	}
}
	
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
		AbsoluteCoordinate(x);			
	}
	return Outflag;
}

inline bool Surface::GoForward(Ray &ray){
	if(ray.Out)			 {copy(ray_pos, pNAN);	return true;}	//	already out of optics
	if(CrossPoint(ray))	 {
		copy(ray_pos, pNAN);	return true;}	//	cross point is in opsite side
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
	if(coef_para>=0.0) coef_para=sqrt(1.0-coef_perp*coef_perp);
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
	if(l<0.0){ // on the convex surface
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
inline void Prism::setSurfaces(){
	surface0.X[0]  = 0.5*size_z*sin(0.5*angle);		surface0.X[1] = 0.0;	surface0.X[2] = 0.0;
	surface0.dX[0] = cos(0.5*angle);				surface0.dX[1] = 0.0;	surface0.dX[2] = sin(0.5*angle);
	surface0.shape = surface0.Rectangular;			surface0.size_y = size_y;	surface0.size_z = size_z;

	surface1.X[0]  =-0.5*size_z*sin(0.5*angle);		surface1.X[1] = 0.0;	surface1.X[2] = 0.0;
	surface1.dX[0] = cos(0.5*angle);				surface1.dX[1] = 0.0;	surface1.dX[2] =-sin(0.5*angle);
	surface1.shape = surface1.Rectangular;			surface1.size_y = size_y;	surface1.size_z = size_z;
}
//---	SphericalSingletLens	---
inline void SphericalSingletLens::setSurfaces(){
	surface0.X[0]=-0.5*fabs(center_thickness);	surface0.X[1]=0.0;	surface0.X[2]=0.0;
	surface0.dX[0]= 1.0;					surface0.dX[1]=0.0;	surface0.dX[2]=0.0;
	surface0.shape = surface0.Circle;		surface0.diameter=diameter;
	surface0.R = radius0;					//surface0.surface_type=surface0.Convex;

	surface1.X[0]= 0.5*fabs(center_thickness);	surface1.X[1]=0.0;	surface1.X[2]=0.0;
	surface1.dX[0]=1.0;					surface1.dX[1]=0.0;	surface1.dX[2]=0.0;
	surface1.shape = surface1.Circle;		surface1.diameter=diameter;
	surface1.R = radius1;					//surface1.surface_type=surface1.Concave;
}

//---	SphericalDoubletLens	---
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
