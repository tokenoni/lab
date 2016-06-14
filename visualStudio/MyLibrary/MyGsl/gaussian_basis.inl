

namespace mygsl{
//-----------------------------------------------------------------------//
//																		 //
//				class for One Gaussian Functions						 //
//																		 //
//-----------------------------------------------------------------------//
	//---	copy constructor	---
	inline OneGaussian::OneGaussian(const OneGaussian& obj){
		xi=obj.l();	yi=obj.m(); zi=obj.n();
		alpha=obj.a();	coef = obj.c();
		calculateCoef();
	}

	//---	setting functions	---
	inline void OneGaussian::set(double exponent, size_t xi_, size_t yi_, size_t zi_){
		xi=xi_;	yi=yi_;	zi=zi_;
		alpha=exponent;
		calculateCoef();
	}
	inline void OneGaussian::set(double exponent, size_t xi_){
		xi=xi_;	yi=0;	zi=0;
		alpha=exponent;
		calculateCoef();
	}

	//---	determines the normalization coefficient	---
	inline void OneGaussian::calculateCoef(){
		double integ = OneIntegral(2*xi, 2.0*alpha)
					 * OneIntegral(2*yi, 2.0*alpha)
					 * OneIntegral(2*zi, 2.0*alpha);
		coef = 1.0/sqrt(integ);
		area = coef
			 * OneIntegral(xi,   alpha)
			 * OneIntegral(yi,   alpha)
			 * OneIntegral(zi,   alpha);
		moment[0] = coef
			 * OneIntegral(xi+1, alpha)
			 * OneIntegral(yi  , alpha)
			 * OneIntegral(zi  , alpha);
		moment[1] = coef
			 * OneIntegral(xi  , alpha)
			 * OneIntegral(yi+1, alpha)
			 * OneIntegral(zi  , alpha);
		moment[2] = coef
			 * OneIntegral(xi  , alpha)
			 * OneIntegral(yi  , alpha)
			 * OneIntegral(zi+1, alpha);
	}

	//---	overloaden operators for mathematical description	---
	//---	returns the overlap	---
	inline double OneGaussian::operator * (const OneGaussian& obj)const {	return OneIntegral(*this, obj, 0, 0, 0);}

	inline OneGaussian& OneGaussian::operator = (const OneGaussian& obj){
		xi=obj.l();	yi=obj.m(); zi=obj.n();
		alpha=obj.a();	coef = obj.c();
		calculateCoef();
		return *this;
	}

	inline bool OneGaussian::operator!= (const OneGaussian& obj)const{
		if( l() != obj.l()) return true;
		if( m() != obj.m()) return true;
		if( n() != obj.n()) return true;
		if( a()*DBL_EPSILON < obj.a() || obj.a()*DBL_EPSILON < a()) return true;	//	machine epsilon is taken into account
		return false;
	}

	inline bool OneGaussian::operator== (const OneGaussian& obj)const{
		if(*this != obj) return false;
		else return true;
	}

	//---	functions returning the values at a certain point	---
	//---	3 dimensional	---
	inline double OneGaussian::get(double x, double y, double z)const{
		return coef*pow(x, 1.0*xi)*pow(y,1.0*yi)*pow(z,1.0*zi)*exp(-alpha*(x*x+y*y+z*z));
	}
	//---	2 dimensional	---
	inline double OneGaussian::get(double x, double y)const{
		return coef*pow(x, 1.0*xi)*pow(y,1.0*yi)*exp(-alpha*(x*x+y*y))*OneIntegral(zi, alpha);
	}
	//---	1 dimensional	---
	inline double OneGaussian::get(double x)const{
		return coef*pow(x, 1.0*xi)*exp(-alpha*(x*x))*OneIntegral(yi, alpha)*OneIntegral(zi, alpha);
	}
	//---	static functions (very basic functions)	---
	//	returns the value of integ[x^i exp(-a x^2)]dx
	inline double OneGaussian::OneIntegral(size_t i, double alpha_){
		if(i%2!=0 || alpha_ <= 0.0) return 0.0;
		else{
			double J = sqrt(M_PI/alpha_);
			for(size_t j=0; j<i; j+=2)	J *= (1.0*j+1.0)/alpha_;
			return J;
		}
	}
	//	returns the value of <phi1| x^i y^j z^k |phi2>
	inline double OneGaussian::OneIntegral(  const OneGaussian& phi1, const OneGaussian& phi2,
									const size_t pow_x, const size_t pow_y, size_t pow_z){
		return phi1.c()*phi2.c()
			*OneIntegral(phi1.l() + phi2.l() + pow_x, phi1.a() + phi2.a())
			*OneIntegral(phi1.m() + phi2.m() + pow_y, phi1.a() + phi2.a())
			*OneIntegral(phi1.n() + phi2.n() + pow_z, phi1.a() + phi2.a());
	}

//-----------------------------------------------------------------------//
//																		 //
//				class for ContractedGaussian Functions					 //
//																		 //
//-----------------------------------------------------------------------//
	//---	setting functions for Non-contracted Gaussian basis	---
	inline void Gaussian::set(double exponent, size_t xi_, size_t yi_, size_t zi_){
		gaussian_set.resize(1);
		ratio.resize(1, 1.0);
		gaussian_set[0].set(exponent, xi_, yi_, zi_);
	}	
	inline void Gaussian::set(double exponent, size_t xi_){
		set(exponent, xi_, 0, 0);
	}
	//---	setting functions for Contracted Gaussian basis	---
	inline void Gaussian::set(const std::vector<OneGaussian>& gaussian_set_, const std::vector<double>& ratio_){
		gaussian_set = gaussian_set_;
		ratio = ratio_;
		double sum = 0.0;
		for(size_t i=0; i<ratio.size(); ++i) sum+=ratio[i];
		for(size_t i=0; i<ratio.size(); ++i) ratio[i]/=sum;
		area=0.0;	moment[0]=0.0;	moment[1]=0.0;	moment[2]=0.0;
		for(size_t i=0; i<ratio.size(); ++i) {
			area +=      ratio[i]*gaussian_set[i].getArea();
			moment[0] += ratio[i]*gaussian_set[i].getMoment_x();
			moment[1] += ratio[i]*gaussian_set[i].getMoment_y();
			moment[2] += ratio[i]*gaussian_set[i].getMoment_z();
		}
	}
	
	inline double Gaussian::get(double x, double y, double z)const{	//	3 dimensional
		double val = 0.0;
		for(size_t i=0; i<gaussian_set.size();++i) val += ratio[i]*gaussian_set[i].get(x, y, z);
		return val;
	}
	inline double Gaussian::get(double x, double y)const{;			//	2 dimensional
		double val = 0.0;
		for(size_t i=0; i<gaussian_set.size();++i) val += ratio[i]*gaussian_set[i].get(x, y);
		return val;
	}
	inline double Gaussian::get(double x)const{;						//	1 dimensional
		double val = 0.0;
		for(size_t i=0; i<gaussian_set.size();++i) val += ratio[i]*gaussian_set[i].get(x);
		return val;
	}

	//---	static functions (very basic functions)	---
	//	returns the value of <phi1|phi2>
	inline double Gaussian::OneIntegral(const Gaussian& phi1, const Gaussian& phi2,
									const size_t pow_x, const size_t pow_y, size_t pow_z){
		double val=0.0;
		for(size_t i=0; i<phi1.size();++i)
			for(size_t j=0; j<phi2.size();++j)
				val+= phi1.getCoef(i)*phi2.getCoef(j)
					* OneGaussian::OneIntegral(phi1[i], phi2[j], pow_x, pow_y, pow_z);
		return val;
	}

	//---	overriden operators	---
	inline double OneGaussian::operator * (const Gaussian& obj)const {
		double val=0.0;
		for(size_t i=0; i<obj.size(); ++i)
			val += *this * obj[i];
		return val;
	}
	inline double Gaussian::operator * (const Gaussian& obj)const {
		double val=0.0;
		for(size_t i=0; i<size(); ++i)
			for(size_t j=0; j<obj.size(); ++j)
				val += (*this)[i] * obj[j];
		return val;
	}

//-----------------------------------------------------------------------//
//																		 //
//					class for Overlap Matricies							 //
//																		 //
//-----------------------------------------------------------------------//


	//---	constructors	---
	inline Overlap::Overlap(void){	allocated = false; Phi.resize(0);}
	inline Overlap::Overlap(const Overlap& obj){				 Copy(obj);}
	inline Overlap& Overlap::operator = (const Overlap& obj){	 Copy(obj); return *this;}
	inline Overlap::~Overlap(void){ Free();}

	inline void Overlap::Copy(const Overlap& obj){
		MemoryControl(obj.size());
		Phi = getPhi();
		gsl_matrix_memcpy(delta, obj());
		gsl_matrix_memcpy(delta_x, obj.x());
		gsl_matrix_memcpy(delta_y, obj.y());
		gsl_matrix_memcpy(delta_z, obj.z());
		gsl_matrix_memcpy(delta_inv, obj.inv());
		gsl_matrix_memcpy(delta_invsqrt, obj.invsqrt());
	}
	//---	memory controling functions	---
	inline void Overlap::MemoryControl(size_t size_){
		if(allocated && size() != size_)	Free();
		if(!allocated) Allocate(size_);
	}
	inline void Overlap::Allocate(size_t size_){
		if(!allocated){
			delta         = gsl_matrix_calloc(size_, size_);
			delta_x       = gsl_matrix_calloc(size_, size_);
			delta_y       = gsl_matrix_calloc(size_, size_);
			delta_z       = gsl_matrix_calloc(size_, size_);
			delta_inv     = gsl_matrix_calloc(size_, size_);
			delta_invsqrt = gsl_matrix_calloc(size_, size_);
		}
		allocated = true;
	}
	inline void Overlap::Free(){
		if(allocated){
			gsl_matrix_free(delta        );
			gsl_matrix_free(delta_x      );
			gsl_matrix_free(delta_y      );
			gsl_matrix_free(delta_z      );
			gsl_matrix_free(delta_inv    );
			gsl_matrix_free(delta_invsqrt);
		}
		allocated = false;
	}

	//---	setting functions	---
	inline void Overlap::set(std::vector<Gaussian>& Phi_){ Phi = Phi_; MemoryControl(size());}
	inline void Overlap::resize(size_t size_){		Phi.resize(size_); MemoryControl(size());}

	//---	core functions for calculating the overlap	---
	inline void Overlap::calculate(){
		MemoryControl(size());
		for(size_t i=0; i<size(); ++i){
			for(size_t j=0; j<size(); ++j){
				gsl_matrix_set(delta  , i, j, Gaussian::OneIntegral(Phi[i], Phi[j], 0, 0, 0));
				gsl_matrix_set(delta_x, i, j, Gaussian::OneIntegral(Phi[i], Phi[j], 1, 0, 0));
				gsl_matrix_set(delta_y, i, j, Gaussian::OneIntegral(Phi[i], Phi[j], 0, 1, 0));
				gsl_matrix_set(delta_z, i, j, Gaussian::OneIntegral(Phi[i], Phi[j], 0, 0, 1));
			}
		}
		//---	deriving the inverse matrix	---
		gsl_matrix_memcpy(delta_inv, delta);
		gsl_linalg_cholesky_decomp(delta_inv);
		gsl_linalg_cholesky_invert(delta_inv);

		//---	deriving the sqrt inverse matrix	---		
		//	(I have no confidence on the following equation. should be checked)
		gsl_matrix_memcpy(delta_invsqrt, delta_inv);
		gsl_linalg_cholesky_decomp(delta_invsqrt);
	}


	//---	core functions for calculating the overlap	---
	inline Overlap Overlap::getTransformedBasisSet(const gsl_matrix* transform_matrix)const{
		Overlap copy(*this);
		copy.BasisTransformation(transform_matrix);
		return copy;
	}

	inline void Overlap::BasisTransformation(const gsl_matrix* transform_matrix){
		std::vector<OneGaussian> one_gaussian_set(0);
		std::vector<std::vector<size_t>> index(Phi.size(), std::vector<size_t>(0));

		//---	once the gaussian components are gathered into one std::vector<OneGaussian>	---
		for(size_t i=0; i<Phi.size();++i){
			index[i].resize(Phi[i].size());
			for(size_t j=0; j<Phi[i].size(); ++j){
				size_t k;
				for(k=0; k<one_gaussian_set.size();++k)
					if(one_gaussian_set[k] == Phi[i][j]) break;
				if(k==one_gaussian_set.size()) one_gaussian_set.push_back(Phi[i][j]);
				index[i][j] = k;
			}
		}

		std::vector<std::vector<double>> coef(Phi.size(), std::vector<double>(one_gaussian_set.size(), 0.0));
		for(size_t i=0; i<Phi.size();++i){
			for(size_t j=0; j<Phi[i].size(); ++j){
				coef[i][index[i][j]] += Phi[i].getCoef(j);
			}
		}

		//---	then the transfered bases are generated	---
		for(size_t i=0; i<Phi.size();++i){
			std::vector<double> coef_trans(Phi.size(), 0.0);
			for(size_t k=0; k<one_gaussian_set.size();++k){
				for(size_t j=0; j<Phi.size();++j){
					coef_trans[k] += gsl_matrix_get(transform_matrix, i, j) * coef[j][k];
				}
			}
			Phi[i].set(one_gaussian_set, coef_trans);
		}
		calculate();
	}
};