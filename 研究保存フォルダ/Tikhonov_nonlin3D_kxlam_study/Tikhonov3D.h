#pragma once
#include <Eigen/eigen>
#include <IGOR/IGORitx.hpp>
#include "Tikhonov.h"

//#define DEBUG

class Tikhonov3D
{
public:
	//	constructor
	Tikhonov3D
		(const size_t size_i_, const size_t size_j_,
		 const double vlight_, const double lambda0_,
		 Eigen::VectorXd wavelength_, Eigen::VectorXd width_);
	// * whether the parameters of Tikhonov except for size_i and j are need or not should be checked *

	//	destructor
	~Tikhonov3D();

	//	set the L matrix
	void set_Lf_Amp();
	void set_Lf_Deriv();
	void set_Lf_2ndDeriv();
	void set_Lv_Amp();
	void set_Lv_Deriv();
	void set_Lv_2ndDeriv();

	//	set the A matrix
	void set_A(const Eigen::MatrixXd& A_){ A = A_; }
	//	set the cos_theta matrix
	void set_cosTheta(const Eigen::MatrixXd& CosTheta_){ CosTheta = CosTheta_; }

	//	set the d vector
	void set_d(const Eigen::MatrixXd& d_){ d = d_; }

	//	set initial values of f0
	void setF0(const Eigen::VectorXd& f0_){ f = f0_; }
	void setV0(const Eigen::VectorXd& v0_){ v = v0_; }
	void setlam0(const double& lam1_, const double& lam2_){ lam1 = lam1_; lam2 = lam2_; }
	void setkandx(const double& k0_, const double& x0_){ k = k0_; x = x0_; }


	//  set sigma
	void setSigma(const Eigen::VectorXd& sigma_){ sigma = sigma_; }

	//	solve 
	bool solve(const double relativeError0, const double vError1, const double lam1Error, const double lam2Error, const double kError, const double xError);

	//	get functions
	Eigen::MatrixXd getA()const{ return A; }
	Eigen::VectorXd getF()const{ return f; }
	Eigen::VectorXd getZ()const{ return f.array().exp().matrix(); }
	Eigen::VectorXd getV()const{ return v; }
	Eigen::VectorXd getLam()const{ return lam.array().exp().matrix(); }
	Eigen::MatrixXd getAf()const;
	Eigen::MatrixXd getCenterAf()const;
	Eigen::MatrixXd getRightAf()const;
	Eigen::MatrixXd getLeftAf()const;
	Eigen::MatrixXd getRes()const{ return d - getAf(); }
	Eigen::MatrixXd getRes_sigma()const{ return ((d - getAf()).array() * sigma.array().inverse()).matrix(); }
	Eigen::VectorXd getDf()const{ return dff.block(0, 0, size_j, 1); }
	Eigen::VectorXd getDv()const{ return dff.block(size_j, 0, size_j, 1); }
	Eigen::VectorXd getKX()const;
	double getLfNorm()const{ return (Lf*f).squaredNorm(); }
	double getLvNorm()const{ return (Lv*v).squaredNorm(); }
	double getResNorm()const{ return (d - getAf()).squaredNorm(); }
	double getResNorm_sigma()const{ return getRes_sigma().squaredNorm(); }

	double getChi2()const{
		return 0.5 * (getResNorm_sigma() + exp(lam1) * (Lf * f).squaredNorm()
			+ exp(lam2) * (Lv * v).squaredNorm() - size_j * (lam1 + lam2));
	}

	//	parameters
	const size_t size_i;
	const size_t size_j;
	const double vlight;	//	the velocity of the light
	const double lambda0;	//  the rest wavelength
	Eigen::VectorXd wavelength;
	Eigen::VectorXd width;

private:
	bool updateFandVandLamandKandX
		(const double relativeError0, const double vError1, const double lam1Error, const double lam2Error, const double kError, const double xError);
	bool updateF
		(const double relativeError0);
	bool updateV
		(const double relativeError1);
	bool updateFLam1
		(const double relativeError0, const double lam1Error);
	bool updatek
		(const double lam1Error);
	bool updateVLam2
		(const double vError1, const double lam2Error);
	bool updatex
		(const double lam2Error);
	bool updateFk
		(const double relativeError0, const double kError);
	bool updateVx
		(const double relativeError0, const double xError);
	bool setJandH();


	//	parameters
	Eigen::VectorXd finf;
	Eigen::MatrixXd A, Aright, Aleft;	//	matrix of los length
	Eigen::MatrixXd Aij, Aijright, Aijleft;
	Eigen::MatrixXd CosTheta;	//	matrix of cos_theta
	Eigen::MatrixXd Lf, Lv;
	Eigen::MatrixXd d;	//	experimental data
	Eigen::VectorXd fv;	//	solution of intensity and velocity
	Eigen::VectorXd sigma; //regidual vector

	Eigen::VectorXd J, Jf, Jv;	        //	Jacobian
	Eigen::MatrixXd H, Hff, Hvv, Hvf;	//	Hessian
	Eigen::VectorXd H1f, H2f, H1v, H2v; //  Hessian
	Eigen::VectorXd Hfk, Hvk, Hfx, Hvx; //  Hessian

	Eigen::VectorXd f;	 //	solution of intensity
	Eigen::VectorXd v;	 //	solution of velocity
	Eigen::VectorXd lam; //	solution of lambda1&2


	Eigen::VectorXd dff;
	Eigen::VectorXd count;

	double lam1, lam2;   // solution of lambda1&2
	double J1, J2, Jk, Jx;            // Jacobian
	double H11, H22, H21, Hkk, Hxx, Hxk, H1k, H1x, H2k, H2x;     // Hessian
	double dlam1, dlam2;
	double k, x, dk, dx;


	//	 small functions
	double A_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adv_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adk_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adkdv_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adx_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adxdv_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adxdx_zee(const size_t i, const size_t j, const double k, const double x)const;
	double Adkdk_zee(const size_t i, const size_t j, const double k, const double x)const;
	double A_(const size_t i, const size_t j, const double velocity)const;
	double Adv_(const size_t i, const size_t j, const double velocity)const;
	double Advdv_(const size_t i, const size_t j, const double velocity)const;
	double Az(const double aa, const double zz);
public:
	static double gauss(double dx, double s);

};

static const double pi = 3.14159265359;

inline double Tikhonov3D::gauss(double dx, double s){
	double s2_inv = 1.0/(2.0*s*s);
	return 1.0 / (sqrt(2.0*pi)*s) * exp(-dx*dx*s2_inv);
}

