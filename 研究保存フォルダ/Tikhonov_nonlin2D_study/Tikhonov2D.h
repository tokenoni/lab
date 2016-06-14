#pragma once
#include <Eigen/eigen>
#include <IGOR/IGORitx.hpp>
#include "Tikhonov.h"

class Tikhonov2D
{
public:
	//	constructor
	Tikhonov2D
		(const size_t size_i_, const size_t size_j_,
		 const double vlight_, const double lambda0_,
		 Eigen::VectorXd wavelength_, Eigen::VectorXd width_);
	// * whether the parameters of Tikhonov except for size_i and j are need or not should be checked *

	//	destructor
	~Tikhonov2D();

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

	//  set sigma
	void setSigma(const Eigen::VectorXd& sigma_){ sigma = sigma_; }

	//	solve 
	bool solve(const double lambda0, const double lambda1,
		const double relativeError0, const double vError1);

	//	get functions
	Eigen::MatrixXd getA()const{ return A; }
	Eigen::VectorXd getF()const{ return f; }
	Eigen::VectorXd getZ()const{ return f.array().exp().matrix(); }
	Eigen::VectorXd getV()const{ return v; }
	Eigen::MatrixXd getAf()const;
	Eigen::MatrixXd getRes()const{ return d - getAf(); }
	Eigen::MatrixXd getRes_sigma()const{ return ((d - getAf()).array() * sigma.array().inverse()).matrix(); }
	Eigen::VectorXd getDf()const{ return dff.block(0, 0, size_j, 1); }
	Eigen::VectorXd getDv()const{ return dff.block(size_j, 0, size_j, 1); }
	double getLfNorm()const{ return (Lf*f).squaredNorm(); }
	double getLvNorm()const{ return (Lv*v).squaredNorm(); }
	double getResNorm()const{ return (d - getAf()).squaredNorm(); }
	double getResNorm_sigma()const{ return getRes_sigma().squaredNorm(); }


	//	parameters
	const size_t size_i;
	const size_t size_j;
	const double vlight;	//	the velocity of the light
	const double lambda0;	//  the rest wavelength
	Eigen::VectorXd wavelength;
	Eigen::VectorXd width;

private:
	bool solveF(const double lambda0, const double relativeError0);

	bool updateFandV
		(const double lambda, const double lambda1, 
		 const double relativeError0, const double vError1);
	bool updateF
		(const double lambda,
		const double relativeError0);
	bool updateV
		(const double lambda1,
		const double relativeError1);
	bool setJandH(const double lambdaSquared0, const double lambdaSquared1);


	//	parameters
	Eigen::VectorXd finf;
	Eigen::MatrixXd A;	//	matrix of los length
	Eigen::MatrixXd Aij;
	Eigen::MatrixXd CosTheta;	//	matrix of cos_theta
	Eigen::MatrixXd Lf, Lv;
	Eigen::MatrixXd d;	//	experimental data
	Eigen::VectorXd fv;	//	solution of intensity and velocity
	Eigen::VectorXd sigma; //regidual vector

	Eigen::VectorXd J, Jf, Jv;	//	Jacobian
	Eigen::MatrixXd H, Hff, Hvv, Hvf;	//	Hessian

	Eigen::VectorXd f;	//	solution of intensity
	Eigen::VectorXd v;	//	solution of velocity

	Eigen::VectorXd dff;

	Eigen::VectorXd count;

	//	 small functions
	double getA_ij(const size_t i, const size_t j)const;
	double getAdv_ij(const size_t i, const size_t j)const;
public:
	static double gauss(double dx, double s);

};

static const double pi = 3.14159265359;

inline double Tikhonov2D::gauss(double dx, double s){
	double s2_inv = 1.0/(2.0*s*s);
	return 1.0 / (sqrt(2.0*pi)*s) * exp(-dx*dx*s2_inv);
}

