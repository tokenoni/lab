#pragma once
#include <Eigen/eigen>
#include <IGOR/IGORitx.hpp>

class Tikhonov
{
public:
	//	constructor
	Tikhonov(const size_t size_i_, const size_t size_j_);
	//	destructor
	~Tikhonov();

	//	set the default value
	void setFinf(const Eigen::VectorXd finf_){ finf = finf_; }
	void setFinf_zero(){ finf.resize(size_j); finf.setZero(); }

	//	set the L matrix
	void set_L_Amp();
	void set_L_Deriv();
	void set_L_2ndDeriv();

	//	set the A matrix
	void set_A(const Eigen::MatrixXd& A_){ A = A_; }

	//	set the d vector
	void set_d(const Eigen::VectorXd& d_){ d = d_; }

	//	set the alpha matrix
	void set_alpha(const Eigen::MatrixXd& alpha_){ alpha = alpha_; }

	//	set initial values of f0
	void setF0(const Eigen::VectorXd& f0_){ f0 = f0_; }

	//  set sigma
	void setSigma(const Eigen::VectorXd& sigma_){ sigma = sigma_; }

	//	solve 
	bool solve(const double lambda, const double relativeError);

	//	get roughness

	double getRough();

	//	get residual norm

	double getNorm();

	//	get functions
	Eigen::MatrixXd getA()const{ return A; }
	Eigen::MatrixXd getLF()const{ return L*(f0 - finf); }
	Eigen::VectorXd getF()const{ return f0; }
	Eigen::VectorXd getAf()const{ return A*f0.array().exp().matrix(); }
	Eigen::VectorXd getRes()const{ return d - A*f0.array().exp().matrix(); }

	//	parameters
	const size_t size_i;
	const size_t size_j;

private:
	bool updateF(const double lambda, const double relativeError);
	bool setJandH(const double lambdaSquared);

	//	parameters
	Eigen::VectorXd finf;
	Eigen::MatrixXd A;
	Eigen::MatrixXd L;
	Eigen::VectorXd d;	//	experimental data
	Eigen::VectorXd f;	//	solution
	Eigen::MatrixXd alpha; //regidual matrix
	Eigen::VectorXd sigma; //regidual vector

	Eigen::VectorXd J;	//	Jacobian
	Eigen::MatrixXd H;	//	Hessian
	Eigen::VectorXd f0;	//	current values of f
};

