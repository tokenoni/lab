#include "Tikhonov.h"


Tikhonov::Tikhonov(const size_t size_i_, const size_t size_j_)
:
size_i(size_i_),	//	size of data point
size_j(size_j_)		//	size of solution
{
	J.resize(size_j);
	H.resize(size_j, size_j);
}

Tikhonov::~Tikhonov(){}

bool Tikhonov::solve(const double lambda, const double relativeError){
//	f = f0;
	for (size_t trial = 0; trial < 20; ++trial){
		if (updateF(lambda, relativeError)) {
//			f = f0;
			return true;
		}
	}
	return false;
}

bool Tikhonov::updateF(const double lambda, const double relativeError){
	double l2 = lambda*lambda;
	setJandH(l2);

	//	solve df
	Eigen::VectorXd df = -H.lu().solve(J);
	f0 += df;
	
	//-------DEBUG----------
	IGORdata::write_itx(f0, "tmp/f0.itx", "f0");
	IGORdata::write_itx(df, "tmp/df.itx", "df");
	//----------------------

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((df.array().abs() < log(relativeError)).all()) 
		return true; 
	//	otherwise
	else return false;
}

bool Tikhonov::setJandH(const double l2){
	//	z0_j = exp(f0_j)
	Eigen::VectorXd z0 = f0.array().exp().matrix();

	//	res = d - Az
	Eigen::VectorXd res = d - A*z0;

	//	second term of J
	J = 2.0 * l2 * L.transpose() * L * f0;
	//	second term of H
	H = 2.0 * l2 * L.transpose() * L;

	//	add 1st term of J and H
	for (size_t j = 0; j < size_j; ++j){
		double j1 = 0.0;
		for (size_t i = 0; i < size_i; ++i)
			j1 += res[i] * A(i, j) * z0(j) / sigma[i] / sigma[i];
		J[j] += -2.0 * j1;

		for (size_t k = 0; k < size_j; ++k){
			double h1 = 0.0;
			for (size_t i = 0; i < size_i; ++i){
				h1 -= A(i, k) * z0[k] * A(i, j) * z0[j] / sigma[i] / sigma[i];
				if (j == k)
					h1 += res[i] * A(i, j)*z0[j] / sigma[i] / sigma[i];
			}
			H(j, k) += -2.0*h1;
		}
	}

	//-------DEBUG----------
//	Eigen::VectorXd Az = A*z0;
//	IGORdata::write_itx(z0, "tmp/z0.itx", "z0");
//	IGORdata::write_itx(Az, "tmp/Az.itx", "Az");
//	IGORdata::write_itx(J, "tmp/JJ.itx", "JJ");
//	IGORdata::write_itx(H, "tmp/HH.itx", "HH");
	//----------------------
	return true;
}

double Tikhonov::getNorm(){
	return getRes().squaredNorm();
}

double Tikhonov::getRough(){
	return getLF().squaredNorm();
}

void Tikhonov::set_L_Amp(){
	L.resize(size_j, size_j);
	L.setIdentity();
}
void Tikhonov::set_L_Deriv(){
	L.resize(size_j-1, size_j);
	L.setZero();
	for (size_t j = 0; j < size_j-1; ++j){
		L(j, j    ) = -1.0;
		L(j, j + 1) =  1.0;
	}
}
void Tikhonov::set_L_2ndDeriv(){
	L.resize(size_j - 1, size_j);
	L.setZero();

	L(0, 0) = 2.0;
	L(0, 1) =-2.0;
	for (size_t j = 0; j < size_j - 3; ++j){
		L(j+1, j)     = 1.0;
		L(j+1, j + 1) =-2.0;
		L(j+1, j + 2) = 1.0;
	}
}
