#include "Tikhonov2D.h"


Tikhonov2D::Tikhonov2D
(const size_t size_i_, const size_t size_j_,
const double vlight_, const double lambda0_,
Eigen::VectorXd wavelength_, Eigen::VectorXd width_)
:
size_i(size_i_),	//	number of wavelength(pixel)
size_j(size_j_),	//	number of temperature + 1(continuous light)
vlight(vlight_),	
lambda0(lambda0_),
wavelength(wavelength_),
width(width_)

{
	J.resize(size_j * 2);	//	first size_j for f and the rest for v
	Jf.resize(size_j);
	Jv.resize(size_j);
	H.resize(size_j * 2, size_j * 2);
	Hff.resize(size_j, size_j);
	Hvf.resize(size_j, size_j);
	Hvv.resize(size_j, size_j);
	Aij.resize(size_i, size_j);
}

Tikhonov2D::~Tikhonov2D()
{
}

bool Tikhonov2D::solve
(const double lambda0, const double lambda1,
 const double relativeError0, const double vError1)
{
//	solveF(lambda0, relativeError0);
	
	count.resize(3);
	count(0) = 100;
	count(1) = 100;
	count(2) = 100;

	for (size_t trial = 0; trial < 20; ++trial){
		if (updateF(lambda0, relativeError0)) { count(0) = trial;  break; }
	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");

	for (size_t trial = 0; trial < 20; ++trial){
		if (updateV(lambda1, vError1)) { count(1) = trial;  break; }
		updateF(lambda0, relativeError0);
	}

	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");


	for (size_t trial = 0; trial < 20; ++trial){
		if (updateFandV(lambda0, lambda1, relativeError0, vError1)){
			count(2) = trial;
			IGORdata::write_itx(count, "output/prac/count.itx", "count");
			IGORdata::write_itx(count, "output/debug/count.itx", "count");
			return true;
		}
	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");
	return false;
}   

bool Tikhonov2D::updateF
(const double lambda, const double relativeError0){
	setJandH(lambda*lambda, 0.0);

	//	solve df

	Eigen::VectorXd df = -Hff.lu().solve(Jf);
	f += df;

	
	IGORdata::write_itx(f, "output/debug/f.itx", "FF");
	IGORdata::write_itx(df, "output/debug/df.itx", "df");
	IGORdata::write_itx(v, "output/debug/v.itx", "vv");
	Eigen::VectorXd Af = A*f.array().exp().matrix();
	IGORdata::write_itx(Af, "output/debug/Af.itx", "Af");
/*	IGORdata::write_itx(df, "output/debug/df.itx", "df");
	IGORdata::write_itx(J, "output/debug/J.itx", "JJ");
	IGORdata::write_itx(H, "output/debug/H.itx", "HH");
	IGORdata::write_itx(Jv, "output/debug/Jv.itx", "Jv");
	IGORdata::write_itx(Hvv, "output/debug/Hvv.itx", "Hvv");
	IGORdata::write_itx(Hvf, "output/debug/Hvf.itx", "Hvf");
	IGORdata::write_itx(Aij, "output/debug/Aij.itx", "Aij");

	*/

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((df.array().abs() < relativeError0).all())
		
		return true;
	//	otherwise
	else return false;
}
    
bool Tikhonov2D::updateV
(const double lambda1, const double vError1){
	setJandH(0.0, lambda1*lambda1);

	//	solve df
	
	Eigen::VectorXd dv = -Hvv.lu().solve(Jv);
	v += dv;
	
	IGORdata::write_itx(f, "output/debug/f.itx", "FF");
	IGORdata::write_itx(v, "output/debug/v.itx", "vv");
	IGORdata::write_itx(dv, "output/debug/dv.itx", "dv");
	Eigen::VectorXd Af = A*f.array().exp().matrix();
	IGORdata::write_itx(Af, "output/debug/Af.itx", "Af");
	/*
	IGORdata::write_itx(J, "output/debug/J.itx", "JJ");
	IGORdata::write_itx(H, "output/debug/H.itx", "HH");
	IGORdata::write_itx(Jv, "output/debug/Jv.itx", "Jv");
	IGORdata::write_itx(Hvv, "output/debug/Hvv.itx", "Hvv");
	IGORdata::write_itx(Hvf, "output/debug/Hvf.itx", "Hvf");
	IGORdata::write_itx(Aij, "output/debug/Aij.itx", "Aij");

	IGORdata::write_itx(dv, "output/debug/dv.itx", "dv");
	*/


	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((dv.array().abs() < vError1).all())
		return true;
	//	otherwise
	else return false;
}


bool Tikhonov2D::updateFandV
(const double lambda0, const double lambda1,
 const double relativeError0, const double vError1){

	setJandH(lambda0*lambda0, lambda1*lambda1);

	//	solve df
	dff = -H.lu().solve(J);
	f += dff.block(0, 0, size_j, 1);
	v += dff.block(size_j, 0, size_j, 1);

	IGORdata::write_itx(f, "output/debug/f.itx", "FF");
	IGORdata::write_itx(v, "output/debug/v.itx", "vv");
	IGORdata::write_itx(dff, "output/debug/dff.itx", "dff");
/*
	IGORdata::write_itx(J, "output/debug/J.itx", "JJ");
	IGORdata::write_itx(H, "output/debug/H.itx", "HH");
	IGORdata::write_itx(Jv, "output/debug/Jv.itx", "Jv");
	IGORdata::write_itx(Hvv, "output/debug/Hvv.itx", "Hvv");
	IGORdata::write_itx(Hvf, "output/debug/Hvf.itx", "Hvf");
	IGORdata::write_itx(Aij, "output/debug/Aij.itx", "Aij");
	Eigen::VectorXd Af = A*f.array().exp().matrix();
	IGORdata::write_itx(Af, "output/debug/Af.itx", "Af");
	*/

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	for (size_t j = 0; j < size_j; ++j){
		if (fabs(dff[j]) > relativeError0) { return false; }
	}
	for (size_t j = 0; j < size_j; ++j){
		if (fabs(dff[j + size_j]) > vError1){ return false; }
	}
	//	otherwise
	return true;
}

bool Tikhonov2D::solveF(const double lambda, const double relativeError0){
	//	get the initial value of f
	Tikhonov tikhonov0(size_i, size_j);
	tikhonov0.set_L_2ndDeriv();
	tikhonov0.set_A(A);
	tikhonov0.setFinf_zero();
	tikhonov0.set_d(d);
	tikhonov0.setF0(f);
	tikhonov0.setSigma(sigma);
	tikhonov0.solve(lambda, relativeError0);
	f = tikhonov0.getF();

	return true;
}

bool Tikhonov2D::setJandH(const double lambdaSquared0, const double lambdaSquared1){
	Jf.setZero();
	Jv.setZero();
	Hff.setZero();
	Hvf.setZero();
	Hvv.setZero();
	J.setZero();
	H.setZero();

	Jf = 2.0 * lambdaSquared0 * Lf.transpose() * Lf * f;
	Jv = 2.0 * lambdaSquared1 * Lv.transpose() * Lv * v;
	Hff = 2.0 * lambdaSquared0 * Lf.transpose() * Lf;
	Hvv = 2.0 * lambdaSquared1 * Lv.transpose() * Lv;

	//	z0_j = exp(f0_j)
	Eigen::VectorXd z = f.array().exp().matrix();

	Aij.setZero();

	for (size_t i = 0; i < size_i; ++i){
		double res_i = d(i);

		for (size_t l = 0; l < size_j; ++l)
			res_i -= A(i, l) * z[l];

		for (size_t j = 0; j < size_j; ++j){
			
			double A_ij   = getA_ij(i, j);
			double Adv_ij = getAdv_ij(i, j);
			Aij(i, j) = A_ij;
			           
			Jf[j] += -2.0 * res_i * A_ij   * z[j] / sigma[i] / sigma[i];
			Jv[j] += -2.0 * res_i * Adv_ij * z[j] / sigma[i] / sigma[i];

			for (size_t k = 0; k < size_j; ++k){
				double A_ik   = getA_ij(i, k);
				double Adv_ik = getAdv_ij(i, k);

				Hff(j, k) += 2.0 * A_ij   * z[j] * A_ik   * z[k] / sigma[i] / sigma[i];
				Hvf(j, k) += 2.0 * Adv_ij * z[j] * A_ik   * z[k] / sigma[i] / sigma[i];
				Hvv(j, k) += 2.0 * Adv_ij * z[j] * Adv_ik * z[k] / sigma[i] / sigma[i];
			}	
		}
	}
	A = Aij;

	J.block(     0, 0, size_j, 1) = Jf;
	J.block(size_j, 0, size_j, 1) = Jv;
	H.block(     0,      0, size_j, size_j) = Hff;
	H.block(     0, size_j, size_j, size_j) = Hvf;
	H.block(size_j,      0, size_j, size_j) = Hvf.transpose();
	H.block(size_j, size_j, size_j, size_j) = Hvv;
	return true;
}

double Tikhonov2D::getA_ij(const size_t i, const size_t j)const{
	
	if (j == size_j - 1)
		return 1.0;
	else
		return (exp(-((wavelength(i) - lambda0 +v[j] * lambda0 / vlight) / width(j)) 
	        	* ((wavelength(i) - lambda0 + v[j] * lambda0 / vlight) / width(j))));
	
}
double Tikhonov2D::getAdv_ij(const size_t i, const size_t j)const{
	if (j == size_j - 1)
		return 0;
	else
		return -(getA_ij(i, j) * 2 * lambda0 * (wavelength(i) - lambda0 + v[j] * lambda0 / vlight)
			/ vlight / width(j) / width(j));
}

Eigen::MatrixXd Tikhonov2D::getAf()const{
	Eigen::VectorXd Af = A*f.array().exp().matrix();
	return Af;
}

void Tikhonov2D::set_Lf_Amp(){
	Lf.resize(size_j, size_j);
	Lf.setIdentity();
}
void Tikhonov2D::set_Lf_Deriv(){
	Lf.resize(size_j, size_j);
	Lf.setZero();
	Lf(0, 0) = 2.0;
	for (size_t j = 0; j < size_j - 1; ++j){
		Lf(j+1, j    ) = -1.0;
		Lf(j+1, j + 1) = 1.0;
	}
}
void Tikhonov2D::set_Lf_2ndDeriv(){
	Lf.resize(size_j - 1, size_j);
	Lf.setZero();

	Lf(0, 0) = 2.0;
	Lf(0, 1) = -2.0;
	for (size_t j = 0; j < size_j - 2; ++j){
		Lf(j + 1, j) = 1.0;
		Lf(j + 1, j + 1) = -2.0;
		Lf(j + 1, j + 2) = 1.0;
	}
}
void Tikhonov2D::set_Lv_Amp(){
	Lv.resize(size_j, size_j);
	Lv.setIdentity();
}
void Tikhonov2D::set_Lv_Deriv(){
	Lv.resize(size_j, size_j);
	Lv.setZero();
	Lv(0, 0) = 2.0;
	for (size_t j = 0; j < size_j - 1; ++j){
		Lv(j+1, j) = -1.0;
		Lv(j+1, j + 1) = 1.0;
	}
}
void Tikhonov2D::set_Lv_2ndDeriv(){
	Lv.resize(size_j - 1, size_j);
	Lv.setZero();

	Lv(0, 0) = 2.0;
	Lv(0, 1) = -2.0;
	for (size_t j = 0; j < size_j - 2; ++j){
		Lv(j + 1, j) = 1.0;
		Lv(j + 1, j + 1) = -2.0;
		Lv(j + 1, j + 2) = 1.0;
	}
}
