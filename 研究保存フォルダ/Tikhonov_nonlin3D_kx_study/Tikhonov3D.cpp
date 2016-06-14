#include "Tikhonov3D.h"


Tikhonov3D::Tikhonov3D
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
	J.resize(size_j * 2 + 2);	//	first size_j for f and the rest for v
	Jf.resize(size_j);
	Jv.resize(size_j);
	H.resize(size_j * 2 + 2, size_j * 2 + 2);
	Hff.resize(size_j, size_j);
	Hvf.resize(size_j, size_j);
	Hvv.resize(size_j, size_j);
	H1f.resize(size_j);
	H2f.resize(size_j);
	H1v.resize(size_j);
	H2v.resize(size_j);
	Aij.resize(size_i, size_j);
	Aijright.resize(size_i, size_j);
	Aijleft.resize(size_i, size_j);
	lam.resize(2);
}

Tikhonov3D::~Tikhonov3D()
{
}

bool Tikhonov3D::solve
(const double relativeError0, const double vError1
, const double lam1Error, const double lam2Error
, const double kError, const double xError)
{

	count.resize(5);
	count.setOnes();
	count *= 100.0;
	size_t trialmax = 20;





	for (size_t trial = 0; trial < trialmax; ++trial){
		if (updateF(relativeError0)) { count(0) = 1.0*trial;  break; }
		std::cout << "trial1 = " << trial << std::endl;
/*
		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		Eigen::VectorXd Af = A*f.array().exp().matrix();
		IGORdata::write_itx(Af, "output/debug/F/Af.itx", "Af");
		
			IGORdata::write_itx(df, "output/debug/F/df.itx", "df");
			IGORdata::write_itx(df, "output/debug/F/df.itx", "df");
			IGORdata::write_itx(J, "output/debug/F/J.itx", "JJ");
			IGORdata::write_itx(H, "output/debug/F/H.itx", "HH");
			IGORdata::write_itx(Jv, "output/debug/F/Jv.itx", "Jv");
			IGORdata::write_itx(Hvv, "output/debug/F/Hvv.itx", "Hvv");
			IGORdata::write_itx(Hvf, "output/debug/F/Hvf.itx", "Hvf");
			IGORdata::write_itx(Aij, "output/debug/F/Aij.itx", "Aij");

	*/		
		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(getAf(), "output/debug/F/Af.itx", "Af");
		IGORdata::write_itx(getRes(), "output/debug/F/res.itx", "res");
		IGORdata::write_itx(getCenterAf(), "output/debug/F/centerAf.itx", "centerAf");
		IGORdata::write_itx(getLeftAf(), "output/debug/F/leftAf.itx", "leftAf");
		IGORdata::write_itx(getRightAf(), "output/debug/F/rightAf.itx", "rightAf");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		IGORdata::write_itx(getZ(), "output/debug/F/z.itx", "ZZ");

	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");

	IGORdata::write_itx(sigma, "output/prac/sigma.itx", "sigma");






	for (size_t trial = 0; trial < trialmax; ++trial){
		if (updateV(vError1)) { count(1) = 1.0*trial;  break; }
		updateF(relativeError0);
		std::cout << "trial2 = " << trial << std::endl;
/*
		IGORdata::write_itx(getLam(), "output/debug/V/lam.itx", "lam");
		IGORdata::write_itx(f, "output/debug/V/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/V/v.itx", "vv");
		Eigen::VectorXd Af = A*f.array().exp().matrix();
		IGORdata::write_itx(Af, "output/debug/V/Af.itx", "Af");
		
		IGORdata::write_itx(dv, "output/debug/V/dv.itx", "dv");
		IGORdata::write_itx(J, "output/debug/V/J.itx", "JJ");
		IGORdata::write_itx(H, "output/debug/V/H.itx", "HH");
		IGORdata::write_itx(Jv, "output/debug/V/Jv.itx", "Jv");
		IGORdata::write_itx(Hvv, "output/debug/V/Hvv.itx", "Hvv");
		IGORdata::write_itx(Hvf, "output/debug/V/Hvf.itx", "Hvf");
		IGORdata::write_itx(Aij, "output/debug/V/Aij.itx", "Aij");

		IGORdata::write_itx(dv, "output/debug/V/dv.itx", "dv");
		*/

		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(getAf(), "output/debug/F/Af.itx", "Af");
		IGORdata::write_itx(getRes(), "output/debug/F/res.itx", "res");
		IGORdata::write_itx(getCenterAf(), "output/debug/F/centerAf.itx", "centerAf");
		IGORdata::write_itx(getLeftAf(), "output/debug/F/leftAf.itx", "leftAf");
		IGORdata::write_itx(getRightAf(), "output/debug/F/rightAf.itx", "rightAf");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		IGORdata::write_itx(getZ(), "output/debug/F/z.itx", "ZZ");
		
	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");


	std::cout << "chi2 = " << getChi2() << std::endl;






	for (size_t trial = 0; trial < trialmax; ++trial){
		if (updateFLam1(relativeError0, lam1Error)) { count(2) = 1.0*trial;  break; }
		//		updateF(relativeError0);
		updateV(vError1);
		std::cout << "trial5 = " << trial << std::endl;
		std::cout << "chi2 = " << getChi2() << std::endl;
		//		std::cout << "1st = " << (Lf * f).squaredNorm() << " 2nd =" << (f.size() - 1.0) << std::endl;
		//		std::cout << "J1 = " << J1 << " H11 =" << H11 << std::endl;
		/*
		IGORdata::write_itx(getLam(), "output/debug/L1/lam.itx", "lam");
		IGORdata::write_itx(f, "output/debug/L1/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/L1/v.itx", "vv");
		Eigen::VectorXd Af = A*f.array().exp().matrix();
		IGORdata::write_itx(Af, "output/debug/L1/Af.itx", "Af");

		IGORdata::write_itx(J, "output/debug/L1/J.itx", "JJ");
		IGORdata::write_itx(H, "output/debug/L1/H.itx", "HH");
		IGORdata::write_itx(Jv, "output/debug/L1/Jv.itx", "Jv");
		IGORdata::write_itx(Hvv, "output/debug/L1/Hvv.itx", "Hvv");
		IGORdata::write_itx(Hvf, "output/debug/L1/Hvf.itx", "Hvf");
		IGORdata::write_itx(Aij, "output/debug/L1/Aij.itx", "Aij");

		IGORdata::write_itx(dv, "output/debug/L1/dv.itx", "dv");
	*/	
		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(getAf(), "output/debug/F/Af.itx", "Af");
		IGORdata::write_itx(getRes(), "output/debug/F/res.itx", "res");
		IGORdata::write_itx(getCenterAf(), "output/debug/F/centerAf.itx", "centerAf");
		IGORdata::write_itx(getLeftAf(), "output/debug/F/leftAf.itx", "leftAf");
		IGORdata::write_itx(getRightAf(), "output/debug/F/rightAf.itx", "rightAf");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		IGORdata::write_itx(getZ(), "output/debug/F/z.itx", "ZZ");

	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");







	for (size_t trial = 0; trial < trialmax; ++trial){
		if (updateVLam2(vError1, lam2Error)) { count(3) = 1.0*trial;  break; }
		updateFLam1(relativeError0, lam1Error);
		//		updateV(vError1);
		//		updateLam1(lam1Error);
		std::cout << "trial6 = " << trial << std::endl;
		std::cout << "chi2 = " << getChi2() << std::endl;
		//		std::cout << "1st = " << (Lv * v).squaredNorm() << " 2nd =" << (v.size() - 1.0) << std::endl;
		//		std::cout << "J2 = " << J2 << " H22 =" << H22 << std::endl;
		/*
		IGORdata::write_itx(getLam(), "output/debug/L2/lam.itx", "lam");
		IGORdata::write_itx(f, "output/debug/L2/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/L2/v.itx", "vv");
		Eigen::VectorXd Af = A*f.array().exp().matrix();
		IGORdata::write_itx(Af, "output/debug/L2/Af.itx", "Af");

		IGORdata::write_itx(J, "output/debug/L2/J.itx", "JJ");
		IGORdata::write_itx(H, "output/debug/L2/H.itx", "HH");
		IGORdata::write_itx(Jv, "output/debug/L2/Jv.itx", "Jv");
		IGORdata::write_itx(Hvv, "output/debug/L2/Hvv.itx", "Hvv");
		IGORdata::write_itx(Hvf, "output/debug/L2/Hvf.itx", "Hvf");
		IGORdata::write_itx(Aij, "output/debug/L2/Aij.itx", "Aij");

		IGORdata::write_itx(dv, "output/debug/L2/dv.itx", "dv");
	*/	

		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(getAf(), "output/debug/F/Af.itx", "Af");
		IGORdata::write_itx(getRes(), "output/debug/F/res.itx", "res");
		IGORdata::write_itx(getCenterAf(), "output/debug/F/centerAf.itx", "centerAf");
		IGORdata::write_itx(getLeftAf(), "output/debug/F/leftAf.itx", "leftAf");
		IGORdata::write_itx(getRightAf(), "output/debug/F/rightAf.itx", "rightAf");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		IGORdata::write_itx(getZ(), "output/debug/F/z.itx", "ZZ");

	}
	IGORdata::write_itx(count, "output/prac/count.itx", "count");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");







	for (size_t trial = 0; trial < trialmax; ++trial){
		if (updateFandVandLamandKandX(relativeError0, vError1, lam1Error, lam2Error, kError, xError)){
			count(4) = 1.0*trial;

#ifdef DEBUG
			IGORdata::write_itx(count, "output/prac/count.itx", "count");
			IGORdata::write_itx(count, "output/debug/count.itx", "count");
#endif
			break;
		}
#ifdef DEBUG
		IGORdata::write_itx(getLam(), "output/debug/FVL/lam.itx", "lam");
		IGORdata::write_itx(f, "output/debug/FVL/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/FVL/v.itx", "vv");
		Eigen::VectorXd Af = A*f.array().exp().matrix();
		IGORdata::write_itx(Af, "output/debug/FVL/Af.itx", "Af");
#endif
		std::cout << "trial7 = " << trial << std::endl;
		/**/
		IGORdata::write_itx(getLam(), "output/debug/F/lam.itx", "lam");
		IGORdata::write_itx(getAf(), "output/debug/F/Af.itx", "Af");
		IGORdata::write_itx(getRes(), "output/debug/F/res.itx", "res");
		IGORdata::write_itx(getCenterAf(), "output/debug/F/centerAf.itx", "centerAf");
		IGORdata::write_itx(getLeftAf(), "output/debug/F/leftAf.itx", "leftAf");
		IGORdata::write_itx(getRightAf(), "output/debug/F/rightAf.itx", "rightAf");
		IGORdata::write_itx(f, "output/debug/F/f.itx", "FF");
		IGORdata::write_itx(v, "output/debug/F/v.itx", "vv");
		IGORdata::write_itx(getZ(), "output/debug/F/z.itx", "ZZ");
	
	}







	Eigen::MatrixXd cov = H.lu().inverse();
	Eigen::VectorXd ferr = cov.diagonal().block(0, 0, size_j, 1).cwiseSqrt();
	Eigen::VectorXd verr = cov.diagonal().block(size_j, 0, size_j, 1).cwiseSqrt();
	Eigen::VectorXd lamerr = cov.diagonal().block(2 * size_j, 0, 2, 1).cwiseSqrt();

	Eigen::VectorXd Af = A*f.array().exp().matrix();
	IGORdata::write_itx(Af, "output/debug/FVL/Af.itx", "Af");

	IGORdata::write_itx(getLam(), "output/debug/FVL/lam.itx", "lam");
	IGORdata::write_itx(f, "output/debug/FVL/f.itx", "FF");
	IGORdata::write_itx(v, "output/debug/FVL/v.itx", "vv");
	IGORdata::write_itx(getZ(), "output/debug/FVL/z.itx", "ZZ");


	IGORdata::write_itx(ferr, "output/debug/FVL/ferr.itx", "ferr");
	IGORdata::write_itx(verr, "output/debug/FVL/verr.itx", "verr");
	IGORdata::write_itx(lamerr, "output/debug/FVL/lamerr.itx", "lamerr");
	IGORdata::write_itx(H, "output/debug/FVL/Hess.itx", "Hess");
	IGORdata::write_itx(count, "output/debug/count.itx", "count");
	return false;
}   









bool Tikhonov3D::updateF
(const double relativeError0){
	setJandH();

	//	solve df

	Eigen::VectorXd df = -Hff.lu().solve(Jf);
	f += df;

	lam(0) = lam1;
	lam(1) = lam2;
	

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((df.array().abs() < relativeError0).all())
		
		return true;
	//	otherwise
	else return false;
}
    



bool Tikhonov3D::updateV
(const double vError1){
	setJandH();

	//	solve df
	
	Eigen::VectorXd dv = -Hvv.lu().solve(Jv);
	v += dv;
	
	lam(0) = lam1;
	lam(1) = lam2;
	


	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((dv.array().abs() < vError1).all())
		return true;
	//	otherwise
	else return false;
}



bool Tikhonov3D::updateFLam1
(const double relativeError0, const double lam1Error)
{
	setJandH();
	Eigen::VectorXd Jac(f.size() + 1);
	Jac.block(0, 0, f.size(), 1) = Jf;
	Jac[f.size()] = J1;

	Eigen::MatrixXd Hes(f.size() + 1, f.size() + 1);
	Hes.setZero();
	Hes.block(0, 0, f.size(), f.size()) = Hff;
	Hes(f.size(), f.size()) = H11;
	Hes.block(0, f.size(), f.size(), 1) = H1f;	
	Hes.block(f.size(), 0, 1, f.size()) = H1f.transpose();

	double chi2 = getChi2();
	Eigen::VectorXd dflam1 = -Hes.lu().solve(Jac);
	f += dflam1.block(0, 0, f.size(), 1);
	lam1 += dflam1[f.size()];
	
	double nu = 0.01;
	size_t count = 0;
	while (!(chi2 > getChi2())){
		f -= dflam1.block(0, 0, f.size(), 1);
		lam1 -= dflam1[f.size()];

		dflam1 = -(Hes + nu * Hes.diagonal().asDiagonal().toDenseMatrix()).lu().solve(Jac);
		f += dflam1.block(0, 0, f.size(), 1);
		lam1 += dflam1[f.size()];

		std::cout << "  LM method chi2 = " << getChi2() << std::endl;

		nu *= 3.0;
		count++;
		if (count > 5)
			break;
	}

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((dflam1.array().abs() < relativeError0).all() && count < 4)

		return true;
	//	otherwise
	else return false;
}




bool Tikhonov3D::updateVLam2
(const double vError1, const double lam2Error)
{
	setJandH();
	Eigen::VectorXd Jac(v.size() + 1);
	Jac.block(0, 0, v.size(), 1) = Jv;
	Jac[v.size()] = J2;
	Eigen::MatrixXd Hes(v.size() + 1, v.size() + 1);
	Hes.setZero();
	Hes.block(0, 0, v.size(), v.size()) = Hvv;
	Hes(v.size(), v.size()) = H22;
	Hes.block(0, v.size(), v.size(), 1) = H2v;
	Hes.block(v.size(), 0, 1, v.size()) = H2v.transpose();

	double chi2 = getChi2();
	Eigen::VectorXd dvlam2 = -Hes.lu().solve(Jac);
	v += dvlam2.block(0, 0, v.size(), 1);
	lam2 += dvlam2[v.size()];

	double nu = 0.01;
	size_t count = 0;
	std::cout << "old chi2 = " << chi2 << ", new chi2 = "<< getChi2() <<  std::endl;
	while (!(chi2 > getChi2())){
		v -= dvlam2.block(0, 0, v.size(), 1);
		lam2 -= dvlam2[v.size()];

		dvlam2 = -(Hes + nu * Hes.diagonal().asDiagonal().toDenseMatrix()).lu().solve(Jac);
		v += dvlam2.block(0, 0, v.size(), 1);
		lam2 += dvlam2[v.size()];

		std::cout << "  LM method chi2 = " << getChi2() << std::endl;

		nu *= 3.0;
		count++;
		if (count > 5)
			break;
	}

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	if ((dvlam2.array().abs() < vError1).all() && fabs(dvlam2(size_j)) < lam2Error && count < 4)

		return true;
	//	otherwise
	else return false;
}




bool Tikhonov3D::updateFandVandLamandKandX
(const double relativeError0, const double vError1
, const double lam1Error, const double lam2Error
, const double kError, const double xError){

	setJandH();
	double chi2 = getChi2();

	//	solve df
	dff.setZero();
	dff = -H.lu().solve(J);
	f += dff.block(0, 0, size_j, 1);
	v += dff.block(size_j, 0, size_j, 1);
	lam1 += dff(2 * size_j);
	lam2 += dff(2 * size_j + 1);
	
	double nu = 0.01;
	size_t count = 0;
	while (!(chi2 > getChi2())){
		f -= dff.block(0, 0, size_j, 1);
		v -= dff.block(size_j, 0, size_j, 1);
		lam1 -= dff(2 * size_j);
		lam2 -= dff(2 * size_j + 1);

		dff = -(H + nu * H.diagonal().asDiagonal().toDenseMatrix()).lu().solve(J);
		f += dff.block(0, 0, size_j, 1);
		v += dff.block(size_j, 0, size_j, 1);
		lam1 += dff(2 * size_j);
		lam2 += dff(2 * size_j + 1);

		std::cout << "  LM method chi2 = " << getChi2() << std::endl;

		nu *= 3.0;
		count++;
		if (count > 5)
			break;
	}

	//	if |df|j / |f|j < relativeError is satisfied for all the elements, finish itaration
	for (size_t j = 0; j < size_j; ++j){
		if (fabs(dff[j]) > relativeError0) { return false; }
	}
	for (size_t j = 0; j < size_j; ++j){
		if (fabs(dff[j + size_j]) > vError1){ return false; }
	}
	if (fabs(dff[2*size_j]) > lam1Error){ return false; }
	if (fabs(dff[2*size_j+1]) > lam2Error){ return false; }

	//	otherwise
	return true;
}





bool Tikhonov3D::setJandH(){
	Jf.setZero();
	Jv.setZero();
	J1 = 0.0;
	J2 = 0.0;
	Jk = 0.0;
	Jx = 0.0;

	Hff.setZero();
	Hvv.setZero();
	H11 = 0.0;
	H22 = 0.0;
	Hvf.setZero();
	H1f.setZero();
	H2f.setZero();
	H1v.setZero();
	H2v.setZero();
	H21 = 0.0;

///Ç±Ç±Ç‹Ç≈Ç®Çã
	//à»â∫ÅCç°Ç©ÇÁkÇ∆xÇÃóvëfÇì¸ÇÍÇƒÇ¢Ç≠

	J.setZero();
	H.setZero();

	Jf = exp(lam1) * Lf.transpose() * Lf * f;
	Jv = exp(lam2) * Lv.transpose() * Lv * v;
	J1 = 0.5 * (exp(lam1) * (Lf * f).squaredNorm() - size_j);
	J2 = 0.5 * (exp(lam2) * v.transpose() * Lv.transpose() * Lv * v - 30.0);

	Hff = exp(lam1) * Lf.transpose() * Lf;
	Hvv = exp(lam2) * Lv.transpose() * Lv;
	H11 = 0.5 * exp(lam1) * (Lf * f).squaredNorm();
	H22 = 0.5 * exp(lam2) * v.transpose() * Lv.transpose() * Lv * v;
	H1f = exp(lam1) * Lf.transpose() * Lf * f;
	H2v = exp(lam2) * Lv.transpose() * Lv * v;


	//	z0_j = exp(f0_j)
	Eigen::VectorXd z = f.array().exp().matrix();

	Aij.setZero();
	Aijright.setZero();
	Aijleft.setZero();

	for (size_t i = 0; i < size_i; ++i){
		double res_i = d(i);

		for (size_t l = 0; l < size_j; ++l){
			res_i -= A_zee(i, l, k, x) * z[l];
		}

		for (size_t l = 0; l < size_j; ++l){
			Jk -= res_i * Adk_zee(i, l, k, x)   * z[l] / sigma[i] / sigma[i];
			Jx -= res_i * Adx_zee(i, l, k, x)   * z[l] / sigma[i] / sigma[i];

		}

		for (size_t j = 0; j < size_j; ++j){

			Aij(i, j) = A_(i, j, v[j]);
			Aijright(i, j) = A_(i, j, v[j] - x);
			Aijleft(i, j) = A_(i, j, v[j] + x);
			           
			Jf[j] -= res_i * A_zee(i, j, k, x)     * z[j] / sigma[i] / sigma[i];
			Jv[j] -= res_i * Adv_zee(i, j, k, x)   * z[j] / sigma[i] / sigma[i];
			
			for (size_t m = 0; m < size_j; ++m){
				Hff(j, m) += A_zee(i, j, k, x)   * z[j] * A_zee(i, m, k, x)   * z[m] / sigma[i] / sigma[i];
				Hvf(j, m) += Adv_zee(i, j, k, x) * z[j] * A_zee(i, m, k, x)   * z[m] / sigma[i] / sigma[i];
				Hvv(j, m) += Adv_zee(i, j, k, x) * z[j] * Adv_zee(i, m, k, x) * z[m] / sigma[i] / sigma[i];
			}	
		}
	}
	A = Aij;
	Aright = Aijright;
	Aleft = Aijleft;

	J.block(     0, 0, size_j, 1) = Jf;
	J.block(size_j, 0, size_j, 1) = Jv;
	J(  2*size_j) = J1;
	J(2*size_j+1) = J2;

	H.block(         0,          0, size_j, size_j) = Hff;
	H.block(    size_j,          0, size_j, size_j) = Hvf.transpose(); //Hfv
	H.block(  2*size_j,          0,      1, size_j) = H1f.transpose(); //Hf1
	H.block(         0,     size_j, size_j, size_j) = Hvf;
	H.block(    size_j,     size_j, size_j, size_j) = Hvv;
	H.block(2*size_j+1,     size_j,      1, size_j) = H2v.transpose(); //Hv2
	H.block(         0,   2*size_j, size_j,      1) = H1f;
	H.block(    size_j, 2*size_j+1, size_j,      1) = H2v;

	H(  2*size_j,   2*size_j) = H11;
	H(2*size_j+1, 2*size_j+1) = H22;

	return true;
}

double Tikhonov3D::A_zee(const size_t i, const size_t j, const double k, const double x)const{
	return 1 / (2 + exp(k))*A_(i, j, v[j]) + (0.5-0.5 / (2 + exp(k)))*A_(i, j, v[j] - x) + (0.5-0.5 / (2 + exp(k)))*A_(i, j, v[j] + x);
}

double Tikhonov3D::Adv_zee(const size_t i, const size_t j, const double k, const double x)const{
	return 1 / (2 + exp(k))*Adv_(i, j, v[j]) + (0.5-0.5 / (2 + exp(k))) * Adv_(i, j, v[j] - x) + (0.5-0.5 / (2 + exp(k))) * Adv_(i, j, v[j] + x);
}

double Tikhonov3D::Adk_zee(const size_t i, const size_t j, const double k, const double x)const{
	return -exp(k) / (2 + exp(k)) / (2 + exp(k)) * A_(i, j, v[j]) + 0.5*exp(k) / (2 + exp(k)) / (2 + exp(k))*A_(i, j, v[j] - x) + 0.5 * exp(k) / (2 + exp(k)) / (2 + exp(k))*A_(i, j, v[j] + x);
}

double Tikhonov3D::Adkdv_zee(const size_t i, const size_t j, const double k, const double x)const{
	return -exp(k) / (2 + exp(k)) / (2 + exp(k)) * Adv_(i, j, v[j]) + 0.5*exp(k) / (2 + exp(k)) / (2 + exp(k))*Adv_(i, j, v[j] - x) + 0.5 * exp(k) / (2 + exp(k)) / (2 + exp(k))*Adv_(i, j, v[j] + x);
}

double Tikhonov3D::Adx_zee(const size_t i, const size_t j, const double k, const double x)const{
	return -(0.5-0.5 / (2 + exp(k))) * Adv_(i, j, v[j] - x) + (0.5-0.5 / (2 + exp(k)))*Adv_(i, j, v[j] + x);
}

double Tikhonov3D::Adxdv_zee(const size_t i, const size_t j, const double k, const double x)const{
	return -(0.5-0.5 / (2 + exp(k))) * Advdv_(i, j, v[j] - x) + (0.5-0.5 / (2 + exp(k)))*Advdv_(i, j, v[j] + x);
}

double Tikhonov3D::Adxdx_zee(const size_t i, const size_t j, const double k, const double x)const{
	return (0.5-0.5 / (2 + exp(k))) * Advdv_(i, j, v[j] - x) + (0.5-0.5 / (2 + exp(k)))*Advdv_(i, j, v[j] + x);
}

double Tikhonov3D::Adkdk_zee(const size_t i, const size_t j, const double k, const double x)const{
	return 0.5 * exp(k) * (2 - exp(k)) / (2 + exp(k)) / (2 + exp(k)) / (2 + exp(k)) * (-2 * A_(i, j, v[j]) + A_(i, j, v[j] - x) + A_(i, j, v[j] + x));
}

double Tikhonov3D::A_(const size_t i, const size_t j, const double velocity)const{
	
	if (j == size_j - 1)
		return 1.0;
	else
		return (exp(-(wavelength(i) - lambda0 + velocity * lambda0 / vlight)
	        	* (wavelength(i) - lambda0 + velocity * lambda0 / vlight)/(width(j)*width(j)+0.0004)));
	
}

double Tikhonov3D::Adv_(const size_t i, const size_t j, const double velocity)const{
	if (j == size_j - 1)
		return 0;
	else
		return -(A_(i, j, v[j]) * 2 * lambda0 * (wavelength(i) - lambda0 + velocity * lambda0 / vlight)
			/ vlight /( width(j) * width(j) + 0.0004));
}

double Tikhonov3D::Advdv_(const size_t i, const size_t j, const double velocity)const{
	if (j == size_j - 1)
		return 0;
	else
		return A_(i, j, v[j]) * 4 * lambda0 * lambda0 * (wavelength(i) - lambda0 + velocity * lambda0 / vlight) * (wavelength(i) - lambda0 + velocity * lambda0 / vlight) / vlight / vlight / width(j) / width(j) / width(j) / width(j) 
		- A_(i, j, v[j]) * 2 * lambda0 * lambda0 / vlight / vlight / width(j) / width(j);
}





Eigen::MatrixXd Tikhonov3D::getCenterAf()const {
	Eigen::VectorXd centerAf = 1 / (2 + exp(k))*A*f.array().exp().matrix();
	return centerAf;
}

Eigen::MatrixXd Tikhonov3D::getLeftAf()const {
	Eigen::VectorXd leftAf = 0.5 * (1 + exp(k)) / (2 + exp(k))*Aleft*f.array().exp().matrix();
	return leftAf;
}

Eigen::MatrixXd Tikhonov3D::getRightAf()const {
	Eigen::VectorXd rightAf = 0.5 * (1 + exp(k)) / (2 + exp(k))*Aright*f.array().exp().matrix();
	return rightAf;
}

Eigen::MatrixXd Tikhonov3D::getAf()const{
	Eigen::VectorXd Af = getCenterAf() + getLeftAf() + getRightAf();
	return Af;
}

Eigen::VectorXd Tikhonov3D::getKX()const{
	Eigen::VectorXd KX;
	KX.resize(2);
	KX(0) = (1+exp(k))/(2+exp(k));
	KX(1) = x;
	return KX;
}

void Tikhonov3D::set_Lf_Amp(){
	Lf.resize(size_j, size_j);
	Lf.setIdentity();
}
void Tikhonov3D::set_Lf_Deriv(){
	Lf.resize(size_j, size_j);
	Lf.setZero();
//	Lf(0, 0) = 2.0;
	for (size_t j = 0; j < size_j - 1; ++j){
		Lf(j+1, j    ) = -1.0;
		Lf(j+1, j + 1) = 1.0;
	}
}
void Tikhonov3D::set_Lf_2ndDeriv(){
	Lf.resize(size_j - 1, size_j);
	Lf.setZero();

//	Lf(0, 0) = 2.0;
//	Lf(0, 1) = -2.0;
	for (size_t j = 0; j < size_j - 2; ++j){
		Lf(j + 1, j) = 1.0;
		Lf(j + 1, j + 1) = -2.0;
		Lf(j + 1, j + 2) = 1.0;
	}
}
void Tikhonov3D::set_Lv_Amp(){
	Lv.resize(size_j, size_j);
	Lv.setIdentity();
}
void Tikhonov3D::set_Lv_Deriv(){
	Lv.resize(size_j, size_j);
	Lv.setZero();
	Lv(0, 0) = 2.0;
	for (size_t j = 0; j < size_j - 1; ++j){
		Lv(j+1, j) = -1.0;
		Lv(j+1, j + 1) = 1.0;
	}
}
void Tikhonov3D::set_Lv_2ndDeriv(){
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
