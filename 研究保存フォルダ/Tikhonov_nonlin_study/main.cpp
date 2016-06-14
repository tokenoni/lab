#include "ExpConfig.h"
#include "Tikhonov.h"

int main(){
	mymisclib::ParameterFileImport params;
	params.readFile("input.txt");

	//	experimental configuration
	ExpConfig expConfig(params);
	IGORdata::write_itx(expConfig.getA(), "output/A.itx", "AA");

	//	experimental data
	Eigen::VectorXd d = IGORdata::itx_VectorXd(params.getValue_string("path_expData"));		//get d observed
	Eigen::VectorXd sigma_original = IGORdata::itx_VectorXd("data/sigma_original.itx");						//get sigma_original
	sigma_original[2047] = sigma_original[2046];
	Eigen::VectorXd sigma = expConfig.getSigma(sigma_original);
//	Eigen::MatrixXd alpha = (1.0 / (sigma.array() * sigma.array())).matrix().asDiagonal();	//cul alpha
	IGORdata::write_itx(sigma, "output/sigma.itx", "sigma");


	//	Tikhonov regularization
	Tikhonov tikhonov(expConfig.size_i, expConfig.size_j);		//get object tikhonov
	tikhonov.set_L_2ndDeriv();									//get L
	tikhonov.set_A(expConfig.getA());							//set A to tikhonov.h
	tikhonov.setFinf_zero();									//set finf = 0
	tikhonov.set_d(d);											//set d to tikhonov.h
//	tikhonov.set_alpha(alpha);									//set alpha to tikhonov.h


	size_t size_lambda = expConfig.numberOfLambda;
	Eigen::VectorXd norm(size_lambda);							//norm about each lambda 
	Eigen::VectorXd rough(size_lambda);							//rough about each lambda

	double lambda_start = expConfig.lambdaStart;									//minimum lambda
	Eigen::VectorXd f0(tikhonov.size_j);
	f0.setOnes();
	f0 *= 1.0;
	tikhonov.setSigma(sigma);
	tikhonov.setF0(f0);
	tikhonov.solve(10.0, 1.0e-3);								//solve
	IGORdata::write_itx(tikhonov.getA(), "A.itx", "AA");		//output A
	IGORdata::write_itx(tikhonov.getF(), "f.itx", "FF");		//output f culculated
	IGORdata::write_itx(tikhonov.getLF(), "Lf.itx", "Lf");		//output Lf culculated
	IGORdata::write_itx(tikhonov.getAf(), "Af.itx", "Af");		//output Af culculated


	for (size_t i = 0; i < size_lambda; i++){					//
		double lambda = lambda_start * pow(expConfig.lambdaRatio, i);				//define lambda
		Eigen::VectorXd f0(tikhonov.size_j);
		f0.setOnes();
		f0 *= 0.1;
		tikhonov.setF0(f0);
		tikhonov.solve(lambda, 1.0e-3);								//solve
		norm[i] = tikhonov.getNorm();							//get resNorm
		rough[i] = tikhonov.getRough();							//get LfNorm

		std::stringstream ss;									//make string ss
		ss << "output/" << i;									//string ss = output/i
		//ss.str() = string ss
		IGORdata::write_itx(tikhonov.getA(), ss.str() + "/A.itx", "AA");		//output A
		IGORdata::write_itx(tikhonov.getF(), ss.str() + "/f.itx", "FF");		//output f culculated
		IGORdata::write_itx(tikhonov.getLF(), ss.str() + "/Lf.itx", "Lf");		//output Lf culculated
		IGORdata::write_itx(tikhonov.getAf(), ss.str() + "/Af.itx", "Af");		//output Af culculated
		IGORdata::write_itx(tikhonov.getRes(), ss.str() + "/res.itx", "res");	//output res of fit
	}

	IGORdata::write_itx(norm, "output/norm.itx", "nnorm");						//output resNorm of fit
	IGORdata::write_itx(rough, "output/rough.itx", "rrough");					//output rough of fit (LfNorm)





}