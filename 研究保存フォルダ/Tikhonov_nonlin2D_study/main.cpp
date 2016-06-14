#include "ExpConfig.h"
#include "Tikhonov2D.h"

int main(){
	mymisclib::ParameterFileImport params;
	params.readFile("15.txt");

	//	experimental configuration
	ExpConfig expConfig(params);
	IGORdata::write_itx(expConfig.getA(), "data/A.itx", "AA");

	//	experimental data
	Eigen::MatrixXd d = IGORdata::itx_MatrixXd(params.getValue_string("path_expData"));
	Eigen::VectorXd sigma = IGORdata::itx_VectorXd(params.getValue_string("sigmaData"));
	IGORdata::write_itx(sigma, "output/sigma.itx", "sigma");

	//	Tikhonov regularization
	Tikhonov2D tikhonov(
		expConfig.size_i, 
		expConfig.size_j,
		expConfig.vlight,
		expConfig.lambda0,
		expConfig.wavelength,
		expConfig.width
		);
	// * whether the parameters of Tikhonov except for size_i and j are need or not should be checked *
	tikhonov.set_Lf_2ndDeriv();
	tikhonov.set_Lv_2ndDeriv();
	tikhonov.set_A(expConfig.getA());
	tikhonov.set_d(d);
	tikhonov.setSigma(sigma);

	// * maybe vector norm and rough should be maked here *

	Eigen::VectorXd f0(tikhonov.size_j), v0(tikhonov.size_j);
	f0.setOnes();
	f0*= 3.0;
	v0.setOnes();
	v0 *= 500;

	// * from here, calcurate with decided lambda0 and 1 
	/**/
	
	tikhonov.setF0(f0);
	tikhonov.setV0(v0);
	tikhonov.solve(1000.0, 0.07, 1.0e-3, 0.1);

	Eigen::VectorXd f = tikhonov.getF();
	Eigen::VectorXd v = tikhonov.getV();

	IGORdata::write_itx(tikhonov.getA(), "output/prac/A.itx", "AA");
	IGORdata::write_itx(tikhonov.getF(), "output/prac/f.itx", "FF");
	IGORdata::write_itx(tikhonov.getV(), "output/prac/v.itx", "vv");
	IGORdata::write_itx(tikhonov.getAf(), "output/prac/Af.itx", "Af");
	IGORdata::write_itx(tikhonov.getRes_sigma(), "output/prac/res_sigma.itx", "res_sigma");
	IGORdata::write_itx(tikhonov.getRes(), "output/prac/res.itx", "res");
	IGORdata::write_itx(tikhonov.getDf(), "output/prac/df.itx", "df");
	IGORdata::write_itx(tikhonov.getDv(), "output/prac/dv.itx", "dv");
	IGORdata::write_itx(tikhonov.getZ(), "output/prac/z.itx", "ZZ");
	
	std::cout << "res = " << tikhonov.getResNorm_sigma()
		<< ", Lf = " << tikhonov.getLfNorm()
		<< ", Lv = " << tikhonov.getLvNorm() << std::endl;

//	getchar();

	
	// * from here, searching the optimal lambda0 and 1 *
/*
	size_t size_f = 50;
	size_t size_v = 50;

	Eigen::MatrixXd lambda_f, lambda_v, resNorm, LfNorm, LvNorm;;
	resNorm.resize(size_f, size_v);
	LfNorm.resize( size_f, size_v);
	LvNorm.resize( size_f, size_v);
	lambda_f.resize(size_f, size_v);
	lambda_v.resize(size_f, size_v);

	for (size_t j = 0; j < size_v; ++j){
		for (size_t i = 0; i < size_f; ++i){
			lambda_f(i, j) = 50.0 * exp(log(2000.0 / 50.0) / (size_f - 1) * i);
			lambda_v(i, j) = 0.007 * exp(log(0.3/0.007) / (size_v - 1) * j);

			tikhonov.set_A(expConfig.getA());
			tikhonov.setF0(f0);
			tikhonov.setV0(v0);
			tikhonov.solve(lambda_f(i ,j), lambda_v(i, j), 1.0e-3, 1.0e-1);
			resNorm(i, j) = tikhonov.getResNorm_sigma();
			LfNorm( i, j)  = tikhonov.getLfNorm();
			LvNorm( i, j)  = tikhonov.getLvNorm();
			std::cout << "i, j = " << i << "," << j << std::endl;


		}
	}
	Eigen::MatrixXd resNormlog = resNorm.array().log().matrix();
	Eigen::MatrixXd LfNormlog = LfNorm.array().log().matrix();
	Eigen::MatrixXd LvNormlog = LvNorm.array().log().matrix();


	IGORdata::write_itx(resNormlog, "output/norm/resNormlog.itx", "resNormlog");
	IGORdata::write_itx(LfNormlog, "output/norm/LfNormlog.itx", "LfNormlog");
	IGORdata::write_itx(LvNormlog, "output/norm/LvNormlog.itx", "LvNormlog");
	IGORdata::write_itx(resNorm, "output/norm/resNorm.itx", "resNorm");
	IGORdata::write_itx(LfNorm, "output/norm/LfNorm.itx", "LfNorm");
	IGORdata::write_itx(LvNorm, "output/norm/LvNorm.itx", "LvNorm");
	IGORdata::write_itx(lambda_f, "output/norm/lambda_f.itx", "lambda_f");
	IGORdata::write_itx(lambda_v, "output/norm/lambda_v.itx", "lambda_v");

	
*/



}