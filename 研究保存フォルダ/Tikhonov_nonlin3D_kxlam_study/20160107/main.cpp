#include "ExpConfig.h"
#include "Tikhonov3D.h"

int main(){
	mymisclib::ParameterFileImport params;
	params.readFile("input/15.txt");

	//	experimental configuration
	ExpConfig expConfig(params);
	IGORdata::write_itx(expConfig.getA(), "data/A.itx", "AA");

	//	experimental data
	Eigen::MatrixXd d = IGORdata::itx_MatrixXd(params.getValue_string("path_expData"));
	Eigen::VectorXd sigma = IGORdata::itx_VectorXd(params.getValue_string("sigmaData"));						//get sigma
	//	IGORdata::write_itx(sigma, "output/sigma.itx", "sigma");

	//	Tikhonov regularization
	Tikhonov3D tikhonov(
		expConfig.size_i, 
		expConfig.size_j,
		expConfig.vlight,
		expConfig.lambda0,
		expConfig.wavelength,
		expConfig.width
		);

	tikhonov.set_Lf_2ndDeriv();
	tikhonov.set_Lv_Deriv();
	tikhonov.set_A(expConfig.getA());
	tikhonov.set_d(d);
	tikhonov.setSigma(expConfig.getSigma(sigma));

	Eigen::VectorXd f0(tikhonov.size_j), v0(tikhonov.size_j);
	double lam1 = 8.3;
	double lam2 = -6.9;               //exp(5.3)=200 exp(-3.9)=0.02
	double x0 = 100;
	double k0 = 0;
	f0.setOnes();
	f0*= 5.0;
	v0.setOnes();
	v0 *= 500;
	
	IGORdata::write_itx(d, "output/debug/FVL/dd.itx", "dd");

	
	tikhonov.setF0(f0);
	tikhonov.setV0(v0);
	tikhonov.setlam0(lam1,lam2);
	tikhonov.setkandx(k0, x0);
	tikhonov.solve(1.0e-3, 0.1, 0.001, 0.001);


	IGORdata::write_itx(tikhonov.getLam(), "output/prac/lam.itx", "lam");
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
/*
 */
    
//	getchar();

/*	
	// * from here, searching the optimal lambda0 and 1 *

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
			lambda_f(i, j) = 5.0 * exp(log(2000.0 / 5.0) / (size_f - 1) * i);
			lambda_v(i, j) = 0.007 * exp(log(3/0.007) / (size_v - 1) * j);

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
/*
	mymisclib::ParameterFileImport params;
	params.readFile("input/change.txt");

	ExpConfig expConfig(params);

	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data1")), "output/change/v1.itx", "v1");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data2")), "output/change/v2.itx", "v2");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data3")), "output/change/v3.itx", "v3");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data4")), "output/change/v4.itx", "v4");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data5")), "output/change/v5.itx", "v5");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data6")), "output/change/v6.itx", "v6");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data7")), "output/change/v7.itx", "v7");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data8")), "output/change/v8.itx", "v8");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data9")), "output/change/v9.itx", "v9");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data10")), "output/change/v10.itx", "v10");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data11")), "output/change/v11.itx", "v11");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data12")), "output/change/v12.itx", "v12");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data13")), "output/change/v13.itx", "v13");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data14")), "output/change/v14.itx", "v14");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data15")), "output/change/v15.itx", "v15");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data16")), "output/change/v16.itx", "v16");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data17")), "output/change/v17.itx", "v17");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data18")), "output/change/v18.itx", "v18");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data19")), "output/change/v19.itx", "v19");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data20")), "output/change/v20.itx", "v20");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data21")), "output/change/v21.itx", "v21");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("v_data22")), "output/change/v22.itx", "v22");

	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data1")), "output/change/verr1.itx", "verr1");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data2")), "output/change/verr2.itx", "verr2");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data3")), "output/change/verr3.itx", "verr3");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data4")), "output/change/verr4.itx", "verr4");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data5")), "output/change/verr5.itx", "verr5");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data6")), "output/change/verr6.itx", "verr6");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data7")), "output/change/verr7.itx", "verr7");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data8")), "output/change/verr8.itx", "verr8");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data9")), "output/change/verr9.itx", "verr9");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data10")), "output/change/verr10.itx", "verr10");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data11")), "output/change/verr11.itx", "verr11");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data12")), "output/change/verr12.itx", "verr12");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data13")), "output/change/verr13.itx", "verr13");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data14")), "output/change/verr14.itx", "verr14");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data15")), "output/change/verr15.itx", "verr15");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data16")), "output/change/verr16.itx", "verr16");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data17")), "output/change/verr17.itx", "verr17");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data18")), "output/change/verr18.itx", "verr18");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data19")), "output/change/verr19.itx", "verr19");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data20")), "output/change/verr20.itx", "verr20");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data21")), "output/change/verr21.itx", "verr21");
	IGORdata::write_itx(IGORdata::itx_MatrixXd(params.getValue_string("verr_data22")), "output/change/verr22.itx", "verr22");
*/
 }