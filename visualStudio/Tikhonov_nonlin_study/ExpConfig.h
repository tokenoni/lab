#pragma once
#include <Eigen/eigen>
#include <MyMiscLib/ParameterFileImport.hpp>
#include <IGOR/IGORitx.hpp>

class ExpConfig
{
public:
	//	constructor
	ExpConfig(mymisclib::ParameterFileImport& params);

	//	destructor
	~ExpConfig();

	//	prepare everything
	void initialize();

	//	get functions
	Eigen::MatrixXd getA()const{ return A; }

	//get sigma
	Eigen::VectorXd getSigma(Eigen::VectorXd sigma_original);

	//	parameters
	const size_t size_i;	//	number of wavelength
	const size_t size_j;	//	number of parameters
	const size_t size_temp;	//	number of temperature
	const double lambda0, vlight, kb, mass;	//	other values
	//  lambda parameters
	const size_t numberOfLambda;
	const double lambdaStart, lambdaRatio;
	//  temperature parameters
	const double tempMin;
	const double tempMax;
	//  sigma change
	const size_t sigmaStart;
	const size_t sigmaEnd;
	const double sigmaValue;

protected:
	std::vector<double> radius;
	//	length matrix
	Eigen::MatrixXd A;

private:
	//	default constructer is not used
	ExpConfig();
	
};

