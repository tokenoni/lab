#include "ExpConfig.h"


ExpConfig::ExpConfig(mymisclib::ParameterFileImport& params_)
:
size_i(params_.getValue_size_t("numberOfWavelength")),
size_temp(params_.getValue_size_t("numberOfTemp")),
size_j(params_.getValue_size_t("numberOfTemp") + 1),
numberOfLambda(params_.getValue_size_t("numberOfLambda")),
lambda0(params_.getValue_double("lambda0")),
vlight(params_.getValue_double("vlight")),
kb(params_.getValue_double("kb")),
mass(params_.getValue_double("mass")),
lambdaStart(params_.getValue_double("lambdaStart")),
lambdaRatio(params_.getValue_double("lambdaRatio")),
tempMin(params_.getValue_double("tempMin")),
tempMax(params_.getValue_double("tempMax")),
sigmaStart(params_.getValue_size_t("sigmaStart")),
sigmaEnd(params_.getValue_size_t("sigmaEnd")),
sigmaValue(params_.getValue_double("sigmaValue"))
{
	initialize();
}


ExpConfig::~ExpConfig()
{
}

void ExpConfig::initialize()
{
//
	Eigen::VectorXd wavelength = IGORdata::itx_VectorXd("data/wavelength_o.itx");
	Eigen::VectorXd Temp;
	Temp.resize(size_j);
	Eigen::VectorXd width;
	width.resize(size_j);

	//	initializing A
	A.resize(size_i, size_j);
	A.setZero();

	for (size_t i = 0; i < size_temp; i++){
		Temp(i) = tempMin * exp(log(tempMax / 1.0) / (size_temp - 1)*i);
		width(i) = lambda0 / vlight*sqrt(2 * kb*Temp(i) / mass);
	}
	IGORdata::write_itx(Temp, "output/temp.itx", "temp");
	for (size_t i = 0; i < size_i; i++){
		for (size_t j = 0; j < size_temp; j++){
			A(i, j) = exp(-((wavelength(i) - lambda0) / width(j))*((wavelength(i) - lambda0) / width(j)));
		}
		A(i, size_temp) = 1.0;
	}
}

Eigen::VectorXd ExpConfig::getSigma(Eigen::VectorXd sigma_original)
{
	Eigen::VectorXd sigma = sigma_original;
		for (size_t i = sigmaStart; i < sigmaEnd; i++){
			sigma[i] = sigmaValue;
		}
		return sigma;
}