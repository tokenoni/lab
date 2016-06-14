#include <IGOR\IGORitx.hpp>
#include <LHDlib\getRho.h>

int main(){
	LHDlib::RhoTable rhoTable;

	
//		rhoTable.read("C:\\Users\\KeisukeFujii\\Dropbox\\visual_studio\\LHDdata\\flx\\lhd-r375q100b000a8020.flx");
	rhoTable.read("C:\\Users\\KeisukeFujii\\Dropbox\\visual_studio\\LHDdata\\flx\\lhd-r390q100g120b000a8020.flx");

	size_t size = 100;
	std::vector<double> R_10(size), Z_10(size);
	std::vector<double> R_11(size), Z_11(size);
	std::vector<double> R_12(size), Z_12(size);
	for (size_t i = 0; i < size; ++i){
		double theta = 2.0*myconst::pi / size*i;
		rhoTable.getRZfromRhoThetaPhi(1.0, theta, myconst::pi/10.0, R_10[i], Z_10[i]);
		rhoTable.getRZfromRhoThetaPhi(1.1, theta, myconst::pi/10.0, R_11[i], Z_11[i]);
		rhoTable.getRZfromRhoThetaPhi(1.2, theta, myconst::pi/10.0, R_12[i], Z_12[i]);
	}
	IGORdata::write_itx(R_10, Z_10, "RZ_10.itx", "RR_10", "ZZ_10");
	IGORdata::write_itx(R_11, Z_11, "RZ_11.itx", "RR_11", "ZZ_11");
	IGORdata::write_itx(R_12, Z_12, "RZ_12.itx", "RR_12", "ZZ_12");

	size = 48;
	std::vector<double> rho(size), Vp(size), dVpdrho(size), dVdrho(size);
	for (size_t i = 0; i < size; ++i){
		rho[i] = 1.2 / size*i;
		Vp[i] = rhoTable.get_Vp(0.0, rho[i]);
		dVpdrho[i] = rhoTable.get_dVpdrho(rho[i]);
		dVdrho[i] = rhoTable.get_dVdrho(rho[i]);
	}

	IGORdata::write_itx(rho, Vp, "rho_Vp.itx", "rho", "Vp");
	IGORdata::write_itx(dVpdrho, dVdrho, "dVdrho.itx", "dVpdrho", "dVdrho");

	size_t xsize = 100, ysize(101);
	std::vector<double> xaxis(xsize), yaxis(ysize);
	std::vector<std::vector<double>> rhoM(xsize), theta(xsize), drhodx(xsize), drhody(xsize), drhodz(xsize);
	for (size_t i = 0; i < xsize; ++i)
		xaxis[i] = 2.5 + 2.5 / xsize * i;
	for (size_t i = 0; i < ysize; ++i)
		yaxis[i] = -1.0 + 2.0 / ysize * i;

	for (size_t i = 0; i < xsize; ++i){
		for (size_t j = 0; j < ysize; ++j){
			double rho_, theta_;
			double drhodx_[3];
			rhoTable.get_rho_drho(xaxis[i], yaxis[j], myconst::pi*0.1, &rho_, &theta_, drhodx_);
			rhoM[i].push_back(rho_);
			theta[i].push_back(theta_);
			drhodx[i].push_back(drhodx_[0]);
			drhody[i].push_back(drhodx_[1]);
			drhodz[i].push_back(drhodx_[2]);
		}
	}
	IGORdata::write_edgeVector(xaxis, "xaxis.itx", "xaxis");
	IGORdata::write_edgeVector(yaxis, "yaxis.itx", "yaxis");
	IGORdata::write_itx(rhoM, "rhoMatrix.itx", "rhoMatrix");
	IGORdata::write_itx(theta, "theta.itx", "theta");
	IGORdata::write_itx(drhodx, "drhodx.itx", "drhodx");
	IGORdata::write_itx(drhody, "drhody.itx", "drhody");
	IGORdata::write_itx(drhodz, "drhodz.itx", "drhodz");
	return 0;
}