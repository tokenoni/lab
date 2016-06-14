#include <IGOR\IGORitx.hpp>
#include <PhysicalData\PlasmaFormulary.hpp>
#include <MyGsl\constants.hpp>

int main(){
	size_t tnum = 200, nnum = 202;
	std::vector<double> te(tnum), ne(nnum);
	std::vector<std::vector<double>> rate_ii(tnum, std::vector<double>(nnum));
	std::vector<std::vector<double>> rate_ie(tnum, std::vector<double>(nnum));
	std::vector<std::vector<double>> rate_ei(tnum, std::vector<double>(nnum));
	std::vector<std::vector<double>> rate_ee(tnum, std::vector<double>(nnum));

	plasmaFormulary::IonProperty proton(myconst::pmass, 1);

	for (size_t i = 0; i < tnum; ++i){
		te[i] = 0.1 * exp(log(1.0e4 / 0.1) / (tnum - 1) * i);
		for (size_t j = 0; j < nnum; ++j){
			ne[j] = 1.0e17 * exp(log(1.0e22 / 1.0e17) / (nnum - 1) * j);
			rate_ii[i][j] = plasmaFormulary::nu_ii_slow(ne[j], te[i], ne[j], te[i], proton, ne[j], te[i], proton)/ne[j];
			rate_ie[i][j] = plasmaFormulary::nu_ie_slow(ne[j], te[i], ne[j], te[i], proton) / ne[j];
			rate_ei[i][j] = plasmaFormulary::nu_ei_slow(ne[j], te[i], ne[j], te[i], proton) / ne[j];
			rate_ee[i][j] = plasmaFormulary::nu_ee_slow(ne[j], te[i]) / ne[j];
		}
	}
	IGORdata::write_itx(rate_ii, "rate_coef_ii.itx", "rate_coef_ii");
	IGORdata::write_itx(rate_ie, "rate_coef_ie.itx", "rate_coef_ie");
	IGORdata::write_itx(rate_ei, "rate_coef_ei.itx", "rate_coef_ei");
	IGORdata::write_itx(rate_ee, "rate_coef_ee.itx", "rate_coef_ee");
	IGORdata::write_edgeVector(te, "te.itx", "te");
	IGORdata::write_edgeVector(ne, "ne.itx", "ne");
}