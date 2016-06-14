#include "../../MolecularConstants.hpp"
#include <vector>
#include <IGOR\IGORitx.hpp>
int main(){
	std::vector<double> N(5), EngXv0(5);
	for(size_t i=0; i<N.size(); ++i){
		N[i] = i;
		EngXv0[i] = myMolecularConst::H2_X1Sigma.getEnergy_in_eV(0, N[i]);
	}
	IGORdata::write_itx(N, EngXv0, "N_EngXv0.itx", "Nrot", "EngXv0");
}