#include <PhysicalData\Janev.hpp>
#include <IGOR\IGORitx.hpp>
#include <vector>
int main(){
	const size_t num = 201;
	std::vector<double> E(num), H2plus_ionization(num), H2plus_dissociation(num);
	std::vector<double> H2plus_dissociation_n3(num),H2plus_dissociation_n4(num),H2plus_dissociation_n5(num);
	for(size_t i=0; i<num; ++i){
		E[i] = 1.0*exp(log(1000.0/1.0)/(num-1)*i);
		H2plus_ionization[i]  = Janev::H2plus_ionization(E[i]);
		H2plus_dissociation[i] = Janev::H2plus_dissociation(E[i]);
		H2plus_dissociation_n3[i] = Janev::H2plus_dissociative_excitation_to_n3(E[i]);
		H2plus_dissociation_n4[i] = Janev::H2plus_dissociative_excitation_to_n4(E[i]);
		H2plus_dissociation_n5[i] = Janev::H2plus_dissociative_excitation_to_n5(E[i]);
	}
	std::string dir = "output/";
	IGORdata::write_itx(E                     , dir + "E.itx", "E_in_eV");
	IGORdata::write_itx(     H2plus_ionization, dir + "H2plus_ionization.itx", "H2plus_ionization");
	IGORdata::write_itx(   H2plus_dissociation, dir + "H2plus_dissociation.itx", "H2plus_dissociation");
	IGORdata::write_itx(H2plus_dissociation_n3, dir + "H2plus_dissociation_n3.itx", "H2plus_dissociation_n3");
	IGORdata::write_itx(H2plus_dissociation_n4, dir + "H2plus_dissociation_n4.itx", "H2plus_dissociation_n4");
	IGORdata::write_itx(H2plus_dissociation_n5, dir + "H2plus_dissociation_n5.itx", "H2plus_dissociation_n5");
};