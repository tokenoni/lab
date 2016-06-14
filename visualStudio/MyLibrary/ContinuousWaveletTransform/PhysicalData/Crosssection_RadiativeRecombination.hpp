#ifndef __CROSSSECTION_RADIATIVE_RECOMBINATION_HPP__
#define __CROSSSECTION_RADIATIVE_RECOMBINATION_HPP__

#include <MyGsl\constants.hpp>

namespace crossSection{
	//	see equation 3.18 p.46, "plasma spectroscopy" by Fujimoto, 
	//	TODO
	//	gaunt factor should be taken into account
	class RadiativeRecombination{
	public:
		RadiativeRecombination(const double Z_, const double p_);
		~RadiativeRecombination(void){};
		
		double get(const double electron_energy_in_eV)const;
	private:
		const double Z, p;
		const double coef0, coef1;
	};

	inline RadiativeRecombination::RadiativeRecombination(const double Z_, const double p_)
		:Z(Z_), p(p_),
		coef0(myconst::alpha * myconst::pi * 64 / (3.0*sqrt(3.0))
			*Z*Z*myconst::Ryd_eV / (myconst::emass * myconst::c * myconst::c)
			*(p*p*myconst::a0 / Z)*(p*p*myconst::a0 / Z)
			/ (p*p*p)),
		coef1(Z*Z*myconst::Ryd_eV/(p/p))
	{
	}

	inline double RadiativeRecombination::get(const double Ee)const{
		//	TODO
		//	gaunt factor should be taken into account
		return coef0 * 1.0 / (1.0 + Ee*coef1) * 1.0 / (Ee*coef1);
	}

};



#endif
