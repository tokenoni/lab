#ifndef __MYGSL_MOLECULAR_CONSTANTS_HPP__
#define __MYGSL_MOLECULAR_CONSTANTS_HPP__

#include <string>
#include <vector>
#include "constants.hpp"

namespace myMolecularConst{
	//----definition of the energy---	
	class molState{
	public:
		enum plus_minus{plus, minus};
		enum gerade_ungerade{gerade, ungerade};
		molState(const std::string name_, const double Lambda_, const double S_, molState::plus_minus pm_, molState::gerade_ungerade gu_,
			const double Te_, const double we_, const double wexe_, const double weye_, 
			const double Be_, const double alpha_e_, const double De_, const double beta_e_, const double gamma_e_)
		{		
			name = name_; Lambda = Lambda_; S = S_; pm = pm_; gu = gu_; 
			Te=Te_; we=we_; wexe=wexe_; weye=weye_;
			Be=Be_; alpha_e=alpha_e_; De=De_; beta_e=beta_e_; gamma_e=gamma_e_;
		}
		~molState(void){};

		double getEnergy_in_cm(const double v, const double J)const{
			double Gv = (v+0.5)*(we + (v+0.5)*(-wexe + weye*(v+0.5)));
			double Bv = Be + (v+0.5)*(-alpha_e + (v+0.5)*(gamma_e));
			double Dv = De + (v+0.5)*( beta_e);
			return Gv + J*(J+1.0)*(Bv + J*(J+1.0)*(-Dv));
		};
		double getEnergy_in_eV(const double v, const double J)const{
			return myconst::cm_to_J(getEnergy_in_cm(v, J));
		}
	private:
		std::string name;
		double Te, we, wexe, weye;
		double Be, alpha_e, De, beta_e, gamma_e;
		plus_minus pm;
		gerade_ungerade gu;
		double Lambda, S;
	};


	//---	practical definitions	---
	//---H2----                 name       Lambda      S       plus_minus          gerade_ungerade
	const molState H2_X1Sigma("X1Sigma",     0.0,     0.0,     molState::plus,     molState::gerade,  
	//                          Te          we           wexe       weye
								0.0,     4401.213,     121.336,     0.0, 
	//                          Be          alpha_e      De_        beta_e_  gamma_e
								60.853,     3.062,     4.81e-2,     0.0,     0.74144);


};
#endif