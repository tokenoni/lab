#ifndef __MYGSL_CONSTANTS_HPP__
#define __MYGSL_CONSTANTS_HPP__

namespace myconst{
//-------------------------------------------//
//											 //
//		elementary physical constants		 //
//											 //
//-------------------------------------------//

	//	speed of light in vacuum in [m/s]
	static const double cvel = 2.99792458e8;
	static const double c    = 2.99792458e8;

	//	Boltzmann's constant for eV [J/eV]
	static const double kb_eV = 1.602176487e-19;
	//	Boltzmann's constant for eV [J/K]
	static const double kb_K  = 1.3806503e-23;	// m2 kg s-2 K-1;

	//	proton mass in [kg]
	static const double pmass = 1.672621637e-27;

	//	electron mass in [kg]
	static const double emass = 9.10938215e-31;

	//	neutral hydrogen mass in [kg]
	static const double Hmass = pmass + emass;

	//	neutral hydrogen molecule mass in [kg]
	static const double H2mass = 2.0*Hmass;

	//	planck's constant
	static const double h = 6.62606896e-34;	//	in [J s]

	//	circular constant
	static const double pi = 3.14159265358979323846264338327950288;
	static const double twopi = 2.0*pi;
	static const double sqrt_twopi = sqrt(2.0*pi);

	//	return eV value from cm^-1
	static const double cm_to_J(const double cm){
		return cm/8065.54445;
	}

	static const double eV_K = 11604.0;	// [K/eV]
	static const double K_eV = 1.0/11604.0;	// [eV/K]

	static const double e = 1.60217657e-19;	//	elementary charge in C


	static const double alpha = 7.29735257e-3;	//	[1] fine structure constant
	static const double a0 = 5.2917721092e-11;	//	[m] first Bohr radius
	static const double Ryd_eV = 13.60569253;	//	[eV] Rydnerg constant


//-------------------------------------------//
//											 //
//     		 atomic units   		         //
//											 //
//-------------------------------------------//
	static const double au_to_m  = 5.2917720859e-11;	//	au/m
	static const double au_to_m2 = au_to_m * au_to_m;
	static const double au_to_kg = 9.10938291e-31;		//	au/kg

//-------------------------------------------//
//											 //
//		  spectroscopic constants   		 //
//											 //
//-------------------------------------------//
	//	return hv in [J]
	static const double hv(const double wl_in_nm){
		return h*c/(wl_in_nm*1.0e-9);
	}

	//	A coefficients in [1/s] and wavelengths in [nm]
	
	//	hydrogen
	static const double Acoef_LymanBeta = 4.410e7;
	static const double Acoef_BalmerAlpha = 5.575e7;
	static const double Acoef_BalmerBeta = 8.419e6;
	static const double Acoef_BalmerGamma = 2.530e6;
	static const double Acoef_BalmerDelta = 9.732e5;

	static const double Lambda_BalmerAlpha = 656.280;
	static const double Lambda_BalmerBeta  = 486.132;
	static const double Lambda_BalmerGamma = 434.046;
	static const double Lambda_BalmerDelta = 410.173;

	// helium 
	static const double Lambda_He21P31S = 728.1350742;
	static const double Lmabda_He21S31P = 501.56783;
	static const double Lambda_He21P31D = 667.8151;
	static const double Lambda_He23P33S = 706.5;
	static const double Lambda_He23S33P = 388.865;
	static const double Lambda_He23P33D = 587.6;


};

#endif