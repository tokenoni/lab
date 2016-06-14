#ifndef __PLASMA_FORMULARY_HPP__
#define __PLASMA_FORMULARY_HPP__

#include <string>

namespace plasmaFormulary{

	//------------------------------------
	//		physical constants in si unit		
	//		NRL plasma formulary p.16
	//------------------------------------
	namespace si{
		static const double k = 1.3807e-23;	//	Boltzmann constant [J/K]
		static const double e = 1.6022e-19;	//	elementary charge [C]
		static const double me = 9.1094e-31;	//	electron mass [kg]
		static const double mp = 1.6726e-27;	//	proton mass [kg]
		static const double h = 6.6261e-34;	//	Planck constant [J s]
		static const double hbar = 1.0546e-34;	//	h/2pi
		static const double c = 2.9979e8;		//	speed of light [m/s]
	};
	//------------------------------------
	//		physical constants in cgs unit		
	//		NRL plasma formulary p.16
	//------------------------------------
	namespace cgs{
		static const double k = 1.3807e-16;	//	Boltzmann constant [erg/deg(K)]
		static const double e = 4.8032e-10;	//	elementary charge [statcoul]
		static const double me = 9.1094e-28;	//	electron mass [g]
		static const double mp = 1.6726e-24;	//	proton mass [g]
		static const double h = 6.6261e-27;	//	Planck constant [erg sec]
		static const double hbar = 1.0546e-27;	//	h/2pi
		static const double c = 2.9979e10;		//	speed of light [cm/s]
	};

	//------------------------------------
	//	properties of the ion
	//------------------------------------
	class IonProperty{
	public:
		IonProperty(const double mi_, const double Z_)
			:mi_cgs(mi_*1.0e3), Z(Z_), mu(mi_cgs / cgs::mp)
		{};
		~IonProperty(void){};
		const double mi_cgs;	//	ion mass in g
		const double Z;		//	ion charge in the unit of the elementary charge
		const double mu;
	};

	static const IonProperty H_plus(si::mp, 1.0);
	static const IonProperty D_plus(2.0*si::mp, 1.0);
	static const IonProperty T_plus(3.0*si::mp, 1.0);
	//------------------------------------
	//		relaxation rate		
	//		NRL plasma formulary p.31
	//		ne [m^-3], te [eV] for field particle, 
	//		eng [eV] is kinetic energy of test particle
	//------------------------------------
	//	electron-electron frequency collision for slow electron
	static double nu_ee_slow(const double ne, const double te);
	//	electron-electron frequency collision for fast electron
	static double nu_ee_fast(const double ne, const double te, const double eng);

	//	electron-ion collision frequency for slow electron
	static double nu_ei_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ionProperty);
	//	electron-ion collision frequency for fast electron
	static double nu_ei_fast(const double ne, const double te, const double ni, const double ti, const IonProperty ionProperty, const double eng);

	//	ion-electron collision frequency for slow ion
	static double nu_ie_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ion);
	//	ion-electron collision frequency for fast ion
	static double nu_ie_fast(const double ne, const double te, const double ni, const double ti, const IonProperty ion, const double eng);

	//	ion-ion collision frequency for slow ion
	static double nu_ii_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ion, 
		const double ni_test, const double ti_test, const IonProperty ion_test);
	//	ion-ion collision frequency for fast ion
	static double nu_ii_fast(const double ne, const double te, const double ni, const double ti, const IonProperty ion, const double eng);


	//------------------------------------
	//		coulomb logarithm		
	//		NRL plasma formulary p.34
	//		ne [m^-3], te [eV], mi [kg]
	//------------------------------------
	//	thermal electron-electron collisions
	static double Lambda_ee(const double ne, const double te);	//	unit less
	
	//	Electron-ion collision
	static double Lambda_ei(const double ne, const double te, const double ti, const IonProperty ion);//	unit less

	//	mixed ion-ion collisions
	static double Lambda_ii(const double ni1, const double ti1, const IonProperty ionProperty1,
		const double ni2, const double ti2, const IonProperty ionProperty2);//	unit less

};
#include "PlasmaFormulary.inl"
#endif
