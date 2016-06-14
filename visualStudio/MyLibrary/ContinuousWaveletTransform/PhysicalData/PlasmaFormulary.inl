namespace plasmaFormulary{

	//------------------------------------
	//		relaxation rate		
	//		NRL plasma formulary p.31
	//		ne [m^-3], te [eV] for field particle, 
	//		eng [eV] is kinetic energy of test particle
	//------------------------------------
	//	electron-electron collision for slow electron
	inline double nu_ee_slow(const double ne, const double te){
		if (te <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		double lambda_ee = Lambda_ee(ne, te);
		return 5.8e-6*pow(te, -1.5) * ne*1.0e-6 * lambda_ee;
	}
	//	electron-electron collision for fast electron
	inline double nu_ee_fast(const double ne, const double te, const double eng){
		if (te <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		double lambda_ee = Lambda_ee(ne, te);
		return 7.7e-6*pow(eng, -1.5) * ne*1.0e-6 * lambda_ee;
	}

	//	electron-ion collision for slow electron
	inline double nu_ei_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ion){
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		if (ni <= 0.0) return 0.0;
		double lambda_ei = Lambda_ei(ne, te, ti, ion);
		return 0.23 * pow(ion.mu / te, 1.5) * ni*1.0e-6 * ion.Z*ion.Z * lambda_ei;
	}
	//	electron-ion collision for fast electron
	inline double nu_ei_fast(const double ne, const double te, const double ni, const double ti, const IonProperty ion, const double eng){
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		if (ni <= 0.0) return 0.0;
		double lambda_ei = Lambda_ei(ne, te, ti, ion);
		return 3.9e-6 * pow(eng, -1.5) * ni*1.0e-6 * ion.Z*ion.Z * lambda_ei;
	}

	//	ion-electron collision frequency for slow ion
	inline double nu_ie_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ion){
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		if (ni <= 0.0) return 0.0;
		double lambda_ie = Lambda_ei(ne, te, ti, ion);
		return 1.6e-9 / ion.mu * pow(ti, -1.5) * ne*1.0e-6 * ion.Z*ion.Z * lambda_ie;
	}
	//	ion-electron collision frequency for fast ion
	inline double nu_ie_fast(const double ne, const double te, const double ni, const double ti, const IonProperty ion, const double eng){
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ni <= 0.0) return 0.0;
		double lambda_ie = Lambda_ei(ne, te, ti, ion);
		return 1.7e-4 * sqrt(ion.mu) * pow(eng, -1.5) * ne*1.0e-6 * ion.Z*ion.Z * lambda_ie;
	}

	//	ion-ion collision frequency for slow ion
	inline double nu_ii_slow(const double ne, const double te, const double ni, const double ti, const IonProperty ion,
		const double ni_test, const double ti_test, const IonProperty ion_test)
	{
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ti_test <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		if (ni <= 0.0) return 0.0;
		double lambda_ii = Lambda_ii(ni, ti, ion, ni_test, ti_test, ion_test);
		return 6.8e-8 * sqrt(ion.mu) / ion_test.mu / sqrt(1.0 + ion.mu / ion_test.mu) * pow(ti, -1.5)
			* ni*1.0e-6* ion.Z*ion.Z * ion_test.Z*ion_test.Z * lambda_ii;
	}

	//------------------------------------
	//		coulomb logarithm		
	//		NRL plasma formulary p.34
	//		ne [m^-3], te [eV], mi [kg]
	//------------------------------------
	//	thermal electron-electron collisions
	inline double Lambda_ee(const double ne, const double te){
		if (te <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		double ne_cgs = ne*1.0e-6;
		return 23.5 - log(sqrt(ne_cgs) * pow(te, -1.25)) - sqrt(1.0e-5 + (log(te) - 2.0)*(log(te) - 2.0)/16.0);
	};	//	unit less

	//	Electron-ion collision
	inline double Lambda_ei(const double ne, const double te, const double ti, const IonProperty ionProperty){
		if (te <= 0.0) return 0.0;
		if (ti <= 0.0) return 0.0;
		if (ne <= 0.0) return 0.0;
		double mi_cgs = ionProperty.mi_cgs;
		double Z = ionProperty.Z;
		double time_mi = ti * cgs:: me / mi_cgs;
		double ne_cgs = ne*1.0e-6;
		if (time_mi < te ){
			if (te < 10.0*Z*Z)
				return 23.0 - log(sqrt(ne_cgs)*Z*pow(te, -1.5));
			else
				return 24.0 - log(sqrt(ne_cgs) / te);
		}
		else{
			double mu = ionProperty.mu;
			return 30.0 - log(sqrt(ne_cgs) * pow(te, -1.5) * Z*Z / mu);
		}
	};//	unit less

	//	mixed ion-ion collisions
	static double Lambda_ii(const double ni1, const double ti1, const IonProperty ionProperty1,
		const double ni2, const double ti2, const IonProperty ionProperty2)
	{
		if (ti1 <= 0.0) return 0.0;
		if (ti2 <= 0.0) return 0.0;
		if (ni1 <= 0.0 && ni2 < 0.0) return 0.0;
		double Z1 = ionProperty1.Z;
		double mu1 = ionProperty1.mu;
		double Z2 = ionProperty2.Z;
		double mu2 = ionProperty2.mu;
		return 23.0 - log(Z1*Z2*(mu1 + mu2) / (mu1*ti2 + mu2*ti1) * sqrt(ni1*1.0e-6*Z1*Z1 / ti1 + ni2*1.0e-6*Z2*Z2 / ti2));
	};//	unit less


};