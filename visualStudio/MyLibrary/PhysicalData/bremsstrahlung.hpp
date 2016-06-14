#ifndef __BREMSSTRAHLUNG_HPP__
#define __BREMSSTRAHLUNG_HPP__

#include <string>
#include <MyGsl/interp.hpp>

namespace physicalData{

	//	te [eV], lambda [nm]
	//	te > 1eV and 100 nm < lambda < 1e5 nm
	//	from 1974 Plasma Phys. 16 1187
	static double getGff(const double te, const double lambda){
		double log10te1 = log10(te);
		double log10lm1 = log10(lambda*0.001);
		double log10lm2 = log10lm1*log10lm1;
		double gff;
		if(lambda <= 1000.0)
			gff =   1.1740 + 0.0662*log10lm1 + 0.0217*log10lm2
			    + ((0.1982 + 0.4819*log10lm1 + 0.2240*log10lm2) 
				+  (0.3230 + 0.0454*log10lm1 - 0.0831*log10lm2)*log10te1)*log10te1;
		else
			gff =   1.1750 + 0.0789*log10lm1 + 0.1872*log10lm2
			    + ((0.1812 + 0.5586*log10lm1 + 0.0304*log10lm2) 
				+  (0.3305 + 0.0203*log10lm1 - 0.0657*log10lm2)*log10te1)*log10te1;
		return gff;
	}

	//	te [eV], lambda [nm], return value of PB / ne^2 [W m^3 nm^-1]
	static double getBremss(const double te, const double lambda){
		return 1.89e-35/(sqrt(te)*lambda*lambda)*exp(-1240/(te*lambda))*getGff(te, lambda);
	//	te [eV], lambda [nm], return value of PB / ne^2 [W m^3 nm^-1 sr^-1]
//		return 1.516e-36/(te*lambda*lambda)*exp(-1240/(te*lambda))*getGff(te, lambda);
	}

/*	class bremss{
	public:
		bremss(void);
		~bremss(void);

	//	te [eV], lambda [nm], return value of PB / ne^2 [W m^3 nm^-1]
		double getBremss(const double te, const double lambda)
		{
			return 0.0;
		}
	};

	class gff{
	public:
		get(const double te, const double lambda){

		}
	};
*/	
};

#endif
