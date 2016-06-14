#ifndef __PHYSICAL_DATA_JANEV_CALCULATION_HPP__
#define __PHYSICAL_DATA_JANEV_CALCULATION_HPP__
#include <cmath>

namespace Janev{
	//---	class which contains the asymptotic coefficients	---
	class JanevCrossSection{
	public:
		//---	constructors	---
		JanevCrossSection(
			const double a0, const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, const double a7, const double a8, 
			const double Eth_in_eV);	//	Eth in [eV]

		JanevCrossSection(const JanevCrossSection& obj);
		JanevCrossSection& operator = (const JanevCrossSection& obj);
		~JanevCrossSection(void){};

		//---	calculation of the cross section	---
		//			(return in [m^2])
		double getCrossSection(const double E_in_eV)const;

		//---	access to the contained data	---
		double getEth()const{return Eth;}

	private:
		double a[9];
		double Eth;
	};

	class JanevRateCoef{
	public:
		//---	constructors	---
		JanevRateCoef(void){};
		JanevRateCoef(
			const double b0, const double b1, const double b2, const double b3, const double b4, const double b5, const double b6, const double b7, const double b8,
			const double Tmin_in_eV=0.0);	//	Tmin in [eV]

		JanevRateCoef(const JanevRateCoef& obj);
		JanevRateCoef& operator = (const JanevRateCoef& obj);
		~JanevRateCoef(void){};

		//---	calculation of the rate coefficient	---
		//			(return in [m^3/s])
		double getRateCoef(const double T_in_eV)const;
		double operator()(const double T_in_eV)const{return getRateCoef(T_in_eV);}
		
		//---	access to the contained data	---
		double getTmin()const{return Tmin;}

		//---	for scaled data	---
		JanevRateCoef operator * (const double coef)const;
		JanevRateCoef& operator *= (const double coef);

	private:
		double b[9];
		double Tmin;
	};
	

//------------	definitions	------
	//---	constructor	---
	inline JanevCrossSection::JanevCrossSection(
		const double a0, const double a1, const double a2, const double a3, const double a4, const double a5, const double a6, const double a7, const double a8, 
		const double Eth_)
	{
		a[0] = a0;	a[1]=a1;	a[2]=a2;	a[3]=a3;	a[4]=a4;	a[5]=a5;	a[6]=a6;	a[7]=a7;	a[8]=a8;
		Eth = Eth_;
	}

	//---	copy constructors	---
	inline JanevCrossSection::JanevCrossSection(const JanevCrossSection& obj){
		for(size_t i=0; i<9; ++i)
			a[i] = obj.a[i];
		
		Eth = obj.Eth;
	}
	
	inline JanevCrossSection& JanevCrossSection::operator = (const JanevCrossSection& obj){
		for(size_t i=0; i<9; ++i)
			a[i] = obj.a[i];	
		Eth = obj.Eth;
		return *this;
	}

	//---	calculation of the cross section	---
	inline double JanevCrossSection::getCrossSection(const double E)const
	{
		if(E<=Eth) return 0.0;
		else{
			double sum=0.0;
			double lnE = log(E);
			double power = 1.0;
			size_t i;
			for(i=0; i<9; ++i){
				sum += a[i] * power;
				power *= lnE;
			}
			return sum * 1.0e-4;		//	returns the value in [m^2]
		}
	}


//------------	definitions	------
	//---	constructor	---
	inline JanevRateCoef::JanevRateCoef(
		const double b0, const double b1, const double b2, const double b3, const double b4, const double b5, const double b6, const double b7, const double b8,
		const double Tmin_)
	{
		b[0] = b0;	b[1]=b1;	b[2]=b2;	b[3]=b3;	b[4]=b4;	b[5]=b5;	b[6]=b6;	b[7]=b7;	b[8]=b8;
		Tmin=Tmin_;
	}

	//---	copy constructors	---
	inline JanevRateCoef::JanevRateCoef(const JanevRateCoef& obj){
		for(size_t i=0; i<9; ++i) b[i]=obj.b[i];
		Tmin = obj.Tmin;
	}
	
	inline JanevRateCoef& JanevRateCoef::operator = (const JanevRateCoef& obj){
		for(size_t i=0; i<9; ++i)
			b[i]=obj.b[i];
		Tmin = obj.Tmin;
		return *this;
	}

	//---	calculation of the rate coefficient	---
	inline double JanevRateCoef::getRateCoef(const double T)const{
		if(T<=Tmin) return 0.0;
		else{
			double sum=0.0;
			double lnT = log(T);
			double power = 1.0;
			size_t i;
			for(i=0; i<9; ++i){
				sum += b[i] * power;
				power *= lnT;
			}
			return exp(sum) * 1.0e-6;		//	returns the value in [m^3 /s]
		}
	}
	
	//---	for scaled data	---
	inline JanevRateCoef& JanevRateCoef::operator *= (const double coef){
		double ln_coef = log(coef);
		b[0] += ln_coef;
		return *this;
	}
	inline JanevRateCoef JanevRateCoef::operator * (const double coef)const{
		JanevRateCoef scaled_coef(*this);
		scaled_coef *= coef;
		return scaled_coef;
	}



};

#endif
