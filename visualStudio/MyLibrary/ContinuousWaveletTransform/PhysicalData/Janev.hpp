#ifndef __PHYSICAL_DATA_JANEV_HPP__
#define __PHYSICAL_DATA_JANEV_HPP__

#include "JanevCalculation.hpp"
namespace Janev{
	//---	rate coefficients	---	

	//	(2.2.11)
	const JanevRateCoef H2plus_ionization
		(
		-3.746192301092e+01,	 1.559355031108e+01,	-6.693238367093e+00,	 
		 1.981700292134e+00,	-4.044820889297e-01,	5.352391623039e-02,	 
		-4.317451841436e-03,	1.918499873454e-04,	-3.591779705419e-06,	
		3.16e+00
		); 

	//	(2.2.12)
	const JanevRateCoef H2plus_dissociation
		(
		-1.781416067709e+01,	2.277799785711e+00,		-1.266868411626e+00, 
		 4.296170447419e-01,	-9.609908013189e-02,	 1.387958040699e-02, 
		-1.231349039470e-03,	 6.042383126281e-05,	-1.247521040900e-06, 
		2.00e-01
		);

	//	(2.2.13) H2+ + e- -> e- + H+ + H(n=2)
	const JanevRateCoef H2plus_dissociative_ionization_to_n2
		(
		-3.408905929046e+01,	 1.573560727511e+01,	-6.992177456733e+00,
		 1.852216261706e+00,	-3.130312806531e-01,	 3.383704123189e-02,
		-2.265770525273e-03,	 8.565603779673e-05,	-1.398131377085e-06, 
		2.00e+00
		);

	// (2.2.14)
	const JanevRateCoef H2plus_dissociative_excitation_to_n
		(
		-1.670435653561e+01,	-6.035644995682e-01,	-1.942745783445e-08,
		-2.005952284492e-07,	 2.962996104431e-08,	 2.134293274971e-08,
		-6.353973401838e-09,	 6.152557460831e-10,	-2.025361858319e-11,
		1.00e-01
		);

	const JanevRateCoef H2plus_dissociative_recombination_to_n2 = H2plus_dissociative_excitation_to_n * 0.10;
	const JanevRateCoef H2plus_dissociative_recombination_to_n3 = H2plus_dissociative_excitation_to_n * 0.45;
	const JanevRateCoef H2plus_dissociative_recombination_to_n4 = H2plus_dissociative_excitation_to_n * 0.22;
	const JanevRateCoef H2plus_dissociative_recombination_to_n5 = H2plus_dissociative_excitation_to_n * 0.12;
	const JanevRateCoef H2plus_dissociative_recombination_to_n6 = H2plus_dissociative_excitation_to_n * 0.069;



};

#endif
