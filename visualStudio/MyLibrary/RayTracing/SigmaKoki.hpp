#ifndef __RAYTRACING_SIGMA_KOKI_HPP__
#define __RAYTRACING_SIGMA_KOKI_HPP__

#include "Optics.hpp"
namespace SigmaKoki{
	static const SphericalDoubletLens AchromaticLens_d100f200( 100.0, 30.0, 7.0, 120.0,  -90.372, -290.929, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens AchromaticLens_d100f300( 100.0, 21.0, 7.0, 181.65, -130.56,   -415.2, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens AchromaticLens_d100f500( 100.0, 14.1, 7.0, 290.929,	-219.8,	-778.5, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens AchromaticLens_d100f800( 100.0, 10.4, 7.0, 442.977,	-350.0,	-1400.0,Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens AchromaticLens_d100f1000(100.0,  9.1, 7.0, 544.2,	-441.2,	-1840.0, Optics::BK7, Optics::SF2);

	static const SphericalSingletLens PlanoConcaveLens_d50f170(50.0, 3.0, -88.23, 0.0, Optics::BK7);

	static const SphericalSingletLens PlanoConvexLens_d100f300(100.0, 12.4, 138.0, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConvexLens_d100f350(100.0, 11.0, 161.0, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConvexLens_d100f400(100.0,  9.9, 184.0, 0.0, Optics::BK7);

	static const CylindricalSingletLens CylindricalConcaveLens_y30_z50_f100(30.0, 50.0, 3.3,-47.47, 0.0, Optics::BK7);
	static const CylindricalSingletLens CylindricalConcaveLens_y30_z50_f130(30.0, 50.0, 2.8,-51.9, 0.0, Optics::BK7);

	static const CylindricalSingletLens CylindricalConvexLens_y30_z50_f50(  30.0, 50.0, 7.0, 25.95, 0.0, Optics::BK7);
};

#endif