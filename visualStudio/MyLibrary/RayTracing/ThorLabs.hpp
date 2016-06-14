#ifndef __RAYTRACING_ThorLabs_HPP__
#define __RAYTRACING_ThorLabs_HPP__

#include "Optics.hpp"
namespace ThorLabs{
	static const CylindricalSingletLens CylindricalConvexLens_y30_z32_f50(30.0, 32.0,  6.8, 31.0, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y50_z53_f50(50.8, 50.8, 21.6, 26.4, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y20_z40_f30(20.0, 40.0,  5.7, 15.5, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y20_z40_f40(20.0, 40.0,  4.6, 20.7, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y30_z60_f50(30.0, 60.0,  6.8, 25.8, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y20_z40_f50(20.0, 40.0,  4.0, 25.8, 0.0, Optics::BK7); 

	static const CylindricalSingletLens CylindricalConcaveLens_y51_z53_f75(50.8, 53.0,  2.0, -38.8, 0.0, Optics::BK7); 


	//	positive meniscus lens
	static const SphericalSingletLens PositiveMeniscusLens_d25_f100(25.4, 3.6,  32.1,  82.2, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d25_f125(25.4, 3.3,  40.6, 106.9, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d25_f150(25.4, 3.1,  49.1, 131.6, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d25_f200(25.4, 2.8,  66.2, 182.2, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d25_f300(25.4, 2.5, 100.9, 288.2, Optics::BK7);

	static const SphericalSingletLens PositiveMeniscusLens_d50_f100(50.8, 9.7,  30.3,  65.8, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d50_f150(50.8, 7.3,  47.9, 119.3, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d50_f200(50.8, 6.2,  65.2, 171.6, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d50_f250(50.8, 5.5,  82.5, 224.7, Optics::BK7);
	static const SphericalSingletLens PositiveMeniscusLens_d50_f300(50.8, 5.1, 100.1, 279.1, Optics::BK7);

//	negative meniscus lens
	static const SphericalSingletLens NegativeMeniscus_dia25_f100( 25.4, 3.0, 100.0, 33.7, Optics::BK7);
	static const SphericalSingletLens NegativeMeniscus_dia25_f150( 25.4, 3.0, 100.0, 43.1, Optics::BK7);
	static const SphericalSingletLens NegativeMeniscus_dia25_f200( 25.4, 3.0, 100.0, 50.2, Optics::BK7);
	static const SphericalSingletLens NegativeMeniscus_dia25_f300( 25.4, 3.0, 250.0, 95.1, Optics::BK7);
	static const SphericalSingletLens NegativeMeniscus_dia25_f400( 25.4, 3.0, 250.0, 112.5, Optics::BK7); 
	static const SphericalSingletLens NegativeMeniscus_dia25_f500( 25.4, 3.0, 250.0, 126.3, Optics::BK7);
	static const SphericalSingletLens NegativeMeniscus_dia25_f1000(25.4, 3.0, 250.0, 253.2, Optics::BK7);

	static const SphericalSingletLens NegativeMeniscus_dia50_no0( 50.0, 4.0, -32.79, -309.45,  Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no1( 50.0, 4.0, -46.95, -446.01,  Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no2( 50.0, 4.0, -73.74, -241.63,  Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no3( 50.0, 4.0, -99.36, -320.20,  Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no4( 50.0, 4.0, -123.82, -402.58, Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no5( 50.0, 4.0, -164.03, -709.83, Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no6( 50.0, 4.0, -174.83, -571.49, Optics::SF10);
	static const SphericalSingletLens NegativeMeniscus_dia50_no7( 50.0, 4.0, -223.20, -663.82, Optics::SF5);


};

#endif

