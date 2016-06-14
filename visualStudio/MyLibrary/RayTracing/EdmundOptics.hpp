#ifndef __RAYTRACING_EDMUND_OPTICS_HPP__
#define __RAYTRACING_EDMUND_OPTICS_HPP__

#include "Optics.hpp"
namespace EdmundOptics{
//	single plano-concave lens
	static const SphericalSingletLens PlanoConcaveLens_dia50_f100(50.0, 5.0, -51.68, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConcaveLens_dia50_f125(50.0, 5.0, -64.60, 0.0, Optics::BK7);

//	single plano-convex lens
	static const SphericalSingletLens PlanoConvexLens_dia75_f500(75.0, 8.0, 258.40, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConvexLens_dia75_f250(75.0, 8.0, 129.20, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConvexLens_dia75_f200(75.0, 8.74,103.36, 0.0, Optics::BK7);
	static const SphericalSingletLens PlanoConvexLens_dia75_f150(75.0, 11.45,77.52, 0.0, Optics::BK7);

//	Spherical Negative lens
//	static const SphericalSingletLens ConcaveLens_dia

//	Spherical Doublet Lens
	static const SphericalDoubletLens AchromaticLens_dia25_f030(25.0,11.04, 3.00,  21.17,  -16.08,-118.66, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia25_f050(25.0, 9.00, 2.50,  34.53,  -21.98,-214.63, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia25_f075(25.0, 7.00, 2.50,  46.44,  -33.77, -95.94, Optics::BK7,   Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia25_f085(25.0, 5.01, 3.00,  55.68,  -38.37, -146.45,Optics::SK11,  Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia25_f100(25.0, 6.00, 2.50,  61.47,  -44.64,-129.94, Optics::BK7,   Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia25_f150(25.0, 5.70, 2.20,  91.37,  -66.21,-197.71, Optics::BK7,   Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia25_f300(25.0, 3.00, 2.50, 184.84, -134.06,-393.98, Optics::BK7,   Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia25_f089(25.4, 6.00, 2.00,  56.98,  -40.06,-162.82, Optics::BaK1,  Optics::SF8);

	static const SphericalDoubletLens AchromaticLens_dia30_f050(30.0, 11.0, 2.20,  34.81,  -22.12,-203.48, Optics::BaF10, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f060(30.0, 8.47, 2.99,  32.60,  -31.81,-799.64, Optics::SSK8, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f075(30.0, 8.40, 3.00,  40.51,  -38.68,-922.04, Optics::SSK8, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f100(30.0, 8.25, 2.80,  61.36,  -44.30,-128.90, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f125(30.0, 4.92, 2.50,  81.53,  -56.59,-220.08, Optics::SK11, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f150(30.0, 8.10, 2.60,  91.31,  -65.57,-195.87, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f200(30.0, 5.00, 2.50, 123.77,  -89.22,-259.43, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens AchromaticLens_dia30_f300(30.0, 4.50, 2.50, 184.40, -133.95,-394.96, Optics::BK7, Optics::SF5);

	static const SphericalDoubletLens AchromaticLens_dia40_f100(40.0, 10.0, 3.00, 46.73, -50.61, -341.70, Optics::BK7, Optics::SF5);

	static const SphericalDoubletLens AchromaticLens_dia50_f075(50.0, 20.00, 4.50,  51.88, -32.79, -309.45, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f100(50.0, 13.92, 4.00,  69.28, -46.95, -446.01, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f150(50.0,  9.50, 4.00,  96.85, -73.74, -241.63, Optics::BaK4,  Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f200(50.0,  9.00, 3.50, 130.48, -99.36, -320.20, Optics::BaK4,  Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f250(50.0,  9.00, 3.50, 130.48, -99.36, -320.20, Optics::BaK4,  Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f300(50.0,  9.00, 3.50, 173.11,-164.03, -709.83, Optics::BaK4,  Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f350(50.0,  8.00, 4.00, 227.16,-174.83, -571.49, Optics::BaK4,  Optics::SF10);
	static const SphericalDoubletLens AchromaticLens_dia50_f500(50.0,  8.00, 4.00, 305.74,-223.20, -663.82, Optics::BK7,   Optics::SF5);

	static const SphericalDoubletLens AchromaticLens_dia75_f200(75.0, 17.94, 6.0, 118.81, -96.37, -288.97, Optics::BK7, Optics::SF5);

//	Large precision achromatic lens
	static const SphericalDoubletLens LargeAchromaticLens_dia64_f355(   63.50, 11.94,  8.86,  203.94, -156.41,  -569.60, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens LargeAchromaticLens_dia64_f486(   63.50, 11.43, 10.16,  290.86, -197.87,  -680.72, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens LargeAchromaticLens_dia77_f850(   76.56, 13.08, 10.16,  497.21, -345.82, -344.20, -1227.86, 0.1, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens LargeAchromaticLens_dia102_f1525(102.31, 12.37,  8.76,  908.13, -626.59, -625.40, -2169.77, 0.1, Optics::BK7, Optics::SF2);
	static const SphericalDoubletLens LargeAchromaticLens_dia128_f1900(128.02, 15.42, 10.92, 1131.72, -780.87, -779.37, -2704.01, 0.1, Optics::BK7, Optics::SF2);

//	Negative Doublet Lens
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f030(25.0, 4.0, 8.00, -21.24, 17.26, 162.24, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f040(25.0, 3.0, 5.50, -27.82, 19.65, 201.68, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f050(25.0, 2.5, 4.06, -35.51, 24.58, 217.99, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f075(25.0, 2.5, 4.06, -46.07, 36.65, 108.00, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f100(25.0, 2.0, 2.57, -70.57, 45.27, 412.67, Optics::BaF10, Optics::SF10);
	static const SphericalDoubletLens NegativeAchromaticLens_dia25_f150(25.0, 2.5, 4.06, -92.30, 68.87, 204.15, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens NegativeAchromaticLens_dia40_f120(40.0, 4.0, 6.5, -75.65, 56.77, 173.40, Optics::BK7, Optics::SF5);
	static const SphericalDoubletLens NegativeAchromaticLens_dia40_f150(40.0, 4.0, 6.5, -91.23, 69.66, 211.28, Optics::BK7, Optics::SF5);


//	Cylindrical Lens
	static const CylindricalSingletLens CylindricalConvexLens_y10_z20_f100(10.0, 20.0,  3.0, 51.63, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y10_z20_f75( 10.0, 20.0,  2.0, 38.72, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y10_z20_f15( 10.0, 20.0,  3.0,  7.75, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y25_z50_f25( 25.0, 50.0,11.03, 12.96, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y25_z50_f50( 25.0, 50.0, 4.24, 25.93, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y25_z50_f75( 25.0, 50.0,10.53, 38.89, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y25_z50_f100(25.0, 50.0, 7.70, 51.85, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_y25_z50_f150(25.0, 50.0, 5.30, 77.78, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_d25_f075(25.0, 3.06, 38.89, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_d25_f100(25.0, 3.16, 64.60, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_d25_f125(25.0, 3.00, 51.85, 0.0, Optics::BK7); 
	static const CylindricalSingletLens CylindricalConvexLens_d25_f150(25.0, 3.50, 77.78, 0.0, Optics::BK7); 

	//	Negative Cylindrical Lens
	static const CylindricalSingletLens NegativeCylindricalLens_y13_z25_f013(12.5, 25.0,  2.5, - 9.81, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y13_z25_f025(12.5, 25.0,  2.5, -12.91, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y13_z25_f050(12.5, 25.0,  2.5, -25.84, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y13_z25_f075(12.5, 25.0,  2.5, -38.76, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y13_z25_f100(12.5, 25.0,  2.5, -51.68, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y25_z50_f050(25.0, 50.0,  4.0, -25.84, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y25_z50_f075(25.0, 50.0,  4.0, -38.76, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y25_z50_f100(25.0, 50.0,  4.0, -51.68, 0.0, Optics::BK7); 
	static const CylindricalSingletLens NegativeCylindricalLens_y25_z50_f150(25.0, 50.0,  4.0, -77.52, 0.0, Optics::BK7); 

	//	AsphericalLens
	static const AsphericalSingletLens AsphericalSingletLens_dia20_f20(20.0, 8.0, 8.107E-002, -6.19614E-1, 0.0, 0.0, -1.292772E-8, -1.932447E-10, 0.0, 0.0, 0.0, 0.0, Optics::BAL35);
	static const AsphericalSingletLens AsphericalSingletLens_dia25_f25(25.0, 7.5, 6.789788158609000300E-002, -1.439137081408, 0.0, 3.416666299952E-5, -6.43804434698E-9, 2.323731121564E-11, -3.519619091892E-14,  0.0, 0.0, 0.0, Optics::BAL35);
	static const AsphericalSingletLens AsphericalSingletLens_dia25_f31(25.0, 6.5, 5.431830526887999900E-002, -1.607912554258, 0.0, 2.063455432593E-5, -7.648976546507E-9, 1.117573177905E-11, -1.010058374283E-14, 0.0, 0.0, 0.0, Optics::BAL35);
	static const AsphericalSingletLens AsphericalSingletLens_dia30_f30(30.0,11.7, 5.402485143165856200E-002, -6.220625E-1,    0.0, 0.0,               -1.772239E-9,       -1.116722E-11,       0.0,                0.0, 0.0, 0.0, Optics::BAL35);
	static const AsphericalSingletLens AsphericalSingletLens_dia30_f22(30.0,14.4, 7.341604874825637100E-002, -6.250402E-1,    0.0, 0.0,                3.055889E-8,       -1.438345E-10,       0.0,				   0.0, 0.0, 0.0, Optics::BAL35);

	//	Correlction plate of spherical abberation
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_plus025lambda(  50.0, 4.0, 0.0, 0.0, 0.0, -1.109024026548E-9, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_minus025lambda( 50.0, 4.0, 0.0, 0.0, 0.0,  1.10902402653E-9,  0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_plus050lambda(  50.0, 4.0, 0.0, 0.0, 0.0, -2.218048058E-9,    0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_minus050lambda( 50.0, 4.0, 0.0, 0.0, 0.0,  2.21804805769E-9,  0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_plus100lambda(  50.0, 4.0, 0.0, 0.0, 0.0, -4.436096151321E-9, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);
	static const AsphericalSingletLens SphericalAbberationCorrectionPlate_d50_minus100lambda( 50.0, 4.0, 0.0, 0.0, 0.0,  4.436096151321E-9, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,Optics::BK7);

	//	concave mirror
	static const ConcaveMirror ConcaveMirror_dia76_f76(76.2, 76.2);
	static const ConcaveMirror ConcaveMirror_dia76_f152(152.4, 76.2);
	static const ConcaveMirror ConcaveMirror_dia76_f203(203.2, 76.2);
};
#endif
