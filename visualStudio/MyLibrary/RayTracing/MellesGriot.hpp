#ifndef __RAYTRACING_MELLES_GRIOT_HPP__
#define __RAYTRACING_MELLES_GRIOT_HPP__

#include "Optics.hpp"
namespace MellesGriot{
	//	Cylindrical Lens
	static const CylindricalSingletLens CylindricalConvexLens_y60_z50_f300(60,50,4.0,155.0,0.0, Optics::BK7);
	static const CylindricalSingletLens CylindricalConvexLens_y60_z50_f150(60,50,6.0, 77.5,0.0, Optics::BK7);
	static const CylindricalSingletLens CylindricalConcaveLens_d51_f100(50.8, 8.0, -50.9, 0.0, Optics::BK7);

};

#endif