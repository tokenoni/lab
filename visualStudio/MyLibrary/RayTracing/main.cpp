#include <RayTracing\Optics.hpp>
#include <RayTracing\EdmundOptics.hpp>
#include <IGOR\IGORitx.hpp>

void testAsphericLens();
void testDoubletLens();

int main(){
//	testAsphericLens();
	testDoubletLens();
	return 0;
};

void testAsphericLens(){
	AsphericalSingletLens lens = EdmundOptics::AsphericalSingletLens_dia20_f20;
	lens.setLocation(10.0,0.0,0.0,1.0,0.0,0.0);

	CreateRay createray_set;
	createray_set.ParallelSource(1000, 17.0);
	createray_set.setLocation(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	
	Screen screen;
	screen.setParameter(10.0, 10.0);
	screen.setLocation(30.1, 0.0, 0.0, 1.0, 0.0, 0.0);

	RayTrajectory ray_trajectory;
	ray_trajectory.clear();

	std::vector<Ray> rayset = createray_set.Go();
	for(size_t i=0; i<rayset.size(); ++i){
		rayset[i].wavelength = 587.6;
		createray_set.getRayTrajectory(ray_trajectory, i);

		lens.Go_Refract_Refract(rayset[i]);
		lens.getRayTrajectory(ray_trajectory);
		
		screen.Go(rayset[i]);
		screen.getRayTrajectory(ray_trajectory);
		ray_trajectory.addSeparator();
	}

	std::string dir = "AsphericalLensTest/";
	IGORdata::write_itx(lens.getShape(), dir+"lens.itx", "lens");
	IGORdata::write_itx(ray_trajectory.x, dir+"ray_trajectory.itx", "ray_trajectory");
	IGORdata::write_itx(screen.getPosition(rayset), dir+"image.itx", "image");
	
}
void testDoubletLens(){
	SphericalDoubletLens lens = EdmundOptics::LargeAchromaticLens_dia128_f1900;
	lens.setLocation(0.0,0.0,0.0,1.0,0.0,0.0);

	CreateRay createray_set;
	createray_set.ParallelSource(1000, 128.0);
	createray_set.setLocation(-30.0, 0.0, 0.0, 1.0, 0.0, 0.0);
	
	Screen screen;
	screen.setParameter(100.0, 100.0);
	screen.setLocation(1903.4, 0.0, 0.0, 1.0, 0.0, 0.0);

	RayTrajectory ray_trajectory;
	ray_trajectory.clear();

	std::vector<Ray> rayset = createray_set.Go();
	for(size_t i=0; i<rayset.size(); ++i){
		rayset[i].wavelength = 587.6;
		createray_set.getRayTrajectory(ray_trajectory, i);

		lens.Go(rayset[i]);
		lens.getRayTrajectory(ray_trajectory);
		
		screen.Go(rayset[i]);
		screen.getRayTrajectory(ray_trajectory);
		ray_trajectory.addSeparator();
	}

	std::string dir = "LargeAsphericalLensTest/";
	IGORdata::write_itx(lens.getShape(), dir+"lens.itx", "lens");
	IGORdata::write_itx(ray_trajectory.x, dir+"ray_trajectory.itx", "ray_trajectory");
	IGORdata::write_itx(screen.getPosition(rayset), dir+"image.itx", "image");
	
}
