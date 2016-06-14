#include <MyGsl\interp.hpp>
#include <MyGsl\log_interp.hpp>
#include <MyGsl\bspline.hpp>
#include <MyGsl\Multifit.hpp>
#include <IGOR\IGORitx.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <MyGsl\voigt.hpp>

void test2dinterp();
void testBspline();
void testCsplineClamped();
void testlog_interp();

int main(){
//	test2dinterp();
//	testBspline();
	testlog_interp();
	testCsplineClamped();
}
void testCsplineClamped(){
	std::string dir = "clamped_spline_dir/";
	std::vector<double> x(200);
	std::vector<double> y(200);
	std::vector<double> y_theory(200);
	double width = 0.1;
	for(size_t i=0; i<x.size();++i){
		x[i] = -1.0 + 2.0/(x.size()-1)*i;
		y[i] = mygsl::lorentz(x[i], 0.0, 1.0*width);
		y_theory[i] = mygsl::lorentz(x[i], 0.0, 3.0*width);
	}
	std::vector<double> x1(100);
	std::vector<double> y1(100);
	for(size_t i=0; i<x1.size();++i){
		x1[i] = -1.0 + 2.0/(x1.size()-1)*i;
		y1[i] = mygsl::lorentz(-x1[i],0.0,2.0*width);
	}

	IGORdata::write_itx(x, y, dir+"data.itx", "xx", "yy");
	mygsl::interp data,data1;
	data.set(  x,  y, mygsl::interp::CsplineFlatEdge);
	data1.set(x1, y1, mygsl::interp::CsplineFlatEdge);
	IGORdata::write_itx(data.get_original_x(), data.get_original_y(), dir+"data0.itx", "x0", "y0");
	IGORdata::write_itx(data1.get_original_x(), data1.get_original_y(), dir+"data1.itx", "x1", "y1");
	IGORdata::write_itx(y_theory, dir+"y_theory.itx", "y_theory");
	
	mygsl::interp conv0 = data.convolute(data1);
	IGORdata::write_itx(conv0.get_original_x(), conv0.get_original_y(), dir+"conv0.itx", "x0_conv", "y0_conv");
	mygsl::interp conv1 = data1.convolute(data);
	IGORdata::write_itx(conv1.get_original_x(), conv1.get_original_y(), dir+"conv1.itx", "x1_conv", "y1_conv");
}

void testlog_interp(){
	
};

void testBspline(){
	const size_t N = 100;
	const size_t NCOEFFS = 12;
	const size_t NBREAK =10;

	size_t i, j;
	gsl_rng *r;
	gsl_rng_env_setup();
	r=gsl_rng_alloc(gsl_rng_default);
		
	gsl_vector *x = gsl_vector_calloc(N);
	gsl_vector *y = gsl_vector_calloc(N);
	
	mygsl::bspline bspline(N, NBREAK);

	for (i = 0; i < N; ++i) {
//		double xi = (15.0 / (N - 1)) * i;
//		double yi = cos(xi) * exp(-0.1 * xi);
		double xi = -5.0+10.0/(N-1)*i;
		double yi = 1.0/(0.25 + xi*xi);
		double sigma = sqrt(yi);
		double dy = gsl_ran_gaussian(r, 0.1*sigma);
		yi += dy;
		gsl_vector_set(x, i, xi);
		gsl_vector_set(y, i, yi);
//		gsl_vector_set(w, i, 1.0 / (sigma * sigma));

		bspline.set(i, xi, yi, sigma);
	}

	std::string dir = "bspline_test/";
	IGORdata::write_itx(x, dir+"testdata_x.itx", "xx");
	IGORdata::write_itx(y, dir+"testdata_y.itx", "yy");
	
	bspline.setInterp(bspline.Automatic);
	
	std::vector<double> smoothed_line(N);
	std::vector<double> deriv(N);
	for (i = 0; i < N; ++i) {
		smoothed_line[i] = bspline.get(gsl_vector_get(x,i));
		deriv[i] = bspline.getderiv(gsl_vector_get(x,i));
	}
	
	IGORdata::write_itx(smoothed_line, dir+"smoothed_line.itx", "smoothed_line");
	IGORdata::write_itx(deriv, dir+"deriv.itx", "deriv");
	
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_rng_free(r);
}

void test2dinterp(){
	size_t sizex=100, sizey=102;
	std::vector<double> x(sizex);
	std::vector<double> y(sizey);
	std::vector<std::vector<double>> data(sizex, std::vector<double>(sizey));

	for(size_t i=0; i<sizex; ++i){
		x[i] = 0.1*i;
		for(size_t j=0; j<sizey; ++j){
			y[j] = 0.1*j*sqrt(1.0*j);
			data[i][j] = x[i]*sqrt(x[i])+y[j];
		}
	}
	std::string dir = "test2dinterp/";
	IGORdata::write_edgeVector(x, dir+"x_test.itx", "xx");
	IGORdata::write_edgeVector(y, dir+"y_test.itx", "yy");
	IGORdata::write_itx(data, dir+"data.itx", "data");
	mygsl::interp2d interp;
	interp.set(x, y, data);

	//	integration parralell to x axis
	std::vector<double> xinteg1(sizex);
	std::vector<double> xinteg2(sizex);
	std::vector<double> xinteg3(sizex);
	for(size_t i=0; i<sizex; ++i){
		xinteg1[i] = interp.getIntegAlong_x(y[10], x[0], x[i]);
		xinteg2[i] = interp.getIntegAlong_x(y[11], x[0], x[i]);
		xinteg3[i] = interp.getIntegAlong_x(0.5*y[10]+0.5*y[11], x[0], x[i]);
	}
	IGORdata::write_itx(xinteg1, dir+"xinteg1.itx", "xinteg1");
	IGORdata::write_itx(xinteg2, dir+"xinteg2.itx", "xinteg2");
	IGORdata::write_itx(xinteg3, dir+"xinteg3.itx", "xinteg3");

	//	integration along an arbitrary line
	size_t sizes=98;
	std::vector<double> s(sizes);
	std::vector<double> s_x(sizes);
	std::vector<double> s_y(sizes);
	std::vector<double> data_along_s(sizes);
	std::vector<double> integ_along_s(sizes);
		integ_along_s[0]= interp.getIntegAlongLine(0.0,4.0, 101.0,120.0);
/**/	for(size_t i=0; i<sizes; ++i){
		s_x[i] = 0.01*i;
		s_y[i] = 2.0*i-1.0;
		s[i] = sqrt((s_x[i]-s_x[0])*(s_x[i]-s_x[0])+(s_y[i]-s_y[0])*(s_y[i]-s_y[0]));
		data_along_s[i] = interp.get(s_x[i], s_y[i]);
		integ_along_s[i]= interp.getIntegAlongLine(s_x[0],s_x[i], s_y[0], s_y[i]);
	}
	IGORdata::write_itx(s  , dir+"s_test.itx", "ss");
	IGORdata::write_itx(s_x, dir+"s_x_test.itx", "s_x");
	IGORdata::write_itx(s_y, dir+"s_y_test.itx", "s_y");
	IGORdata::write_itx(data_along_s, dir+"data_along_s.itx", "data_along_s");
	IGORdata::write_itx(integ_along_s, dir+"integ_along_s.itx", "integ_along_s");
}