#ifndef __MEANSHIFT_SMOOTHING_HPPP__
#define __MEANSHIFT_SMOOTHING_HPPP__

#include <vector>

namespace myfunc{
	class MeanShiftSmoothing1D{
	public:
		enum Kernel{
			gaussian
		} kernel;

		void set(const std::vector<double>& x, const std::vector<double>& x_sigma, const std::vector<double>& y, const std::vector<double>& y_sigma);
		void set(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& y_sigma);
		double run(const double threashold_, Kernel kernel_);
	private:
		std::vector<double> x, dx, y, dy;
		double threashold;
		size_t data_num;
	};

};

#include "meanShiftSmoothing.inl"
#endif