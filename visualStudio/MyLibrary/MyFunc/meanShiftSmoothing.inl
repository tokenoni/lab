namespace myfunc{
	void MeanShiftSmoothing1D::set(const std::vector<double>& x_, const std::vector<double>& x_sigma, const std::vector<double>& y_, const std::vector<double>& y_sigma){
		x = x_;
		dx = x_sigma;
		y = y_;
		dy = y_sigma;
		data_num = x_.size();
	};

	double MeanShiftSmoothing1D::run(const double threashold_, Kernel kernel_){
		threashold = threashold_;
				
	}


};