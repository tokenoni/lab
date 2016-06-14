
inline void LHDretrieve_MultiChPMT::getDataFromRawData(const void* params){
	Allocate(raw_data.size());

	double ScalingCoef0 = raw_data.getParameter("ScalingCoef0");
	double ScalingCoef1 = raw_data.getParameter("ScalingCoef1");
	double ScalingCoef2 = raw_data.getParameter("ScalingCoef2");
	double ScalingCoef3 = raw_data.getParameter("ScalingCoef3");
	rate = raw_data.getParameter("rate");
	start_t = 0.0;

	for(size_t i=0; i<size(); ++i){
		double x = raw_data.getAsShort_with_BigEndianness(i);
		data[i] = ScalingCoef0 + ScalingCoef1*(1.0*x)
					+ ScalingCoef2*(1.0*x)*(1.0*x)
					+ ScalingCoef3*(1.0*x)*(1.0*x)*(1.0*x);
	}
	end_t = 1.0*size()/rate;
	isPrepared_t = true;
}
