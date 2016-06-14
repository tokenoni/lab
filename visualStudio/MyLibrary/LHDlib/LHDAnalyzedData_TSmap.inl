
inline void LHDAnalyzedData_TSmap::getNeCalibAllTimings(){
	std::vector<double> nl_thomson_3669 =ReduceTo1dimension(Xaxis, "nl_thomson_3669");
	std::vector<double> nl_thomson_3579 =ReduceTo1dimension(Xaxis, "nl_thomson_3579");
	std::vector<double> nl_thomson_3759 =ReduceTo1dimension(Xaxis, "nl_thomson_3759");
	std::vector<double> nl_thomson_3849 =ReduceTo1dimension(Xaxis, "nl_thomson_3849");
	std::vector<double> nl_thomson_3939 =ReduceTo1dimension(Xaxis, "nl_thomson_3939");
	std::vector<double> nl_fir_3669     =ReduceTo1dimension(Xaxis, "nl_fir_3669");
	std::vector<double> nl_fir_3579     =ReduceTo1dimension(Xaxis, "nl_fir_3579");
	std::vector<double> nl_fir_3759     =ReduceTo1dimension(Xaxis, "nl_fir_3759");
	std::vector<double> nl_fir_3849     =ReduceTo1dimension(Xaxis, "nl_fir_3849");
	std::vector<double> nl_fir_3939     =ReduceTo1dimension(Xaxis, "nl_fir_3939");
	
	ne.resize( size_x(), std::vector<double>(size_y()));
	dne.resize(size_x(), std::vector<double>(size_y()));
	int ne_index  = SearchDiagname("n_e");
	int dne_index = SearchDiagname("dn_e");
	for(size_t i=0; i<nl_thomson_3579.size();++i){
		double fir_sum= nl_fir_3669[i]+nl_fir_3579[i]+nl_fir_3759[i]+nl_fir_3849[i]+nl_fir_3939[i];
		double thom_sum= nl_thomson_3669[i]+nl_thomson_3579[i]+nl_thomson_3759[i]+nl_thomson_3849[i]+nl_thomson_3939[i];

		for(size_t j=0; j<size_y(); ++j){
			ne[i][j] = fir_sum/thom_sum  * getMatrix( ne_index)[i][j];
			dne[i][j] = fir_sum/thom_sum * getMatrix(dne_index)[i][j];
		}
	}
}					
