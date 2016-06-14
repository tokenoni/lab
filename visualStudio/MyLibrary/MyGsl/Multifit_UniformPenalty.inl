namespace mygsl{




	bool Multifit_UniformPenalty_params::initialize(
		const std::vector<double>& data_, const std::vector<double>& sigma_, 
		const std::vector<std::vector<std::vector<double>>>& phi_ijk_,
		const std::vector<std::vector<double>>& bik_, const std::vector<double>& delta2_, const std::vector<double>& lambda_)
	{
		num_i = bik_.size();
		if(num_i==0) return false;
		num_k = bik_[0].size();

		num_j = data_.size();

		if(sigma_.size() != num_j) return false;
		if(phi_ijk_.size() != num_i) return false;
		if(phi_ijk_[0].size() != num_j) return false;
		if(phi_ijk_[0][0].size() != num_k) return false;

		//---	allocating memory space	---
		//
		//	and copy the input
		//
		//

		//	data, sigma
		data = new double [num_j];
		sigma= new double [num_j];
		for(size_t j=0; j<num_j; ++j){
			data[j] = data_[j];
			sigma[j] = sigma_[j];
		}
		
		//	bik, phi_ijk
		bik = new double*[num_i];
		phi_ijk = new double** [num_i];
		for(size_t i=0; i<num_i; ++i){
			bik[i] = new double[num_k];
			phi_ijk[i] = new double* [num_j];
			for(size_t j=0; j<num_j; ++j)
				phi_ijk[i][j] = new double [num_k];
		}

		for(size_t i=0; i<num_i; ++i){
			for(size_t k=0; k<num_k; ++k)
				bik[i][k] = bik_[i][k];
			for(size_t j=0; j<num_j; ++j){
				for(size_t k=0; k<num_k; ++k)
					phi_ijk[i][j][k] = phi_ijk_[i][j][k];
			}
		}

		delta2 = new double [num_k];
		lambda = new double [num_k];
	} 

	void Multifit_UniformPenalty_params::free(){
		delete [] data;
		delete [] sigma;
		
		for(size_t i=0; i<num_i; ++i){
			for(size_t j=0; j<num_j; ++j)
				delete [] phi_ijk[i][j];
			delete [] phi_ijk[i];
			delete [] bik[i];
		}
		delete [] phi_ijk;
		delete [] bik;
		
		delete [] delta2;
		delete [] lambda;
	}
};