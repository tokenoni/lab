namespace mygsl{
	//---	constructor no.1	---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, size_t NumData_){
		std::vector<std::string> params_name_(0);
		std::vector<Fixed_NotFixed> FixedParams_(params.size());
		for(size_t i=0; i<params.size(); i++) FixedParams_[i] = NotFixed;
		Process(function, params, params_name_, NumData_, FixedParams_);
	}
	//---	constructor no.2	---
	template< class TyF_ >
	inline Multifit::	Multifit( TyF_& function, std::vector<double>& params, std::vector<double>& sigma_){
		std::vector<std::string> params_name_(0);
		std::vector<Fixed_NotFixed> FixedParams_(params.size());
		for(size_t i=0; i<params.size(); i++) FixedParams_[i] = NotFixed;
		Process(function, params, params_name_, sigma_, FixedParams_);
	}
	//---	constructor no.3	---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData_){
		std::vector<Fixed_NotFixed> FixedParams_(params.size());
		for(size_t i=0; i<params.size(); i++) FixedParams_[i] = NotFixed;
		Process(function, params, params_name_, NumData_, FixedParams_);
	}
	//---	constructor no.4---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_){
		std::vector<Fixed_NotFixed> FixedParams_(params.size());
		for(size_t i=0; i<params.size(); i++) FixedParams_[i] = NotFixed;
		Process(function, params, params_name_, sigma_, FixedParams_);
	}
	//---	constructor no.5---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, size_t NumData_, std::vector<Fixed_NotFixed> FixedParams_){
		std::vector<std::string> params_name_(0);
		Process(function, params, params_name_, NumData_, FixedParams_);
	}
	//---	constructor no.6---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_){
		std::vector<std::string> params_name_(0);
		Process(function, params, params_name_, sigma_, FixedParams_);
	}
	//---	constructor no.7---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData_, std::vector<Fixed_NotFixed> FixedParams_){
		Process(function, params, params_name_, NumData_, FixedParams_);
	};
	//---	constructor no.8---
	template< class TyF_ >
	inline Multifit::Multifit( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_){
		Process(function, params, params_name_, sigma_, FixedParams_);
	}

	template< class TyF_ >
	void inline Multifit::Process( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, size_t NumData_, std::vector<Fixed_NotFixed> FixedParams_){
		params_name = params_name_;
		NumData = NumData_;
		sigma.resize(NumData);
		for(size_t i=0; i<sigma.size(); i++) sigma[i] = 1.0;
		MAXerr = 1.0e-25; 	MAXderr= 1.0e-4; 	MAXiter=100;
		FixedParams = FixedParams_;
		params_err.resize(params.size());
		params_result.resize(params.size());
		CoreMultifit(function, params);
	};
	template< class TyF_ >
	void inline Multifit::Process( TyF_& function, std::vector<double>& params, std::vector<std::string>& params_name_, std::vector<double>& sigma_, std::vector<Fixed_NotFixed> FixedParams_){
		params_name = params_name_;
		NumData = sigma_.size();
		sigma = sigma_;
		MAXerr = 1.0e-25; 	MAXderr= 1.0e-4; 	MAXiter=100;
		FixedParams = FixedParams_;
		params_err.resize(params.size());
		params_result.resize(params.size());
		CoreMultifit(function, params);
	}


	template< class TyF_ >
	inline void Multifit::CoreMultifit(TyF_& function, std::vector<double>& params){

		std::vector<double> data(NumData);
		FuncParamData<TyF_> FuncParamData(function, params, sigma, FixedParams);

		gsl_multifit_function_fdf f;		 //******	definition the classure for multifitting
		f.f   = &Multifit::calc_function_f< TyF_ >;
		f.df  = &Multifit::calc_function_df< TyF_ >;
		f.fdf = &Multifit::calc_function_fdf< TyF_ >;
		f.n =  data.size();	//	number of data points to be fitted
		f.p =  params.size();
		f.params = (void *)& FuncParamData;	

		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;
		int status=1;
		gsl_vector* params_gsl = std2gsl_vector(params);
		T = gsl_multifit_fdfsolver_lmsder;							// fitting method
		s = gsl_multifit_fdfsolver_alloc(T, f.n, f.p); 	// 計算に用いるインスタンスの定義
		gsl_multifit_fdfsolver_set(s, &f, params_gsl);

		size_t iter=0;
		do{
			iter++;
			status = gsl_multifit_fdfsolver_iterate(s);
			std::cerr<<"status = "<<gsl_strerror (status)<<std::endl;
			print_state(iter, s);
			if (status) break;
			status = gsl_multifit_test_delta(s->dx, s->x, MAXerr, MAXderr); // deriving parameters are kept as s->x
			//--		OUTPUT for check		----
			if(iter%10==0) 
				gsl2std_vector(s->x, params_result);
			//--------------------------------------
		}while(status == GSL_CONTINUE && iter < MAXiter);

		//----	error estimation	---
		gsl_matrix* params_cov = gsl_matrix_calloc(params.size(), params.size());
		gsl_multifit_covar(s->J, 0.0, params_cov);
		chi = gsl_blas_dnrm2(s->f);
		double dof = f.n - f.p;
		double c   = GSL_MAX_DBL(1, chi/sqrt(dof));

		chisq = chi/dof;
		for(size_t i=0; i<params.size(); i++)
			params_err[i] = c*sqrt(gsl_matrix_get(params_cov, i, i));
		//-----------------------------

		params_result = gsl2std_vector(s->x);
		ScaledCovMatrix.resize(params.size());
		for(size_t i=0; i<params.size(); i++){
			ScaledCovMatrix[i].resize(params.size());
			params[i] = params_result[i];
			for(size_t j=0; j<params.size(); j++)
				ScaledCovMatrix[i][j] = c*c*gsl_matrix_get(params_cov, i, j)/(params_err[i]*params_err[j]);
		}

		gsl_multifit_fdfsolver_free(s);
		gsl_vector_free(params_gsl);
		gsl_matrix_free(params_cov);

		if(status==GSL_SUCCESS) std::cerr << "\n---\t calculation finished successfully \n \t\t\twith "<<iter << " iterations \t---\n" <<std::endl;
	};

	template<class TyF_ >
	inline int Multifit::calc_function_f (const gsl_vector * param_gsl, void *pFunc, gsl_vector * f){
	//--------------------------------------------------------------
	//	calculate the values of the funciton with parameter params_tmp
	//--------------------------------------------------------------
		FuncParamData< TyF_>* fpd = (FuncParamData< TyF_> *)pFunc;
		Multifit::gsl2std_vector(param_gsl, fpd->param);
		(*(fpd->function))(fpd->param, fpd->data);
	
		for (size_t i = 0; i < fpd->data.size(); i++) gsl_vector_set(f, i, fpd->data[i] / fpd->sigma[i]);
		return GSL_SUCCESS;
	}

	template<class TyF_ >
	inline int Multifit::calc_function_df (const gsl_vector * param_gsl, void *pFunc, gsl_matrix * J){	
	//--------------------------------------------------------------
	//	calculate the differential matrix (Jacobian matrix)
	//	Jij = defined by d/dpj f(i)
	//--------------------------------------------------------------
		FuncParamData< TyF_>* fpd = (FuncParamData< TyF_> *)pFunc;
		Multifit::gsl2std_vector(param_gsl, fpd->param);
		(*(fpd->function))(fpd->param, fpd->data);
	
		for(size_t j=0;j<fpd->param.size();j++)
		{
			if(fpd->FixedParams[j] == Multifit::Fixed){
				for (size_t i = 0; i < fpd->data.size(); i++) gsl_matrix_set(J,i,j,0.0);
			}
			else{
				for(size_t k=0; k<fpd->param.size(); k++)	fpd->param_dif[k] = fpd->param[k];

				double dp =  1.0e-5*fpd->param[j];		//	differential is calculated by {F(p+dp) - F(p-dp)}/(2dp)
				if(dp == 0.0 || dp != dp) dp=1.0e-7;
				fpd->param_dif[j] = fpd->param[j] + dp;				
				(*(fpd->function))(fpd->param_dif, fpd->data_plus);

				fpd->param_dif[j] = fpd->param[j] - dp;				
				(*(fpd->function))(fpd->param_dif, fpd->data_minus);

				for (size_t i = 0; i < fpd->data.size(); i++) {
					double Jij = (fpd->data_plus[i] - fpd->data_minus[i])/(2.0*dp*fpd->sigma[i]);
					gsl_matrix_set(J,i,j,Jij);
				}
			}
		}

		return GSL_SUCCESS;
	}

	template<class TyF_ >
	inline int Multifit::calc_function_fdf (const gsl_vector * params, void *pFunc, gsl_vector * f, gsl_matrix * J) {
		Multifit::calc_function_f< TyF_ >( params, pFunc, f);
		Multifit::calc_function_df< TyF_ >(params, pFunc, J);

		return GSL_SUCCESS;
	}


	inline void Multifit::print_state(size_t iter, gsl_multifit_fdfsolver * s){
		std::cerr <<"iter: "<<iter<< ",\t |f(x)|= "<< gsl_blas_dnrm2(s->f) <<std::endl;
		if(params_name.size()==s->x->size)
			for(size_t i=0;i<s->x->size;i++)
				std::cerr << "\t"<< params_name[i] + " \t= " <<gsl_vector_get(s -> x , i)<<std::endl;
		else 
			for(size_t i=0;i<s->x->size;i++)	
				std::cerr << "\tp("<< i<<")\t= " <<gsl_vector_get(s -> x , i)<<std::endl;
		return;
	}
};