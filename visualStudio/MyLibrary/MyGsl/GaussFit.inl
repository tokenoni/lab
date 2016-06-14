namespace mygsl{
	inline GaussFit::GaussFit(const size_t dataSize)
		:size(dataSize)
	{
		params = gsl_vector_alloc(4);	//	y0, A, x0, sigma	
		cov = gsl_matrix_alloc(4, 4);
		T = gsl_multifit_fdfsolver_lmsder;
		s = gsl_multifit_fdfsolver_alloc(T, size, 4);
		result = std::vector<double>(4);
		result_err = std::vector<double>(4);
		isSigmaProvided = true;
		fixedFlag = std::vector<int>(4, 0);
		func.p = params->size;
		func.n = size;
		func.f = &GaussFit::f;
		func.df = &GaussFit::df;
		func.fdf = &GaussFit::fdf;
		func.params = (void*)this;
	}

	inline GaussFit::~GaussFit(void){
		gsl_vector_free(params);
		gsl_matrix_free(cov);
		gsl_multifit_fdfsolver_free(s);
	}

	inline int GaussFit::fit(const std::vector<double>& x_, const std::vector<double>& data_,
		const std::vector<double>& sigma, const std::vector<double>& initial_guess,
		double epsabs, double epsrel)
	{
		x = x_;
		data = data_;
		sigma_inv = std::vector<double>(size);
		for (size_t i = 0; i < size; ++i)
			sigma_inv[i] = 1.0 / sigma[i];
		for (size_t i = 0; i < params->size; ++i)
			gsl_vector_set(params, i, initial_guess[i]);

		isSigmaProvided = true;
		return fitCore(epsabs, epsrel);
	}

	inline int GaussFit::fit(const std::vector<double>& x_, const std::vector<double>& data_,
		const std::vector<double>& initial_guess,
		double epsabs, double epsrel)
	{
		x = x_;
		data = data_;
		sigma_inv = std::vector<double>(size);
		for (size_t i = 0; i < size; ++i)
			sigma_inv[i] = 1.0;
		for (size_t i = 0; i < params->size; ++i)
			gsl_vector_set(params, i, initial_guess[i]);

		isSigmaProvided = false;
		return fitCore(epsabs, epsrel);
	}
	inline void GaussFit::setFixedParameters(const std::vector<int>& fixedFlag_){
		fixedFlag = fixedFlag_;
	}

	inline std::vector<double> GaussFit::getFitResult(){
		gsl_vector *rslt_gslv = gsl_vector_alloc(size);
		std::vector<double> rslt(data.size(), 0.0);
		f(params, (void*)this, rslt_gslv);
		for (size_t i = 0; i < size; ++i){
			rslt[i] = (gsl_vector_get(rslt_gslv, i) + data[i]) / sigma_inv[i];
		}
		return rslt;
	}

	inline int GaussFit::fitCore(double epsabs, double epsrel){
		gsl_multifit_fdfsolver_set(s, &func, params);

		size_t iter = 0, status;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate(s);
	//		printf("status = %s\n", gsl_strerror(status));
			if (status) break;
			status = gsl_multifit_test_delta(s->dx, s->x, epsabs, epsrel);
		} while (status == GSL_CONTINUE && iter < 1500);
		
		gsl_multifit_covar(s->J, 0.0, cov);

		//	TODO
		//	must be implemented if sigma is not provided
		chi = gsl_blas_dnrm2(s->f);
		double dof = size - params->size;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));
		for (size_t i = 0; i < params->size; ++i){
			if (fixedFlag[i] == 0){
				result[i] = gsl_vector_get(s->x, i);
				result_err[i] = c*sqrt(gsl_matrix_get(cov, i, i));
				gsl_vector_set(params, i, gsl_vector_get(s->x, i));
			}
			else{
				result[i] = gsl_vector_get(params, i);
				result_err[i] = 0.0;
			}
		}
		return status;
	}

	//	----------	static functions	----------------

	inline int GaussFit::f(const gsl_vector*v, void* params, gsl_vector* f){
		GaussFit* p = (GaussFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A  = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_two_width_sq = 1.0 / (2.0*width*width);
		double coef = A / (myconst::sqrt_twopi*width);

		for (size_t i = 0; i < p->size; ++i){
			double dx = (p->x)[i] - x0;
			double y = y0 + coef*exp(-dx*dx*one_over_two_width_sq);
			gsl_vector_set(f, i, (y-(p->data)[i])*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	}

	inline int GaussFit::df(const gsl_vector*v, void* params, gsl_matrix* J){
		GaussFit* p = (GaussFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_width = 1.0 / width;								//	1/width
		double one_over_two_width_sq = 0.5*one_over_width*one_over_width;	//	1/(2*width*width)
		double coef_without_A = 1.0 / (myconst::sqrt_twopi*width);			//	1/(sqrt(2pi)*width)
		double coef = A*coef_without_A;					//	A/(sqrt(2pi)*width)

		for (size_t i = 0; i < p->size; ++i){
			double dx = (p->x)[i] - x0;
			double exp_part = exp(-dx*dx*one_over_two_width_sq*(p->sigma_inv)[i]);
			//	dG/dy0
			gsl_matrix_set(J, i, 0, 1.0*(p->sigma_inv)[i]);
			//	dG/dA
			gsl_matrix_set(J, i, 1, coef_without_A*exp_part*(p->sigma_inv)[i]);
			//	dG/dx0
			gsl_matrix_set(J, i, 2, coef*exp_part*dx*one_over_two_width_sq*2.0*(p->sigma_inv)[i]);
			//	dG/dwidth
			gsl_matrix_set(J, i, 3, -coef*exp_part*one_over_width                 *(p->sigma_inv)[i]
				+ coef*exp_part*dx*dx*one_over_width*one_over_width*one_over_width*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	}

	inline int GaussFit::fdf(const gsl_vector*v, void* params, gsl_vector* f, gsl_matrix* J){
		GaussFit* p = (GaussFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_width = 1.0 / width;								//	1/width
		double one_over_two_width_sq = 0.5*one_over_width*one_over_width;	//	1/(2*width*width)
		double coef_without_A = 1.0 / (myconst::sqrt_twopi*width);			//	1/(sqrt(2pi)*width)
		double coef = A*coef_without_A;					//	A/(sqrt(2pi)*width)

		for (size_t i = 0; i < p->size; ++i){
			double dx = (p->x)[i] - x0;
			double exp_part = exp(-dx*dx*one_over_two_width_sq*(p->sigma_inv)[i]);
			double y = y0 + coef*exp_part;
			//	G
			gsl_vector_set(f, i, (y - (p->data)[i])*(p->sigma_inv)[i]);
			//	dG/dy0
			gsl_matrix_set(J, i, 0, 1.0*(p->sigma_inv)[i]);
			//	dG/dA
			gsl_matrix_set(J, i, 1, coef_without_A*exp_part*(p->sigma_inv)[i]);
			//	dG/dx0
			gsl_matrix_set(J, i, 2, coef*exp_part*dx*one_over_two_width_sq*2.0*(p->sigma_inv)[i]);
			//	dG/dwidth
			gsl_matrix_set(J, i, 3, -coef*exp_part*one_over_width                                    *(p->sigma_inv)[i]
				+ coef*exp_part*dx*dx*one_over_width*one_over_width*one_over_width*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	};

	//-------------------------
	//	fit by a Airy function with
	//	y0 + A/(pi * sigma) * ( sin((x-x0)/width)/((x-x0)/width) )^2
	//
	//--------------------------
	inline AiryFit::AiryFit(const size_t dataSize) : GaussFit(dataSize){
		func.f = &AiryFit::f;
		func.df = &AiryFit::df;
		func.fdf = &AiryFit::fdf;
	}
	
	inline std::vector<double> AiryFit::getFitResult(){
		gsl_vector *rslt_gslv = gsl_vector_alloc(size);
		std::vector<double> rslt(data.size(), 0.0);
		f(params, (void*)this, rslt_gslv);
		for (size_t i = 0; i < size; ++i){
			rslt[i] = (gsl_vector_get(rslt_gslv, i) + data[i]) / sigma_inv[i];
		}
		return rslt;
	}

	inline int AiryFit::f(const gsl_vector*v, void* params, gsl_vector* f){
		AiryFit* p = (AiryFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_width = 1.0 / width;	//	1/width
		double coef = A / (myconst::pi*width);

		for (size_t i = 0; i < p->size; ++i){
			double dx_over_width = ((p->x)[i] - x0) * one_over_width;
			double sinc = 1.0;
			if (dx_over_width != 0.0){
				sinc = sin(dx_over_width) / dx_over_width;
			}
			double y = y0 + coef*sinc*sinc;
			gsl_vector_set(f, i, (y - (p->data)[i])*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	}

	inline int AiryFit::df(const gsl_vector*v, void* params, gsl_matrix* J){
		AiryFit* p = (AiryFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_width = 1.0 / width;	//	1/width
		double coef_without_A = 1.0 / (myconst::pi*width);
		double coef = A * coef_without_A;

		for (size_t i = 0; i < p->size; ++i){
			double dx_over_width = ((p->x)[i] - x0) * one_over_width;
			double sinc = 1.0;
			double dsinc_ddx = 0.0;
			if (dx_over_width != 0.0){
				sinc = sin(dx_over_width) / dx_over_width;
				dsinc_ddx = (-sinc + cos(dx_over_width))/dx_over_width;
			}
			double y = y0 + coef*sinc*sinc;
			//	df/fy
			gsl_matrix_set(J, i, 0, 1.0*(p->sigma_inv)[i]);
			//	df/dA
			gsl_matrix_set(J, i, 1, coef_without_A*sinc*sinc*(p->sigma_inv)[i]);
			
			//	df/dx0
			gsl_matrix_set(J, i, 2, -2.0*coef*sinc*dsinc_ddx*one_over_width*(p->sigma_inv)[i]);
			//	df/dwidth
			double df_dwidth = -1.0*coef*one_over_width*sinc*sinc - 2.0*coef*sinc*dsinc_ddx*dx_over_width*one_over_width;
			gsl_matrix_set(J, i, 3, df_dwidth*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	}

	inline int AiryFit::fdf(const gsl_vector*v, void* params, gsl_vector* f, gsl_matrix* J){
		AiryFit* p = (AiryFit*)params;
		double y0 = gsl_vector_get(v, 0);
		double A = gsl_vector_get(v, 1);
		double x0 = gsl_vector_get(v, 2);
		double width = gsl_vector_get(v, 3);
		if (p->fixedFlag[0]) y0 = gsl_vector_get(p->params, 0);
		if (p->fixedFlag[1]) A = gsl_vector_get(p->params, 1);
		if (p->fixedFlag[2]) x0 = gsl_vector_get(p->params, 2);
		if (p->fixedFlag[3]) width = gsl_vector_get(p->params, 3);

		double one_over_width = 1.0 / width;	//	1/width
		double coef_without_A = 1.0 / (myconst::pi*width);
		double coef = A * coef_without_A;

		for (size_t i = 0; i < p->size; ++i){
			double dx_over_width = ((p->x)[i] - x0) * one_over_width;
			double sinc = 1.0;
			double dsinc_ddx = 0.0;
			if (dx_over_width != 0.0){
				sinc = sin(dx_over_width) / dx_over_width;
				dsinc_ddx = (-sinc + cos(dx_over_width)) / dx_over_width;
			}
			double y = y0 + coef*sinc*sinc;
			//	f
			gsl_vector_set(f, i, (y - (p->data)[i])*(p->sigma_inv)[i]);
			//	df/fy
			gsl_matrix_set(J, i, 0, 1.0*(p->sigma_inv)[i]);
			//	df/dA
			gsl_matrix_set(J, i, 1, coef_without_A*sinc*sinc*(p->sigma_inv)[i]);

			//	df/dx0
			gsl_matrix_set(J, i, 2, -2.0*coef*sinc*dsinc_ddx*one_over_width*(p->sigma_inv)[i]);
			//	df/dwidth
			double df_dwidth = -1.0*coef*one_over_width*sinc*sinc - 2.0*coef*sinc*dsinc_ddx*dx_over_width*one_over_width;
			gsl_matrix_set(J, i, 3, df_dwidth*(p->sigma_inv)[i]);
		}
		return GSL_SUCCESS;
	}
};