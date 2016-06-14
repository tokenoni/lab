namespace myEigen{
	
	//---------------------------------------------------//
	//													 //
	//				basic class for linear fit			 //
	//													 //
	//---------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	template<typename _Scalar >
	inline SparseLinearFit< _Scalar >::SparseLinearFit
		(const size_t size_time_, const size_t size_d_, const size_t size_f_)
		:
		size_time(size_time_), size_d(size_d_), size_f(size_f_),
		max_size(2000)
	{
		A = Eigen::SparseMatrix< _Scalar >(size_time*size_d, size_time*size_f);
		A.setZero();	
		f.resize(size_time * size_f);	f.setZero();
		d.resize(size_time * size_d);	d.setZero();
		res.resize(size_time * size_f);	res.setZero();
		Hes = Eigen::SparseMatrix< _Scalar >(size_time*size_f, size_time*size_f);
		Hes.setZero();
	}


	template<typename _Scalar >
	inline bool SparseLinearFit< _Scalar >::solve(){
		set_Jac_and_Hes();
		update_f();
		calculate_cov();
		return true;
	}

	template<typename _Scalar >
	inline double SparseLinearFit< _Scalar >::update_f(){
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> fprev = f;
		if (f.size() < max_size){
			Eigen::SimplicialCholesky<Eigen::SparseMatrix< _Scalar > > solver(Hes);
			f = solver.solve(Jac);
		}
		else{
			Eigen::ConjugateGradient<Eigen::SparseMatrix< _Scalar > > solver;
			solver.compute(Hes);
			f = solver.solveWithGuess(Jac, f);
		}

		return (f - fprev).norm() / f.norm();
	}

	template<typename _Scalar >
	inline bool SparseLinearFit< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * d;
		Hes = A.transpose() * A;
		return true;
	}

	template<typename _Scalar >
	inline bool SparseLinearFit< _Scalar >::calculate_cov(){
		if (f.size() < max_size){
			Eigen::SimplicialCholesky<Eigen::SparseMatrix< _Scalar > > solver(Hes);
			Eigen::SparseMatrix< _Scalar > I(size_time*size_f, size_time*size_f);
			I.setIdentity();
			cov = solver.solve(I);
		}
		else{
			Eigen::ConjugateGradient<Eigen::SparseMatrix< _Scalar > > solver;
			Eigen::SparseMatrix< _Scalar > I(size_time*size_f, size_time*size_f);
			solver.compute(Hes);
			I.setIdentity();
			cov = solver.solve(I);
		}
		return true;
	}

	template<typename _Scalar >
	inline void SparseLinearFit< _Scalar >::initialize_by_copy(const SparseLinearFit< _Scalar >& obj){
		f = obj.f;
	}
	template<typename _Scalar >
	inline Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> SparseLinearFit< _Scalar >::getError()const{
		return cov.diagonal().array().sqrt().matrix();
	}
	template<typename _Scalar >
	inline std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > SparseLinearFit< _Scalar >::get_f()const{
		std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > rslt(size_time);
		for (size_t t = 0; t < size_time; ++t){
			rslt[t] = f.block(t * size_f, 0, size_f, 1);
		}
		return rslt;
	}
	template<typename _Scalar >
	inline std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > SparseLinearFit< _Scalar >::get_f_sqerr()const{
		std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > rslt(size_time);
		for (size_t t = 0; t < size_time; ++t){
			rslt[t] = cov.diagonal().block(t * size_f, 0, size_f, 1);
		}
		return rslt;
	}

	//--------------------------------------------------//
	//													//
	//			linear fit with specified sigma			//
	//													//
	//--------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	//	the width of the Gauss distribution should be set as alpha, 
	//	where alpha = 1.0 / (sigma * sigma)
	template<typename _Scalar >
	inline SparseLinearFit_Weighted< _Scalar >::SparseLinearFit_Weighted
		(const size_t size_time_, const size_t size_d_, const size_t size_f_)
		:
		SparseLinearFit(size_time_, size_d_, size_f_)
	{
		alpha.resize(size_time * size_d);
		alpha.setOnes();
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_Weighted< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * alpha.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * A;
		return true;
	}

	//--------------------------------------------------//
	//													//
	//			linear fit with Regularization			//
	//													//
	//--------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	//	the width of the Gauss distribution should be set as alpha, 
	//	where alpha = 1.0 / (sigma * sigma)
	//	the regularization matrix L should be specified before "run"
	template<typename _Scalar >
	inline SparseLinearFit_Regularized_2D< _Scalar >::SparseLinearFit_Regularized_2D
		(const size_t size_time_, const size_t size_d_, const size_t size_f_,
		const size_t size_lx_, const size_t size_ly_)
		:
		SparseLinearFit_Weighted< _Scalar >(size_time_, size_d_, size_f_),
		size_lx(size_lx_), size_ly(size_ly_)
	{
		Lx.resize(size_time*size_lx, size_time*size_f);
		Lx.setZero();
		Ly.resize(size_time*size_ly, size_time*size_f);
		Ly.setZero();
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_Regularized_2D< _Scalar >::set_Jac_and_Hes(){
		setL();
		Jac = A.transpose() * alpha.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * A
			+ lambda_x * LLx + lambda_y * LLy;
		return true;
	}
	template<typename _Scalar >
	inline bool SparseLinearFit_Regularized_2D< _Scalar >::solve
		(const double lambda_x_, const double lambda_y_)
	{
		lambda_x = lambda_x_;
		lambda_y = lambda_y_;
		set_Jac_and_Hes();
		update_f();
		calculate_cov();
		return true;
	}

	//--------------------------------------------------//
	//													//
	//			linear fit with L1-Regularization		//
	//													//
	//--------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	//	the width of the Gauss distribution should be set as alpha, 
	//	where alpha = 1.0 / (sigma * sigma)
	//	the regularization matrix L should be specified before "run"
	template<typename _Scalar >
	inline SparseLinearFit_L1Regularized_2D< _Scalar >::SparseLinearFit_L1Regularized_2D
		(const size_t size_time_, const size_t size_d_, const size_t size_f_,
		const size_t size_lx_, const size_t size_ly_)
		:
		SparseLinearFit_Regularized_2D< _Scalar >(size_time_, size_d_, size_f_, size_lx_, size_ly_)
	{
		Lx_weight.resize(size_time * size_lx);
		Lx_weight.setOnes();
		Ly_weight.resize(size_time * size_ly);
		Ly_weight.setOnes();
		Lx_scale.resize(size_time);
		Lx_scale.setOnes();
		Ly_scale.resize(size_f);
		Ly_scale.setOnes();
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_L1Regularized_2D< _Scalar >::set_Jac_and_Hes(){
		setL();
		Jac = A.transpose() * alpha.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * A
			+ lambda_x * LLx + lambda_y * LLy;
		return true;
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_L1Regularized_2D< _Scalar >::solve
		(const double lambda_x_, const double lambda_y_,
		 const size_t iterMax, const double epsrel)
	{
		lambda_x = lambda_x_;
		lambda_y = lambda_y_;
		setL();		//	calculate L * weight * L
		set_Jac_and_Hes();
		update_f();
		for (size_t t = 0; t < iterMax; ++t){
			update_L_weight();
			set_Jac_and_Hes();
			double norm = update_f();
			if (norm < epsrel){
				calculate_cov();
				return true;
			}
		}
		calculate_cov();
		return false;
	}

	template<typename _Scalar >
	inline void SparseLinearFit_L1Regularized_2D< _Scalar >::setL(){
		LLx = Lx.transpose() * Lx_weight.asDiagonal() * Lx; 
		LLy = Ly.transpose() * Ly_weight.asDiagonal() * Ly;
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_L1Regularized_2D< _Scalar >::update_L_weight()
	{
		//	evaluate MAD for Lx
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > Lxf = Lx * f;
		Lxf.array() = Lxf.array().abs();
		for (size_t t = 0; t < size_time; ++t){
			Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > Lxf1 = Lxf.block(t * size_lx, 0, size_lx, 1);
			std::sort(Lxf1.data(), Lxf1.data() + Lxf1.size());
			double MAD = Lxf1[Lxf1.size() / 2];
			Lx_scale[t] = 1.4826 * MAD;
		}

		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > Lyf = Ly * f;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > Lyf1(size_time);
		Lyf.array() = Lyf.array().abs();
		for (size_t f = 0; f < size_f; ++f){
			for (size_t t = 0; t < size_time; ++t)	Lyf1[t] = Lyf[t * size_f + f];
			std::sort(Lyf1.data(), Lyf1.data() + Lyf1.size());
			double MAD = Lyf1[Lyf1.size() / 2];
			Ly_scale[f] = 1.4826 * MAD;
		}

		//	evaluating MAD
		return update_L_weight_from_scale(Lxf, Lyf);
	}

	template<typename _Scalar >
	inline bool SparseLinearFit_L1Regularized_2D< _Scalar >::
		update_L_weight_from_scale
		(const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& Lxf, 
		 const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& Lyf)
	{	//	evaluate the weight matrix
		for (size_t t = 0; t < size_time; ++t){
			for (size_t f = 0; f < size_f; ++f){
				Lx_weight[t * size_f + f] = Lx_scale[t] / sqrt(Lxf[t * size_f + f] * Lxf[t * size_f + f] + Lx_scale[t] * Lx_scale[t]);
				Ly_weight[t * size_f + f] = Ly_scale[f] / sqrt(Lyf[t * size_f + f] * Lyf[t * size_f + f] + Ly_scale[f] * Ly_scale[f]);
			}
		}
		return true;
	}




	//--------------------------------------------------//
	//													//
	//			robust fit with Regularization			//
	//													//
	//--------------------------------------------------//
	//	the regularization matrix L should be specified before "run"
	//	using turkey's biweight.
	template<typename _Scalar >
	inline SparseRobustFit_Regularized_2D< _Scalar >::SparseRobustFit_Regularized_2D
		(const size_t size_time_, const size_t size_d_, const size_t size_f_,
		const size_t size_lx_, const size_t size_ly_)
		:
		SparseLinearFit_Regularized_2D(size_time_, size_d_, size_f_, size_lx_, size_ly_)
	{
		weight.resize(size_time * size_d);
		weight.setOnes();
		a_robust = 2.560;
	}


	template<typename _Scalar >
	inline bool SparseRobustFit_Regularized_2D< _Scalar >::solve
		(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel)
	{
		lambda_x = lambda_x_;
		lambda_y = lambda_y_;
		setL();
		set_Jac_and_Hes();
		update_f();
		for (size_t t = 0; t < iterMax; ++t){
			update_weight();
			set_Jac_and_Hes();
			double norm = update_f();
			if (norm < epsrel){
				calculate_cov();
				return true;
			}
		}
		calculate_cov();
		return false;
	}

	template<typename _Scalar >
	inline bool SparseRobustFit_Regularized_2D< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * A 
			+ lambda_x * LLx + lambda_y * LLy;
		return true;
	}

	template<typename _Scalar >
	inline bool SparseRobustFit_Regularized_2D< _Scalar >::update_weight(){
		res = d - A * f;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> weighted_res = res;
		//	evaluating MAD
		weighted_res.array() = weighted_res.array().abs();
		for (size_t t = 0; t < size_time; ++t){
			Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> wres1 = weighted_res.block(t * size_lx, 0, size_lx, 1);
			std::sort(wres1.data(), wres1.data() + wres1.size());
			double MAD = wres1[wres1.size() / 2];
			double scale_inv = 1.0/(1.4826 * MAD);

			for (int i = 0; i < size_d; ++i)
				weight[t*size_d + i] = get_weight_from_u_over_a(res[t*size_d + i] * scale_inv / a_robust);
		}
		return true;
	}

	template<typename _Scalar >
	inline double SparseRobustFit_Regularized_2D< _Scalar >::get_weight_from_u_over_a(const double u_over_a)const{
		double one_minus_u_over_a2 = 1.0 - u_over_a*u_over_a;
		if (one_minus_u_over_a2 > 0.0)
			return one_minus_u_over_a2 * one_minus_u_over_a2;
		else return 0.0;
	}

	//--------------------------------------------------//
	//													//
	//			robust fit with L1-Regularization		//
	//													//
	//--------------------------------------------------//
	//	the regularization matrix L should be specified before "run"
	//	using turkey's biweight.
	//	a coefficient is 2.560
	template<typename _Scalar >
	inline SparseRobustFit_L1Regularized_2D< _Scalar >::SparseRobustFit_L1Regularized_2D
		(const size_t size_time_, const size_t size_d_, const size_t size_f_,
		const size_t size_lx_, const size_t size_ly_)
		:
		SparseLinearFit_L1Regularized_2D(size_time_, size_d_, size_f_, size_lx_, size_ly_)
	{
		weight.resize(size_time * size_d);
		weight.setOnes();
		a_robust = 2.560;
	}


	template<typename _Scalar >
	inline bool SparseRobustFit_L1Regularized_2D< _Scalar >::solve
		(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel)
	{
		lambda_x = lambda_x_;
		lambda_y = lambda_y_;
		setL();		//	calculate L * weight * L
		set_Jac_and_Hes();
		update_f();
		for (size_t t0 = 0; t0 < iterMax; ++t0){
			for (size_t t = 0; t < iterMax; ++t){
				update_weight();
				set_Jac_and_Hes();
				double norm = update_f();
				if (norm < epsrel)
					break;
			}
			update_weight();
			update_L_weight();
			setL();
			set_Jac_and_Hes();
			double norm = update_f();
			if (norm < epsrel){
				calculate_cov();
				return true;
			}
		}
		calculate_cov();
		return false;
	}

	template<typename _Scalar >
	inline bool SparseRobustFit_L1Regularized_2D< _Scalar >::solve_updateLfirst
		(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel)
	{
		lambda_x = lambda_x_;
		lambda_y = lambda_y_;
		setL();		//	calculate L * weight * L
		set_Jac_and_Hes();
		update_f();
		for (size_t t0 = 0; t0 < iterMax; ++t0){
			for (size_t t = 0; t < iterMax; ++t){
				update_L_weight();
				set_Jac_and_Hes();
				double norm = update_f();
				if (norm < epsrel)
					break;
			}
			update_L_weight();
			update_weight();
			setL();
			set_Jac_and_Hes();
			double norm = update_f();
			if (norm < epsrel){
				calculate_cov();
				return true;
			}
		}
		calculate_cov();
		return false;
	}


	template<typename _Scalar >
	inline bool SparseRobustFit_L1Regularized_2D< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * A
			+ lambda_x * LLx + lambda_y * LLy;
		return true;
	}

	template<typename _Scalar >
	inline bool SparseRobustFit_L1Regularized_2D< _Scalar >::update_weight(){
		res = d - A * f;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> weighted_res = res;
		//	evaluating MAD
		weighted_res.array() = weighted_res.array().abs();
		for (size_t t = 0; t < size_time; ++t){
			Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> wres1 = weighted_res.block(t * size_lx, 0, size_lx, 1);
			std::sort(wres1.data(), wres1.data() + wres1.size());
			double MAD = wres1[wres1.size() / 2];
			double scale_inv = 1.0 / (1.4826 * MAD);

			for (int i = 0; i < size_d; ++i)
				weight[t*size_d + i] = get_weight_from_u_over_a(res[t*size_d + i] * scale_inv / a_robust);
		}
		return true;
	}

	template<typename _Scalar >
	inline double SparseRobustFit_L1Regularized_2D< _Scalar >::get_weight_from_u_over_a(const double u_over_a)const{
		double one_minus_u_over_a2 = 1.0 - u_over_a*u_over_a;
		if (one_minus_u_over_a2 > 0.0)
			return one_minus_u_over_a2 * one_minus_u_over_a2;
		else return 0.0;
	}

	/*
	---------------------	under construction  -------------------------
	*/
};