namespace myEigen{
	
	//---------------------------------------------------//
	//													 //
	//				basic class for linear fit			 //
	//													 //
	//---------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	template<typename _Scalar >
	inline LinearFit< _Scalar >::LinearFit(const size_t size_d_, const size_t size_f_)
		:size_f(size_f_), size_d(size_d_)
	{
		A.resize(size_d, size_f);	A.setZero();
		f.resize(size_f);			f.setZero();
		d.resize(size_d);			d.setZero();
		res.resize(size_d);			res.setZero();
		Jac.resize(size_f);			Jac.setZero();
		Hes.resize(size_f, size_f);	Hes.setZero();
	}

	template<typename _Scalar >
	inline bool LinearFit< _Scalar >::solve(){
		set_Jac_and_Hes();
		update_f();
		calculate_cov();
		return true;
	}

	template<typename _Scalar >
	inline double LinearFit< _Scalar >::update_f(){
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> fprev = f;
		f = Hes.lu().solve(Jac);
		return (f - fprev).norm() / f.norm();
	}

	template<typename _Scalar >
	inline bool LinearFit< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * d;
		Hes = A.transpose() * A;
		return true;
	}

	template<typename _Scalar >
	inline bool LinearFit< _Scalar >::calculate_cov(){
		set_Jac_and_Hes();
		cov = Hes.lu().inverse();
		return true;
	}

	template<typename _Scalar >
	inline void LinearFit< _Scalar >::initialize_by_copy(const LinearFit< _Scalar >& obj){
		f = obj.f;
	}

	template<typename _Scalar >
	inline Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> LinearFit< _Scalar >::getError()const
	{
		return cov.diagonal().array().sqrt().matrix();
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
	inline LinearFit_Weighted< _Scalar >::LinearFit_Weighted(const size_t size_d_, const size_t size_f_)
		: LinearFit< _Scalar >(size_d_, size_f_)
	{
		alpha.resize(size_d);
		alpha.setOnes();
	}

	template<typename _Scalar >
	inline bool LinearFit_Weighted< _Scalar >::set_Jac_and_Hes(){
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
	inline LinearFit_Regularized< _Scalar >::LinearFit_Regularized(const size_t size_d_, const size_t size_f_, const size_t size_l_)
		: LinearFit_Weighted< _Scalar >(size_d_, size_f_),
		size_l(size_l_)
	{
		L.resize(size_l, size_f);	L.setZero();
		LL.resize(size_f, size_f);	LL.setZero();
	}

	template<typename _Scalar >
	inline bool LinearFit_Regularized< _Scalar >::set_Jac_and_Hes(){
		setL();
		Jac = A.transpose() * alpha.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * A + lambda * LL;
		return true;
	}
	template<typename _Scalar >
	inline bool LinearFit_Regularized< _Scalar >::solve(const double lambda_){
		lambda = lambda_;
		set_Jac_and_Hes();
		update_f();
		calculate_cov();
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
	inline RobustFit_Regularized< _Scalar >::RobustFit_Regularized
		(const size_t size_d_, const size_t size_f_, const size_t size_l_)
		: LinearFit_Regularized< _Scalar >(size_d_, size_f_, size_l_)
	{
		weight.resize(size_d);
		weight.setOnes();
		a_robust = 2.560;
	};

	template<typename _Scalar >
	inline bool RobustFit_Regularized< _Scalar >::solve
		(const double lambda_, const size_t iterMax, const double epsrel)
	{
		lambda = lambda_;
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
	inline bool RobustFit_Regularized< _Scalar >::set_Jac_and_Hes(){
		Jac = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * d;
		Hes = A.transpose() * alpha.asDiagonal() * weight.asDiagonal() * A + lambda * LL;
		return true;
	}

	template<typename _Scalar >
	inline bool RobustFit_Regularized< _Scalar >::update_weight(){
		res = d - A * f;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > weighted_res = res;
		//	evaluating MAD
		weighted_res.array() = weighted_res.array().abs();
		std::sort(weighted_res.data(), weighted_res.data() + weighted_res.size());
		double MAD = weighted_res[weighted_res.size() / 2];
		double scale_inv = 1.0 / (1.4826 * MAD);
		//	evaluate the weight matrix
		weight = alpha;
		for (int i = 0; i < d.size(); ++i){
			weight[i] = get_weight_from_u_over_a(res[i] * scale_inv / a_robust);
		}
		return true;
	}

	template<typename _Scalar >
	inline double RobustFit_Regularized< _Scalar >::get_weight_from_u_over_a(const double u_over_a)const{
		double one_minus_u_over_a2 = 1.0 - u_over_a*u_over_a;
		if (one_minus_u_over_a2 > 0.0)
			return one_minus_u_over_a2 * one_minus_u_over_a2;
		else return 0.0;
	}
	//	initializing
	template<typename _Scalar >
	inline void RobustFit_Regularized< _Scalar >::initialize_by_copy(const RobustFit_Regularized< _Scalar >& obj){
		LinearFit_Weighted< _Scalar >::initialize_by_copy(obj);
		weight = obj.weight;
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
	inline RobustFit_L1Regularized< _Scalar >::RobustFit_L1Regularized
		(const size_t size_d_, const size_t size_f_, const size_t size_l_)
		: RobustFit_Regularized< _Scalar >(size_d_, size_f_, size_l_)
	{
		L_weight.resize(size_l);
		L_weight.setOnes();
	};

	template<typename _Scalar >
	inline bool RobustFit_L1Regularized< _Scalar >::solve
		(const double lambda_, const size_t iterMax, const double epsrel)
	{
		lambda = lambda_;
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
	inline bool RobustFit_L1Regularized< _Scalar >::update_L_weight()
	{
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1 > Lf = L * f;
		L_weight = Lf;

		//	evaluating MAD
		Lf.array() = Lf.array().abs();
		std::sort(Lf.data(), Lf.data() + Lf.size());
		double MAD = Lf[Lf.size() / 2];
		double scale = 1.4826 * MAD;
		return update_L_weight_from_scale(scale);
	}

	template<typename _Scalar >
	inline bool RobustFit_L1Regularized< _Scalar >::update_L_weight_from_scale(const double scale)
	{	//	evaluate the weight matrix
		for (int i = 0; i < L_weight.size(); ++i)
			L_weight[i] = scale / sqrt(L_weight[i] * L_weight[i] + scale*scale);
		return true;
	}

	template<typename _Scalar >
	inline void RobustFit_L1Regularized< _Scalar >::setL(){
		if (L_weight.size() != L.cols()) LL = L.transpose() * L;
		else LL = L.transpose() * L_weight.asDiagonal() * L;
	}

	//	initializing
	template<typename _Scalar >
	inline void RobustFit_L1Regularized< _Scalar >::initialize_by_copy(const RobustFit_L1Regularized< _Scalar >& obj)
	{
		LinearFit_Regularized< _Scalar >::initialize_by_copy(obj);
		L_weight = obj.L_weight;
	}
};