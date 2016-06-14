#pragma once
#include <eigen/Sparse>
#include <algorithm>

namespace myEigen{
	//---------------------------------------------------//
	//													 //
	//				basic class for linear fit			 //
	//													 //
	//---------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	template<typename _Scalar >
	class SparseLinearFit{
	public:
		const size_t max_size;
		SparseLinearFit(const size_t size_time_, const size_t size_d_, const size_t size_f_);
		virtual ~SparseLinearFit(void){};
		
		//	initializing
		virtual void initialize_by_copy(const SparseLinearFit< _Scalar >& obj);

		//	solving
		virtual bool solve();
		virtual double update_f();

		//	get functions
		virtual Eigen::VectorXd getResult()const{ return A * f; };

		//	const variables
		const size_t size_time;	//	number of frames
		const size_t size_d;	//	number of data point in one frame
		const size_t size_f;	//	number of parameters in one frame
	protected:
		//	member variables
		//	d = A * f + e  (e: noise)
		Eigen::SparseMatrix< _Scalar > A;		//	evaluation matrix
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> f;		//	parameters to be derived
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> d;		//	experimental data

		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> res;		//	residual
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> Jac;		//	Jacobian	
		Eigen::SparseMatrix< _Scalar > Hes;	//	minus of Hessian
		Eigen::SparseMatrix< _Scalar > cov;	//	minus of Hessian

		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool calculate_cov();

	public:
		virtual Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> getError()const;
		virtual std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > get_f()const;
		virtual std::vector< Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> > get_f_sqerr()const;
	};

	//--------------------------------------------------//
	//													//
	//			linear fit with specified sigma			//
	//													//
	//--------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	//	the width of the Gauss distribution should be set as alpha, 
	//	where alpha = 1.0 / (sigma * sigma)
	template<typename _Scalar >
	class SparseLinearFit_Weighted : public SparseLinearFit< _Scalar >{
	public:
		SparseLinearFit_Weighted(const size_t size_time_, const size_t size_d_, const size_t size_f_);
	protected:
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> alpha;
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
	};

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
	class SparseLinearFit_Regularized_2D : public SparseLinearFit_Weighted< _Scalar >{
	public:
		SparseLinearFit_Regularized_2D
			(const size_t size_time_, const size_t size_d_, const size_t size_f_, 
			 const size_t size_lx_, const size_t size_ly_);
		//	solving
		virtual bool solve(const double lambda_x, const double lambda_y);

		const size_t size_lx, size_ly;
	protected:
		double lambda_x, lambda_y;
		Eigen::SparseMatrix< _Scalar > Lx, Ly;
		Eigen::SparseMatrix< _Scalar > LLx, LLy;
		virtual void setL(){ LLx = Lx.transpose() * Lx; LLy = Ly.transpose() * Ly; }
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();

	private:
		//	solving
		virtual bool solve(){return true;};
	};

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
	class SparseLinearFit_L1Regularized_2D : public SparseLinearFit_Regularized_2D< _Scalar >{
	public:
		SparseLinearFit_L1Regularized_2D
			(const size_t size_time_, const size_t size_d_, const size_t size_f_,
			 const size_t size_lx_, const size_t size_ly_);

		//	solving
		virtual bool solve
			(const double lambda_x, const double lambda_y, 
			 const size_t iterMax, const double epsrel);

	protected:
		virtual void setL();
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool update_L_weight();
		virtual bool update_L_weight_from_scale
			(const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& Lxf, 
			 const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& Lyf);
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> Lx_weight, Ly_weight;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> Lx_scale, Ly_scale;
	private:
		//	solving
		virtual bool solve(){ return true; };
	};

	//--------------------------------------------------//
	//													//
	//			robust fit with Regularization			//
	//													//
	//--------------------------------------------------//
	//	the regularization matrix L should be specified before "run"
	//	using turkey's biweight.
	//	a coefficient is 2.560
	template<typename _Scalar >
	class SparseRobustFit_Regularized_2D : public SparseLinearFit_Regularized_2D< _Scalar >{
	public:
		SparseRobustFit_Regularized_2D
			(const size_t size_time_, const size_t size_d_, const size_t size_f_,
			 const size_t size_lx_, const size_t size_ly_);

		//	solving
		virtual bool solve(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel);

	protected:
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> weight;	//	weight := alpha * weight. nxn matrix
		double a_robust;	//	a coefficient for the Turkey's biweight
		void set_a(const double a){ a_robust = a; }
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool update_weight();
		virtual double get_weight_from_u_over_a(const double u_over_a)const;
	private:
		//	prohobited
		virtual bool solve(){ return true; };

	public:
	};

	//--------------------------------------------------//
	//													//
	//			robust fit with L1-Regularization		//
	//													//
	//--------------------------------------------------//
	//	the regularization matrix L should be specified before "run"
	//	using turkey's biweight.
	//	a coefficient is 2.560
	template<typename _Scalar >
	class SparseRobustFit_L1Regularized_2D : public SparseLinearFit_L1Regularized_2D< _Scalar >{
	public:
		SparseRobustFit_L1Regularized_2D
			(const size_t size_time_, const size_t size_d_, const size_t size_f_,
			 const size_t size_lx_, const size_t size_ly_);
		
		//	solving
		virtual bool solve(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel);
		virtual bool solve_updateLfirst(const double lambda_x_, const double lambda_y_, const size_t iterMax, const double epsrel);

	protected:
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> weight;	//	weight := alpha * weight. nxn matrix
		double a_robust;	//	a coefficient for the Turkey's biweight
		void set_a(const double a){ a_robust = a; }
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool update_weight();
		virtual double get_weight_from_u_over_a(const double u_over_a)const;
	private:
		//	prohobited
		virtual bool solve(){ return true; };
	};

	/*
	---------------------	under construction  -------------------------

	*/
};

#include "SparseLinearFit.inl"