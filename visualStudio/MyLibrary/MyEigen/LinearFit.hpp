#pragma once
#include <eigen/eigen>
#include <algorithm>

namespace myEigen{
	//---------------------------------------------------//
	//													 //
	//				basic class for linear fit			 //
	//													 //
	//---------------------------------------------------//
	//	the noise distribution of the data is assumed as a Gauss distribution
	template<typename _Scalar >
	class LinearFit{
	public:
		LinearFit(const size_t size_f_, const size_t size_d_);
		virtual ~LinearFit(void){};
		
		//	initializing
		virtual void initialize_by_copy(const LinearFit< _Scalar >& obj);

		//	solving
		virtual bool solve();
		virtual double update_f();

		
		const size_t size_f, size_d;
	protected:
		//	member variables
		//	d = A * f + e  (e: noise)
		Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic> A;		//	evaluation matrix
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> f;		//	parameters to be derived
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> d;		//	experimental data
		
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> res;	//	residual
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> Jac;					//	Jacobian	
		Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic> Hes;	//	minus of Hessian
		Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic> cov;	//	minus of Hessian

		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool calculate_cov();
		
	public:	//	get functions
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& get_d()const{ return d; }
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& get_f()const{ return f; }
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic>& get_A()const{ return A; }

		virtual Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> getResult()const{ return A * f; };
		virtual Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> getResult(const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& fgiven)const{	return A * fgiven;	};

		virtual Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> getError()const;
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
	class LinearFit_Weighted : public LinearFit< _Scalar >{
	public:
		LinearFit_Weighted(const size_t size_f_, const size_t size_d_);

	protected:
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> alpha;
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
	
	public://	get functions
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& get_alpha()const{ return alpha; }
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
	class LinearFit_Regularized : public LinearFit_Weighted< _Scalar >{
	public:
		LinearFit_Regularized(const size_t size_f_, const size_t size_d_, const size_t size_l_);
		//	solving
		virtual bool solve(const double lambda_);

		const size_t size_l;	//	size of regularization matrix
	protected:
		double lambda;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic> L;
		Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic> LL;
		virtual void setL(){ LL = L.transpose() * L; }
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();

	private:
		//	solving
		virtual bool solve(){return true;};

	public:
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic>& get_L()const{ return L; }
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, Eigen::Dynamic>& get_LL()const{ return LL; }
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
	class RobustFit_Regularized : public LinearFit_Regularized< _Scalar >{
	public:
		RobustFit_Regularized(const size_t size_f_, const size_t size_d_, const size_t size_l_);
		//	solving
		virtual bool solve(const double lambda_, const size_t iterMax, const double epsrel);
		//	initializing
		virtual void initialize_by_copy(const RobustFit_Regularized< _Scalar >& obj);

	protected:
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> weight;	//	weight := alpha * weight. nxn matrix
		double a_robust;	//	a coefficient for the Turkey's biweight
		void set_a(const double a){ a_robust = a; }
		//	calculate Jacobian and Hessian
		virtual bool set_Jac_and_Hes();
		virtual bool update_weight();
		virtual double get_weight_from_u_over_a(const double u_over_a)const;
	private:
		//	solving
		virtual bool solve(){ return true; };

	public://	get functions
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& get_weight()const{ return weight; }
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
	class RobustFit_L1Regularized : public RobustFit_Regularized< _Scalar >{
	public:
		RobustFit_L1Regularized(const size_t size_f_, const size_t size_d_, const size_t size_l_);
		//	initialize
		virtual void initialize_by_copy(const RobustFit_L1Regularized< _Scalar >& obj);
		//	solving
		virtual bool solve(const double lambda_, const size_t iterMax, const double epsrel);
	protected:
		//	update L weight
		virtual bool update_L_weight();
		virtual bool update_L_weight_from_scale(const double scale);
		//	set L : LL = Lt * reg_weight * L;
		virtual void setL();
		Eigen::Matrix< _Scalar, Eigen::Dynamic, 1> L_weight;	//	weight := alpha * weight. nxn matrix

	public://	get functions
		const Eigen::Matrix< _Scalar, Eigen::Dynamic, 1>& get_L_weight()const{ return L_weight; }
	};

};

#include "LinearFit.inl"