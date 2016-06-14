#ifndef __MultiscaleMULTIMIN_HPP__
#define __MultiscaleMULTIMIN_HPP__
#include "Multimin.hpp"
namespace mygsl{

	class MultiscaleMultimin{
	public:
	//---	Constructor	---
		template<class TyF_ >
		MultiscaleMultimin
			( TyF_& function, std::vector<double>& params, std::vector<size_t>& params_order, double step_ratio = 1.0e-4, double relative_err= 1.0e-4, bool isprint_=true);

/*		template<class TyF_ >
		MultiscaleMultimin
			( TyF_& function, std::vector<double>& params, std::vector<size_t>& params_order, std::vector<double>& params_step, double relative_err, IsPrint isprint_);

		template<class TyF_ >
		MultiscaleMultimin
			( TyF_& function, std::vector<double>& params, std::vector<size_t>& params_order, std::vector<double>& params_step, double relative_err, std::vector<std::string> params_name_, IsPrint isprint_);

		template<class TyF_ >
		MultiscaleMultimin
			( TyF_& function, std::vector<double>& params, std::vector<size_t>& params_order, double step_ratio, double relative_err, std::vector<std::string> params_name_, IsPrint isprint_);
*/


		~MultiscaleMultimin(void){};

	};

	template < class TyF_ >
	class OnescaleMultimin{
	public:
		OnescaleMultimin(void){};
		~OnescaleMultimin(void){};

		void setOriginal( TyF_ * OriginalObject_, std::vector<double>* initial_all_params, const std::vector<size_t>& free_params_index_){
			OriginalObject = OriginalObject_;
			LowerscaleObject = NULL;
			all_params = initial_all_params;
			free_params_index = free_params_index_;
			IsOriginal = true;
		}
		void setScale( OnescaleMultimin< TyF_ > * LowerscaleObject_, std::vector<double>* initial_all_params, const std::vector<size_t>& free_params_index_){
			OriginalObject = NULL;
			LowerscaleObject = LowerscaleObject_;
			all_params = initial_all_params;
			free_params_index = free_params_index_;
			IsOriginal = false;
		}

		double operator() (std::vector<double>& free_params){
			//	first, the free_params is copied to the all_params
			for(size_t i=0; i<free_params_index.size(); ++i)	
				(*all_params)[free_params_index[i]] = free_params[i];

			//	if ObjectToBeMinimized is the functors to be minimized (if it is the most bottom-scale class)
			if(IsOriginal){
				return ( *OriginalObject)( * all_params);
			}
			else{
				//	prepare the free_parameters for minimizing the "ObjectToBeMinimized "
				std::vector<double> free_params_previous(LowerscaleObject -> getFreeParamsIndex().size());
				for(size_t i=0; i<free_params_previous.size(); ++i)
					free_params_previous[i] = ( *all_params)[LowerscaleObject -> getFreeParamsIndex()[i]];

				//---	derived the parameters	----
				Multimin multimin( * LowerscaleObject, free_params_previous, 1.0e-4, 1.04-4, Multimin::NoPrint);
			
				//---	copy to the *all_params from the derived parameters	---
				for(size_t i=0; i<free_params_previous.size(); ++i)
					(*all_params)[LowerscaleObject-> getFreeParamsIndex()[i]] = free_params_previous[i];

				return (*LowerscaleObject)(free_params_previous);
			}
		}
		const std::vector<size_t>& getFreeParamsIndex()const{return free_params_index;}
	private:
		bool IsOriginal;		//	if the TyF_ is OnescaleMultimin, IsOriginal is false
		TyF_ * OriginalObject;
		OnescaleMultimin< TyF_ > * LowerscaleObject;
		std::vector<double>* all_params;			//	pointing the all parameters which are always common among every scale.
		std::vector<size_t> free_params_index;		//	stores the indices of the free_parameters in the all_parameter:
	};


};
#include "MultiscaleMultimin.inl"
#endif