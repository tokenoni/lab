#ifndef __MYGSL_INTERP_VECTOR_HPP__
#define __MYGSL_INTERP_VECTOR_HPP__

#include "interp.hpp"

namespace mygsl{

	class interp_vector{
	public:
	//--	constructors	---
		interp_vector(void){};
		interp_vector(size_t size){data.resize(size);}
		~interp_vector(void){};
		interp_vector(const interp_vector& obj);
		interp_vector& operator = (const interp_vector& obj);

	//---	functions	---
		void resize(size_t size){data.resize(size);}
		size_t size()const{return data.size();}
		mygsl::interp& operator[] (size_t i) {return data[i];}

	private:
		std::vector<mygsl::interp> data;

	//---	function for gsl::vector	---
#ifdef CCGSL_VECTOR_HPP
	public:	
		const gsl::vector& operator () (double x);
	private: 
		gsl::vector vdata;
#endif
	};
};

#include "interp_vector.inl"
#endif
