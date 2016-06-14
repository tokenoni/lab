#ifndef __MYGSL_INTERP_MATRIX_HPP__
#define __MYGSL_INTERP_MATRIX_HPP__

#include "interp.hpp"
#include "interp_vector.hpp"

namespace mygsl{

	class interp_matrix{
	public:
	//--	constructors	---
		interp_matrix(void){};
		interp_matrix(size_t size){ data.resize(size);}
		interp_matrix(size_t size1, size_t size2){data.resize(size1, mygsl::interp_vector(size2));}
		~interp_matrix(void){};
		interp_matrix(const interp_matrix& obj);
		interp_matrix& operator = (const interp_matrix& obj);

	//---	functions	---
		void resize(const size_t size1, const size_t size2);
		size_t size()const{return data.size();}
		mygsl::interp_vector& operator[] (size_t i) {return data[i];}

	private:
		std::vector<mygsl::interp_vector> data;
	//---	function for gsl::vector	---
#ifdef CCGSL_VECTOR_HPP
	public:	
		const gsl::matrix& operator () (double x);
		bool set(const std::vector<double>& xsrc, const std::vector<gsl::matrix>& msrc);
	private: 
		gsl::matrix mdata;
#endif
	};
};

#include "interp_matrix.inl"
#endif
