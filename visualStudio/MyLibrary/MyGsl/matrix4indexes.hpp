#ifndef __MYGSL_MATRIX4INDEXES_HPP__
#define __MYGSL_MATRIX4INDEXES_HPP__
#include <ccgsl\matrix.hpp>

namespace mygsl{
	class matrix4indexes{
	public:
		void resize(const size_t i0num_, const size_t i1num_, const size_t j0num_, const size_t j1num_){
			i0num = i0num_;	j0num = j0num_;
			data = gsl::matrix(i0num_*i1num_, j0num_*j1num_);
			data.set_zero();
		}

		void set_zero(){	data.set_zero();}

		gsl::matrix& operator ()() {return data;}

		void set( size_t const i0, size_t const i1, size_t const j0, size_t const j1, double x){
			data.set(i0+i0num*i1, j0+j0num*j1, x);
		}
		
		double get( size_t const i0, size_t const i1, size_t const j0, size_t const j1)const	{
			return data.get(i0+i0num*i1, j0+j0num*j1);
		}

		gsl_matrix* get(){
			return data.get();
		}

	private:
		size_t i0num, j0num;
		gsl::matrix data;
	};

};
#endif
