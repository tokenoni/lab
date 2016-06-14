#ifndef __MYGSL_VECTOR2INDEXES_HPP__
#define __MYGSL_VECTOR2INDEXES_HPP__
#include <ccgsl/vector.hpp>

namespace mygsl{
	class vector2indexes {
	public:
		vector2indexes(void){i0num = 1;}
		vector2indexes(const size_t i0num_, const size_t i1num_){
			i0num = i0num_;
			data = gsl::vector(i0num_* i1num_);
		}

		void setSize(const size_t i0num_){
			i0num = i0num_;
		}

		void resize(size_t const i0num_, size_t const i1num_){
			i0num = i0num_;
			data = gsl::vector(i0num*i1num_);
		}
	
		void set( size_t const i0, size_t const i1, double x)	{
			data.set(i0+i0num*i1, x);
		}
	
		double get( size_t const i0, size_t const i1)const	{
			return data.get(i0+i0num*i1);
		}

		void set_zero(){ data.set_zero();}

		gsl_vector* get(){ return data.get();}
		gsl::vector& operator()(){return data;}

/*		gsl::matrix getMatrix()const{
			return gsl::matrix::const_view_array_with_tda(data.get()->data, i0num, data.size()/i0num, i0num*data.get()->stride);
		}
		*/
	private:
		size_t i0num;
		gsl::vector data;
	};
};
#endif
