#ifndef __MYGSL_matrix4indexes__
#define __MYGSL_matrix4indexes__

#include <gsl\gsl_vector.h>
#include <gsl\gsl_matrix.h>

//---	vector having two index	---
//	the element can be accessed as m[i*jmax + j]
//	where index i (spatial), j (velocity) 
namespace mygsl{
	class matrix4indexes{
	public:
		matrix4indexes(void);
		matrix4indexes(size_t size_i, size_t size_j, size_t size_k, size_t size_l);
		void resize(size_t size_i, size_t size_j, size_t size_k, size_t l);
		
		double& operator () (size_t i, size_t j, size_t k, size_t l);
		size_t size1()const{return data->size1;}
		size_t size2()const{return data->size2;}
		
		void setzero(){if(allocated) gsl_matrix_set_zero(data);}
		gsl_matrix* get(){return data;}

	private:
		size_t jmax, lmax;
		gsl_matrix* data;
		bool allocated;
		void allocate(size_t size1, size_t size2);
		void free();
	};
	
	//---	definition	---
	inline matrix4indexes::matrix4indexes(void){	allocated = false;}
	inline matrix4indexes::matrix4indexes(size_t size_i, size_t size_j, size_t size_k, size_t size_l){
		allocated = false; allocate(size_i*size_j, size_k*size_l);
		jmax = size_j;	lmax = size_l;
	}
	inline void matrix4indexes::resize(size_t size_i, size_t size_j, size_t size_k, size_t size_l){
		//	if the size_i and size_j is same 
		if(allocated && 
			(size_i * size_j == data->size1 && size_k*size_l == data->size2)) {
				jmax = size_j, lmax = size_l;	return;}
		free();
		jmax = size_j;	lmax = size_l;
		allocate(size_i * size_j, size_k * size_l);
	}
	
	inline double& matrix4indexes::operator() (size_t i, size_t j, size_t k, size_t l){
		return data->data[(i*jmax+j)*data->tda + (k*lmax +l)];
	}
	inline void matrix4indexes::allocate(size_t size1, size_t size2){
		if(!allocated)
			data = gsl_matrix_calloc(size1, size2);
		allocated = true;
	}
	inline void matrix4indexes::free(){
		if(allocated)
			gsl_matrix_free(data);
		allocated = false;
	}
};

#endif