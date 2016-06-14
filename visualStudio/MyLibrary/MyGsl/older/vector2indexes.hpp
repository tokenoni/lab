#ifndef __MYGSL_VECTOR2INDEXES__
#define __MYGSL_VECTOR2INDEXES__

#include <gsl\gsl_vector.h>

//---	vector having two index	---
//	the element can be accessed as m[i*jmax + j]
//	where index i (spatial), j (velocity) 
namespace mygsl{
	class vector2indexes{
	public:
		vector2indexes(void);
		vector2indexes(size_t size_i, size_t size_j);
		void resize(size_t size_i, size_t size_j);
		
		double& operator () (size_t i, size_t j);
		size_t size()const{return data->size;}
		
		void set_zero(){	if(allocated) gsl_vector_set_zero(data);};
		gsl_vector* get(){ return data;}

	private:
		size_t jmax;
		gsl_vector* data;
		bool allocated;
		void allocate(size_t size_);
		void free();
	};
	
	//---	definition	---
	inline vector2indexes::vector2indexes(void){	allocated = false;}
	inline vector2indexes::vector2indexes(size_t size_i, size_t size_j){
		allocated = false; allocate(size_i*size_j);
		jmax = size_j;
	}
	inline void vector2indexes::resize(size_t size_i, size_t size_j){
		//	if the size_i and size_j is same 
		if(allocated && size_j == jmax && size_i * size_j == data->size) return;
		free();
		jmax = size_j;
		allocate(size_i * size_j);
	}
	
	inline double& vector2indexes::operator() (size_t i, size_t j){
		return data->data[(i*jmax+j)*data->stride];
	}
	inline void vector2indexes::allocate(size_t size_){
		if(!allocated)
			data = gsl_vector_calloc(size_);
		allocated = true;
	}
	inline void vector2indexes::free(){
		if(allocated)
			gsl_vector_free(data);
		allocated = false;
	}
};

#endif