#ifndef CCGSL_FUJII_MATRIX_NOTATION
#define CCGSL_FUJII_MATRIX_NOTATION

#include "matrix.hpp"
#include "vector.hpp"

namespace gsl{
	void matrix::mul(const matrix &a, const matrix& b);
	matrix matrix::operator * (const matrix &a);

	inline void matrix::mul(const matrix &a, const matrix& b){
		size_t size1 = a.size1();
		size_t size2 = a.size2();
		size_t size3 = b.size2();
		if(size2 != b.size1())
		gsl_error( "size doesn't match between two matrices ", __FILE__, __LINE__, exception::GSL_EFAILED );
		else{
			setzero();
			for(size_t i=0; i<size1; i++){
				for(size_t k=0; k<size3; k++){
					for(size_t j=0; j<size2; j++){
						set(i,k,get(i,k) + a.get(i,j)*b.get(j,k));
					}}}
		}return;
	}

	inline matrix matrix::operator * (const matrix& a)const{
		matrix product(size1(), a.size2());
		product.mul(*this, a);
		return product;
	}
	

	/*class vector{
	public:
		vector operator* (vector& a){
			size_t size1 = get()->size1;
			size_t size2 = get()->size2;
			vector product(size1);
			if(size2 != a.get()->size)
				gsl_error( "size doesn't match between two matrices ", __FILE__, __LINE__, exception::GSL_EFAILED );
			else{
				for(size_t i=0; i<size1; i++){
					product[i] = 0.0;
					for(size_t j=0; j<size2; j++){
						product[i] += this[i][j]*a[j];
					}}
			}	return product;
		}
	};*/
};

#endif
