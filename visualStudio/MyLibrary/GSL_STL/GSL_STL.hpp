#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_permutation.h>
#include <vector>
#include <complex>

namespace GSL_STL{
//---	get stl contener from gsl_vector	--
	static std::vector<double> vector(gsl_vector *p){
		std::vector<double> v(p->size);
		for(size_t i=0; i<v.size(); i++) v[i] = gsl_vector_get(p, i);
		return v;}

	static std::vector<std::complex<double>> vector(gsl_vector_complex *p){
		std::vector<std::complex<double>> v(p->size);
		for(size_t i=0; i<v.size(); i++) v[i] = std::complex<double>(GSL_VECTOR_REAL(p, i), GSL_VECTOR_IMAG(p, i));
		return v;}

	static std::vector<size_t> vector(gsl_permutation *p){
		std::vector<size_t> v(p->size);
		for(size_t i=0; i<v.size(); i++) v[i] = gsl_permutation_get(p, i);
		return v;}
		
//---	get stl contener matrix from gsl_matrix	---
	static std::vector<std::vector<double>> (gsl_matrix *p){
		std::vector<std::vector<double>> v(p.size2);
		for(size_t i=0; i<v.size(); i++) v[i] = gsl_permutation_get(p, i);
		return v;}
		

};