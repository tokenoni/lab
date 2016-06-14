#pragma once
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include < vector >
#include <iostream>

// class for calculating the Singular Value Decomposition Methods;
class SVD
{
public:
	gsl_matrix* A;
	gsl_vector* S;
	gsl_matrix* V;
	gsl_matrix* X;
	gsl_vector* work;
	void solve(std::vector<std::vector<double>>& src);
		//	src[i][j] = Aji		<--		order will be changed! caution!
	void getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv);
	//	application functions
	std::vector<double> getSingularValues();
	std::vector<std::vector<double>> getEigenFunctions();
	std::vector<std::vector<double>> getChannelDistribution();

	SVD(void){Allocated=false;m=0; n=0;}
	~SVD(void);
private:
	bool Allocated;
	size_t m, n;
	void Allocation(size_t m_, size_t n_);
	void gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src);
	void gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src);
};

inline void SVD::getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv){
	gslv_to_stdv(s_stdv, S);
	gslm_to_stdv(u_stdv, A);
	gslm_to_stdv(v_stdv, V);
}

inline std::vector<double> SVD::getSingularValues(){
	std::vector<double> tmp;
	gslv_to_stdv(tmp, S);
	return tmp;
}
inline std::vector<std::vector<double>> SVD::getEigenFunctions(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, A);
	return tmp;
}
inline std::vector<std::vector<double>> SVD::getChannelDistribution(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, V);
	return tmp;
}

inline void SVD::gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src){
	dest.resize(src->size);
	for(size_t i=0; i<dest.size(); i++)	dest[i] = gsl_vector_get(src,i);
}

inline void SVD::gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src){
	dest.resize(src->size1);
	for(size_t i=0; i<dest.size(); i++)	{
		dest[i].resize(src->size2);
		for(size_t j=0; j<src->size2; ++j)	dest[i][j] = gsl_matrix_get(src, i, j);
	}
}
inline void SVD::solve(std::vector<std::vector<double> > &src){
	if(src.size()<=0)
	{ std::cerr << "invalid matrix" <<std::endl; return;}
	//---	matrix allocation	---
	Allocation(src[0].size(), src.size());
	//	copy into gsl_matrix* A
	for(size_t i=0; i<m; ++i) for(size_t j=0; j<n; ++j) gsl_matrix_set(A, i, j, src[j][i]);

	if(m<n)	{ std::cerr << "invalid matrix" <<std::endl; return;}

	//---	the condition of m >> n is assumed as m > 100*n
	if(m < 100*n)	//	Golub Reinsch method;
		gsl_linalg_SV_decomp(A, V, S, work);
	else	//	modified Golub Reinsch method; 
		gsl_linalg_SV_decomp_mod(A, X, V, S, work);
}

//---	function for memory allocation	---
inline void SVD::Allocation(size_t m_, size_t n_){
	//	if memory is already allocated, it should be free
	if(Allocated && (m != m_ || n != n_))	{
		gsl_matrix_free(A);	gsl_vector_free(S);	gsl_matrix_free(V);	
		gsl_matrix_free(X);	gsl_vector_free(work);
		Allocated = false;
	}
	m=m_;	n=n_;
	//	if memory is not yet allocated
	if(!Allocated){
		A = gsl_matrix_calloc(m, n);
		S = gsl_vector_calloc(n);
		V = gsl_matrix_calloc(n, n);
		X = gsl_matrix_calloc(n, n);
		work = gsl_vector_calloc(n);
		Allocated = true;
	}
}

inline SVD::~SVD(void){
	if(Allocated){
		gsl_matrix_free(A);	gsl_vector_free(S);	gsl_matrix_free(V);	
		gsl_matrix_free(X);	gsl_vector_free(work);
		Allocated = false;
	}
}