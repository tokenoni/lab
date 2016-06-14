#ifndef __LHD_SVD_HPP__
#define __LHD_SVD_HPP__


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include < vector >
#include <iostream>
#include "LHDretrieve_dataset.hpp"

// class for calculating the Singular Value Decomposition Methods;
class LHD_SVD:public LHDretrieve_dataset
{
public:
	gsl_matrix* A;
	gsl_vector* S;
	gsl_matrix* V;
	gsl_matrix* X;
	gsl_vector* work;

//	void solve(const LHDretrieve_dataset& dataset);
	bool solve(const LHDretrieve_dataset& dataset,const size_t start_p=0, const size_t end_p=0, const size_t increment=1);
	void setUncommonData();	//	all data is de

	void solve(std::vector<std::vector<double>>& src);
		//	src[i][j] = Aji		<--		order will be changed! caution!
	void getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv);
	//	application functions
	std::vector<double> getSingularValues();
	std::vector<std::vector<double>> getEigenFunctions();
	std::vector<std::vector<double>> getChannelDistribution();

	LHD_SVD(void){Allocated=false;m=0; n=0;}
	~LHD_SVD(void);
private:
	bool Allocated;
	size_t m, n;
	void Allocation(size_t m_, size_t n_);
	void gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src);
	void gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src);
};

inline bool LHD_SVD::solve(const LHDretrieve_dataset& dataset,const size_t start_p, const size_t end_p_, const size_t increment){
	size_t size = 0;
	for(size=0; size<dataset.size() && dataset.IsAllocated(size); ++size){}
	if(size==0) return false;
	
	size_t end_p = end_p_;
	if(end_p==0) end_p = dataset.t_size();
	if(end_p<=start_p) return false;
	size_t num = (size_t)floor(1.0*(end_p-start_p)/increment);

	data_v.resize(size, LHDretrieve_data(num));
	setRate(   dataset.getRate()/increment);
	setStart_t(dataset.getStart_t()+start_p/dataset.getRate());
	setEnd_t(  getStart_t() + t_size()/getRate());

	//---	matrix allocation	---
	Allocation(t_size(), size);
	if(m<n) return false;
	//	copy into gsl_matrix* A
	for(size_t i=0; i<m; ++i){
		for(size_t j=0; j<n; ++j) gsl_matrix_set(A, i, j, dataset[j][start_p+i*increment]);
	};

	//---	the condition of m >> n is assumed as m > 100*n
	if(m < 100*n)	//	Golub Reinsch method;
		gsl_linalg_SV_decomp(A, V, S, work);
	else	//	modified Golub Reinsch method; 
		gsl_linalg_SV_decomp_mod(A, X, V, S, work);
	return true;
}

inline void LHD_SVD::setUncommonData(){
	if(!Allocated) return;
	for(size_t i=0; i<m; ++i){
		for(size_t k=0; k<n; ++k){
			data_v[k].get(i) = 0.0;
			for(size_t j=0; j<n; ++j){
				if(j!=0)
					data_v[j].get(i) += gsl_matrix_get(A, i, j)*gsl_vector_get(S, j)*gsl_matrix_get(V, k,j);
			}
		}
	}
}



inline void LHD_SVD::getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv){
	gslv_to_stdv(s_stdv, S);
	gslm_to_stdv(u_stdv, A);
	gslm_to_stdv(v_stdv, V);
}

inline std::vector<double> LHD_SVD::getSingularValues(){
	std::vector<double> tmp;
	gslv_to_stdv(tmp, S);
	return tmp;
}
inline std::vector<std::vector<double>> LHD_SVD::getEigenFunctions(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, A);
	return tmp;
}
inline std::vector<std::vector<double>> LHD_SVD::getChannelDistribution(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, V);
	return tmp;
}

inline void LHD_SVD::gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src){
	dest.resize(src->size);
	for(size_t i=0; i<dest.size(); i++)	dest[i] = gsl_vector_get(src,i);
}

inline void LHD_SVD::gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src){
	dest.resize(src->size1);
	for(size_t i=0; i<dest.size(); i++)	{
		dest[i].resize(src->size2);
		for(size_t j=0; j<src->size2; ++j)	dest[i][j] = gsl_matrix_get(src, i, j);
	}
}
inline void LHD_SVD::solve(std::vector<std::vector<double> > &src){
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
inline void LHD_SVD::Allocation(size_t m_, size_t n_){
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

inline LHD_SVD::~LHD_SVD(void){
	if(Allocated){
		gsl_matrix_free(A);	gsl_vector_free(S);	gsl_matrix_free(V);	
		gsl_matrix_free(X);	gsl_vector_free(work);
		Allocated = false;
	}
}

#endif

