#ifndef __TIMEDOMAINDATASET_SVD_HPP__
#define __TIMEDOMAINDATASET_SVD_HPP__


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include < vector >
#include <iostream>
#include "TimeDomainDataSet.hpp"

// class for calculating the Singular Value Decomposition Methods;
class TimeDomainDataSet_SVD:public TimeDomainDataSet
{
public:
	gsl_matrix* A;
	gsl_vector* S;
	gsl_matrix* V;
	gsl_matrix* X;
	gsl_vector* work;

	bool solve(const TimeDomainDataSet& dataset);
	bool solve(const TimeDomainDataSet& dataset, const double start_t, const double delta_t, const double avg_t);
	void setUncommonData();	//	all data is de

	void solve(std::vector<std::vector<double>>& src);
		//	src[i][j] = Aji		<--		order will be changed! caution!
	void getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv);
	//	application functions
	std::vector<double> getSingularValues();
	std::vector<std::vector<double>> getEigenFunctions();
	std::vector<std::vector<double>> getChannelDistribution();

	TimeDomainDataSet_SVD(void){Allocated=false;m=0; n=0;}
	~TimeDomainDataSet_SVD(void);
	TimeDomainData common_component;
private:
	bool Allocated;
	size_t m, n;
	void Allocation(size_t m_, size_t n_);
	void gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src);
	void gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src);
};

//	dataset should be already excuted "eraseUnallocatedData"
inline bool TimeDomainDataSet_SVD::solve(const TimeDomainDataSet& dataset, const double start_t, const double delta_t, const double avg_t){
	size_t size = 0;
	for(size=0; size<dataset.size() && dataset.get(size).IsAllocated(); ++size){}
	if(size==0) return false;

	// prepare the timing 
	size_t avgnum = (size_t)floor(avg_t*dataset[0].getRate());
	if(avgnum<1) avgnum = 1;
	size_t start_p = dataset.getRate()*start_t;

	//---	matrix allocation	---
	size_t delta_p = (size_t)floor(delta_t*dataset.getRate()/avgnum);
	Allocation(delta_p , size);
	common_component.Allocate(delta_p);
	data_v.resize(size, TimeDomainData(delta_p));
	setTiming(start_t+0.5*avg_t, dataset.get_t(start_t, avgnum));
	if(m<n) return false;
	

	//	copy into gsl_matrix* A
	for(size_t i=0; i<m; ++i){
		for(size_t j=0; j<n; ++j){
			gsl_matrix_set(A, i, j, dataset[j].getAvg(start_p + i*avgnum, avgnum));
		}
	}

	//---	the condition of m >> n is assumed as m > 100*n
	if(m < 100*n)	//	Golub Reinsch method;
		gsl_linalg_SV_decomp(A, V, S, work);
	else	//	modified Golub Reinsch method; 
		gsl_linalg_SV_decomp_mod(A, X, V, S, work);
	return true;
}
inline bool TimeDomainDataSet_SVD::solve(const TimeDomainDataSet& dataset){
	size_t size = 0;
	for(size=0; size<dataset.size() && dataset.get(size).IsAllocated(); ++size){}
	if(size==0) return false;
	data_v.resize(size, TimeDomainData(dataset.get(0).size()));
	//---	matrix allocation	---
	Allocation(data_v[0].size(), size);
	if(m<n) return false;
	//	copy into gsl_matrix* A
	for(size_t i=0; i<m; ++i) for(size_t j=0; j<n; ++j) gsl_matrix_set(A, i, j, dataset[j][i]);

	//---	the condition of m >> n is assumed as m > 100*n
	if(m < 100*n)	//	Golub Reinsch method;
		gsl_linalg_SV_decomp(A, V, S, work);
	else	//	modified Golub Reinsch method; 
		gsl_linalg_SV_decomp_mod(A, X, V, S, work);
	return true;
}

inline void TimeDomainDataSet_SVD::setUncommonData(){
	if(!Allocated) return;
	for(size_t i=0; i<m; ++i){
		common_component[i] = 0.0;
		for(size_t k=0; k<n; ++k){
			data_v[k][i] = 0.0;
			for(size_t j=0; j<n; ++j){
				if(j!=0)
					data_v[k][i] += gsl_matrix_get(A, i, j)*gsl_vector_get(S, j)*gsl_matrix_get(V, k,j);
				if(j==0)
					common_component[i] += gsl_matrix_get(A, i, j)*gsl_vector_get(S, j)*gsl_matrix_get(V, k,j);
			}
		}
	}
}



inline void TimeDomainDataSet_SVD::getResult(std::vector<double>& s_stdv, std::vector<std::vector<double>>& u_stdv, std::vector<std::vector<double>>& v_stdv){
	gslv_to_stdv(s_stdv, S);
	gslm_to_stdv(u_stdv, A);
	gslm_to_stdv(v_stdv, V);
}

inline std::vector<double> TimeDomainDataSet_SVD::getSingularValues(){
	std::vector<double> tmp;
	gslv_to_stdv(tmp, S);
	return tmp;
}
inline std::vector<std::vector<double>> TimeDomainDataSet_SVD::getEigenFunctions(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, A);
	return tmp;
}
inline std::vector<std::vector<double>> TimeDomainDataSet_SVD::getChannelDistribution(){
	std::vector<std::vector<double>> tmp;
	gslm_to_stdv(tmp, V);
	return tmp;
}

inline void TimeDomainDataSet_SVD::gslv_to_stdv(std::vector<double>& dest, const gsl_vector* src){
	dest.resize(src->size);
	for(size_t i=0; i<dest.size(); i++)	dest[i] = gsl_vector_get(src,i);
}

inline void TimeDomainDataSet_SVD::gslm_to_stdv(std::vector<std::vector<double>>& dest, const gsl_matrix* src){
	dest.resize(src->size1);
	for(size_t i=0; i<dest.size(); i++)	{
		dest[i].resize(src->size2);
		for(size_t j=0; j<src->size2; ++j)	dest[i][j] = gsl_matrix_get(src, i, j);
	}
}
inline void TimeDomainDataSet_SVD::solve(std::vector<std::vector<double> > &src){
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
inline void TimeDomainDataSet_SVD::Allocation(size_t m_, size_t n_){
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

inline TimeDomainDataSet_SVD::~TimeDomainDataSet_SVD(void){
	if(Allocated){
		gsl_matrix_free(A);	gsl_vector_free(S);	gsl_matrix_free(V);	
		gsl_matrix_free(X);	gsl_vector_free(work);
		Allocated = false;
	}
}

#endif

