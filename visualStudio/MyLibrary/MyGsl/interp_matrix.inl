namespace mygsl{
	//---	constructors	---
	inline interp_matrix::interp_matrix(const interp_matrix& obj){
		data.resize(obj.size());
		for(size_t i=0; i<size();++i)
			data[i] = obj.data[i];
#ifdef CCGSL_MATRIX_HPP
		mdata = gsl::matrix(size(), obj.data[0].size());
#endif
	}

	inline interp_matrix& interp_matrix::operator = (const interp_matrix& obj){
		data.resize(obj.size());
		for(size_t i=0; i<size();++i)
			data[i] = obj.data[i];
#ifdef CCGSL_MATRIX_HPP
		mdata = gsl::matrix(size(), obj.data[0].size());
		return *this;
#endif
	}
	//---	functions	---
	inline void interp_matrix::resize(const size_t size1, const size_t size2){
		data.resize(size1, interp_vector(size2));
#ifdef CCGSL_MATRIX_HPP
		mdata = gsl::matrix(size1, size2);
#endif
	}

#ifdef CCGSL_MATRIX_HPP
	inline const gsl::matrix& interp_matrix::operator () (double x){
		for(size_t i=0; i<size();++i){
			for(size_t j=0; j<data[i].size();++j)
				mdata.set(i, j, data[i][j].get(x));
		}
		return mdata;
	}

	inline bool interp_matrix::set(const std::vector<double>& xsrc, const std::vector<gsl::matrix>& msrc){
		if(xsrc.size() != msrc.size()) return false;
		std::vector<double> mtmp(xsrc.size());
		
		resize(msrc[0].size1(),msrc[0].size2());
		for(size_t i=0; i<msrc[0].size1(); ++i){
			for(size_t j=0; j<msrc[0].size2(); ++j){
				for(size_t k=0; k<xsrc.size(); ++k)		mtmp[k] = msrc[k].get(i,j);
				interp mij_interp;
				data[i][j].set(xsrc, mtmp);
			}
		}
		return true;
	}

#endif



};