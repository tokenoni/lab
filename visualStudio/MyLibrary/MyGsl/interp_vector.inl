namespace mygsl{
	//---	constructors	---
	inline interp_vector::interp_vector(const interp_vector& obj){
		data.resize(obj.size());
		for(size_t i=0; i<size();++i)
			data[i] = obj.data[i];
#ifdef CCGSL_VECTOR_HPP
		vdata = gsl::vector(size());
#endif
	}

	inline interp_vector& interp_vector::operator = (const interp_vector&obj){
		data.resize(obj.size());
		for(size_t i=0; i<size();++i)
			data[i] = obj.data[i];
#ifdef CCGSL_VECTOR_HPP
		vdata = gsl::vector(size());
		return *this;
#endif
	}

	//---	functions	---
#ifdef CCGSL_VECTOR_HPP
	inline const gsl::vector& interp_vector::operator () (double x){
		for(size_t i=0; i<size();++i)
			vdata.set(i, data[i].get(x));
		return vdata;
	}
#endif
};
