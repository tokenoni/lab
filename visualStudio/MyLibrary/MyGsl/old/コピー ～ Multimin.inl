namespace mygsl{
	//---	constructor no.1	---
	template< class TyF_ >
	inline Multimin::Multimin
		( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err)
	{
		params_name.resize(0);
		MAXiter = 1000;
		process( function, params, params_step, relative_err);
	}
	//---	constructor no.2	---
	template< class TyF_ >
	inline Multimin::Multimin
		( TyF_& function, std::vector<double>& params, double step_ratio, double relative_err)
	{
		params_name.resize(0);
		std::vector<double> params_step(params.size());
		for(size_t i=0; i<params.size(); i++){
			params_step[i] = params[i]*step_ratio;
			if(params_step[i] == 0.0) params_step[i]=step_ratio;}
		MAXiter = 1000;
		process( function, params, params_step, relative_err);
	}
	//---	constructor no.3	---
	template<class TyF_ >
	inline Multimin::Multimin
		( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err, std::vector<std::string> params_name_)
	{
		params_name=params_name_;
		MAXiter = 1000;
		process( function, params, params_step, relative_err);
	}
	//---	constructor no.4	---
	template<class TyF_ >
	inline Multimin::Multimin
		( TyF_& function, std::vector<double>& params, double step_ratio, double relative_err, std::vector<std::string> params_name_)
	{
		std::vector<double> params_step(params.size());
		MAXiter = 1000;
		params_name=params_name_;
		for(size_t i=0; i<params.size(); i++){
			params_step[i] = params[i]*step_ratio;
			if(params_step[i] == 0.0) params_step[i]=step_ratio;}
		process( function, params, params_step, relative_err);
	}
	
	//---	core	---
	template< class TyF_ >
	int inline Multimin::process
	( TyF_& function, std::vector<double>& params, std::vector<double>& params_step, double relative_err)
	{
		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
		gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, params.size());
		gsl_multimin_function minex_func;

		gsl_vector *x  = std2gsl_vector(params);
		gsl_vector *ss = std2gsl_vector(params_step);	//	step size

	//	minex_func.f = reinterpret_cast< double(*)(const gsl_vector*, void*) >(& (func));
		Cfunction< TyF_ > cfunc(function, params.size());	//	functor is reinterpreted as c function
		minex_func.f = cfunc.adopt();

		minex_func.n = params.size();
		minex_func.params = NULL;

		gsl_multimin_fminimizer_set(s,  &minex_func,  x,  ss);
		size_t iter=0, status;
		do{
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if(status) break;

			double size = gsl_multimin_fminimizer_size(s);
			status = gsl_multimin_test_size(size, relative_err);	//	<<	this value should be optimized	

			if(status == GSL_SUCCESS) std::cerr << "converged to minimum at" <<std::endl;
			print_state(iter, s, size);
		}while(status == GSL_CONTINUE && iter<MAXiter);

		params = gsl2std_vector(s->x);
		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free(s);

		return status;
	}

	void inline Multimin::print_state(size_t iter, gsl_multimin_fminimizer *s, double size){
		std::cerr <<"iter: "<<iter<< ",\t size= "<< size <<std::endl;
		if(params_name.size()==s->x->size)
			for(size_t i=0;i<s->x->size;i++)
				std::cerr << "\t"<< params_name[i] + " \t= " <<gsl_vector_get(s -> x , i)<<std::endl;
		else 
			for(size_t i=0;i<s->x->size;i++)	
				std::cerr << "\tp("<< i<<")\t= " <<gsl_vector_get(s -> x , i)<<std::endl;
		return;
	}
};