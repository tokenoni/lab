namespace mygsl{

	template<class TyF_ >
	inline MultiscaleMultimin::MultiscaleMultimin
		( TyF_& function, std::vector<double>& params, std::vector<size_t>& params_order, double step_ratio, double relative_err, bool isprint_)
	{

		std::vector<std::vector<size_t>> free_params_index(1,std::vector<size_t>(0));
		std::vector<size_t> index_vs_order(params.size());
		size_t order=0;
		for(size_t i=0; i<params.size(); ++i){
			for(size_t j=0; j<params_order.size(); ++j){
				if(params_order[j]==i){
					free_params_index[order].push_back(j);
					index_vs_order[j] = order;
				}
			}
			if(free_params_index[order].size()>0){
				++order;
				free_params_index.push_back(std::vector<size_t>(0));
			}
		}

		std::vector< OnescaleMultimin < TyF_ >> multiscale_multimin(order);
		multiscale_multimin[0].setOriginal( &function, &params, free_params_index[0]);

		for(size_t i=1; i<multiscale_multimin.size(); ++i)
			multiscale_multimin[i].setScale( & multiscale_multimin[i-1], &params, free_params_index[i]);
	
		std::vector<double> free_params_top;
		for(size_t i=0; i<params.size(); ++i){
			if(index_vs_order[i] == order-1)
				free_params_top.push_back(params[i]);
		}

		Multimin multimin(multiscale_multimin[order-2], free_params_top);
	}
}