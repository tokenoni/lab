namespace mygsl{
//-----------------------------------------------------//
//													   //
//		class for basis functions interpolation		   //
//													   //
//-----------------------------------------------------//
	inline bspline::bspline(void){allocated = false;	allocated_workspace = false;	}
	
	inline bspline::bspline(const size_t size_, const size_t nbreak_, const size_t k_)
	{
		allocated = false;	allocated_workspace = false;
		MemoryControl(size_, nbreak_, k_);
	}
	
	inline bspline::bspline(const bspline& obj){
		allocated = false;	allocated_workspace = false;	
		if(obj.allocated){
			Allocate(obj.num);
			copyData(obj);
			if(obj.allocated_workspace){
				AllocateWorkspace(obj.nbreak, obj.k);
				copyWorkspace(obj);
			}
		}

	}
	inline bspline& bspline::operator = (const bspline& obj){
		if(allocated) FreeMemory();
		if(allocated_workspace) FreeMemoryWorkspace();
		if(obj.allocated){
			Allocate(obj.num);
			copyData(obj);
			if(obj.allocated_workspace){
				AllocateWorkspace(obj.nbreak, obj.k);
				copyWorkspace(obj);
			}
		}
		return *this;
	}
	//----------	memory controling functions		------------
	inline void bspline::MemoryControl(const size_t size_, const size_t nbreak_, const size_t k_ ){
		if(!allocated || num != size_){
			FreeMemory();
			Allocate(size_);
		}
		if(!allocated_workspace || (nbreak !=nbreak_ || k!=k_)){
			FreeMemoryWorkspace();
			AllocateWorkspace(nbreak_, k_);
		}
	}
		
	inline void bspline::FreeMemory(){
		if(allocated){
			delete [] pos;	delete [] data;	delete [] w;
		}
		allocated = false;
	}

	inline void bspline::Allocate(const size_t size_){
		num = size_;
		if(!allocated){
			pos = new double [size_];
			data = new double [size_];
			w = new double [size_];
		}
		allocated = true;
	}
	inline void bspline::FreeMemoryWorkspace(){
		if(allocated_workspace){
			gsl_bspline_free(bw);
			gsl_bspline_deriv_free(dbw);
			gsl_vector_free(Bk);
			gsl_matrix_free(dBk);
			gsl_vector_free(c);
		}
		allocated = false;
	}
		
	inline void bspline::AllocateWorkspace(const size_t nbreak_, const size_t k_){
		k = k_;
		nbreak = nbreak_;
		if(!allocated_workspace){
			bw  = gsl_bspline_alloc(k, nbreak);
			dbw = gsl_bspline_deriv_alloc(k);
			Bk  = gsl_vector_calloc(k);
			dBk = gsl_matrix_calloc(k,k);
			c   = gsl_vector_calloc(nbreak - 2 + k);
		}
		allocated_workspace = true;
	}

	inline void bspline::copyData(const bspline& obj){
		for(size_t i=0; i<num; ++i){
			pos[i] = obj.pos[i];
			data[i] = obj.data[i];
			w[i] = obj.w[i];
		}
	}
	inline void bspline::copyWorkspace(const bspline& obj){
		//	fow gsl_bspline_workspace
		bw->k      = obj.bw->k;
		bw->km1    = obj.bw->km1;
		bw->l      = obj.bw->l;
		bw->nbreak = obj.bw->nbreak;
		bw->n      = obj.bw->n;
		gsl_vector_memcpy(bw->knots  ,  obj.bw->knots  );
		gsl_vector_memcpy(bw->deltal ,  obj.bw->deltal );
		gsl_vector_memcpy(bw->deltar ,  obj.bw->deltar );
		gsl_vector_memcpy(bw->B      ,  obj.bw->B      );
		
		//	fow gsl_bspline_derive_workspace
		dbw->k  = obj.dbw->k;
		gsl_matrix_memcpy(dbw->A , obj.dbw->A);
		gsl_matrix_memcpy(dbw->dB, obj.dbw->dB);

		gsl_vector_memcpy(Bk, obj.Bk);
		gsl_matrix_memcpy(dBk, obj.dBk);
		gsl_vector_memcpy(c, obj.c);
	}

	//---	get bsplined value at x	---
	inline double bspline::get(const double x)const{
		size_t istart, iend;
		gsl_bspline_eval_nonzero(x, Bk, &istart, &iend, bw);
		double val=0.0;
		for(size_t i=0; i<k; ++i)
			val += gsl_vector_get(c, i+istart)*gsl_vector_get(Bk,i);
		return val;
	}

	inline double bspline::getderiv(const double x)const{
		size_t istart, iend;
		gsl_bspline_deriv_eval_nonzero(x, 1, dBk, &istart, &iend, bw, dbw);
		double val=0.0;
		for(size_t i=0; i<k; ++i)
			val += gsl_vector_get(c, i+istart)*gsl_matrix_get(dBk,i,1);
		return val;
	}

	//---	setting functions	---
	inline void bspline::set(const size_t i, const double xi, const double yi, const double sigma_){
		pos[i] = xi;	data[i] = yi;	
		if(sigma_ !=0) w[i] = 1.0/sigma_;
		else w[i]=1.0;
	}

	template < class Tyv_ >
	inline void bspline::set(const Tyv_& x, const Tyv_& y, const size_t nbreak_, const size_t k_, KnotsGettingMethod knots_method_){
		MemoryControl(x.size(), nbreak_, k_);
		for(size_t i=0; i<size(); ++i)
			set(i, x[i], y[i]);
		setInterp(knots_method_);
	}
	inline void bspline::set_y_zero(){
		for(size_t i=0; i<size();++i) data[i]=0.0;
	}

	inline void bspline::setInterp(KnotsGettingMethod knots_method_, HoldingMethod holdingMethod_  ){
		//	sort at first
		knots_method = knots_method_;
		holdingMethod = holdingMethod_;
		sort();
		switch (knots_method_){
		case Uniform:
			getKnotsUniform();
			setBspline(holdingMethod_  );
			break;
		case Greville:
			getKnotsGreville();
			setBspline(holdingMethod_  );
			break;
		case Automatic:
			setBsplineAutomatic();
			break;
		};
	}
	inline void bspline::setInterp(std::vector<double> breakpoints_, HoldingMethod holdingMethod_ ){
		holdingMethod = holdingMethod_;
		gsl_vector *breakpoints;
		breakpoints = gsl_vector_calloc(breakpoints_.size());
		for(size_t i=0; i<breakpoints_.size();++i)
			gsl_vector_set(breakpoints, i, breakpoints_[i]);
		gsl_bspline_knots(breakpoints,bw);
		gsl_vector_free(breakpoints);
		setBspline(holdingMethod_  );
	}

	inline void bspline::getKnotsUniform(){
		gsl_bspline_knots_uniform(pos[0], pos[size()-1], bw);
//		gsl_bspline_knots_uniform(pos[10], pos[size()-11], bw);
	}
	inline void bspline::getKnotsGreville(){
	//	---	under constructions	---
		gsl_bspline_knots_uniform(pos[0], pos[size()-1], bw);
	}

	inline void bspline::sort(){
		//	sort
		size_t *p = new size_t[size()];
		double *x = new double[size()];
		double *y = new double[size()];
		double *dy = new double[size()];

		for(size_t i=0; i<num; ++i){	x[i] = pos[i];	y[i] = data[i];	dy[i] = w[i];}

		gsl_sort_index(p, x, 1, num);
		for(size_t i=0; i<num; ++i){	pos[i] = x[p[i]];	data[i] = y[p[i]];	w[i] = dy[p[i]];}

		delete [] p;
		delete [] x;
		delete [] y;
		delete [] dy;
	}

	inline double bspline::setBspline(HoldingMethod holdingMethod_ ){
		switch(holdingMethod_){
		case NoHoldings:					return setBspline_NoHoldings();				break;
		case Flat_at_FirstEdge:				return setBspline_Flat_at_FirstEdge();			break;
		case Flat_at_LastEdge:				return setBspline_Flat_at_LastEdge();			break;
		case Flat_at_BothEdge:				return setBspline_Flat_at_BothEdge();			break;
		case Zero_at_FirstEdge:				return setBspline_Zero_at_FirstEdge();			break;
		case Zero_at_LastEdge:				return setBspline_Zero_at_LastEdge();			break;
		case Zero_at_BothEdge:				return setBspline_Zero_at_BothEdge();			break;
		case ZeroAndFlat_at_FirstEdge:		return setBspline_ZeroAndFlat_at_FirstEdge();	break;
		case ZeroAndFlat_at_LastEdge:		return setBspline_ZeroAndFlat_at_LastEdge();	break;
		case ZeroAndFlat_at_BothEdge:		return setBspline_ZeroAndFlat_at_BothEdge();	break;
		case ZeroAndFlat_at_FirstEdge_Zero_at_LastEdge:	return setBspline_ZeroAndFlat_at_FirstEdge_Zero_at_LastEdge();	break;
		case ZeroAndFlat_at_FirstEdge_Flat_at_LastEdge:	return setBspline_ZeroAndFlat_at_FirstEdge_Flat_at_LastEdge();	break;
		case ZeroAndFlat_at_LastEdge_Zero_at_FirstEdge:	return setBspline_ZeroAndFlat_at_LastEdge_Zero_at_FirstEdge();	break;
		case ZeroAndFlat_at_LastEdge_Flat_at_FirstEdge:	return setBspline_ZeroAndFlat_at_LastEdge_Flat_at_FirstEdge();	break;

		}
		return 0.0;
	}
	inline double bspline::setBspline_NoHoldings(){
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size());
		gsl_matrix *cov = gsl_matrix_calloc(coef_size(), coef_size());
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size());
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j)
				gsl_matrix_set(X, i, j, gsl_vector_get(Bk,j-jstart));
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c, cov, &chisq, mw);
		
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBsplineAutomatic(){
		BsplineAutomaticKnotsDetermine automatic
			(size(), coef_size(), nbreak, pos, data, w, Bk, c,bw);
		
		std::vector<double> breakpts_angle(nbreak-2);
		for(size_t i=0; i<nbreak-2; ++i)
			breakpts_angle[i] = M_PI/(nbreak-3)*(i+1);
		Multimin multimin(automatic, breakpts_angle);

		automatic(breakpts_angle);
	}


	inline std::vector<double> bspline::getOriginal_x()const{
		std::vector<double> rslt(size());
		for(size_t i=0; i<size(); ++i)
			rslt[i] = pos[i];
		return rslt;
	}
	inline std::vector<double> bspline::getOriginal_y()const{
		std::vector<double> rslt(size());
		for(size_t i=0; i<size(); ++i)
			rslt[i] = data[i];
		return rslt;
	}



	//---	class only for used in the automatically knots dentemining in Bspline class	---
	inline BsplineAutomaticKnotsDetermine::BsplineAutomaticKnotsDetermine
		(const size_t data_size_, const size_t coef_size_, const size_t break_size_, 
		const double *x_, const double *y_, const double *w_, 
		gsl_vector *Bk_, gsl_vector *c_, gsl_bspline_workspace *bw_)
		: data_size(data_size_), coef_size(coef_size_), nbreak(break_size_)
	{
		x = x_;	y=y_; w=w_; bw=bw_;
		X = gsl_matrix_calloc(data_size, coef_size);
		cov = gsl_matrix_calloc(coef_size, coef_size);
		mw = gsl_multifit_linear_alloc(data_size,coef_size);
		breakpts = gsl_vector_calloc(nbreak);
		Bk = Bk_;
		c = c_;
		p = gsl_permutation_calloc(nbreak-2);
	}

	inline double BsplineAutomaticKnotsDetermine::operator () (std::vector<double>& breakpts_angle){
		for(size_t i=0; i<nbreak-2; ++i){
			double ratio_tmp = 0.5*(cos(breakpts_angle[i])+1.0);
			gsl_vector_set(breakpts, i, ratio_tmp);
		}
		//---	sorting	---
		gsl_sort_index(p->data, breakpts->data, breakpts->stride, nbreak-2);
		gsl_vector_set(breakpts, 0, x[0]);
		gsl_vector_set(breakpts, nbreak-1, x[data_size-1]);
		for(size_t i=0; i<nbreak-2; ++i){
			double ratio_tmp = 0.5*(cos(breakpts_angle[gsl_permutation_get(p,i)])+1.0);
			gsl_vector_set(breakpts, i+1, x[0] + ratio_tmp*(x[data_size-1] - x[0]));
		}
		//---	set to the bspline workspace	---
		gsl_bspline_knots(breakpts, bw);

		size_t jstart, jend;
		for(size_t i=0; i<data_size; ++i){
			gsl_bspline_eval_nonzero(x[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=0; j<jend-jstart+1; ++j)
				gsl_matrix_set(X, i, j+jstart, gsl_vector_get(Bk,j));
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, data_size);
		gsl_vector_const_view y_view = gsl_vector_const_view_array(y,  data_size);

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c, cov, &chisq, mw);
		return chisq;
	}

	inline BsplineAutomaticKnotsDetermine::~BsplineAutomaticKnotsDetermine(void){
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);	
		gsl_vector_free(breakpts);
		gsl_permutation_free(p);
	}




	//----------------------------------------------------------------------//
	//																		//
	//				setBspline methods for edge holding 					//
	//																		//
	//----------------------------------------------------------------------//

	inline double bspline::setBspline_Flat_at_FirstEdge(){
		size_t coef_size_new = coef_size()-1;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=0 && j!=1)
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
				else if (j==0)
					gsl_matrix_set(X, i, 0, gsl_vector_get(Bk,0-jstart) + gsl_vector_get(Bk,1-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set(c, 0, gsl_vector_get(c_tmp, 0));
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));

		gsl_vector_free(c_tmp);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_Zero_at_FirstEdge(){
		size_t coef_size_new = coef_size()-1;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=0)
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set(c, 0, 0.0);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));

		gsl_vector_free(c_tmp);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_ZeroAndFlat_at_FirstEdge(){
		size_t coef_size_new = coef_size()-2;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=0 && j!=1)
					gsl_matrix_set(X, i, j-2, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set(c, 0, 0.0);
		gsl_vector_set(c, 1, 0.0);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+2, gsl_vector_get(c_tmp, i));

		gsl_vector_free(c_tmp);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_Flat_at_LastEdge(){
		size_t coef_size_new = coef_size()-1;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=coef_size_new && j!=coef_size_new-1)
					gsl_matrix_set(X, i, j, gsl_vector_get(Bk,j-jstart));
				else if (j==coef_size_new )
					gsl_matrix_set(X, i, coef_size_new-1, gsl_vector_get(Bk,coef_size_new-1-jstart) + gsl_vector_get(Bk,coef_size_new-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i, gsl_vector_get(c_tmp, i));
		gsl_vector_set(c, coef_size_new, gsl_vector_get(c_tmp, coef_size_new-1));

		gsl_vector_free(c_tmp);
		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_Zero_at_LastEdge(){
		size_t coef_size_new = coef_size()-1;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=coef_size_new)
					gsl_matrix_set(X, i, j, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set(c, coef_size_new, 0.0);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_ZeroAndFlat_at_LastEdge(){
		size_t coef_size_new = coef_size()-2;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if(j!=coef_size_new && j!=coef_size_new+1)
					gsl_matrix_set(X, i, j, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_Flat_at_BothEdge(){
		size_t coef_size_new = coef_size()-2;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 && j!=1) && (j!=coef_size_new +1&& j!=coef_size_new))
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
				else if (j==coef_size_new)
					gsl_matrix_set(X, i, coef_size_new-1, gsl_vector_get(Bk,coef_size()-1-jstart) + gsl_vector_get(Bk,coef_size()-2-jstart));
				else if (j==0)
					gsl_matrix_set(X, i, 0, gsl_vector_get(Bk,0-jstart) + gsl_vector_get(Bk,1-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));
		gsl_vector_set(c, coef_size()-1, gsl_vector_get(c_tmp, coef_size_new-1));
		gsl_vector_set(c, 0, gsl_vector_get(c_tmp, 0));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_Zero_at_BothEdge(){
		size_t coef_size_new = coef_size()-2;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 ) && (j!=coef_size_new +1))
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_ZeroAndFlat_at_BothEdge(){
		size_t coef_size_new = coef_size()-4;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 && j!=1) && (j!=coef_size_new +3&& j!=coef_size_new+2))
					gsl_matrix_set(X, i, j-2, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+2, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}

	inline double bspline::setBspline_ZeroAndFlat_at_FirstEdge_Zero_at_LastEdge(){
		size_t coef_size_new = coef_size()-3;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 && j!=1) && (j!=coef_size_new+2))
					gsl_matrix_set(X, i, j-2, gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+2, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	}
	inline double bspline::setBspline_ZeroAndFlat_at_FirstEdge_Flat_at_LastEdge(){
		size_t coef_size_new = coef_size()-3;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 && j!=1) && (j!=coef_size_new +3&& j!=coef_size_new+2))
					gsl_matrix_set(X, i, j-2, gsl_vector_get(Bk,j-jstart));
				else if (j==coef_size_new+2 && coef_size()-1-jstart>=0)
					gsl_matrix_set(X, i, j-3, gsl_vector_get(Bk,coef_size()-1-jstart) + gsl_vector_get(Bk,coef_size()-2-jstart));
			}
		}
//		IGORdata::write_itx(X, "X.itx","XX");
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+2, gsl_vector_get(c_tmp, i));
		gsl_vector_set(c, coef_size_new+1, gsl_vector_get(c_tmp, coef_size_new-1));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	};
	inline double bspline::setBspline_ZeroAndFlat_at_LastEdge_Zero_at_FirstEdge(){
		size_t coef_size_new = coef_size()-3;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
		
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0) && (j!=coef_size_new +1&& j!=coef_size_new+2))
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
			}
		}
//		IGORdata::write_itx(X, "X.itx","XX");
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	};
	inline double bspline::setBspline_ZeroAndFlat_at_LastEdge_Flat_at_FirstEdge(){
		size_t coef_size_new = coef_size()-3;
		gsl_matrix *X = gsl_matrix_calloc(size(), coef_size_new);
		gsl_matrix *cov = gsl_matrix_calloc(coef_size_new, coef_size_new);
		gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(size(),coef_size_new);
		gsl_vector *c_tmp = gsl_vector_calloc(coef_size_new);
		size_t jstart, jend;
		for(size_t i=0; i<size(); ++i){
			gsl_bspline_eval_nonzero(pos[i], Bk, &jstart, &jend, bw);
					
			for(size_t j=jstart; j<jend+1; ++j){
				if((j!=0 && j!=1) && (j!=coef_size_new +1&& j!=coef_size_new+2))
					gsl_matrix_set(X, i, j-1, gsl_vector_get(Bk,j-jstart));
				else if (j==1 && j>jstart)
					gsl_matrix_set(X, i, 0, gsl_vector_get(Bk,j-1-jstart) + gsl_vector_get(Bk,j-jstart));
			}
		}
		gsl_vector_const_view w_view = gsl_vector_const_view_array(w, size());
		gsl_vector_const_view y_view = gsl_vector_const_view_array(data,  size());

		double chisq;
		gsl_multifit_wlinear(X, &w_view.vector, &y_view.vector,c_tmp, cov, &chisq, mw);
		
		gsl_vector_set_zero(c);
		gsl_vector_set(c, 0, gsl_vector_get(c_tmp, 0));
		for(size_t i=0; i<coef_size_new; ++i)
			gsl_vector_set(c, i+1, gsl_vector_get(c_tmp, i));

		gsl_matrix_free(X);
		gsl_vector_free(c_tmp);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
		return chisq;
	};

};