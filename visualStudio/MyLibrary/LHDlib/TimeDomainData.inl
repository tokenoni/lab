//-------------------------------------------------//
//												   //
//		constructors and related functions 		   //
//												   //
//-------------------------------------------------//

inline TimeDomainData::TimeDomainData(void){
	isAllocated = false;
	isPrepared_t = false;
}
inline TimeDomainData::TimeDomainData(const size_t num){
	isAllocated = false;
	isPrepared_t = false;
	Allocate(num);
}
inline TimeDomainData::~TimeDomainData(void){clear();}

inline TimeDomainData::TimeDomainData(const TimeDomainData& obj){
	isAllocated = false;
	isPrepared_t = false;
	copy(obj);
}
inline TimeDomainData& TimeDomainData::operator = (const TimeDomainData& obj){
	return copy(obj);
}
inline TimeDomainData& TimeDomainData::copy(const TimeDomainData& obj){
	if(obj.IsAllocated()){
		Allocate(obj.size());
		for(size_t i=0; i<size(); ++i) data[i] = obj.get(i);
		if(obj.IsPrepared_t()) rate    = obj.getRate();		else rate = 0.0;
		if(obj.IsPrepared_t()) start_t = obj.getStart_t();	else start_t = 0.0;
		if(obj.IsPrepared_t()) end_t   = obj.getEnd_t();	else end_t   = 0.0;
		if(obj.IsPrepared_t()) isPrepared_t = obj.IsPrepared_t();
	}
	else{
		clear();
	}
	return *this;
};

inline void TimeDomainData::Allocate(const size_t num){
	if(isAllocated && num==size()) return;
	else{
		if(isAllocated) delete[] data;
		data = new double [ num ];
		isAllocated = true;
		data_num =num;
		return;
	}
}

inline void TimeDomainData::clear(){
	if(isAllocated) delete [] data;
	isAllocated = false;
	isPrepared_t= false;
	return;
}


//--- smart copy function. the data in the copied object cannot be used any more. ---
//--- this function will be used in the inheriting object ---
inline void TimeDomainData::SmartCopy(TimeDomainData& obj){
	clear();
	if(!obj.IsAllocated()) return;
	else{	
		data_num = obj.size();
		isAllocated = true;
		data = obj._SmartCopied();

		if(obj.IsPrepared_t()) rate    = obj.getRate();		else rate = 0.0;
		if(obj.IsPrepared_t()) start_t = obj.getStart_t();	else start_t = 0.0;
		if(obj.IsPrepared_t()) end_t   = obj.getEnd_t();	else end_t   = 0.0;
		if(obj.IsPrepared_t()) isPrepared_t = obj.IsPrepared_t();
	}
}

//--- function for smart copy. Don't use this function!! ---
inline double* TimeDomainData::_SmartCopied(){
	if(!IsAllocated())return NULL;
	else{
		isAllocated = false;
		double* data_address = data;
		data = NULL;
		return data_address;
	}
}


//-------------------------------------------------//
//												   //
//				Tip functions   				   //
//												   //
//-------------------------------------------------//

inline size_t TimeDomainData::getIndexFromTime(const double time_)const{
	int index=0;
	if(isPrepared_t) index = (time_ - start_t)*rate;
	else             index = 0;

	if(index<0) return 0;
	else if(index>= size()) return size()-1;
	else return (size_t)index;
};

inline void TimeDomainData::setTiming(const double rate_, const double start_t_){
	rate = rate_;
	start_t = start_t_;
	end_t = start_t + 1.0*data_num / rate;
	isPrepared_t = true;
}

inline bool TimeDomainData::getAvgCopy(TimeDomainData& dest, const size_t avgnum)const{
	size_t num = (size_t)floor(1.0*size()/avgnum);
	TimeDomainData tmp(num);
	for(size_t i=0; i<num && i*avgnum<size(); ++i){
		tmp[i] = getAvg(i*avgnum, avgnum);
	}
	tmp.setTiming(rate/avgnum, start_t + 0.5*avgnum/rate);
	dest.SmartCopy(tmp);
}

//-------------------------------------------------//
//												   //
//				data operation  				   //
//												   //
//-------------------------------------------------//
inline double TimeDomainData::getAvg(const int start_index_, const size_t avgnum_)const{
	double val = 0.0;
	for(int i=start_index_; i<start_index_+avgnum_; ++i)
		val += data[i];
	return val/=avgnum_;
}

inline double TimeDomainData::getAvg_t(const int start_index, const size_t avgnum)const{
	return start_t + (start_index+0.5*avgnum)/rate;
}

inline std::vector<double> TimeDomainData::getAvg(const size_t avgnum)const{
	size_t num = (size_t)floor(1.0*size()/avgnum);
	std::vector<double> rslt(num, 0.0);
	for(size_t i=0; i<num; ++i)
		rslt[i] = getAvg(i*avgnum, avgnum);
	return rslt;
}

inline std::vector<double> TimeDomainData::getAvg_t(const size_t avgnum)const{
	size_t num = (size_t)floor(1.0*size()/avgnum);
	std::vector<double> rslt(num, 0.0);
	for(size_t i=0; i<num; ++i)
		rslt[i] = getAvg_t(i*avgnum, avgnum);
	return rslt;
}