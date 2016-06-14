
//-----------		constructors			-----------
inline TimeDomainDataSet::TimeDomainDataSet(void){
	data_v.resize(0, TimeDomainData());
}
inline TimeDomainDataSet::TimeDomainDataSet(const size_t chnum){
	data_v.resize(chnum);
}
inline TimeDomainDataSet::TimeDomainDataSet(const TimeDomainDataSet& obj){
	copy(obj);
}
inline TimeDomainDataSet& TimeDomainDataSet::copy(const TimeDomainDataSet& obj){
	data_v = obj.data_v;
	return *this;
}
inline TimeDomainDataSet& TimeDomainDataSet::operator =(const TimeDomainDataSet& obj){
	return copy(obj);
}
inline void TimeDomainDataSet::eraseUnallocatedData(){
	for(size_t ch=0; ch<size(); ++ch){
		if(!data_v[ch].IsAllocated()){
			data_v.erase(data_v.begin()+ch);
			--ch;
		}
	}
}

inline void TimeDomainDataSet::push_back(const TimeDomainData& data){
	if(data.IsAllocated())
		data_v.push_back(data);
}

inline double TimeDomainDataSet::getAvg(const int t_index, const int ch_index, const size_t t_avgnum, const size_t ch_avgnum)const{
	double val = 0.0;
	size_t d_ch_index = (size_t)floor(1.0*(ch_avgnum-1)/2);
	size_t count=0;
	for(int ch = ch_index - d_ch_index; ch<ch_index+ch_avgnum && ch<size(); ++ch){
		if(data_v[ch].IsAllocated()){
			val += data_v[ch].getAvg(t_index, t_avgnum); 
			++count;
		}
	}
	return val/=count;
}


//---	tip functions ---
inline size_t TimeDomainDataSet::t_size()const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsAllocated()) return data_v[ch].size();
	}
	return 0;
}

inline double TimeDomainDataSet::getRate()const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsPrepared_t()) return data_v[ch].getRate();
	}
	return 0.0;
}
inline double TimeDomainDataSet::getStart_t()const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsPrepared_t()) return data_v[ch].getStart_t();
	}
	return 0.0;
}
inline double TimeDomainDataSet::getEnd_t()const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsPrepared_t()) return data_v[ch].getEnd_t();
	}
	return 0.0;
}

inline double TimeDomainDataSet::get_t(const size_t index, const size_t avgnum)const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsPrepared_t()) return data_v[ch].getAvg_t(index, avgnum);
	}
	return 0.0;
}

inline size_t TimeDomainDataSet::getIndexFromTime(const double time)const{
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsPrepared_t()) return data_v[ch].getIndexFromTime(time);
	}
	return 0;
}


inline void TimeDomainDataSet::setTiming(const double rate, const double start_t){
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsAllocated()) data_v[ch].setTiming(rate, start_t);
	}
	return;
}

//---	outputfunctions ---
inline std::vector<std::vector<double>> TimeDomainDataSet::getAvg_stdvv(const size_t tsmth, const size_t chsmth_)const{
	std::vector<std::vector<double>> rslt(0);
	size_t chsmth = chsmth_;
	if(chsmth==0) ++chsmth;
	for(size_t ch=0; ch<size(); ch+=chsmth){
		if(data_v[ch].IsAllocated()){
			rslt.push_back(std::vector<double>((size_t)floor(1.0*data_v[ch].size()/tsmth), 0.0));
			for(size_t i=0; i<rslt[rslt.size()-1].size();++i)
				rslt[rslt.size()-1][i] = getAvg(i*tsmth, ch, tsmth, chsmth);
		}
	}
	return rslt;
}

inline std::vector<std::vector<double>> TimeDomainDataSet::get_stdvv()const{
	std::vector<std::vector<double>> rslt(0);
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsAllocated()){
			rslt.push_back(std::vector<double>(data_v[ch].size()));
			for(size_t i=0; i<rslt[rslt.size()-1].size();++i)
				rslt[rslt.size()-1][i] = get(ch).get(i);
		}
	}
	return rslt;
}

inline std::vector<std::vector<double>> TimeDomainDataSet::get_stdvv(const double start_t, const double end_t)const{
	std::vector<std::vector<double>> rslt(0);
	for(size_t ch=0; ch<size(); ++ch){
		if(data_v[ch].IsAllocated()){
			size_t start_p = data_v[ch].getIndexFromTime(start_t);
			size_t end_p   = data_v[ch].getIndexFromTime(end_t);
			rslt.push_back(std::vector<double>(end_p - start_p + 1));
			for(size_t i=0; i<end_p-start_p+1;++i)
				rslt[rslt.size()-1][i] = get(ch).get(i+start_p);
		}
	}
	return rslt;
}

inline std::vector<double> TimeDomainDataSet::get_t_stdv()const{
	std::vector<double> rslt(t_size());
	for(size_t i=0; i<t_size(); ++i)
		rslt[i] = get_t(i);
	return rslt;
}
inline std::vector<double> TimeDomainDataSet::get_t_stdv(const double start_t, const double end_t)const{
	int start_p = getIndexFromTime(start_t);
	int end_p   = getIndexFromTime(end_t);
	std::vector<double> rslt(abs(end_p - start_p+1));
	for(size_t i=0; i<rslt.size();++i)
		rslt[i] = get_t(i+start_p);
	return rslt;
}
