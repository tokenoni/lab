inline LHDretrieve_data::LHDretrieve_data(void) : TimeDomainData(){	channel = 0;}
inline LHDretrieve_data::LHDretrieve_data(const size_t num): TimeDomainData(num){	channel = 0;}

inline LHDretrieve_data::LHDretrieve_data(const LHDretrieve_data& obj){	copy(obj);	channel = 0;}

inline LHDretrieve_data& LHDretrieve_data::copy(const LHDretrieve_data& obj){
	copy(obj);
	if(obj.IsAllocated()){
		raw_data = obj.getRawData();
		channel  = obj.getCh();
	}
	return *this;
}

inline LHDretrieve_data& LHDretrieve_data::operator=(const LHDretrieve_data& obj){	return copy(obj);}

//-------------------------------------------------//
//												   //
//				data preparation 				   //
//												   //
//-------------------------------------------------//
inline bool LHDretrieve_data::retrieve(const std::string& diagname, const size_t shotnum, const size_t chnum, const double bgtime){
	bool isExist = raw_data.run(diagname, shotnum, chnum);
	if(!isExist) return false;
	
	getDataFromRawData();
	raw_data.clear();
	channel = chnum;

	getBackground(bgtime);
	return true;
}

inline void LHDretrieve_data::getDataFromRawData(const void* params){
	Allocate(raw_data.size());
	//-- data copy --
	for(size_t i=0; i<size();++i)
		data[i] = 1.0*raw_data[i];
	
	//-- get rate and start_t --
	if(raw_data.IsAllocated_t()){
		size_t one_point     = raw_data.size()/4-1;
		size_t another_point = raw_data.size()/2-1;
		
		size_t t_point_num = another_point - one_point;
		rate = (1.0*t_point_num)/(raw_data.get_t(another_point) - raw_data.get_t(one_point));
		start_t = raw_data.get_t(one_point) - 1.0*one_point/rate;
		end_t = start_t + 1.0*size()/rate;
		isPrepared_t = true;
	}
	else isPrepared_t = false;
}

inline void LHDretrieve_data::getBackground(const double bgtime){
	size_t bg_num = (size_t)floor(bgtime*rate+0.5);
	double bg = 0.0;
	for(size_t i=0; i<bg_num && i<size(); ++i) bg += data[i];
	if(bg_num !=0) bg /= bg_num;
	for(size_t i=0; i<size(); ++i)	data[i] -= bg;
}
