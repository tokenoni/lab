

//------- get retrieve --------
inline bool LHDretrieve_dataset::getRetrieve(const std::string diagname, const size_t shotnum, const size_t chstart, const size_t chnum, const double bgtime){
	size_t count=0;
	resize(chnum);
	for(size_t ch=chstart; ch< chstart + chnum; ++ch){
			if(diagname == "MultiChPMT" || diagname == "MultiChPMT2"){
			LHDretrieve_MultiChPMT data;
			data.retrieve(diagname, shotnum, ch, bgtime);
			data_v[count].SmartCopy(data);
		}else{
			LHDretrieve_data data;
			data.retrieve(diagname, shotnum, ch, bgtime);
			data_v[count].SmartCopy(data);
		}
		if(data_v[count].IsAllocated()) ++count;
	}
	resize(count);
	if(count ==0) return false;
	else return true;
}

inline bool LHDretrieve_dataset::getRetrieveOne(const std::string diagname, const size_t shotnum, const size_t ch, const double bgtime){
	data_v.resize(data_v.size()+1);
	if(diagname == "MultiChPMT" || diagname == "MultiChPMT2"){
		LHDretrieve_MultiChPMT data;
		data.retrieve(diagname, shotnum, ch, bgtime);
		data_v[data_v.size()-1].SmartCopy(data);
	}else{
		LHDretrieve_data data;
		data.retrieve(diagname, shotnum, ch, bgtime);
		data_v[data_v.size()-1].SmartCopy(data);
	}
	if(!data_v[data_v.size()-1].IsAllocated()){
		data_v.resize(data_v.size()-1);
		return false;
	}else
		return true;
}
