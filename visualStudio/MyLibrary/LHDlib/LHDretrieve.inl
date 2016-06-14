//---	constructors	---

inline LHDretrieve::LHDretrieve(void){
	data_length = 0;
	comment_size = 256;
	retrieve_index = -1;
	isAllocated = false;
	isAllocated_t = false;
}

inline LHDretrieve::LHDretrieve(const LHDretrieve& obj){
	copy(obj);
}

inline void LHDretrieve::copy(const LHDretrieve& obj){
	if(isAllocated) clear();

	if(obj.isAllocated){
		//---	copy the parameters		---
		retrieve_index = obj.retrieve_index;
		data_length        = obj.data_length;
		comp_length        = obj.comp_length;
		data_length_actual = obj.comp_length;
		param_count        = obj.param_count;
		value_len          = obj.value_len;
		data_type          = obj.data_type;
		image_type         = obj.image_type;
		is_nframe          = obj.is_nframe;
		strcpy(management  , obj.management);
		strcpy(server      , obj.server);
		strcpy(comment     , obj.comment);
	
		//---	copy the data		---
		data = new char [data_length];
		isAllocated = true;
		for(size_t i=0; i<data_length; ++i)	data[i] = obj.data[i];
		if(obj.isAllocated_t){
			data_t = new double [data_length];
			for(size_t i=0; i<data_length; ++i)	data_t[i] = obj.data_t[i];
		}
	}

}

inline LHDretrieve& LHDretrieve::operator = (const LHDretrieve& obj){
	copy(obj);
	return *this;
}


//---	core functions	---

inline bool LHDretrieve::run(const std::string diagname, const size_t shot_number, const size_t channel_number){
	clear();
	
	retrieve_index = retrieveOpen(diagname.c_str(), NULL, shot_number, 1);
	if( retrieve_index <0 ){
		retrieveClose(retrieve_index);
		return false;
	};

	int int_isExist = retrieveChInfo(retrieve_index, 1, &data_length, &comp_length, &param_count, &data_type, &image_type, &value_len, &is_nframe, management, comment, comment_size);
	if( int_isExist <0 ){
		retrieveClose(retrieve_index);
		clear();
		return false;
	}
	data = new char [data_length];
	isAllocated = true;
	int_isExist = retrieveChData(retrieve_index, channel_number, data, data_length, &data_length_actual);
	if( int_isExist <0 ){
		retrieveClose(retrieve_index);
		clear();
		return false;
	}

	getParams(channel_number);
	retrieveClose(retrieve_index);

	isTimingGot = getTiming(diagname, shot_number, channel_number);

	return true;
}

inline void LHDretrieve::clear(){
	if(isAllocated) delete [] data;
	isAllocated = false;
	if(isAllocated_t) delete [] data_t;
	isAllocated_t = false;
}


inline bool LHDretrieve::getParams(const size_t channel_number){
	if(param_count > 0){
		char ** param_namep  = new char* [param_count];
		char ** param_valuep = new char* [param_count];
		int *   param_typep  = new int   [param_count];
		for(size_t i=0; i<param_count; ++i){
			param_namep[ i] = new char[128];
			param_valuep[i] = new char[128];
		}
		int int_isExist = retrieveChParams(retrieve_index, channel_number, param_namep, param_valuep, param_typep);

		param_names.resize( param_count );
		param_values.resize(param_count );
		for(size_t i=0; i<param_count; ++i){
			std::stringstream ss;
			param_names[i] = std::string(param_namep[i]);
			ss << std::string(param_valuep[i]);
			ss >> param_values[i];
			
			delete[] param_namep[i];
			delete[] param_valuep[i];
		}
		delete[] param_namep;
		delete[] param_valuep;
		delete[] param_typep;
		if(int_isExist<0) return false;
		else return true;
	}
	else return true;
}

inline bool LHDretrieve::getTiming(const std::string diagname, const size_t shot_number, const size_t channel_number){
	unsigned short sub_shot = 1;
	unsigned int retShot;						// ショット番号
	unsigned short retSubShot;					// サブショット番号
	unsigned short retChannel;					// 含まれるチャネルデータ数
	short retYear;									// 年
	short retMonth;									// 月
	short retDay;									// 日
	short retHour;									// 時
	short retMinute;									// 分
	short retSecond;									// 秒
	
	char	DTSsource[64], DTShostID[64], DTSmoduleID[64], clockSource[64], clockHostID[64], clockModuleID[64];
	short	DTStriggerChannel, DTSclockChannel;
	char	strDTStriggerChannel[64], strDTSclockChannel[64], ExtOrInt[64], InternalClock[64], SamplingInterval[64], PreSampling[64];
	int		DTSuserDefine, DTStimeArraySize;
	
	char *IndexServer=getenv("INDEXSERVERNAME");					// インデックスサーバ
//	short Fastch=channel_number;								// 最始チャネル番号
//	short Lastch=channel_number;								// 最終チャネル番号
	unsigned int ch=channel_number;
	
	unsigned int	preSamplingNum;
	unsigned int	triggerTiming_ns;
	unsigned int	cycleTime_ns;
	unsigned int	cycleTime_Hz;
	short bSilent=0;
	double start_t, cycle_t;
	int ret =0;
	int TimeoutSec = 10;

	ret = retrieveGetDTSLastChannel(diagname.c_str(), shot_number, sub_shot, IndexServer, &ch, TimeoutSec);
	if(ret) return false;

	if(isAllocated_t)delete[] data_t;
	isAllocated_t = false;

	//---	 get link information from parameters ----
	ch = channel_number;
	int param_count;
	char strParameter[2048];
	double	dClockCycle, dStartTiming;
	short ch_short = channel_number;
	int rd = retrieveOpenWait(diagname.c_str(),IndexServer,shot_number,sub_shot,TimeoutSec);
	ret = retrieveGetParameterString(rd, ch, &param_count, strParameter);
	if(ret) return false;

	ret = retrieveGetDTSInfoFromParams(param_count, strParameter, &retShot, &retSubShot, 
								DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel,  
								clockSource, clockHostID, clockModuleID, strDTSclockChannel,
								&DTSuserDefine, &DTStimeArraySize, IndexServer, ExtOrInt,
								InternalClock, SamplingInterval, PreSampling);

	if(ret) return false;
//	data_t = new double [size()]; isAllocated_t =true;
	data_t = new double [DTStimeArraySize]; isAllocated_t =true;
	ret=retrieveGetDTSDataDBL2(IndexServer, &ch_short, (double*)data_t,
													DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel, 
													clockSource, clockHostID, clockModuleID, strDTSclockChannel,
													&DTSuserDefine, DTStimeArraySize, retShot, retSubShot,
													ExtOrInt, InternalClock, SamplingInterval, PreSampling, bSilent,
													&dClockCycle, &dStartTiming);
/*	for(size_t i=0; i<DTStimeArraySize; ++i)
		data_t[i] = dStartTiming + dClockCycle*i;
*/	if(ret) return false;
	else return true;
}

inline double LHDretrieve::getParameter(const std::string param_name)const{
	for(size_t i=0; i<param_names.size();++i)
		if(param_names[i] == param_name) return param_values[i];
	return 0.0;
}

inline unsigned short LHDretrieve::getAsUnsignedShort_with_LittleEndianness( const size_t i)const{
	return (unsigned short)(((data[2*i+1] & 0x00FF)<<8) | (data[2*i] & 0x00FF));
};
inline unsigned short LHDretrieve::getAsUnsignedShort_with_BigEndianness( const size_t i)const{
	return (unsigned short)(((data[2*i] & 0x00FF)<<8) | (data[2*i+1] & 0x00FF));
};
inline short LHDretrieve::getAsShort_with_LittleEndianness( const size_t i)const{
	return (short)(((data[2*i+1] & 0x00FF)<<8) | (data[2*i] & 0x00FF));
};
inline short LHDretrieve::getAsShort_with_BigEndianness( const size_t i)const{
	return (short)(((data[2*i] & 0x00FF)<<8) | (data[2*i+1] & 0x00FF));
};
