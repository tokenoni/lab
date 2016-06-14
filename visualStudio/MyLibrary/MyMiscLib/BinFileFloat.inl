namespace mymisclib{
	inline BinFileFloat::BinFileFloat(void)
	{
	}


	inline BinFileFloat::~BinFileFloat(void)
	{
	}

	inline void BinFileFloat::add(const std::string param_name, const double param_value){
		param_names.push_back(param_name);
		param_values.push_back(param_value);
	}

	inline bool BinFileFloat::open_for_write(const std::string filename){
		ofs.open(filename.c_str(), std::ofstream::out);
		if(!ofs.good()) return false;
	
		ofs << "number of rows for header = " << param_names.size() + 2 << std::endl << std::endl;
		for(size_t i=0; i<param_names.size(); ++i)
			ofs << param_names[i] << " = " << param_values[i] << std::endl;

		ofs.close();
		ofs.open(filename.c_str(), std::ofstream::app | std::ofstream::binary);

		if(!ofs.good()) return false;
		return true;
	}

	inline bool BinFileFloat::add_bindata(float data){
		ofs.write(reinterpret_cast<char*>( &data ), sizeof data);
		return true;
	}
	inline void BinFileFloat::finish_writing(){
		if(ofs.good()) ofs.close();
	}

	inline bool BinFileFloat::open_for_read(const std::string filename){
		ifs.open(filename, std::ifstream::in | std::ifstream::binary);
		if(!ifs.good()) return false;
		
		std::string sdata;
		getline(ifs,sdata);
		if(sdata.find("number of rows for header = ")== std::string::npos)
			return false;
		
		size_t num_row = atoi(sdata.substr(sdata.find("number of rows for header = ")+28).c_str());
		param_names = std::vector<std::string>(num_row);
		param_values = std::vector<double>(num_row);
		
		param_names[0] = "number of rows for header";
		param_values[0] = num_row;
		for(size_t i=1; i<num_row; ++i){
			getline(ifs,sdata);
			param_values[i] = getValueAndParamsFrom1Row(sdata, param_names[i]);
		}
		return true;
	}
	inline float BinFileFloat::read1data(){
		float data;
		ifs.read((char*)&data, sizeof data);
		return data;
	}

	inline double BinFileFloat::getValueFrom1Row(const std::string &sdata)const{
		return (float)atof(sdata.substr(sdata.find(" = ")+3).c_str());
	}
	inline double BinFileFloat::getValueAndParamsFrom1Row(const std::string& sdata, std::string& param_name)const{
		param_name = sdata.substr(0, sdata.find(" = ")).c_str();
		if(sdata.find(" = ") == std::string::npos) return sqrt(-1.0);
		return atof(sdata.substr(sdata.find(" = ")+3).c_str());
	}

	inline double BinFileFloat::getParams(const std::string para_name)const{
		for(size_t i=0; i<param_values.size();++i){
			if((param_names.size()>0 && param_names[i][0] != '#') &&
				param_names[i].find(para_name) != std::string::npos) 
				return param_values[i];
		}
		return sqrt(-1.0);
	}


};