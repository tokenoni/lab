namespace mygsl{

	inline generalBinfile::generalBinfile(datatype datatype_){dataType = datatype_;}

	inline generalBinfile::~generalBinfile(void){}

	inline void generalBinfile::addParams(const std::string param_name, const double param_value){
		param_names.push_back(param_name);
		param_values.push_back(param_value);
	}

	inline void generalBinfile::clearParams(){
		param_names.clear();
		param_values.clear();
	}

	inline bool generalBinfile::Parafile_write(const std::string para_file_){
		ofs_para.open(para_file_.c_str(), std::ofstream::out);
		if(ofs_para.good()){
			//---	first row is the number of the parameters (number of rows until the data started)	---
			ofs_para<< "number_of_parameters = " << param_values.size() << std::endl;

			for(size_t i=0;i<param_values.size(); ++i)
				ofs_para << param_names[i].c_str() << " = " << param_values[i] << std::endl;

			//---	data types	---
			ofs_para<< "data_type = " << getDataType_string(dataType) << std::endl;

			ofs_para.close();
			return true;
		}
		else return false;
	}

	inline bool generalBinfile::Datafile_write_open(const std::string data_file_){
		ofs_data.open(data_file_.c_str(), std::ofstream::binary, std::ofstream::app);
		if(ofs_data.good()){
			return true;
		}
		else return false;
	}
	inline bool generalBinfile::Datafile_write_add( void* const data){
		if(ofs_data.good()){
		//	append to file in appropriate type
			switch(dataType){
			case intiger32bit:
				ofs_data << (int*) data;
				break;
			case unsigned_integer32bit:
				ofs_data << (unsigned int*) data;
				break;
			case short_intiger16bit:
				ofs_data << (short*) data;
				break;
			case unsigned_short_integer16bit:
				ofs_data << (unsigned short*) data;
				break;
			case double64bit:
				ofs_data << (double*) data;
				break;
			case float32bit:
				ofs_data << (float*) data;
				break;
			};
			return true;
		}else return false;	
		
	}
	
	inline bool generalBinfile::Datafile_write_close(){
		if(ofs_data.good()){
			ofs_data.close();
			return true;
		}else return false;	
	}

	inline std::string generalBinfile::getDataType_string(const datatype datatype_)const{
		switch(dataType){
			case intiger32bit:
				return "integer32bit\n";
				break;
			case unsigned_integer32bit:
				return "unsigned_integer32bit\n";
				break;
			case short_intiger16bit:
				return "short_intiger16bit\n";
				break;
			case unsigned_short_integer16bit:
				return "short_intiger16bit\n";
				break;
			case double64bit:
				return "double64bit";
				break;
			case float32bit:
				return "float32bit";
			};
	}
};