class LHDAnalyzedData;

inline bool LHDAnalyzedData::downloadFile(const std::string igetfilePath, const size_t shotnum, const std::string diagname, const std::string outputDir){
	std::stringstream command;
	command << igetfilePath << " -s " << shotnum << " -d " << diagname << " -o " << outputDir;
	system(command.str().c_str());
	if (errno == -1) {
		std::string errormessage = "error in excuting " + command.str();
		throw errormessage;
	};
	return true;
}

inline bool LHDAnalyzedData::read(const std::string filename){

	std::ifstream ifs(filename.c_str(), std::ifstream::in);

	if (!ifs.is_open()){
		std::string errormessage = filename + " is not found";
		throw errormessage;
	}

	//---	reading	---
	std::string sdata;
	std::vector<unsigned short> data_1line_us(8);	//	this vector is used to convert from double to unsigned short
	while(!ifs.eof()){
		getline(ifs,sdata);

		//	reading the header	
		if( sdata.size() > 0 && sdata[0] == '#'){
			//	dimension
			if(sdata.find("DimNo") != std::string::npos){
				std::string dimnum = sdata.substr(sdata.find("=")+1);
				dimension = (size_t)floor(atof(dimnum.c_str())+0.5);
			}

			//	x, y axis
			else if(sdata.find("DimName") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				x_name = substr.substr(0, substr.find_first_of(","));
				y_name = substr.substr(substr.find_first_of(",")+1);
			}

			//	x, y number
			else if(sdata.find("DimSize") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				x_num = (size_t)floor(atof(substr.substr(0, substr.find_first_of(",")).c_str())+0.5);
				y_num = (size_t)floor(atof(substr.substr(substr.find_first_of(",")+1).c_str())+0.5);
				x_data.resize(0);
				y_data.resize(0);
			}

			//	x, y unit
			else if(sdata.find("DimUnit") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				x_unit = substr.substr(0, substr.find_first_of(","));
				y_unit = substr.substr(substr.find_first_of(",")+1);
			}

			//	x, y unit
			else if(sdata.find("DimUnit") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				x_unit = substr.substr(0, substr.find_first_of(","));
				y_unit = substr.substr(substr.find_first_of(",")+1);
			}

			//	number of variables
			else if(sdata.find("ValNo") != std::string::npos){
				val_num = (size_t) floor(atof(sdata.substr(sdata.find("=")+1).c_str()) + 0.5);
				data_table.resize(val_num);
				data_name.resize(val_num);
				data_unit.resize(val_num);
				data_table.resize(0);
			}

			//	variables
			else if(sdata.find("ValName") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				for(size_t i=0; i<val_num; ++i){
					data_name[i] = substr.substr(0, substr.find_first_of(","));
					substr = substr.substr(substr.find_first_of(",")+1);
				}
			}

			//	variables unit
			else if(sdata.find("ValUnit") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=")+1);
				for(size_t i=0; i<val_num; ++i){
					data_unit[i] = substr.substr(0, substr.find_first_of(","));
					substr = substr.substr(substr.find_first_of(",")+1);
				}
			}
		}

		else if(sdata.size() > 0 ){	//	reading the data
			std::string substr = sdata.substr(sdata.find("=")+1);
			data_table.push_back(std::vector<double> (val_num, 0.0));
			
			x_data.push_back(atof(substr.substr(0, substr.find_first_of(",")).c_str()));
			substr = substr.substr(substr.find_first_of(",")+1);
			if(dimension==2){
				y_data.push_back(atof(substr.substr(0, substr.find_first_of(",")).c_str()));
				substr = substr.substr(substr.find_first_of(",")+1);
			}

			for(size_t i=0; i<val_num; ++i){
				data_table[data_table.size()-1][i]= atof(substr.substr(0, substr.find_first_of(",")).c_str());
				substr = substr.substr(substr.find_first_of(",")+1);
			}
		}
	}
	ifs.close();

	//	set into table if the dimension is 2
	if(dimension ==2)
		return SetIntoMatrix();
	else if(dimension ==1)
		return SetIntoVector();
	else return false;
	return true;
}

inline bool LHDAnalyzedData::SetIntoMatrix(){
	data_table_vector.resize(val_num);
	for(size_t i=0; i<val_num; ++i){
		data_table_vector[i].resize(x_num, std::vector<double> (y_num));
		for(size_t j=0; j<x_num; ++j){
			for(size_t k=0; k<y_num; ++k){
				data_table_vector[i][j][k] = data_table[j*y_num+k][i];
			}
		}
	}
	x_data_1dim.resize(x_num);
	for(size_t i=0; i<x_num; ++i) x_data_1dim[i] = x_data[i*y_num];
	y_data_1dim.resize(y_num);
	for(size_t j=0; j<y_num; ++j) y_data_1dim[j] = y_data[j];

//	clear the older data
	data_table.clear();
	x_data.clear(); y_data.clear();
	return true;
}

inline bool LHDAnalyzedData::SetIntoVector(){
	if(data_table.size()==0) return false;
	std::vector<std::vector<double>> data_table_new(data_table[0].size(), std::vector<double>(data_table.size()));
	for(size_t i=0; i<data_table_new.size(); ++i){
		for(size_t j=0; j<data_table_new[i].size(); ++j){
			data_table_new[i][j]=data_table[j][i];
		}
	}
	data_table = data_table_new;
	return true;
}

inline int LHDAnalyzedData::SearchDiagname(const std::string diagname)const{
	for(size_t i=0; i<val_num; ++i){
		std::string diagname_tmp = data_name[i].substr(data_name[i].find_first_of('\'')+1);
		diagname_tmp = diagname_tmp.substr(0,diagname_tmp.find_first_of('\''));
		if(diagname == diagname_tmp)
			return i;
	}
	std::string errorMessage = diagname + "is not found in the LHDAnalyzedData";
	return -1;
}

inline const std::vector<double>& LHDAnalyzedData::getVector(const size_t index)const{
	return data_table[index];
}
inline const std::vector<double>& LHDAnalyzedData::getVector(const std::string diagname)const{
	int diag_number = SearchDiagname(diagname);
//	if(diag_number >=0 )
		return getVector(diag_number);
//	else return std::vector<std::vector<double>>(0, std::vector<double>(0));
}
inline const std::vector<std::vector<double>>& LHDAnalyzedData::getMatrix(const size_t index)const{
	return data_table_vector[index];
}
inline const std::vector<std::vector<double>>& LHDAnalyzedData::getMatrix(const std::string diagname)const{
	int diag_number = SearchDiagname(diagname);
//	if(diag_number >=0 )
		return getMatrix(diag_number);
//	else return std::vector<std::vector<double>>(0, std::vector<double>(0));
}

inline std::vector<double> LHDAnalyzedData::ReduceTo1dimension(const Dimension axis, const size_t index)const{
	if(axis == Yaxis){
		return data_table_vector[index][0];
	}
	else{
		std::vector<double> tmp(x_num);
		for(size_t i=0; i<x_num; ++i)
			tmp[i] = data_table_vector[index][i][0];
		return tmp;
	}
}

inline std::vector<double> LHDAnalyzedData::ReduceTo1dimension(const Dimension axis, const std::string diagname)const{
	int diag_number = SearchDiagname(diagname);
	if(diag_number >=0 )
		return ReduceTo1dimension(axis, diag_number);
	else return std::vector<double>(0);
}
inline double LHDAnalyzedData::ReduceTo1Value(const std::string diagname)const{
	int diag_number = SearchDiagname(diagname);
	if(diag_number >=0 )
		return ReduceTo1Value(diag_number);
	else return 0.0;
}

inline size_t LHDAnalyzedData::getNearestIndex(const double timing)const{
	if(x_data.size()==0 && x_data_1dim.size()==0 ) return 0;
	
	if(x_data.size()!=0){
		double tmp=fabs(x_data[0]-timing);
		size_t tmin=0;
		for(size_t t=0; t<x_data.size(); ++t){
			if(fabs(x_data[t] - timing) < tmp){
				tmin = t; tmp = fabs(x_data[t] - timing) ;
			}
		}
		return tmin;
	}else{
		double tmp=fabs(x_data_1dim[0]-timing);
		size_t tmin=0;
		for(size_t t=0; t<x_data_1dim.size(); ++t){
			if(fabs(x_data_1dim[t] - timing) < tmp){
				tmin = t; tmp = fabs(x_data_1dim[t] - timing) ;
			}
		}
		return tmin;
	}
}