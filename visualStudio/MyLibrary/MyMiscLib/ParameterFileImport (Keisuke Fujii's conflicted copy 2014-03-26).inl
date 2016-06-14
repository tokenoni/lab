namespace mymisclib{

inline ParameterFileImport::ParameterFileImport(const ParameterFileImport& obj){
	filename = obj.filename;
	parameter_names = obj.parameter_names;
	parameter_values = obj.parameter_values;
	comment = obj.comment;
}

inline ParameterFileImport& ParameterFileImport::operator= (const ParameterFileImport&obj){
	filename = obj.filename;
	parameter_names = obj.parameter_names;
	parameter_values = obj.parameter_values;
	comment = obj.comment;
	return *this;
}

inline bool ParameterFileImport::readFile(const std::string filename_){
	filename = filename_;
	std::fstream ifs;
	ifs.open(filename, std::ifstream::in);
	if(!ifs.good()) return false;

	//--- clear the memory ---
	comment.clear();
	parameter_names.clear();
	parameter_values.clear();
	//------------------------

	std::string sdata;
	do{
		getline(ifs,sdata);
		if(sdata[0] == '#')
			comment.append(sdata);
		else
			addNameValue(sdata);
	}while(!ifs.eof());

	return true;
}

inline void ParameterFileImport::addNameValue(const std::string sdata){
	if(sdata.find("=") == std::string::npos)
		return;
	
	size_t i = parameter_names.size();
	parameter_names.push_back(sdata.substr(0, sdata.find("=")));
	//	removing spaces
	parameter_names[i].erase(std::remove(parameter_names[i].begin(), parameter_names[i].end(), ' '), parameter_names[i].end());

	parameter_values.push_back(sdata.substr(sdata.find("=")+1));
	//	removing spaces
	parameter_values[i].erase(std::remove(parameter_values[i].begin(), parameter_values[i].end(), ' '), parameter_values[i].end());
}

inline size_t ParameterFileImport::searchIndex(const std::string parameter_name)const
{
	size_t i=0;
	for(i=0; i<parameter_names.size();++i){
		if(parameter_names[i].find(parameter_name) != std::string::npos)
			return i;
	}
	if(i==parameter_names.size()){
		std::cerr << "could not find the parameter name of " << parameter_name << "in the parameter file" << std::endl;
	}
	return i;
}

inline bool ParameterFileImport::isParameterExist(const std::string paraname)const{
	size_t i=0;
	for(i=0; i<parameter_names.size();++i){
		if(parameter_names[i].find(paraname) != std::string::npos)
			return true;
	}
	if(i<parameter_names.size())	return true;
	else return false;
}

inline std::string ParameterFileImport::getValue_string(const std::string paraname)const{
	return parameter_values[searchIndex(paraname)];
}
inline double ParameterFileImport::getValue_double(const std::string paraname)const{
	return atof(parameter_values[searchIndex(paraname)].c_str());
}
inline int ParameterFileImport::getValue_int(const std::string paraname)const{
	return atoi(parameter_values[searchIndex(paraname)].c_str());
}
inline size_t ParameterFileImport::getValue_size_t(const std::string paraname)const{
	return abs(atoi(parameter_values[searchIndex(paraname)].c_str()));
}
inline bool ParameterFileImport::getValue_bool(const std::string paraname)const{
	if(parameter_values[searchIndex(paraname)].find("yes") != std::string::npos)
		return true;
	if(parameter_values[searchIndex(paraname)].find("true") != std::string::npos)
		return true;
	else return false;
}

inline std::vector<std::string> ParameterFileImport::getValueList_string(const std::string paraname)const{
	std::string sdata = parameter_values[searchIndex(paraname)];
	std::vector<std::string> rslt(0);
	while(sdata.find(",")!= std::string::npos){
		std::string a_string = sdata.substr(0, sdata.find(","));
		
		if(a_string.size()!=0)
			rslt.push_back(a_string);
		sdata = sdata.substr(sdata.find(",")+1);
	}
	if(sdata.size()!=0)
		rslt.push_back(sdata);
	return rslt;
}

inline std::vector<double> ParameterFileImport::getValueList_double(const std::string paraname)const{
	std::vector<std::string> rslt_s = getValueList_string(paraname);
	std::vector<double> rslt(rslt_s.size());
	for(size_t i=0; i<rslt.size();++i)
		rslt[i] = atof(rslt_s[i].c_str());
	return rslt;
}
inline std::vector<int> ParameterFileImport::getValueList_int(const std::string paraname)const{
	std::vector<std::string> rslt_s = getValueList_string(paraname);
	std::vector<int> rslt(rslt_s.size());
	for(size_t i=0; i<rslt.size();++i)
		rslt[i] = atoi(rslt_s[i].c_str());
	return rslt;
}
inline std::vector<size_t> ParameterFileImport::getValueList_size_t(const std::string paraname)const{
	std::vector<std::string> rslt_s = getValueList_string(paraname);
	std::vector<size_t> rslt(rslt_s.size());
	for(size_t i=0; i<rslt.size();++i)
		rslt[i] = abs(atoi(rslt_s[i].c_str()));
	return rslt;
}

};