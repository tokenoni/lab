inline bool BWSpec::read(std::string filename){
	data.resize(0);
	std::ifstream ifs;
	//---	for treating the file pass in japanese	---
	std::locale default_locale;
	std::locale::global(std::locale("japanese"));
	ifs.open(filename.c_str(), std::ifstream::in );
	std::locale::global(default_locale);

	if(!ifs.is_open()){	return false;	}

	std::string sdata;
	getline(ifs,sdata);
	if(sdata!="File Version;BWSpec 3.24u_57"){	return false;	}

	//---	integration time	---
	while(sdata.find("intigration") == std::string::npos)	getline(ifs,sdata);
	integ_time = 0.001*atof(sdata.substr(18).c_str());

	//---	reading the diag name	---
	while(sdata[0]!='P') getline(ifs, sdata);
	data_name = ssplit(sdata, ";");
	
	while(sdata[0]!='0') getline(ifs,sdata);
		
	while(1){
		std::vector<double>row_data = split(sdata, ";\n");
		if(row_data.size()>0) data.push_back(row_data);
		getline(ifs,sdata);
		if(ifs.eof()) break;
	}
	ifs.close();
	return true;
}

//	if the BWSpec data contains the name of "diag", return true;
//	the data is stored in data_v
inline bool BWSpec::getv(std::string diag, std::vector<double>& data_v){
	size_t index=data_name.size();
	for(size_t i=0; i<data_name.size(); i++){
		if(data_name[i].find(diag)!= std::string::npos){
			index = i;
			break;
		}
	}
	if(index==data_name.size()) return false;
	
	data_v.resize(data.size());
	for(size_t i=0; i<data.size(); ++i) data_v[i] = data[i][index];
	for(size_t i=0; i<data.size(); ++i)
		if(data_v[i] !=0) return true;
	return false;	//	if the data contain only "0", it is not data
}

