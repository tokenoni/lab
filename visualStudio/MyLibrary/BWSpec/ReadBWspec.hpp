#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace BWSpec{
	static bool IsNumeric(std::string p)
	{
		double val=atof(p.c_str());
		if (val != 0.0) return true;
		else if(val == 0.0 && p[0] != '0') return false;
		else return true;
	}
//---	BWSpec::split function	---
	static std::vector<double> split(const std::string &s, const std::string d) {
		std::vector<double> elements;						
		size_t i_start=0;							
		size_t i_end = s.find_first_of(d.c_str(), i_start);	
		while(i_end!=std::string::npos){					
			elements.push_back(atof((s.substr(i_start, i_end - i_start +1)).c_str()));
			i_start=i_end+1;
			i_end = s.find_first_of(d.c_str(), i_start);	
		}
		elements.push_back(atof((s.substr(i_start, i_end - i_start +1)).c_str()));
		return elements;
	}
//---	split string and return std::vector<std::string>	---
	static std::vector<std::string> ssplit(const std::string &s, const std::string d) {
		std::vector<std::string> elements;						
		size_t i_start=0;							
		size_t i_end = s.find_first_of(d.c_str(), i_start);	
		while(i_end!=std::string::npos){					
			elements.push_back(s.substr(i_start, i_end - i_start +1));
			i_start=i_end+1;
			i_end = s.find_first_of(d.c_str(), i_start);	
		}
		return elements;
	}

//---	check the file existance	---
	static bool open_check(std::string filename){
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			return false;
		}else return true;
	}


//----	read an BWSpec file and put values into std::vector<std::vector<double>>	---
	static std::vector<std::vector<double>> read(std::string filename){
		std::vector<std::vector<double>> all_data;
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			getchar();exit(1);
		}

		std::string sdata;
		getline(ifs,sdata);
		if(sdata!="File Version;BWSpec 3.24u_57"){
			std::cerr << "not a file for BWSpec 3.24u_57"<< std::endl;
			getchar();exit(1);
		}
		
		while(sdata.find("Wavelength")==std::string::npos)
			getline(ifs,sdata);
			
		getline(ifs,sdata);
		while(1){
			std::vector<double>row_data = split(sdata, ";\n");
			if(row_data.size()>0) all_data.push_back(row_data);
			getline(ifs,sdata);
			if(ifs.eof()) break;
		}
		ifs.close();
		return all_data;
	}

//----	read an BWSpec file and put values into std::vector<std::vector<double>>	---
	static std::vector<std::vector<double>> read(std::string filename, const std::vector< std::string> data_type_v){
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			getchar();exit(1);
		}

		std::string sdata;
		getline(ifs,sdata);
		if(sdata!="File Version;BWSpec 3.24u_57"){
			std::cerr << "not a file for BWSpec 3.24u_57"<< std::endl;
			getchar();exit(1);
		}
		
		while(sdata.find("Wavelength")==std::string::npos)
			getline(ifs,sdata);

		std::vector<std::string> data_type_src = ssplit(sdata, ";");
		std::vector<size_t> num_data_v = std::vector<size_t> (data_type_v.size(), sqrt(-1.0));

		getline(ifs,sdata);

		for(size_t i=0; i<num_data_v.size(); ++i){
			for(size_t j=0; j<num_data_v.size(); ++j){
				if(data_type_src[j].find(data_type_v[i]) != std::string::npos)
					num_data_v[j] = i;
			}
		}
			
		getline(ifs,sdata);
		std::vector<std::vector<double>> all_data = std::vector< std::vector<double>>(0);
		while(1){
			std::vector<double>row_data = split(sdata, ";\n");
			std::vector<double>row_data_dest(num_data_v.size(),0.0);
			for(size_t i=0; i<num_data_v.size(); ++i){
				if(num_data_v[i] == num_data_v[i] )
					 row_data_dest[i] = row_data[num_data_v[i]];
				else row_data_dest[i] = sqrt(-1.0);
			}
			all_data.push_back(row_data_dest);
			getline(ifs,sdata);
			if(ifs.eof()) break;
		}
		ifs.close();
		return all_data;
	}
	
//---	read data type	---
	std::vector<std::string> read_type(std::string filename){
		std::vector<std::string> data;
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			getchar();exit(1);
		}

		std::string sdata;
		getline(ifs,sdata);
		if(sdata!="File Version;BWSpec 3.24u_57"){
			std::cerr << "not a file for BWSpec 3.24u_57"<< std::endl;
			getchar();exit(1);
		}

		while(sdata[0]!='P')	getline(ifs,sdata);
		if(ifs.eof()){
			std::cerr << "!!!   file: " + filename << " doesn't contain intigration time data" << std::endl;
			getchar();exit(1);
		}else return ssplit(sdata, ";");
	}

//---	read integration time	---
	double read_integtime(std::string filename){
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			getchar();exit(1);
		}

		std::string sdata;
		getline(ifs,sdata);
		if(sdata!="File Version;BWSpec 3.24u_57"){
			std::cerr << "not a file for BWSpec 3.24u_57"<< std::endl;
			getchar();exit(1);
		}

		while(sdata.find("intigration") == std::string::npos)	getline(ifs,sdata);
		return atof(sdata.substr(18).c_str());
	}

//---	read experimental time	---
	std::vector<double> read_exptime(std::string filename){
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			getchar();exit(1);
		}

		std::string sdata;
		getline(ifs,sdata);
		if(sdata!="File Version;BWSpec 3.24u_57"){
			std::cerr << "not a file for BWSpec 3.24u_57"<< std::endl;
			getchar();exit(1);
		}
	//Date;2010-07-09 10:43:10

		while(sdata.find("Date") == std::string::npos)	getline(ifs,sdata);
		std::vector<double> date = split(sdata, ";- :\n");
		date.erase(date.begin());
		return date;
	}
		
};