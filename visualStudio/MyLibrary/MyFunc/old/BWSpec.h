#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <locale>

class BWSpec
{
public:
	std::vector<std::vector<double>> data;
	std::vector<std::string> data_name;
	double integ_time;
	std::vector<double> exp_time;
	BWSpec(void){};
	~BWSpec(void){};
	bool read(std::string filename);
	bool getv(std::string diag, std::vector<double>& data_v);
	static bool open_check(std::string filename);
private:
	static bool IsNumeric(std::string p);
	static std::vector<double> split(const std::string &s, const std::string d);
	static std::vector<std::string> ssplit(const std::string &s, const std::string d);
};

//--------------------------------------------------------------
//
//		static functions 
//
//--------------------------------------------------------------

inline bool BWSpec::IsNumeric(std::string p)
{
	double val=atof(p.c_str());
	if (val != 0.0) return true;
	else if(val == 0.0 && p[0] != '0') return false;
	else return true;
}

//---	BWSpec::split function	---
inline std::vector<double> BWSpec::split(const std::string &s, const std::string d) {
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
inline std::vector<std::string> BWSpec::ssplit(const std::string &s, const std::string d) {
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
inline bool BWSpec::open_check(std::string filename){
	std::ifstream ifs(filename.c_str(), std::ifstream::in );
	if(!ifs.is_open()) return false;
	std::string sdata;
	getline(ifs,sdata);
	if(sdata!="File Version;BWSpec 3.24u_57")return false;
	else return true;
}

#include "BWSpec.inl"
