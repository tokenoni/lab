#ifndef __PARAMETER_FILE_IMPORT_HPP__
#define __PARAMETER_FILE_IMPORT_HPP__

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm> 

namespace mymisclib{
class ParameterFileImport{
public:
	ParameterFileImport(const std::string filename_){readFile(filename_);};
	ParameterFileImport(void){};
	ParameterFileImport(const ParameterFileImport& obj);
	ParameterFileImport& operator= (const ParameterFileImport&obj);

	~ParameterFileImport(void){};
	bool readFile(const std::string filename_);
	std::string getValue_string(const std::string paraname)const;
	double getValue_double(const std::string paraname)const;
	int getValue_int(const std::string paraname)const;
	size_t getValue_size_t(const std::string paraname)const;
	bool getValue_bool(const std::string paraname)const;
	bool isParameterExist(const std::string paraname)const;
	
	std::vector<double> getValueList_double(const std::string paraname)const;
	std::vector<int> getValueList_int(const std::string paraname)const;
	std::vector<size_t> getValueList_size_t(const std::string paraname)const;
	std::vector<bool> getValueList_bool(const std::string paraname)const;
	std::vector<std::string> getValueList_string(const std::string paraname)const;

private:
	std::string filename;
	std::vector<std::string> parameter_names, parameter_values;
	size_t searchIndex(const std::string parameter_name)const;

	void addNameValue(const std::string sdata);
	std::string comment;
};
};
#include "ParameterFileImport.inl"

#endif
