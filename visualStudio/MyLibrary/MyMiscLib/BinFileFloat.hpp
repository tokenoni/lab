#ifndef __MYGSL_BINFILEFLOAT_HPP__
#define __MYGSL_BINFILEFLOAT_HPP__

#include <vector>
#include <string>
#include <fstream>

namespace mymisclib{
class BinFileFloat
{
public:
	BinFileFloat(void);
	~BinFileFloat(void);

	void clear(){param_names.clear(); param_values.clear();}
	void add(const std::string, const double);

	//	for writing
	bool open_for_write(const std::string filename);
	bool add_bindata(float data);
	void finish_writing();

	//	for reading
	bool open_for_read(const std::string filename);
	float read1data();
	bool eof()const{return ifs.eof();}
	void finish_reading(){ifs.close();};

	double getParams(const std::string para_name)const;
protected:
	std::ofstream ofs;
	std::ifstream ifs;
	std::vector< std::string > param_names;
	std::vector< double > param_values;

	double getValueFrom1Row(const std::string& sdata)const;
	double getValueAndParamsFrom1Row(const std::string& sdata, std::string& param_name)const;
};

};
#include "BinFileFloat.inl"

#endif