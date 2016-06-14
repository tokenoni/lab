#pragma once
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

static const std::string RetrievePath = "C:/LABCOM/Retrieve/bin/Retrieve.exe";
static const std::string Retrieve_tPath = "C:/LABCOM/Retrieve/bin/Retrieve_t.exe";

template< class Ty_ >
class RawData
{
public:
//	size_t data_num;
	std::vector< Ty_ > data;
	std::vector<double> timing;
	bool Exist;
	bool Exist_t;

	enum BitOrder{
		BigEndianness,
		LittleEndianess
	};

	bool ReadBinary(std::string filename, BitOrder bit_order, size_t bit_size);
	bool ReadBinary(BitOrder bit_order, size_t bit_size);
	bool ReadTiming(std::string filename);
	bool ReadTiming();
	bool DeleteBinary(std::string filename){
		if(remove(filename.c_str())==0) return true;
		else return false;};
	bool DeleteBinary(){
		if((remove(RetrievedDataPath.c_str()) ==0 && remove(RetrievedParaPath.c_str())==0)){
			remove(Retrieved_tDataPath.c_str());
			remove(Retrieved_tParaPath.c_str());
			return true;
		}else return false;	}
	bool Retrieve(std::string Expname, size_t shot_num, size_t subshot_num, size_t ch_num);
	bool Retrieve_t(std::string Expname, size_t shot_num, size_t subshot_num, size_t ch_num);


	RawData(void){Exist = false; Exist_t = false;};
	~RawData(void){};
	void Free(){data.clear(); timing.clear();}

	Ty_& operator [](size_t i){return data[i];}
	size_t size(){return data.size();}
	std::string RetrievedDataPath;
	std::string RetrievedParaPath;
	std::string Retrieved_tDataPath;
	std::string Retrieved_tParaPath;
private:
};


template< class Ty_ >
inline bool RawData< Ty_ >::Retrieve(std::string Expname, size_t shot_num, size_t subshot_num, size_t ch_num){
	std::stringstream command;
	command << "cmd /c " << RetrievePath << " " << Expname << " " <<shot_num << " " << subshot_num << " " << ch_num << " ";
	if(system(command.str().c_str())) Exist = true; else Exist = false;
	
	std::stringstream path;
	path << Expname <<"-"<< shot_num <<"-"<<subshot_num<<"-"<<ch_num;
	RetrievedDataPath = path.str()+".dat";
	RetrievedParaPath = path.str()+".prm";
	return true;
}
template< class Ty_ >
inline bool RawData< Ty_ >::Retrieve_t(std::string Expname, size_t shot_num, size_t subshot_num, size_t ch_num){
	std::stringstream command;
	command << "cmd /c " << Retrieve_tPath << " " << Expname << " " <<shot_num << " " << subshot_num << " " << ch_num << " ";
	if(system(command.str().c_str())) Exist_t = true; else Exist_t=false;
	
	std::stringstream path;
	path << Expname <<"-"<< shot_num <<"-"<<subshot_num<<"-"<<ch_num;
	Retrieved_tDataPath = path.str()+".time";
	Retrieved_tParaPath = path.str()+".tprm";
	return true;
}

template< class Ty_ >
inline bool RawData< Ty_ >::ReadBinary(BitOrder bit_order, size_t bit_size){
	return ReadBinary( RetrievedDataPath, bit_order, bit_size);
}

template< class Ty_ >
inline bool RawData< Ty_ >::ReadBinary(std::string filename, BitOrder bit_order, size_t bit_size)
{
	std::ifstream fin;
	fin.open( filename.c_str(), std::ios::in | std::ios::binary);
	if(fin.fail()){
		std::cerr << "file " << filename <<" could not open!" << std::endl;
		Exist=false; return false;
	}
	
	char c;
	int one_byte_data[8];	//	max byte size is 64 bit

	Ty_ one_data;
	data.resize(0);
	
//--	variables for memory allocation	--
	size_t data_num=0, data_num_tmp = 0;
	size_t data_num_increment = 10000;
//----------------------------------------
	while(1){
		one_data=0;
		for(size_t i=0; i<bit_size/8; i++){
			fin.read((char*)&c, 1);
			one_byte_data[i] = c;
		}
		if(fin.eof()) break;
		//	expand the vector data number
		if(data_num >= data_num_tmp){
			data_num_tmp += data_num_increment;
			data.resize(data_num_tmp);
		}
		//---	LittleEndianness	---
		//---	BigEndianness	---
		if(bit_order == BigEndianness){
			for(size_t i=0; i<bit_size/8-1; i++){
				one_data = (one_data | one_byte_data[i]);
				one_data = one_data << 8;
			}
			one_data = (one_data | (one_byte_data[bit_size/8-1] & 255));
		}
		//---	LittleEndianness	---
		else{
			for(int i=bit_size/8-1; i>0; i--){
				one_data = (one_data | (one_byte_data[i] & 255));
				one_data = one_data << 8;
			}
			one_data = (one_data | (one_byte_data[0] & 255));
		}
		data[data_num] = one_data;
		++data_num;
//		data.push_back(one_data);
	}
	fin.close();
	data.resize(data_num);
	return true;
}

template< class Ty_ >
inline bool RawData< Ty_ >::ReadTiming(){
	return ReadTiming(Retrieved_tParaPath);
}

template< class Ty_ >
inline bool RawData< Ty_ >::ReadTiming(std::string filename){
	std::vector<std::string> param_list(2);
	param_list[0] =	"ClockCycle,";
	param_list[1] =	"StartTiming,";

	std::fstream ifs(filename.c_str(), std::ios::in);
	if(ifs.fail()){
		std::cerr << "file " << filename <<" could not open!" << std::endl;
		return false;
	}
	
	//	read parameter
	std::vector<std::string> sdata_list(0);
	std::string sdata;
	while(!ifs.eof()){
		getline(ifs,sdata);
		sdata_list.push_back(sdata);
	}

	//	sort parameter
	std::vector<std::string> param_value(param_list.size());
	for(size_t i=0; i<param_list.size();i++){
		for(size_t j=0; j<sdata_list.size(); j++){
			if(sdata_list[j].find(param_list[i]) != std::string::npos){
				size_t i_colon =sdata_list[j].find(param_list[i])+param_list[i].length();
//				size_t i_num =sdata_list[j].find("\n", i_colon);
				param_value[i] = sdata_list[j].substr(i_colon);
				break;
			}
		}

	}
	double rate  =	1.0/atof(param_value[0].c_str());
	double start_t =   atof(param_value[1].c_str());
	timing.resize(size());
	for(size_t i=0; i<timing.size(); ++i) timing[i] = start_t + 1.0/rate*i;
	ifs.close();

	return true;
}
