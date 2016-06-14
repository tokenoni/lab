#pragma once
#ifndef _IGORITX_HPP_
#define _IGORITX_HPP_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <limits>
#include <cstdlib>
#include <ios>


namespace IGORdata{
	static enum LinefeedCode{ CR, LF, CRLF};
	static inline std::istream& getline(std::istream& is, std::string&str, LinefeedCode lc){
		char charCR = 0x0D; char charLF = 0x0A;
		char onechar;
		switch (lc){
		case CR:
			getline(is, str, charCR);
			break;
		case LF:
			getline(is, str, charLF);
			break;
		case CRLF:
			getline(is, str, charCR);
			is.get(onechar);
			break;
		};
		return is;
	}
	static inline LinefeedCode checkLinefeedCode(std::istream& ifs){
		char charCR = 0x0D; char charLF = 0x0A;
		std::streamoff length = ifs.tellg();
		LinefeedCode lc;
		std::string sdata; char onechar;
		getline(ifs,sdata, charCR);
		if(ifs.eof()) 
			lc = LF;
		else{
			ifs.get(onechar);
			if(onechar == charLF)
				lc = CRLF;
			else
				lc = CR;
		}	
		ifs.seekg(length);//	get back to the original location of the file
		return lc;
	}
//---	check whether numeric or string ---
	static inline bool IsNumeric(std::string p)
	{
		double val=atof(p.c_str());
		if (val != 0.0) return true;
		else if(val == 0.0 && p[p.find_first_not_of(" \t\n")]== '0') return true;
		else if(val == 0.0 && p[p.find_first_not_of(" \t\n")]== '-') return true;
		else if(val == 0.0 && p[p.find_first_not_of(" \t\n")]== 'i') return true;
		else if(val == 0.0 && p[p.find_first_not_of(" \t\n")]== 'N') return true;
		else return false;
	}
//---	IGORdata::split function	---
	static inline std::vector<double> split(const std::string &s, const std::string d) 
	{
		std::vector<double> elements;							//	IGORdata::split std::string s and convert to double, 	
		size_t i_start=0, i_end=0;								//	and set to vector element
		while(i_end!=std::string::npos && i_start < s.size()){						//
			i_end = s.find_first_of(d.c_str(), i_start);		//
			if(i_end > i_start){								//
				if(IGORdata::IsNumeric(s.substr(i_start))){		//
					elements.push_back(atof((s.substr(i_start)).c_str()));
					i_start=i_end+1;
				}
			}else{	i_start = i_end+1;}
		}
		return elements;
	}

	enum OpenCondition{
		WaitFileOpenUntilReturnPressed,
		Ignore,
		TryJapanese
	};

//---	wait for the file opening	---
	static bool open(std::ifstream &ifs, const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		ifs.open(filename.c_str(), std::ios::in);
		if(open_condition == WaitFileOpenUntilReturnPressed){
			while(!ifs.is_open()){
				std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
				std::cerr << "!!!   ---  press return for retry ---		 !!!" << std::endl;
				getchar();
				ifs.close();ifs.clear();
				ifs.open(filename.c_str(), std::ifstream::in);
			}
			return true;
		}
		else if(open_condition == Ignore){
			return true;
		}
		else if(open_condition == TryJapanese){
			std::locale default_locale;
			std::locale::global(std::locale("japanese"));
			ifs.open(filename.c_str(), std::ifstream::in );
			std::locale::global(default_locale);
			return ifs.is_open();
		}
		return false;
	}
	static bool open(std::ofstream &ofs, const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		ofs.open(filename.c_str(), std::ios::out);
		if(open_condition == WaitFileOpenUntilReturnPressed){
			while(!ofs.is_open()){
				std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
				std::cerr << "!!!   ---  press return for retry ---		 !!!" << std::endl;
				getchar();
				ofs.close();ofs.clear();
				ofs.open(filename.c_str());
			}
			return true;
		}
		else if(open_condition == Ignore){
			return true;
		}
		else if(open_condition == TryJapanese){
			std::locale default_locale;
			std::locale::global(std::locale("japanese"));
			ofs.open(filename.c_str());
			std::locale::global(default_locale);
			return ofs.is_open();
		}
		return false;
	}

//---	check the file existance	---
	static bool open_check(const std::string filename, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			return false;
		}else return true;
	}


//---	operation for infinity or NAN data	---
	inline std::string WriteData(double data, 
			const int precision, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::stringstream ss;
		ss<<std::setprecision(precision);
/*		if(data.IsPositiveInfinity(data))		return "inf";
		else if(data.IsNegativeInfinity(data))	return "-inf";
		else if(data.IsNaN(data))				return "NAN";
*/		if(data != data) return "NAN";
		else{
			ss << data;
			return ss.str();
		}
	}
	
//----	read an itx file and put values into std::vector<std::vector<double>>	---
	static std::vector<std::vector<std::vector<double>>> itx_read(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> all_data;
		std::ifstream ifs;
		ifs.open( filename, std::ios_base::binary | std::ios_base::in | open_condition);
		LinefeedCode lc = checkLinefeedCode(ifs);

		std::string sdata;
		getline(ifs,sdata,lc);
		if(sdata.find("IGOR")==std::string::npos){
			std::cerr << "missing IGOR keyword"<< std::endl;
			getchar();exit(1);
		}
		
		while(!ifs.eof()){
			getline(ifs,sdata,lc);
			while(sdata.find("BEGIN") != std::string::npos){
				std::vector<std::vector<double>> block_data;
				getline(ifs,sdata,lc);
				while(sdata.find("END")==std::string::npos){
					std::vector<double>row_data = IGORdata::split(sdata, ",\t \n\r");
					if(row_data.size()>0) block_data.push_back(row_data);
					getline(ifs, sdata,lc);
					if(ifs.eof()) break;
				}
				for(size_t i=0; i<block_data.size(); i++){
					if(block_data[0].size() != block_data[i].size()){
						std::cerr << "wrong format data (each row must have same dimensions)"<< std::endl;
						getchar();exit(1);
					}
				}
				all_data.push_back(block_data);
			}
		}
		ifs.close();
		return all_data;
	}
	//	if the file is known to contain only 1 matrix
	static std::vector<std::vector<double>> itx_1stdvv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		if(datas.size()!=1){
			std::cerr << "function in IGORdata::itx_stdvv\n the file "<<filename<<" contains more than 1 matrix"<< std::endl;
			getchar();exit(1);
		}
		return datas[0];
	}
	//	if the file is known to contain only 1 matrix
	static std::vector < std::vector < std::vector< double >>> itx_1stdvvv(const std::string filename,
		OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> all_data;
		std::ifstream ifs;
		ifs.open(filename, std::ios_base::binary | std::ios_base::in | open_condition);
		LinefeedCode lc = checkLinefeedCode(ifs);

		std::string sdata;
		getline(ifs, sdata, lc);
		if (sdata.find("IGOR") == std::string::npos){
			std::cerr << "missing IGOR keyword" << std::endl;
			getchar(); exit(1);
		}
		do{
			getline(ifs, sdata, lc);
		} while (sdata.find("BEGIN") == std::string::npos);

		std::vector<std::vector<double>> block_data;
		while (!ifs.eof()){
			getline(ifs, sdata, lc);
			if (sdata.find("END") != std::string::npos) break;
			std::vector<double>row_data = IGORdata::split(sdata, ",\t \n\r");
			if (row_data.size() != 0){
				block_data.push_back(row_data);
			}
			else{
				all_data.push_back(block_data);
				block_data.clear();
			}
		}
		size_t size3 = all_data.size();
		size_t size1 = all_data[0].size();
		size_t size2 = all_data[0][0].size();
		std::vector<std::vector<std::vector<double>>> data(size1, std::vector < std::vector<double>> (size2, std::vector<double>(size3)));
		for (size_t i = 0; i < size1; ++i){
			for (size_t j = 0; j < size2; ++j){
				for (size_t k = 0; k < size3; ++k){
					data[i][j][k] = all_data[k][i][j];
				}
			}
		}
		ifs.close();
		return data;
	}
	//	if the file is known to contain only 2 vectors
	static std::vector<std::vector<double>> itx_multiple_stdv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		if(datas.size()!=1){
			std::cerr << "function in IGORdata::itx_stdvv\n the file contains more than 1 matrix"<< std::endl;
			getchar();exit(1);
		}
		std::vector<std::vector<double>> multi_v(datas[0][0].size(), std::vector<double> (datas[0].size(),0.0));
		for(size_t i=0; i<multi_v.size(); ++i){
			for(size_t j=0; j<multi_v[i].size(); ++j){
				multi_v[i][j] = datas[0][j][i];
			}
		}
		return multi_v;
	}
	//	if the file is known to contain only 1 vector
	static std::vector<double> itx_1stdv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		if(datas.size()<=2){
			std::vector<double> stdv1(datas[0].size());
			for(size_t i=0; i<datas[0].size(); ++i){
				if(datas[0][i].size() !=1){
					std::cerr << "function in IGORdata::itx_1stdv\n the file "<<filename<<" contains more than 1 vector"<< std::endl;
					getchar();exit(1);
				}
				stdv1[i] = datas[0][i][0];
			}
			return stdv1;
		}
		else
			return std::vector<double>(0);
	}

	//	if the file is known to contain only 1 string vector
	static std::vector<std::string> itx_1stdv_string(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ifstream ifs;
		ifs.open( filename, std::ios_base::binary | std::ios_base::in | open_condition);
		LinefeedCode lc = checkLinefeedCode(ifs);

		std::string sdata;
		getline(ifs,sdata,lc);
		if(sdata.find("IGOR")==std::string::npos){
			std::cerr << "missing IGOR keyword"<< std::endl;
			getchar();exit(1);
		}
		
		std::vector<std::string> data;
		while(!ifs.eof()){
			getline(ifs,sdata,lc);
			while(sdata.find("BEGIN") != std::string::npos){
				getline(ifs,sdata,lc);
				while(sdata.find("END")==std::string::npos){
					std::string sdata2 = sdata.substr(sdata.find("\"")+1);
					data.push_back(sdata2.substr(0, sdata2.find("\"")));
					getline(ifs, sdata,lc);
					if(ifs.eof()) break;
				}
			}
		}
		ifs.close();
		return data;
	}

/*		std::vector<std::vector<double>> stdv = itx_1stdvv(filename, open_condition);
		if(stdv.size()<=1){
			return stdv[0];
		}else{
			std::vector<double> stdv1(stdv.size());
			for(size_t i=0; i<stdv.size(); ++i){
				if(stdv[i].size() !=1){
					std::cerr << "function in IGORdata::itx_1stdv\n the file "<<filename<<" contains more than 1 vector"<< std::endl;
					getchar();exit(1);
				}
				stdv1[i] = stdv[i][0];
			}
			return stdv1;
		}
	
	}
	*/


//----	writ to itx file from std::vector<size_t>  ----
	static bool write_itx(const std::vector<size_t> data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(data[i], precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<double>  ----
	template <class Ty_ >
	static bool write_itx(const std::vector< Ty_ >& data, const std::string filename,const  std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(data[i], precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<double>  ----
	template <class Ty_ >
	static bool write_itx_log(const std::vector< Ty_ >& data, const std::string filename,const  std::string wavename, 
			const size_t precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(log(data[i]), precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<std::string>  ----
	static bool write_itx(const std::vector<std::string>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/T/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t\""<<data[i]<<"\""<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}

//----	writ to itx file from std::vector<std::vector<>>  ----
	template <class Ty_ >
	static bool write_itx(const std::vector<std::vector< Ty_ >>& data, const std::string filename, const std::string wavename, 
				const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		for(size_t i=0; i<data.size(); i++){
			if(data[0].size()!=data[i].size()){
				std::cerr << "!!!   data should be same size   !!!" << std::endl;
				getchar();getchar();exit(1);
			}
		}
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<", "<<data[0].size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++){
			for(size_t j=0; j<data[0].size(); j++){
				ofs <<"\t"<<WriteData(data[i][j], precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<std::vector<std::vector<>>>  ----
	template <class Ty_ >
	static bool write_itx(const std::vector<std::vector<std::vector< Ty_ >>> & data, const std::string filename, const std::string wavename, 
				const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		for(size_t i=0; i<data.size(); i++){
			if(data[0].size()!=data[i].size()){
				std::cerr << "!!!   data should be same size   !!!" << std::endl;
				getchar();getchar();exit(1);
			}
		}
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<", "<<data[0].size()<<", "<<data[0][0].size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t k=0; k<data[0][0].size(); k++){
			for(size_t i=0; i<data.size(); i++){
				for(size_t j=0; j<data[i].size(); j++){
					ofs <<"\t"<<WriteData(data[i][j][k], precision);
				}
				ofs << std::endl;
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<std::vector<>>  ----
	static bool write_itx_log(const std::vector<std::vector< double >>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		for(size_t i=0; i<data.size(); i++){
			if(data[0].size()!=data[i].size()){
				std::cerr << "!!!   data should be same size   !!!" << std::endl;
				getchar();getchar();exit(1);
			}
		}
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<", "<<data[0].size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++){
			for(size_t j=0; j<data[0].size(); j++){
				ofs <<"\t"<<WriteData(log(data[i][j]), precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from two std::vector ----
	template <class Ty_ >
	static bool write_itx(const std::vector< Ty_ >& datax, const std::vector< Ty_ >& datay, const std::string filename, const std::string wavenamex, const std::string wavenamey, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(datax.size()!=datay.size()){
			std::cerr << "!!!  same number of datapoints are required  !!!" << std::endl;
			getchar();exit(1);
		}
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<datax.size()<< ")\t"<<wavenamex <<"\t"<< wavenamey <<"\nBEGIN\n"; 
		for(size_t i=0; i<datax.size();i++)
			ofs<<"\t"<< WriteData(datax[i], precision) << "\t" << WriteData(datay[i], precision) << std::endl;
		ofs<<"END"<<std::endl;
		ofs.close();
		return true;
	}

//----	writing the edge vector to the file	----
	template< class TyI_ >
	static bool write_edgeVector(const std::vector< TyI_ >& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<< ")\t"<<wavename <<"\nBEGIN\n"; 
		for(size_t i=0; i<data.size();i++)
			ofs<<"\t"<< WriteData(data[i], precision) << std::endl;
		ofs<<"END"<<std::endl;
		
		//---	edge vector	---
		std::vector<double> edge(data.size()+1);
		if(data.size()<=1) return false;
		edge[0] = 1.5*data[0]-0.5*data[1];
		for(size_t i=1; i<data.size(); i++)	edge[i] = 0.5*data[i-1]+0.5*data[i];
		edge[data.size()] = 2.0*data[data.size()-1] - 1.0*data[data.size()-2];

		ofs << "WAVES/D/O/N=("<<data.size()+1<< ")\t"<<wavename+"_edge" <<"\nBEGIN\n"; 
		for(size_t i=0; i<data.size()+1;i++)
			ofs<<"\t"<< WriteData(edge[i], precision) << std::endl;
		ofs<<"END"<<std::endl;

		ofs.close();
		return true;
	}

//----	add comment to existing itx file ----
	static bool itx_add_comment(const std::string filename, const std::string comment, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs<<"x //"<<comment <<std::endl;
		ofs.close();
		return true;
	}

//----	open the igor file on windows	---
	static bool itx_open(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
//		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		std::string command = "cmd.exe /c" + std::string("\"") + filename + std::string("\"");
		system(command.c_str());
		return true;
	}

	static bool write_itx(const double start, const double end, const size_t size, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<size<< ")\t"<<wavename <<"\nBEGIN\n"; 

		for(size_t i=0; i<size;i++){
			double data = start + (end-start)/(size-1)*i;		
			ofs<<"\t"<< data << std::endl;
		}
		ofs<<"END"<<std::endl;
		
		ofs.close();
		return true;
	}
	
	//----	write from class ----
	template <class Ty_ >
	static bool write_itx(const Ty_& data, const size_t size, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<size<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<size; i++)
			ofs <<"\t"<<WriteData(data[i], precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}

	//----	write from pointer	----
	template <class Ty_ >
	static bool write_itx(const Ty_* data, const size_t size, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<size<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<size; i++)
			ofs <<"\t"<<WriteData(data[i], precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}

//----	for Eigen matrices	-----
#ifdef EIGEN_DEFAULT_DENSE_INDEX_TYPE
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	static bool write_itx
		(const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& data, 
		 const std::string filename, const std::string wavename,
		 const int precision = 6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;

		if (!open(ofs, filename, open_condition)) return false;
		
		//	vector case
		if (data.cols() == 1)
			ofs << "IGOR\nWAVES/D/O/N=(" << data.rows()  << ")\t" << wavename << "\nBEGIN\n";
		//	matrix case
		else
			ofs << "IGOR\nWAVES/D/O/N=(" << data.rows() << ", " << data.cols() << ")\t" << wavename << "\nBEGIN\n";

		ofs << data << std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}

	static Eigen::VectorXd itx_VectorXd(const std::string filename,
		OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<double> src = itx_1stdv(filename, open_condition);
		Eigen::VectorXd rslt;
		rslt.resize(src.size());
		for (size_t i = 0; i < src.size(); ++i){
			rslt[i] = src[i];
		}
		return rslt;
	}

	static Eigen::MatrixXd itx_MatrixXd(const std::string filename,
		OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<double> > src = itx_1stdvv(filename, open_condition);
		Eigen::MatrixXd rslt;
		rslt.resize(src.size(), src[0].size());
		rslt.setZero();
		for (size_t i = 0; i < src.size(); ++i){
			for (size_t j = 0; j < src[i].size(); ++j){
				rslt(i, j) = src[i][j];
			}
		}
		return rslt;
	}

#endif

//----	read an itx file and put values into std::vector<gsl_vector>	---
#ifdef __GSL_VECTOR_H__
	static std::vector<gsl_vector *> itx_gslv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		std::vector<gsl_vector *> gslv;
		for(size_t i=0; i<datas.size(); i++){
			for(size_t k=0; k<datas[i][0].size(); k++){
				gsl_vector *gsldata=gsl_vector_calloc(datas[i].size());
				for(size_t j=0; j<datas[i].size(); j++){
					gsl_vector_set(gsldata, j, datas[i][j][k]);
				}
				gslv.push_back(gsldata);
			}
		}
		return gslv;
	}
//----	writ to itx file from gsl_vector ----
	static bool write_itx(const gsl_vector *data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data->size<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data->size; i++){
			ofs <<"\t" <<WriteData(gsl_vector_get(data,i), precision)<<std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from gsl_matrix ----
	static bool write_itx(const gsl_matrix *data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data->size1<<", "<<data->size2<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data->size1; i++){
			for(size_t j=0; j<data->size2; j++){
				ofs <<"\t"<<WriteData(gsl_matrix_get(data,i,j), precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
#endif

//----	writ to itx file from std::vector<gsl_complex>  ----
#ifdef __GSL_COMPLEX_H__
	static bool write_itx(const std::vector<gsl_complex>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/C/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(GSL_REAL(data[i]), precision) << "\t"<<WriteData(GSL_IMAG(data[i]), precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}


	static bool write_itx(const std::vector<std::vector< gsl_complex >>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		for(size_t i=0; i<data.size(); i++){
			if(data[0].size()!=data[i].size()){
				std::cerr << "!!!   data should be same size   !!!" << std::endl;
				getchar();getchar();exit(1);
			}
		}
		ofs << "IGOR\nWAVES/D/O/C/N=("<<data.size()<<", "<<data[0].size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++){
			for(size_t j=0; j<data[0].size(); j++){
				ofs <<"\t"<<WriteData(GSL_REAL(data[i][j]), precision)<<"\t"<<WriteData(GSL_REAL(data[i][j]), precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}


#endif
//----	writ to itx file from std::vector<mygsl::complex>  ----
#ifdef _MYGSL_COMPLEX_HPP_
	static bool write_itx(const std::vector<mygsl::complex>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		if(!open(ofs, filename, WaitFileOpenUntilReturnPressed)) return false;
		ofs << "IGOR\nWAVES/D/O/C/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(data[i].real(), precision) << "\t"<<WriteData(data[i].imag(), precision)<<std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}
#endif
//----	read an itx file and put values into vector<gsl_matrix>	---
#ifdef __GSL_MATRIX_H__
	static std::vector<gsl_matrix *> itx_gslm(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		std::vector<gsl_matrix *> gslm;
		for(size_t i=0; i<datas.size(); i++){
			gsl_matrix *gsldata=gsl_matrix_calloc(datas[i].size(), datas[i][0].size());
			for(size_t j=0; j<datas[i].size(); j++){
				for(size_t k=0; k<datas[i][j].size(); k++){
					gsl_matrix_set(gsldata, j, k, datas[i][j][k]);
				}
			}
			gslm.push_back(gsldata);
		}
		return gslm;
	}
#endif

//----	read an itx file and put values into vector<gsl::vector>	---
#ifdef CCGSL_VECTOR_HPP
	static std::vector<gsl::vector> itx_ccgslv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		std::vector<gsl::vector> ccgslv;
		for(size_t i=0; i<datas.size(); i++){
			for(size_t k=0; k<datas[i][0].size(); k++){
				gsl::vector v_tmp(datas[i].size());
				for(size_t j=0; j<datas[i].size(); j++){
					v_tmp[j]=datas[i][j][k];
				}
				ccgslv.push_back(v_tmp);
			}
		}
		return ccgslv;
	}
	//	if the file is known to contain only 1 vector
	static gsl::vector itx_1ccgslv(const std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<gsl::vector> ccgslv = itx_ccgslv(filename, open_condition);
		if(ccgslv.size()!=1){
			std::cerr << "function in IGORdata::itx_ccgslv1\n the file contains more than 1 vector"<< std::endl;
			getchar();exit(1);
		}
		return ccgslv[0];
	}
//----	writ to itx file from gsl::vector ----
	static bool write_itx(const gsl::vector data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++){
			ofs <<"\t" <<WriteData(data[i], precision)<<std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
//----	writ to itx file from std::vector<gsl::vector>  ----
	static bool write_itx(const std::vector<gsl::vector>& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		for(size_t i=0; i<data.size(); i++){
			if(data[0].size()!=data[i].size()){
				std::cerr << "!!!   data should be same size   !!!" << std::endl;
				getchar();getchar();exit(1);
			}
		}
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<<", "<<data[0].size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++){
			for(size_t j=0; j<data[0].size(); j++){
				ofs <<"\t"<<WriteData(data[i][j], precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
		
//----	writ to itx file from set of gsl::vector ----
	static bool write_itx(const std::vector<gsl::vector>& datav, const std::string filename, const std::vector<std::string> wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(datav.size()!=wavename.size()){
			std::cerr << "!!!  same number of wavenames are required  !!!" << std::endl;
			getchar();exit(1);
		}
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR" <<std::endl;
		size_t i = 0;
		while(i<datav.size()){
			ofs<<"WAVES/D/O"; 
			size_t j=i;
			for(j=i; j<datav.size();j++){			//	datav[i] ~ datav[j] are same size
				if(datav[i].size()!=datav[j].size()) break;
			}
			for(size_t ii=i; ii<j; ii++)
				ofs<<"\t"<<wavename[ii];
			ofs<<std::endl<<"BEGIN"<<std::endl;
			for(size_t k=0; k<datav[i].size();k++){
				for(size_t ii=i; ii<j; ii++)
					ofs<<"\t"<< WriteData(datav[ii][k], precision);
				ofs<<std::endl;
			}
			ofs<<"END"<<std::endl;
			i=j;
		}
		ofs.close();
		return true;
	}

//----	writ to itx file from two gsl::vector ----
	static bool write_itx(const gsl::vector& datax, const gsl::vector& datay, std::string filename, std::string wavenamex, std::string wavenamey, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(datax.size()!=datay.size()){
			std::cerr << "!!!  same number of datapoints are required  !!!" << std::endl;
			getchar();exit(1);
		}
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<datax.size()<< ")\t"<<wavenamex <<"\t"<< wavenamey <<"\nBEGIN\n"; 
		for(size_t i=0; i<datax.size();i++)
			ofs<<"\t"<< WriteData(datax[i], precision) << "\t" << WriteData(datay[i], precision) << std::endl;
		ofs<<"END"<<std::endl;
		ofs.close();
		return true;
	}

//----	writing the edge vector to the file	----
	static bool write_edgeVector(const gsl::vector& data, const std::string filename, const std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size()<< ")\t"<<wavename <<"\nBEGIN\n"; 
		for(size_t i=0; i<data.size();i++)
			ofs<<"\t"<< WriteData(data[i], precision) << std::endl;
		ofs<<"END"<<std::endl;
		
		//---	edge vector	---
		std::vector<double> edge(data.size()+1);
		if(data.size()<=1) return false;
		edge[0] = 1.5*data[0]-0.5*data[1];
		for(size_t i=1; i<data.size(); i++)	edge[i] = 0.5*data[i-1]+0.5*data[i];
		edge[data.size()] = 2.0*data[data.size()-1] - 1.0*data[data.size()-2];

		ofs << "WAVES/D/O/N=("<<data.size()+1<< ")\t"<<wavename+"_edge" <<"\nBEGIN\n"; 
		for(size_t i=0; i<data.size()+1;i++)
			ofs<<"\t"<< WriteData(edge[i], precision) << std::endl;
		ofs<<"END"<<std::endl;

		ofs.close();
		return true;
	}

#endif
#ifdef CCGSL_COMPLEX_HPP
//----	writ to itx file from std::vector<gsl::complex>  ----
	static bool write_itx(const std::vector<gsl::complex>& data, std::string filename, std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/C/N=("<<data.size()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size(); i++)
			ofs <<"\t"<<WriteData(data[i].real(), precision)<<"\t"<<WriteData(data[i].imag(), precision) << std::endl;
		ofs << "END\n";
		ofs.close();
		return true;
	}
#endif
//----	read an itx file and put values into std::vector<gsl::matrix>	---
#ifdef CCGSL_MATRIX_HPP
	static std::vector<gsl::matrix> itx_ccgslm(std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<std::vector<double>>> datas = IGORdata::itx_read(filename, open_condition);
		std::vector<gsl::matrix> ccgslm;
		for(size_t i=0; i<datas.size(); i++){
			gsl::matrix m_tmp(datas[i].size(),datas[i][0].size());
			for(size_t k=0; k<datas[i][0].size(); k++){
				for(size_t j=0; j<datas[i].size(); j++){
					m_tmp[j][k]=datas[i][j][k];
				}
			}
			ccgslm.push_back(m_tmp);
		}
		return ccgslm;
	}
	//	if the file is known to contain only 1 matrix
	static gsl::matrix itx_1ccgslm(std::string filename, 
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<gsl::matrix> ccgslm = itx_ccgslm(filename, open_condition);
		if(ccgslm.size()!=1){
			std::cerr << "function in IGORdata::itx_ccgslm1\n the file contains more than 1 matrix"<< std::endl;
			getchar();exit(1);
		}
		return ccgslm[0];
	}

//----	writ to itx file from gsl::matrix ----
	static bool write_itx(const gsl::matrix& data, std::string filename, std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/N=("<<data.size1()<<", "<<data.size2()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size1(); i++){
			for(size_t j=0; j<data.size2(); j++){
				ofs <<"\t"<<WriteData(data[i][j], precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
#endif
#ifdef CCGSL_MATRIX_COMPLEX_HPP
//----	writ to itx file from gsl::matrix_complex ----
	static bool write_itx(const gsl::matrix_complex& data, std::string filename, std::string wavename, 
			const int precision=6, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ofstream ofs;
		 
		if(!open(ofs, filename, open_condition)) return false;
		ofs << "IGOR\nWAVES/D/O/C/N=("<<data.size1()<<", "<<data.size2()<<")\t"<<wavename<<"\nBEGIN\n";
		for(size_t i=0; i<data.size1(); i++){
			for(size_t j=0; j<data.size2(); j++){
				ofs <<"\t"<<WriteData(data[i][j].real(), precision)<<"\t"<<WriteData(data[i][j].imag(), precision);
			}
			ofs << std::endl;
		}
		ofs << "END\n";
		ofs.close();
		return true;
	}
#endif

}
#endif
