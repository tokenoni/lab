#pragma once
#ifndef _CVS_HPP_
#define _CVS_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <limits>
#include <cstdlib>

namespace CVSdata{
	static bool IsNumeric(std::string p)
	{
		double val=atof(p.c_str());
		if (val != 0.0) return true;
		else if(val == 0.0 && p[p.find_first_not_of(" \t")]!= '0') return false;
		else return true;
	}

//---	CVS::split function	---
	static std::vector<double> split(const std::string &s, const std::string d) 
	{
		std::vector<double> elements;							//	CVS::split std::string s and convert to double, 	
		size_t i_start=0, i_end=0;								//	and set to vector element
		while(i_end!=std::string::npos){						//
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

	static enum OpenCondition{
		WaitFileOpenUntilReturnPressed,
		Ignore,
		TryJapanese
	};

//---	wait for the file opening	---
	static bool open(std::ifstream &ifs, std::string filename, 
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
	static bool open(std::ofstream &ofs, std::string filename, 
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
	static bool open_check(std::string filename, OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::ifstream ifs(filename.c_str(), std::ifstream::in );
		if(!ifs.is_open()){
			return false;
		}else return true;
	}


//----	read an itx file and put values into std::vector<std::vector<double>>	---
	static std::vector<std::vector<double>> read(std::string filename, std::string comment_indicator, std::string deliminator,
			OpenCondition open_condition = WaitFileOpenUntilReturnPressed)
	{
		std::vector<std::vector<double>> all_data;
		std::ifstream ifs;
		open(ifs, filename, open_condition);

		std::string sdata;
		while(!ifs.eof()){
			getline(ifs,sdata);
			//---	comment	---
			size_t comment_found = sdata.find(comment_indicator);
			if(comment_found != std::string::npos)
				sdata = sdata.substr(comment_found-1);

			if(sdata.size() == 0) break;

			std::vector<double>row_data = split(sdata, deliminator);
			if(row_data.size()>0) all_data.push_back(row_data);
		}
		ifs.close();
		return all_data;
	}
}

#endif