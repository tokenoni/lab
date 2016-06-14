#ifndef __MYGSL_GENERAL_BINFILE_HPP__
#define __MYGSL_GENERAL_BINFILE_HPP__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

namespace mygsl{
//-------------------------------------------------------//
//   	class for 										 //
//		saving and reading general binary files			 //
//														 //
//-------------------------------------------------------//
	class generalBinfile{
	public:
		enum datatype{
			intiger32bit, unsigned_integer32bit,
			short_intiger16bit, unsigned_short_integer16bit,
			double64bit, float32bit
		} ;

		generalBinfile(datatype);
		~generalBinfile(void);

	//---	for separated para/bin files	---
		bool Parafile_write(const std::string para_file_);
		bool Parafile_read(const std::string para_file_);

		bool Datafile_write_open(const std::string data_file_);		
		bool Datafile_write_add( void* const data);
		bool Datafile_write_close();

	//---	for convined para/bin files	---
/*		bool ParaDatafile_open_write(const std::string para_file_)const;
		bool Datafile_write_add(const void* data)const;
		bool Datafile_write_close()const;
		bool ParaDatafile_read(const std::string para_file_);
*/		
	//---	for data manipuration	---
		void addParams(const std::string param_name, const double param_value);
		void clearParams();
		void addParams(const std::string comment_);

	private:
		std::ifstream ifs_para, ifs_data, ifs_para_data;
		std::ofstream ofs_para, ofs_data, ofs_para_data;
		void addComment();
		datatype dataType;
		std::string data_filename, para_filename;
		std::vector< std::string > param_names;
		std::vector< double > param_values;
		std::string comment;
		
		std::string getDataType_string(const datatype datatype_)const;
		//---	prohivited	---
		generalBinfile(void){};
	};
};
#include "general_binfile.inl"
#endif