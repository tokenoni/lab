#ifndef ANDORSIF_H
#define ANDORSIF_H
#include "AndorSIF2.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#define MAXPATH _MAX_PATH

class AndorSIF
{
public:
	AndorSIF(void){Allocated = false;};
	~AndorSIF(void){FreeMemory();};
	//	read SIF file 
	bool read(std::string filename);
	//	return the exposure time 
	double ExposureTime()const;
	//	return the number of images
	size_t SequenceSize()const;
	//	return the number of vertical data point
	size_t VerticalSize()const;
	//	return the number of horizontal data point
	size_t HorizontalSize()const;
	//	return 

	//	substract the background data
	bool subBackground(const size_t bg_start, const size_t bg_num=1);

	double get(size_t sequence_index, size_t vertical_index, size_t horizontal_index)const;
	double operator () (size_t sequence_index, size_t vertical_index, size_t horizontal_index)const;
	std::vector<double> get(size_t sequence_index, size_t vertical_index)const;
	std::vector<std::vector<double>> get(size_t sequence_index)const;
	std::vector<std::vector<double>> getHorizontalDevelopment(size_t vertical_index)const;
	std::vector<std::vector<double>> getVerticalDevelopment(size_t horizontal_index)const;
	std::vector<std::vector<std::vector<double>>> getAllFrames()const;

private:
//---	main data	---
	float* data;
	size_t sequence_size, vertical_size, horizontal_size, one_image_size;
//---	SIF data	---
	struct TImage Image;
	struct TInstaImage InstaImage;
	struct TCalibImage CalibImage;

	std::string sif_name;	//	filename
	std::ifstream fin;		//	fstream of SIF
	void read_instaimage();
	void read_calibimage();
	void read_image_structure();
	void read_image();
//---	variables and functions for memory allocation	---
	bool Allocated;
	size_t size;
	void MemoryControl();
	void FreeMemory();
};

//------------------------------------------------//
//												  //
//			small functions for SIF reading		  //
//												  //
//------------------------------------------------//

namespace SIFread{
//---	read 1 string before "terminator" and store in s	---
	inline static void read_string(std::ifstream& fin, std::string& s, char terminator){
		s.clear();
		char ch ='a';
		while (ch != terminator){
	 		fin.read(&ch, 1);
			s.append(&ch,1);
		}
	}
//---	read 1 string before "terminator" and store in s	---
	inline static std::string read_string(std::ifstream& fin, char terminator){
		std::string s="";
		char ch ='a';
		while (ch != terminator){
	 		fin.read(&ch, 1);
			s.append(&ch, 1);
		}
		return s;
	}
//---	read 1 long-integer before "terminator" and return	---
	inline static long read_int(std::ifstream& fin, char terminator){
		std::string s ="";
		char ch='a';
		//	skip spaces
		do{	fin.read(&ch,1); }while(isspace(ch));
		s.append(&ch, 1);
		while(ch!=terminator){
			fin.read(&ch, 1);
			//	skip any letters and spaces
			if(!isspace(ch) && !isalpha(ch)) s.append(&ch, 1);
		}
		return atol(s.c_str());
	}
//---	read 1 long-integer before "terminator" and return	---
	inline static float read_float(std::ifstream& fin, char terminator){
		std::string s ="";
		char ch='a';
		//	skip spaces
		do{	fin.read(&ch,1); }while(isspace(ch));
		s.append(&ch,1);
		while(ch!=terminator){
			fin.read(&ch, 1);
			//	skip any letters and spaces
			if(!isspace(ch) && !isalpha(ch)) s.append(&ch, 1);
		}
		if(s.find("+INF") != std::string::npos) s = "0";
		return (float)atof(s.c_str());
	}
	/*****************************************************************************/
	//Function: read_byte
	//Outputs: A number of type 'int'
	//The purpose of this function is to do a single read of a byte
	//It does this by reading a charachter and outputting its integer value
	//i.e. the Integer equivalent of the ascii value.
	/*****************************************************************************/
	inline static int read_byte_and_skip_terminator(std::ifstream& fin){
		char ch,termin_ch;
		int i;
		ch ='a';
		fin.read(&ch, 1);
		fin.read(&termin_ch, 1);
		i = (int)(ch);                          // gives integer value of ascii code
		return i;
	}
	//---	return the string with src2	---
	template< class TyS_ >
	inline static std::string add_to_string(std::string& src, TyS_& src2){
		std::stringstream s;
		s<< src << src2;
		return s.str();
	}
	/*****************************************************************************/
	//Function: read_len_chars
	//Inputs: The length of the string
	//Outputs: The string
	//The purpose of this function is to read in a predefined number
	//of charachers and return them as a string
	/*****************************************************************************/
	inline static void read_len_chars(std::ifstream& fin, int string_length, std::string& len_chars_buffer)
	{
		char ch;
		ch ='a';
		len_chars_buffer.clear();
		if(string_length >= MAXPATH){
			std::cout<<"file might be broken"<<std::endl;
		}
		for(int i=0;i<string_length;i++){
			fin.read(&ch, 1);         // Reads in a string of length string_length
			if(i<MAXPATH) len_chars_buffer.append(&ch, 1);
		}
	}
	inline static std::string read_len_chars(std::ifstream& fin, int string_length){
		std::string s;
		read_len_chars(fin, string_length, s);
		return s;
	}

};

#include "AndorSIF.inl"
#endif