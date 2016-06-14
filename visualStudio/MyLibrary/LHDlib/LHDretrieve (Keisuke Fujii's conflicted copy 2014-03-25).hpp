#ifndef __LHDRETRIEVE_HPP__
#define __LHDRETRIEVE_HPP__

#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>
#include "include/retrieve.h"
#include <IGOR\IGORitx.hpp>


class LHDretrieve
{
public:
	//---	core functions	---
	LHDretrieve(void);
	LHDretrieve(const LHDretrieve& obj);
	LHDretrieve& operator = (const LHDretrieve& obj);
	void copy(const LHDretrieve& obj);
	~LHDretrieve(void){clear();};

	//---	core functions	---
	bool run(  const std::string diagname, const size_t shot_number, const size_t channel_number);
	void clear();
	
	//--------------------		getting functions	------------------------------
	//
	//	all data is stored as char*, while usually data is short-type
	//	get(i) function returns short value if the data is stored in Little Endiannness
	//	if any other type value is nescessary,
	//	use other getting functions 
	//
	//----------------------------------------------------------------------------
	short get(const size_t i)const{return getAsUnsignedShort_with_LittleEndianness(i);}
	short operator [] (const size_t i)const{return get(i);};
	double get_t(const size_t i)const{return data_t[i];}

	unsigned short getAsUnsignedShort_with_LittleEndianness( const size_t i)const;
	unsigned short getAsUnsignedShort_with_BigEndianness( const size_t i)const;
	short getAsShort_with_LittleEndianness( const size_t i)const;
	short getAsShort_with_BigEndianness( const size_t i)const;

	//---	compareing function ---
	bool operator == (const LHDretrieve& src);

	//---	tip functions ---
	size_t size()const       {if(isAllocated) return data_length_actual/2; else return 0;}
	bool IsAllocated()const  {return isAllocated;}
	bool IsAllocated_t()const{return isAllocated_t;}
	double getStart_t()const {if(!isAllocated) return 0; if(isAllocated_t) return data_t[0];else return 0.0;}
	double getEnd_t()const   {if(!isAllocated) return 0; if(isAllocated_t) return data_t[data_length_actual-1];else return data_length_actual-1;}

	//---	getParameterValue	---
	double getParameter(const std::string param_name)const;

private:
	
	bool isAllocated;
	char *data;
	
//---	parameters for retrieve	---
	int retrieve_index;
	size_t data_length, comp_length, data_length_actual;
	unsigned short param_count, value_len;
	short data_type;
	char image_type;
	int is_nframe;
	char management[32*8], server[32*8];
	char comment[256]; 
	int comment_size;


//--- parameters for the data	---
	bool getParams(const size_t channel_number);
	std::vector<std::string> param_names;
	std::vector<double> param_values;
	std::vector<int> param_types;

	bool getTiming(const std::string diagname, const size_t shot_number, const size_t channel_number);
	bool isTimingGot;	//	true if the timing data is available
	double *data_t;
	bool isAllocated_t;
};

#include "LHDretrieve.inl"
#endif
