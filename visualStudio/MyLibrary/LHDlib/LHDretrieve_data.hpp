#ifndef __LHDRETRIEVE_DATA_HPP__
#define __LHDRETRIEVE_DATA_HPP__

#include "LHDretrieve.hpp"
#include "TimeDomainData.hpp"

class LHDretrieve_data: public TimeDomainData{
public:
	//---	constructor - destroctors	---
	LHDretrieve_data(void);
	LHDretrieve_data(const size_t num);
	LHDretrieve_data& operator=(const LHDretrieve_data& obj);
	LHDretrieve_data(const LHDretrieve_data& obj);
	LHDretrieve_data& copy(const LHDretrieve_data& obj);
	~LHDretrieve_data(void){};

//--- retrieve and prepare the data ---
	bool retrieve(const std::string& diagname, const size_t shotnum, const size_t chnum, const double bgtime = 0.0);

//--- convert raw_data to data, rate, start_t ---
	virtual void getDataFromRawData(const void* params = NULL);
	void getBackground(const double bgtime);
	size_t getCh()const     {if(IsAllocated()) return channel; else return 0;}
	

//--- stored data ---
protected:
	size_t channel;
	LHDretrieve raw_data;
	const LHDretrieve& getRawData()const{return raw_data;}


};

#include "LHDretrieve_data.inl"

#endif
