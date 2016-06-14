#ifndef __LHDRETRIEVE_MULTICHPMT_HPP__
#define __LHDRETRIEVE_MULTICHPMT_HPP__

#include "LHDretrieve_data.hpp"


class LHDretrieve_MultiChPMT: public LHDretrieve_data{
public:
protected:
	//	override functions
	void getDataFromRawData(const void* params = NULL);
private:

};

#include "LHDretrieve_MultiChPMT.inl"
#endif