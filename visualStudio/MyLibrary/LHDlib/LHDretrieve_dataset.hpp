#ifndef __LHDRETRIEVE_DATASET_HPP__
#define __LHDRETRIEVE_DATASET_HPP__

#include <vector>
#include "TimeDomainDataSet.hpp"
#include "LHDretrieve_data.hpp"
#include "LHDretrieve_MultiChPMT.hpp"

class LHDretrieve_dataset:public TimeDomainDataSet{
public:
	//--- retrieve ---
	bool getRetrieve(const std::string diagname, const size_t shotnum, const size_t chstart=1, const size_t chnum=32, const double bgtime=0.0);
	bool getRetrieveOne(const std::string diagname, const size_t shotnum, const size_t channel, const double bgtime=0.0);
private:
};

#include "LHDretrieve_dataset.inl"
#endif
