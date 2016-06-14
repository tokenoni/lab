#ifndef __TIME_DOMAIN_DATA_SET_HPP__
#define __TIME_DOMAIN_DATA_SET_HPP__

#include "TimeDomainData.hpp"
class TimeDomainDataSet{
public:
//--- constructors ---
	TimeDomainDataSet(void);
	TimeDomainDataSet(const size_t chnum);
	TimeDomainDataSet(const TimeDomainDataSet& obj);
	TimeDomainDataSet& operator = (const TimeDomainDataSet& obj);
	TimeDomainDataSet& copy(const TimeDomainDataSet& obj);
	~TimeDomainDataSet(void){};

	//--- data allocation ---
	void resize(const size_t size_){data_v.resize(size_,TimeDomainData());}
	void push_back(const TimeDomainData& data);
	void clear(){data_v.clear();}
	void eraseUnallocatedData();	// erase the unretrieved data. The size of data_v will change.
	void SmartCopy(TimeDomainDataSet& obj);
	void _SmartCopied();

	//--- access to the data ---
	const TimeDomainData& get(const size_t i)const{ return data_v[i];}
	TimeDomainData& get(const size_t i){ return data_v[i];}
	const TimeDomainData& operator [] (const size_t i)const{ return data_v[i];}
	TimeDomainData& operator [] (const size_t i){ return data_v[i];}
	
	//--- no error saving functions for time and channel index ---
	double getAvg(const int t_index, const int ch_index, const size_t t_avgnum=1, const size_t ch_avgnum=1)const;

	//--- tip functions ---
	size_t size() const{return data_v.size();}
	size_t t_size() const;
	double get_t(const size_t index, const size_t avgnum=1)const;
	double getRate()const;
	double getStart_t()const;
	double getEnd_t()const;
	size_t getIndexFromTime(const double time)const;

	void setTiming(const double rate, const double start_t);

	//--- output functions ---
	std::vector<std::vector<double> > getAvg_stdvv(const size_t tsmth, const size_t chsmth=0)const;
	std::vector<std::vector<double> > get_stdvv()const;
	std::vector<std::vector<double> > get_stdvv(const double start_t, const double end_t)const;
	std::vector< double > getAvg_t(const size_t avgnum=0)const{return data_v[0].getAvg_t(avgnum);}
	std::vector<double> get_t_stdv()const;
	std::vector<double> get_t_stdv(const double start_t, const double end_t)const;

protected:
	std::vector<TimeDomainData> data_v;


};

#include "TimeDomainDataSet.inl"

#endif
