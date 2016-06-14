#ifndef __TIME_DOMAIN_DATA_HPP__
#define __TIME_DOMAIN_DATA_HPP__


class TimeDomainData{
public:
	TimeDomainData(void);
	TimeDomainData(const size_t num);
	TimeDomainData(const TimeDomainData& obj);
	~TimeDomainData(void);
	TimeDomainData& operator = (const TimeDomainData& obj);
	TimeDomainData& copy(const TimeDomainData& obj);
	

protected:
	void getBackground(const double bgtime);
//--- access to the data ---
public:
	bool IsAllocated ()const{return isAllocated;}
	bool IsPrepared_t()const{return isPrepared_t;}
	size_t size()const{return data_num;}
	const double* get()const{return data;}		// returns data as a pointer

	// the following functions have no error saving functions
	double& operator [](const size_t i){return data[i];}
	double operator [](const size_t i)const{return data[i];}
	double& get(const size_t i){return data[i];};		
	double get(const size_t i)const{return data[i];};		
	const double* get_p(const size_t i)const{return data+i;}
	double* get_p(const size_t i){return data+i;}
	double getAvg(const int start_index, const size_t avgnum)const;
	bool getAvgCopy(TimeDomainData& dest, const size_t avgnum)const;
	std::vector<double> getAvg(const size_t avgnum=1)const;
	std::vector<double> getAvg_t(const size_t avgnum=1)const;

//--- tip functions ---
	double get_t(const size_t i)const{return start_t + 1.0*i/rate;}
	double getAvg_t(const int start_index, const size_t avgnum)const;
	double getRate()const   {if(isPrepared_t)return rate;else return 0.0;}
	double getStart_t()const{if(isPrepared_t)return start_t;else return 0.0;}
	double getEnd_t()const  {if(isPrepared_t)return end_t;else return 0.0;}
	size_t getTimeFromIndex(const double index_)const;
	size_t getIndexFromTime(const double time_)const;
	void setTiming(const double rate_, const double start_t);

//--- memory controlling functions ----
public:
	void clear();
	void Allocate(const size_t num);

//--- smart copy function. the data in the copied object cannot be used any more. ---
//--- this function will be used in the inheriting object ---
	void SmartCopy(TimeDomainData& obj);
//--- function for smart copy. Don't use this function except in the function of "SmartCopy"!! ---
	double* _SmartCopied();

//--- stored data ---
protected:
	double *data;
	double rate, start_t, end_t;
	bool isPrepared_t;		// true if rate and start_t are preapared
private:
	size_t data_num;
	bool isAllocated;		// true if memory for the data is allocated 

};

#include "TimeDomainData.inl"
#endif
