#ifndef __LHD_RETRIEVE_CCD_TIMING_HPP__
#define __LHD_RETRIEVE_CCD_TIMING_HPP__

#include "LHDretrieve_data.hpp"
#include <vector>

class LHDretrieve_CCDtiming:public LHDretrieve_data{
public:
	LHDretrieve_CCDtiming(void);
	enum camera_name{
		Andor_133VIS,
		Andor_Kyoto,
		Andor_50cmVIS2,
		Hamamatsu_Flash4
	};
	const size_t channel_133mVIS;
	const size_t channel_KyotoDV435;
	const size_t channel_50cmVIS2;
	const size_t channel_Flash4;
	std::string getCameraName()const;
	
	void getDataFromRawData(const void* params = NULL);

	const std::vector<double>& getExposureTiming()const{return exposureTiming;}
	const std::vector<double>& getExposureTime()const{return exposureTime;}
	double getAvgExposureTime()const{return avgExposureTime;}
	bool retrieve(enum camera_name name, const size_t shotnum, const double bgtime);
	bool retrieve(const std::string diagname, const size_t shotnum, const size_t ch, double bgtime){
		return LHDretrieve_data::retrieve(diagname, shotnum, ch, bgtime);
	}
		
private:
	std::vector<double> exposureTiming;
	std::vector<double> exposureTime;
	double avgExposureTime;
	enum camera_name cameraName;
};

inline LHDretrieve_CCDtiming::LHDretrieve_CCDtiming(void):LHDretrieve_data(), 
//--- write the channel information here ---
	channel_133mVIS(124),
	channel_KyotoDV435(123),
	channel_50cmVIS2(125),
	channel_Flash4(125)
{
}

inline bool LHDretrieve_CCDtiming::retrieve(enum camera_name name, const size_t shotnum, const double bgtime){
	cameraName = name;
	switch(cameraName){
	case Andor_133VIS:
		return LHDretrieve_data::retrieve("Halpha", shotnum, channel_133mVIS, bgtime);
		break;
	case Andor_Kyoto:
		return LHDretrieve_data::retrieve("Halpha", shotnum, channel_KyotoDV435, bgtime);
		break;
	case Andor_50cmVIS2:
		return LHDretrieve_data::retrieve("Halpha", shotnum, channel_50cmVIS2, bgtime);
		break;
	case Hamamatsu_Flash4:
		//	16th cycle (2013)
		if (shotnum < 125519){
			return LHDretrieve_data::retrieve("Brems", shotnum, channel_Flash4, bgtime);
		}
		//	17th cycle (2014)
		else{
			return LHDretrieve_data::retrieve("Halpha", shotnum, channel_KyotoDV435, bgtime);
		}
		break;
	}
	return false;
}

inline void LHDretrieve_CCDtiming::getDataFromRawData(const void* params){
	LHDretrieve_data::getDataFromRawData();

	exposureTiming.resize(0);
	exposureTime.resize(0);
	
	if(!IsAllocated()) return;
	double threashold = 2500.0;

	double expStart_t=0;
	double expEnd_t=0;
	for(size_t i=0; i<size()-1; ++i){
		if((1.0*data[i] - threashold) * (1.0*data[i+1] - threashold) <= 0.0){
			// rising edge, exposure_start
			if((1.0*data[i] - threashold)<0.0)
				expStart_t = 0.5*get_t(i)+0.5*get_t(i+1);
			else{
				expEnd_t = 0.5*get_t(i)+0.5*get_t(i+1);
				exposureTiming.push_back(0.5*expStart_t+0.5*expEnd_t);
				exposureTime.push_back(expEnd_t-expStart_t);
			}
		}
	}

	avgExposureTime=0.0;
	for(size_t i=0; i<exposureTime.size(); ++i)
		avgExposureTime += exposureTime[i];
	avgExposureTime/=exposureTime.size();
}

inline 	std::string LHDretrieve_CCDtiming::getCameraName()const{
switch(cameraName){
	case Andor_133VIS:
		return "Andor_133VIS";
		break;
	case Andor_Kyoto:
		return "Andor_Kyoto";
		break;
	case Andor_50cmVIS2:
		return "Andor_50cmVIS2";
		break;
	case Hamamatsu_Flash4:
		return "Hamamatsu_Flash4";
		break;
	}
	return "";

}

#endif