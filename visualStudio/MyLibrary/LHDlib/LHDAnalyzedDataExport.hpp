#ifndef __LHD_ANALYZED_DATA_EXPORT__
#define __LHD_ANALYZED_DATA_EXPORT__

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip>

class LHDAnalyzedDataExport
{
public:
	LHDAnalyzedDataExport(void){};
	~LHDAnalyzedDataExport(void){};
};

class LHDAnalyzedDataExport2D{
public:
	LHDAnalyzedDataExport2D(void){};
	~LHDAnalyzedDataExport2D(void){};

	static bool fileExport(const std::string exportFileName, const std::string analyzedDataName, const size_t analyzedShotnum, 
		const std::vector<double> axis1data, const std::string axis1name, const std::string axis1unit, 
		const std::vector<double> axis2data, const std::string axis2name, const std::string axis2unit,
		const std::vector<std::vector<std::vector<double>>> exportData, 
		const std::vector<std::string> exportDataNames, 
		const std::vector<std::string> exortDataUnits, const std::string comment);


};
#include "LHDAnalyzedDataExport.inl"
#endif