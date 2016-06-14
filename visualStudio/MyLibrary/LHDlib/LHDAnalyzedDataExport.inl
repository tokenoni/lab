
//------------------------------------------------------------
//			LHDAnalyzedDataExport2D
//------------------------------------------------------------

inline bool LHDAnalyzedDataExport2D::fileExport(
	const std::string exportFileName, const std::string analyzedDataName, const size_t analyzedShotnum,
	const std::vector<double> axis1data, const std::string axis1name, const std::string axis1unit,
	const std::vector<double> axis2data, const std::string axis2name, const std::string axis2unit,
	const std::vector<std::vector<std::vector<double>>> exportData, const std::vector<std::string> exportDataNames, const std::vector<std::string> exortDataUnits, 
	const std::string comment)
{
	std::ofstream ofs;
	ofs.open(exportFileName);
	if (!ofs.is_open()){
		std::cerr << "error in opening file of" << exportFileName << std::endl;
		return false;
	}

	size_t valNo = exportData.size(); 
	size_t numAxis1 = axis1data.size();
	size_t numAxis2 = axis2data.size();

	//	data dim check
	//	for valNo
	if (exportDataNames.size() != valNo){
		std::cerr << "size of exportDataNames is invalid" << std::endl;
		return false;
	}
	if (exortDataUnits.size() != valNo){
		std::cerr << "size of exortDataUnits is invalid" << std::endl;
		return false;
	}
	// for numAxis1
	for (size_t k = 0; k < valNo; ++k){
		if (exportData[k].size() != numAxis1){
			std::cerr << "size1 of " << exportDataNames[k] << " is invalid" << std::endl;
			return false;
		}
		if (exportData[k][0].size() != numAxis2){
			std::cerr << "size2 of " << exportDataNames[k] << " is invalid" << std::endl;
			return false;
		}
	}

	ofs << "# [Parameters]" << std::endl;
	ofs << "# Name = '" << analyzedDataName << "'"<< std::endl;
	ofs << "# ShotNo = " << analyzedShotnum << std::endl;
	
	//	obtaining the current time
	time_t timer = time(NULL);
	struct tm date;
	localtime_s(&date, &timer);

	ofs << "# Date = '" 
		<< std::setw(2) << std::setfill('0') << date.tm_mon+1 << "/"
		<< std::setw(2) << std::setfill('0') << date.tm_mday << "/"
		<< std::setw(4) << std::setfill('0') << date.tm_year + 1900 << " "
		<< std::setw(2) << std::setfill('0') << date.tm_hour << ":" 
		<< std::setw(2) << std::setfill('0') << date.tm_min << "'" << std::endl;
	ofs << "# " << std::endl;
	
	ofs << "# DimNo = 2" << std::endl;
	ofs << "# DimName = '" << axis1name << "', '" << axis2name << "'" << std::endl;
	ofs << "# DimSize = " << axis1data.size() << ", " << axis2data.size() << std::endl; 
	ofs << "# DimUnit = '" << axis1unit << "', '" << axis2unit << "'" << std::endl;
	ofs << "# " << std::endl;

	ofs << "# ValNo = " << valNo << std::endl;
	
	ofs << "# ValName = '" << exportDataNames[0];
	for (size_t k = 1; k < valNo; ++k){
		ofs << "', '" << exportDataNames[k];
	}
	ofs << "'"<< std::endl;
	ofs << "# ValUnit = '" << exortDataUnits[0];
	for (size_t k = 1; k < valNo; ++k){
		ofs << "', '" << exortDataUnits[k];
	}
	ofs << "'" << std::endl;
	ofs << "# " << std::endl;

	//	comment
	ofs << "# [Comments]" << comment << std::endl;


	//	data export
	ofs << "# [data]" << std::endl;
	for (size_t i = 0; i < numAxis1; ++i){
		for (size_t j = 0; j < numAxis2; ++j){
			ofs << "\t" << axis1data[i] << ", \t" << axis2data[j];
			for (size_t k = 0; k < valNo; ++k){
				ofs << ", \t" << exportData[k][i][j];
			}
			ofs << std::endl;
		}
	}
	return true;
}
