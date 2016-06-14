#ifndef __LHD_ANALYZED_DATA__
#define __LHD_ANALYZED_DATA__

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

class LHDAnalyzedData
{
public:
	enum Dimension{
		Xaxis,
		Yaxis
	};

	LHDAnalyzedData(const std::string filename){ read(filename); };
	LHDAnalyzedData(void){};
	~LHDAnalyzedData(void){};
	bool downloadFile(const std::string igetfilePath, const size_t shotnum, const std::string diagname, const std::string outputDir);
	bool read(const std::string filename);

	//	returns the data as matrix (xdata vs ydata, diagnostics is specified), used for 2-dim file
	const std::vector<std::vector<double>>& getMatrix(const size_t index)const;	
	const std::vector<std::vector<double>>& getMatrix(const std::string diagname)const;
	//	returns the data as vector (vs ydata, xdata and diagnostics is specified), used for 
	std::vector<double> getMatrixAtNearestTime(const double timing, const size_t index)const{ return getMatrix(index)[getNearestIndex(timing)];}
	std::vector<double> getMatrixAtNearestTime(const double timing, const std::string diagname)const{ return getMatrix(diagname)[getNearestIndex(timing)];};
	
	//	returns the data as vector, used for 1-dim file
	const std::vector<double>& getVector(const size_t index)const;
	const std::vector<double>& getVector(const std::string diagname)const;
	

	std::vector<double> ReduceTo1dimension(const Dimension axis, const size_t index)const;
	std::vector<double> ReduceTo1dimension(const Dimension axis, const std::string diagname)const;

	double ReduceTo1Value(const size_t index)const{return data_table_vector[index][0][0];}
	double ReduceTo1Value(const std::string diagname)const;

	// returns the index for the "diagname"
	int SearchDiagname(const std::string diagname)const;
	std::vector<double> getTiming()const{return getXdata1dim();};

	// used for 1-dim file
	const std::vector<double>& getXdata()const{return x_data;}
	// used for 2-dim file
	const std::vector<double>& getXdata1dim()const{return x_data_1dim;}
	const std::vector<double>& getYdata1dim()const{return y_data_1dim;}

	size_t size()const{return x_num;}
	size_t size_x()const{return x_num;}
	size_t size_y()const{return y_num;}
	//--- get nearest timing ---
	size_t getNearestIndex(const double timing)const;
protected:
	size_t dimension, val_number;
	std::string x_name, y_name;
	std::string x_unit, y_unit;
	size_t x_num, y_num;

	size_t val_num;
	std::vector<double> x_data, y_data, x_data_1dim, y_data_1dim;	//	_1dim is used for 2dimensional file
	std::vector<std::vector<double>> data_table;
	std::vector<std::string> data_name, data_unit;

	bool SetIntoMatrix();
	bool SetIntoVector();
	std::vector<std::vector<std::vector<double>>> data_table_vector;	// used for 2 dimensional file

};


#include "LHDAnalyzedData.inl"
#endif