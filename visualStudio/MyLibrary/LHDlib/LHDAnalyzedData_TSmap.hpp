#ifndef __LHD_ANALYZED_DATA_TSMAP_HPP__
#define __LHD_ANALYZED_DATA_TSMAP_HPP__
#include "LHDAnalyzedData.hpp"

class LHDAnalyzedData_TSmap: public LHDAnalyzedData{
public:
	void getNeCalibAllTimings();
	const std::vector<std::vector<double>>& getTe() const{return getMatrix("Te");};
	const std::vector<std::vector<double>>& getdTe()const{return getMatrix("dTe");};
	const std::vector<std::vector<double>>& getNe() const{return ne;};
	const std::vector<std::vector<double>>& getdNe()const{return dne;};
	std::vector<double> getTiming()const{return getXdata1dim();};
	std::vector<double> getRadius()const{return getYdata1dim();};
	std::vector<double> getLaser() const{return ReduceTo1dimension(Xaxis, "laser number");};
	const std::vector<std::vector<double>>& getReffective()const{return getMatrix("reff");};
	const std::vector<std::vector<double>> getAbsReffective();
	std::vector<double> getMatrixAtNearestTime(const double timing, const std::string diagname)const;

private:
	std::vector<std::vector<double>> ne;
	std::vector<std::vector<double>> dne;

};
#include "LHDAnalyzedData_TSmap.inl"

#endif
