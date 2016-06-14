#ifndef __LHD_ANALYZED_DATA_TSMESH__
#define __LHD_ANALYZED_DATA_TSMESH__

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <gsl/gsl_interp.h>

class LHDAnalyzedData_TSmesh_Handler{
public:
	LHDAnalyzedData_TSmesh_Handler
		(const double* const t_axis_,  const size_t t_num_,
		 const double* const r_axis_,  const size_t r_num_,
		 const double* const z_axis_,  const size_t z_num_,
		 const double* const phi_axis_,const size_t phi_num_,
 		 const double**** const reff_,
 		 const double**** const theta_);
	~LHDAnalyzedData_TSmesh_Handler();

	double get_reff(const double t, const double r, const double z, const double phi);
	double get_theta(const double t, const double r, const double z, const double phi);

	double getTmin()const{ return t_axis[0]; }
	double getTmax()const{ return t_axis[t_num - 2]; }
	double getRmin()const{ return r_axis[0]; }
	double getRmax()const{ return r_axis[r_num - 2]; }
	double getZmin()const{ return z_axis[0]; }
	double getZmax()const{ return z_axis[z_num - 2]; }
	double getPhimin()const{ return phi_axis[0]; }
	double getPhimax()const{ return phi_axis[phi_num - 2]; }
private:
	double get(const double t, const double r, const double z, const double phi, const double **** const data);

	const size_t t_num, r_num, z_num, phi_num;

	const double* const t_axis;
	const double* const r_axis;
	const double* const z_axis;
	const double* const phi_axis;
	const double**** const reff;
	const double**** const theta;

	gsl_interp_accel* const acc_t;
	gsl_interp_accel* const acc_r;
	gsl_interp_accel* const acc_z;
	gsl_interp_accel* const acc_phi;
};

class LHDAnalyzedData_TSmesh
{
public:
	LHDAnalyzedData_TSmesh(const std::string filename);
	~LHDAnalyzedData_TSmesh(void){};
	bool free();
//	bool downloadFile(const std::string igetfilePath, const size_t shotnum, const std::string diagname, const std::string outputDir);
	bool read(const std::string filename);

	LHDAnalyzedData_TSmesh_Handler getHandler()const;
protected:
	size_t dimension, val_number;
	std::string t_name, r_name, z_name, phi_name;
	std::string t_unit, r_unit, z_unit, phi_unit;
	size_t t_num, r_num, z_num, phi_num;

	size_t val_num;
	bool isAllocated;
	double* t_axis;
	double* r_axis;
	double* z_axis;
	double* phi_axis;
	double**** reff;
	double**** theta;

};


#include "LHDAnalyzedData_TSmesh.inl"
#endif