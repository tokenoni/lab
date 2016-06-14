class LHDAnalyzedData_TSmesh;
/*
inline bool LHDAnalyzedData_TSmesh::downloadFile(const std::string igetfilePath, const size_t shotnum, const std::string diagname, const std::string outputDir){
	std::stringstream command;
	command << igetfilePath << " -s " << shotnum << " -d " << diagname << " -o " << outputDir;
	system(command.str().c_str());
	if (errno == -1) {
		std::string errormessage = "error in excuting " + command.str();
		throw errormessage;
	};
	return true;
}
*/
inline LHDAnalyzedData_TSmesh::LHDAnalyzedData_TSmesh(const std::string filename){
	isAllocated = false;
	if (!read(filename)){
		std::string errormessage = "error in reading " + filename + "\n";
		throw errormessage;
	}
}

inline bool LHDAnalyzedData_TSmesh::free(){
	if (isAllocated){
		for (size_t t = 0; t < t_num; ++t){
			for (size_t r = 0; r<r_num; ++r){
				for (size_t z = 0; z<z_num; ++z){
					delete[] reff[t][r][z];
					delete[] theta[t][r][z];
				}
				delete[] reff[t][r];
				delete[] theta[t][r];
			}
			delete[] reff[t];
			delete[] theta[t];
		}
		delete[] reff;
		delete[] theta;

		delete[] t_axis, r_axis, z_axis, phi_axis;
	}
	isAllocated = false;
}

inline bool LHDAnalyzedData_TSmesh::read(const std::string filename){

	std::ifstream ifs(filename.c_str(), std::ifstream::in);

	if (!ifs.is_open()){
		std::string errormessage = filename + " is not found";
		throw errormessage;
	}

	//---	reading	---
	std::string sdata;
	std::vector<unsigned short> data_1line_us(8);	//	this vector is used to convert from double to unsigned short
	while (!ifs.eof()){
		getline(ifs, sdata);

		//	reading the header	
		if (sdata.size() > 0 && sdata[0] == '#'){
			//	dimension
			if (sdata.find("DimNo") != std::string::npos){
				std::string dimnum = sdata.substr(sdata.find("=") + 1);
				dimension = (size_t)floor(atof(dimnum.c_str()) + 0.5);
				if (dimension != 4) return false;
			}

			//	x, y, z, w axis
			else if (sdata.find("DimName") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=") + 1);
				t_name = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				r_name = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				z_name = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				phi_name = substr;
			}

			//	x, y, z, w number
			else if (sdata.find("DimSize") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=") + 1);
				t_num = (size_t)floor(atof(substr.substr(0, substr.find_first_of(",")).c_str()) + 0.5);
				substr = substr.substr(substr.find_first_of(",") + 1);
				r_num = (size_t)floor(atof(substr.substr(0, substr.find_first_of(",")).c_str()) + 0.5);
				substr = substr.substr(substr.find_first_of(",") + 1);
				z_num = (size_t)floor(atof(substr.substr(0, substr.find_first_of(",")).c_str()) + 0.5);
				substr = substr.substr(substr.find_first_of(",") + 1);
				phi_num = (size_t)floor(atof(substr.c_str()) + 0.5);

			}

			//	x, y, z, w unit
			else if (sdata.find("DimUnit") != std::string::npos){
				std::string substr = sdata.substr(sdata.find("=") + 1);
				t_unit = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				r_unit = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				z_unit = substr.substr(0, substr.find_first_of(","));
				substr = substr.substr(substr.find_first_of(",") + 1);
				phi_unit = substr;
			}
		}
		if (sdata[0] != '#') break;
	}

//--------	allocation	----------
	isAllocated = true;

	t_axis = new double[t_num];
	r_axis = new double[r_num];
	z_axis = new double[z_num];
	phi_axis = new double[phi_num];

	reff = new double***[t_num];
	theta = new double***[t_num];
	for (size_t t = 0; t < t_num; ++t){
		reff[t] = new double**[r_num];
		theta[t] = new double**[r_num];
		for (size_t r = 0; r<r_num; ++r){
			reff[t][r] = new double*[z_num];
			theta[t][r] = new double*[z_num];
			for (size_t z = 0; z<z_num; ++z){
				reff[t][r][z] = new double[phi_num];
				theta[t][r][z] = new double[phi_num];
			}
		}
	}
//--------	end of allocation	----------


//--------	reading	----------

	size_t count = 0;
	for (size_t t = 0; t < t_num; ++t){
		for (size_t phi = 0; phi < phi_num; ++phi){
			for (size_t z = 0; z < z_num; ++z){
				for (size_t r = 0; r < r_num; ++r){
					count++;
					std::string substr = sdata;
					t_axis[t] = atof(substr.substr(0, substr.find_first_of(",")).c_str());
					substr = substr.substr(substr.find_first_of(",") + 1);
					r_axis[r] = atof(substr.substr(0, substr.find_first_of(",")).c_str());
					substr = substr.substr(substr.find_first_of(",") + 1);
					z_axis[z] = atof(substr.substr(0, substr.find_first_of(",")).c_str());
					substr = substr.substr(substr.find_first_of(",") + 1);
					phi_axis[phi] = atof(substr.substr(0, substr.find_first_of(",")).c_str());
					substr = substr.substr(substr.find_first_of(",") + 1);

					reff[t][r][z][phi] = fabs(atof(substr.substr(0, substr.find_first_of(",")).c_str()));
					substr = substr.substr(substr.find_first_of(",") + 1);
					theta[t][r][z][phi] = atof(substr.substr(0, substr.find_first_of(",")).c_str());
					substr = substr.substr(substr.find_first_of(",") + 1);

					if (ifs.eof() || sdata.size() == 0){
						ifs.close();
						if (count == t_num * r_num * z_num * phi_num)
							return true;
						else return false;
					}
					getline(ifs, sdata);
				}
			}
		}
	}
	return true;
}

inline LHDAnalyzedData_TSmesh_Handler LHDAnalyzedData_TSmesh::getHandler()const{
	return LHDAnalyzedData_TSmesh_Handler
		( (const double* const) t_axis,   t_num, 
		  (const double* const) r_axis,   r_num, 
		  (const double* const) z_axis,   z_num, 
		 phi_axis, phi_num, 
		 (const double**** const)reff, (const double**** const)theta);
};

//-----------------------------------------------------
//
//			LHDAnalyzedData_TSmesh_Handler
//
//-----------------------------------------------------

inline LHDAnalyzedData_TSmesh_Handler::LHDAnalyzedData_TSmesh_Handler
(const double* const t_axis_, const size_t t_num_,
 const double* const r_axis_, const size_t r_num_,
 const double* const z_axis_, const size_t z_num_,
 const double* const phi_axis_, const size_t phi_num_,
 const double**** const reff_,
 const double**** const theta_)
:
t_axis(t_axis_),
r_axis(r_axis_),
z_axis(z_axis_),
phi_axis(phi_axis_),
t_num(t_num_),
r_num(r_num_),
z_num(z_num_),
phi_num(phi_num_),
reff(reff_),
theta(theta_),
acc_t(gsl_interp_accel_alloc()),
acc_r(gsl_interp_accel_alloc()),
acc_z(gsl_interp_accel_alloc()),
acc_phi(gsl_interp_accel_alloc())
{

}

inline LHDAnalyzedData_TSmesh_Handler::~LHDAnalyzedData_TSmesh_Handler(){
	gsl_interp_accel_free(acc_t);
	gsl_interp_accel_free(acc_r);
	gsl_interp_accel_free(acc_z);
	gsl_interp_accel_free(acc_phi);
}

inline double LHDAnalyzedData_TSmesh_Handler::get(const double t, const double r, const double z, const double phi, const double **** const data){
	const size_t t_index   = gsl_interp_accel_find(acc_t  ,   t_axis,   t_num,   t);
	const size_t r_index   = gsl_interp_accel_find(acc_r  ,   r_axis,   r_num,   r);
	const size_t z_index   = gsl_interp_accel_find(acc_z  ,   z_axis,   z_num,   z);
	const size_t phi_index = gsl_interp_accel_find(acc_phi, phi_axis, phi_num, phi);

	double dt_inv   = 1.0/(  t_axis[  t_index + 1] -   t_axis[  t_index]);
	double dr_inv   = 1.0/(  r_axis[  r_index + 1] -   r_axis[  r_index]);
	double dz_inv   = 1.0/(  z_axis[  z_index + 1] -   z_axis[  z_index]);
	double dphi_inv = 1.0/(phi_axis[phi_index + 1] - phi_axis[phi_index]);

	double ratio_floor_t   = (  t_axis[  t_index + 1] -   t) * dt_inv  ;
	double ratio_floor_r   = (  r_axis[  r_index + 1] -   r) * dr_inv  ;
	double ratio_floor_z   = (  z_axis[  z_index + 1] -   z) * dz_inv  ;
	double ratio_floor_phi = (phi_axis[phi_index + 1] - phi) * dphi_inv;

	double ratio_ceil_t   = -(  t_axis[  t_index    ] -   t) * dt_inv  ;
	double ratio_ceil_r   = -(  r_axis[  r_index    ] -   r) * dr_inv  ;
	double ratio_ceil_z   = -(  z_axis[  z_index    ] -   z) * dz_inv  ;
	double ratio_ceil_phi = -(phi_axis[phi_index    ] - phi) * dphi_inv;

/*
	std::cerr << "  t_index = " <<   t_index << std::endl;
	std::cerr << "  r_index = " <<   r_index << std::endl;
	std::cerr << "  z_index = " <<   z_index << std::endl;
	std::cerr << "phi_index = " << phi_index << std::endl;
*/
	return
		  ratio_floor_t * (ratio_floor_r * (ratio_floor_z * (ratio_floor_phi * data[t_index  ][r_index  ][z_index  ][phi_index  ]
		                                                   +  ratio_ceil_phi * data[t_index  ][r_index  ][z_index  ][phi_index+1])
		                                  +  ratio_ceil_z * (ratio_floor_phi * data[t_index  ][r_index  ][z_index+1][phi_index  ]
		                                                    + ratio_ceil_phi * data[t_index  ][r_index  ][z_index+1][phi_index+1]))
		                  + ratio_ceil_r * (ratio_floor_z * (ratio_floor_phi * data[t_index  ][r_index+1][z_index  ][phi_index  ]
		                                                   +  ratio_ceil_phi * data[t_index  ][r_index+1][z_index  ][phi_index+1])
		                                  +  ratio_ceil_z * (ratio_floor_phi * data[t_index  ][r_index+1][z_index+1][phi_index  ]
		                                                    + ratio_ceil_phi * data[t_index  ][r_index+1][z_index+1][phi_index+1])))
		+  ratio_ceil_t * (ratio_floor_r * (ratio_floor_z * (ratio_floor_phi * data[t_index+1][r_index  ][z_index  ][phi_index  ]
		                                                   +  ratio_ceil_phi * data[t_index+1][r_index  ][z_index  ][phi_index+1])
		                                  +  ratio_ceil_z * (ratio_floor_phi * data[t_index+1][r_index  ][z_index+1][phi_index  ]
		                                                    + ratio_ceil_phi * data[t_index+1][r_index  ][z_index+1][phi_index+1]))
		                  + ratio_ceil_r * (ratio_floor_z * (ratio_floor_phi * data[t_index+1][r_index+1][z_index  ][phi_index  ]
		                                                   +  ratio_ceil_phi * data[t_index+1][r_index+1][z_index  ][phi_index+1])
		                                  +  ratio_ceil_z * (ratio_floor_phi * data[t_index+1][r_index+1][z_index+1][phi_index  ]
		                                                    + ratio_ceil_phi * data[t_index+1][r_index+1][z_index+1][phi_index+1])));
}

inline double LHDAnalyzedData_TSmesh_Handler::get_reff(const double t, const double r, const double z, const double phi){
	return get(t, r, z, phi, reff);
}
inline double LHDAnalyzedData_TSmesh_Handler::get_theta(const double t, const double r, const double z, const double phi){
	return get(t, r, z, phi, theta);
}