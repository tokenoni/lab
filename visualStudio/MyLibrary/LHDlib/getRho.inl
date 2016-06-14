namespace LHDlib{
	inline RhoTable::RhoTable(void){
		//	variables for multimin 
		//	which is used to derive (phi, theta) from (R, Z)
		rho_theta = gsl_vector_alloc(2);
		my_func.n =2;
		my_func.f = &(RhoTable::Psi);
		my_func.df = &(RhoTable::dPsi);
		my_func.fdf = &(RhoTable::Psi_dPsi);
		my_func.params = &rmnZmnSet;
		T = gsl_multimin_fdfminimizer_conjugate_fr;
		s = gsl_multimin_fdfminimizer_alloc(T, 2);
		rmnZmnSet.size=0;
	}

	inline RhoTable::~RhoTable(void){
		gsl_multimin_fdfminimizer_free(s);
		gsl_vector_free(rho_theta);
	}
	inline void RhoTable::copy(const RhoTable& obj){
		rmnZmnSet = obj.rmnZmnSet;
		m_size = obj.m_size;
		rho_size = obj.rho_size;
		Rmn = obj.Rmn;
		Zmn = obj.Rmn;
		m_list = obj.m_list;
		n_list = obj.n_list;
		rho_list = obj.rho_list;
		dV_drho = obj.dV_drho;
		Vp = obj.Vp;
		e_helical = obj.e_helical;
		e_toroidal = obj.e_helical; 
		q = obj.q;
		my_func.params = &rmnZmnSet;

		rho_vs_R2 = obj.rho_vs_R2;
		Rax_avg = obj.Rax_avg;

	}
	inline RhoTable& RhoTable::operator=(const RhoTable& obj){
		copy(obj);
		return *this;
	}
	inline bool RhoTable::read(const std::string filename){
		if (!readFLXfile(filename))
			return false;

		std::string ergoFilename = filename.substr(0, filename.find_last_of('b'));
		readErgofile(ergoFilename + "b000_RmnErgodic.itx", ergoFilename + "b000_ZmnErgodic.itx");
		return prepare_interp();
	}
	
	inline bool RhoTable::readFLXfile(const std::string filename){
			std::ifstream ifs(filename.c_str(), std::ifstream::in);
			//	---	loop for weiting until file open ---
		if(!ifs.is_open()){
			std::cerr << "!!!   file: " + filename << " cannot be opened" << std::endl;
			std::cerr << "!!!   ---  press return for retry ---		 !!!" << std::endl;
			getchar();
			ifs.close();ifs.clear();
			return false;
		}


		//---	reading	---
		std::string sdata;

		sdata = filename.substr(filename.find("lhd-r")+5,3); 
		Rax_avg = atof(sdata.c_str())*1.0e-2;

		//	allocation 
		Rmn.clear();	Zmn.clear();
		m_list.clear();	n_list.clear();
		rho_list.clear();
		dV_drho.clear();
		e_helical.clear();
		e_toroidal.clear();
		q.clear();
		while(!ifs.eof()){
			getline(ifs,sdata);

			//	reading the header	
			if( sdata.size() > 0 && sdata[0] != '#'){
				if(sdata.find("rmaj0 = ") != std::string::npos){
					sdata = sdata.substr(sdata.find("rmaj0 = ")+8);
					sdata = sdata.substr(0, sdata.find(","));
					Rax_avg = atof(sdata.c_str());
				}

				//  for one section
				if(sdata.find("jr") != std::string::npos){
					std::vector<double> rmn_tmp(0), zmn_tmp(0);
					std::vector<int> m_tmp(0), n_tmp(0);

					getline(ifs,sdata);
					//--- for each rho ---
					while(sdata.size() > 0 && sdata.find("m=") != std::string::npos){
						std::string sdata_part=sdata, sub_str;
						// m=
						sdata_part = sdata_part.substr(sdata_part.find("m=")+2);
						sub_str = sdata_part.substr(0, sdata_part.find("n=")-1);
						m_tmp.push_back((int)atof(sub_str.c_str()));
						// n=
						sdata_part = sdata_part.substr(sdata_part.find("n=")+2);
						sub_str = sdata_part.substr(0, sdata_part.find("rmn="));
						n_tmp.push_back((int)atof(sub_str.c_str()));
						// rmn = 
						sdata_part = sdata_part.substr(sdata_part.find("rmn=")+4);
						sub_str = sdata_part.substr(0, sdata_part.find("zmn="));
						rmn_tmp.push_back(atof(sub_str.c_str()));
						// zmn = 
						sdata_part = sdata_part.substr(sdata_part.find("zmn=")+4);
						sub_str = sdata_part;
						zmn_tmp.push_back(atof(sub_str.c_str()));

						getline(ifs,sdata);
					}
					m_list.push_back(m_tmp);
					n_list.push_back(n_tmp);
					Rmn.push_back(rmn_tmp);
					Zmn.push_back(zmn_tmp);

					//--- for rho ---
					while(sdata.find("rho=") == std::string::npos)	getline(ifs,sdata);
					std::string sub_str = sdata.substr(sdata.find("rho=")+4, sdata.find("1/q")-sdata.find("rho=")-4);
					rho_list.push_back(atof(sub_str.c_str()));

					//--- for dV_drho ---
					while (sdata.find("1/q=") == std::string::npos)	getline(ifs, sdata);
					std::string sub_str2 = sdata.substr(sdata.find("1/q=") + 4);
					q.push_back(1.0/atof(sub_str2.c_str()));

					//--- for dV_drho ---
					while (sdata.find("Vp=") == std::string::npos)	getline(ifs, sdata);
					sub_str2 = sdata.substr(sdata.find("Vp=") + 3);
					dV_drho.push_back(atof(sub_str2.c_str()));
					Vp.push_back(8.0*myconst::pi*myconst::pi * dV_drho[dV_drho.size() - 1] * rho_list[rho_list.size() - 1]);

					//--- for eh, et ---
					while (sdata.find("et=") == std::string::npos) getline(ifs, sdata);
					sub_str2 = sdata.substr(sdata.find("et=") + 3);
					e_toroidal.push_back(atof(sub_str2.c_str()));
					while (sdata.find("eh=") == std::string::npos) getline(ifs, sdata);
					sub_str2 = sdata.substr(sdata.find("eh=") + 3);
					e_helical.push_back(atof(sub_str2.c_str()));
				}
			}
		}
		rmnZmnSet.e_helical.set(rho_list, e_helical, mygsl::interp::Linear);
		rmnZmnSet.e_toroidal.set(rho_list, e_helical, mygsl::interp::Linear);
		rmnZmnSet.q.set(rho_list, q, mygsl::interp::Linear);
		return true;
	}

	inline bool RhoTable::readErgofile(const std::string ergoRmnFilename, const std::string ergoZmnFilename){
		std::ifstream ifsR(ergoRmnFilename), ifsZ(ergoZmnFilename);
		if (!ifsR.is_open() || !ifsZ.is_open())	
			return false;
		std::vector<double> rmnErgodic = IGORdata::itx_1stdv(ergoRmnFilename);
		std::vector<double> zmnErgodic = IGORdata::itx_1stdv(ergoZmnFilename);
		
		m_size = m_list[0].size();
		std::vector<int> m_tmp = m_list[0];
		std::vector<int> n_tmp = n_list[0];
		if (rmnErgodic.size() != m_size || zmnErgodic.size() != m_size)
			return false;
		
		double VpErgo = getVpOne(rmnErgodic, zmnErgodic, m_tmp, n_tmp);
		double VpLCFS = getVpOne(Rmn[m_list.size()-1], Zmn[m_list.size()-1], m_tmp, n_tmp);
		double rhoLCFS = rho_list[rho_list.size() - 1];

		mygsl::interp rhovsVp(m_list.size(), mygsl::interp::Linear);
		for (size_t i = 0; i < m_list.size(); ++i){
			double Vptmp = getVpOne(Rmn[i], Zmn[i], m_tmp, n_tmp);
			rhovsVp.set(i, Vptmp, rho_list[i]);
		}
		rhovsVp.setInterp();
		double rhoErgo = rhovsVp.get(VpErgo);

		//		double rhoErgo = sqrt(VpErgo / VpLCFS);

/*		double a = (Vp[Vp.size() - 1] - Vp[Vp.size() - 2]) / (rho_list[rho_list.size() - 1] - rho_list[rho_list.size() - 2]);
		double b = Vp[Vp.size() - 1] - a * rhoLCFS;
		double c = VpLCFS - VpErgo - Vp[Vp.size() - 1]*rhoLCFS;
		double rhoErgo = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
*/
		



//		double dV_drho_ergo = ((VpErgo - VpLCFS) / (4.0 / 3.0*myconst::pi*myconst::pi) - dV_drho_LCFS*(rhoErgo*rhoErgo + rhoErgo*rhoLCFS - 2.0*rhoLCFS*rhoLCFS)) 
//			/ (rhoLCFS*rhoLCFS + rhoErgo*rhoLCFS - 2.0*rhoErgo*rhoErgo);
		double dV_drho_LCFS = dV_drho[dV_drho.size() - 1];
		double dVpdrhoLCFS = Vp[Vp.size() - 1];
		double dVpdrhoErgo = (VpErgo - VpLCFS + rhoLCFS*dVpdrhoLCFS) / rhoErgo;

		Vp.push_back(dVpdrhoErgo);
		dV_drho.push_back(dVpdrhoErgo / (8.0*myconst::pi*myconst::pi * rhoErgo));

		m_list.push_back(m_tmp);
		n_list.push_back(n_tmp);
		Rmn.push_back(rmnErgodic);
		Zmn.push_back(zmnErgodic);
		rho_list.push_back(rhoErgo);
		return true;
	}

	inline bool RhoTable::prepare_interp(){
	//	---	check if all the sections have same number of data set
		rho_size = m_list.size();
		if(rho_size == 0) return false;
		if(rho_size != n_list.size()) return false;
		if(rho_size != Rmn.size()) return false;
		if(rho_size != Zmn.size()) return false;
		if(rho_size != rho_list.size()) return false;
		if(rho_size != dV_drho.size()) return false;
		if(rho_size != Vp.size()) return false;

		m_size = m_list[0].size();
		for (size_t i = 0; i<m_list.size(); ++i) {
			if(m_list[i].size() != m_size) return false;
			if(n_list[i].size() != m_size) return false;
			if(Rmn[i].size() != m_size) return false;
			if(Zmn[i].size() != m_size) return false;
		}

		//	stack into std::vector<mygsl::interp>
		rmnZmnSet.Rmn_interp.setSize(m_size, rho_size, mygsl::interp::Linear);
		rmnZmnSet.Zmn_interp.setSize(m_size, rho_size, mygsl::interp::Linear);
		rmnZmnSet.m.setSize(m_size, mygsl::interp::Linear);
		rmnZmnSet.n.setSize(m_size, mygsl::interp::Linear);
		for(size_t i=0; i<m_size; ++i){
			for(size_t j=0; j<rho_size; ++j){
				rmnZmnSet.Rmn_interp.set(i, j, rho_list[j], Rmn[j][i]);
				rmnZmnSet.Zmn_interp.set(i, j, rho_list[j], Zmn[j][i]);
			}
			rmnZmnSet.m[i] = m_list[0][i];
			rmnZmnSet.n[i] = n_list[0][i];
		}
		rmnZmnSet.Rmn_interp.setInterp();
		rmnZmnSet.Zmn_interp.setInterp();
		rmnZmnSet.size = m_list[0].size();

		rmnZmnSet.dVdrho.set(rho_list,dV_drho, mygsl::interp::Linear);
		rmnZmnSet.Vp.set(    rho_list,Vp, mygsl::interp::Linear);

		std::vector<double> reff_tmp(rho_size);
		for (size_t j = 0; j < rho_size; ++j){
			reff_tmp[j] = sqrt(get_Vp(0.0, rho_list[j]) / (2.0*myconst::pi*myconst::pi*Rax_avg));
		}
		rmnZmnSet.rhoVsReff.set(reff_tmp, rho_list, mygsl::interp::Linear);
		rmnZmnSet.reffVsRho.set(rho_list, reff_tmp, mygsl::interp::Linear);
		return true;
	}

	inline void RhoTable::get_rho_drho_xyz(const double *x, double* rho, double* theta, double *drho_xyz){
		get_rho_drho_xyz(x[0], x[1], x[2], rho, theta, drho_xyz);
	}
	inline void RhoTable::get_rho_drho_xyz(const double x, const double y, const double z, double* rho, double* theta, double* drho_xyz){
		double r = sqrt(x*x + y*y);
		double phi = atan2(y, x) + myconst::pi/10.0;
		double drho_rzphi[3];
		
		get_rho_drho(r, z, phi, rho, theta, drho_rzphi);
		drho_xyz[0] = drho_rzphi[0]*x/r + drho_rzphi[2]*y/r/r; 
		drho_xyz[1] = drho_rzphi[0]*y/r - drho_rzphi[2]*x/r/r; 
		drho_xyz[2] = drho_rzphi[1];
	}

	inline void RhoTable::get_rho_drho(const double *r, double* rho, double* theta, double *drho){
		get_rho_drho(r[0], r[1],r[2], rho, theta, drho);
	}
	inline void RhoTable::get_rho_drho(const double r, const double z, 
								const double phi, double* rho, double* theta, double* drho)
	{
		rmnZmnSet.R0 = r;
		rmnZmnSet.Z0 = z;
		rmnZmnSet.phi= phi;

		double Raxis=0.0, Zaxis=0.0;
		for(size_t i=0; i<rmnZmnSet.size; ++i){
			Raxis += rmnZmnSet.Rmn_interp.get(i, 0.0) * cos(rmnZmnSet.n[i]*phi);
			Zaxis += rmnZmnSet.Zmn_interp.get(i, 0.0) * sin(rmnZmnSet.n[i]*phi);
		}

		// initialize
		double rho_init = sqrt((r-Raxis)*(r-Raxis) + (z-Zaxis)*(z-Zaxis))/0.6;
		gsl_vector_set(rho_theta, 0, rho_init);
		gsl_vector_set(rho_theta, 1, -atan2(-(z-Zaxis),r-Raxis));

		size_t iter=0; int status;
		double err = rho_init / 1e3;
		if(err < 1.0e-4) err = 1.0e-4;
		gsl_multimin_fdfminimizer_set (s, &my_func, rho_theta, err, 1e-6);
		do{
			iter ++;
			status = gsl_multimin_fdfminimizer_iterate(s);
			if(status) break;
			status = gsl_multimin_test_gradient(s->gradient, 1e-4);

			//---	constran : rho > 0 ---
			if(gsl_vector_get(s->x,0) < 0.0){
				gsl_vector_set(s->x,0, -gsl_vector_get(s->x, 0));
				gsl_vector_set(s->x,1,  gsl_vector_get(s->x, 1)+myconst::pi);
			}
		}while (status == GSL_CONTINUE && iter < 100);

		*rho = gsl_vector_get(s->x,0)*gsl_vector_get(s->x,0);
		*theta = gsl_vector_get(s->x,1);
		get_drho(*rho, *theta, phi, drho);
	}

	inline void RhoTable::get_drho(const double rho, const double theta, const double phi, double* drho){
		size_t size = rmnZmnSet.size;
		double R=0.0, Z=0.0;
		double dR_drho=0.0, dZ_drho=0.0, dR_dtheta=0.0, dZ_dtheta=0.0, dR_dphi=0.0, dZ_dphi=0.0;
		
		for(size_t i=0; i<size; ++i){
			double cos_mtheta_nphi = cos(rmnZmnSet.m[i]*theta + rmnZmnSet.n[i]*rmnZmnSet.phi);
			double sin_mtheta_nphi = sin(rmnZmnSet.m[i]*theta + rmnZmnSet.n[i]*rmnZmnSet.phi);
			double rmn_val = rmnZmnSet.Rmn_interp.get(i, rho);
			double zmn_val = rmnZmnSet.Zmn_interp.get(i, rho);
			double dRmn_drho = rmnZmnSet.Rmn_interp.getderiv(i, rho) ;
			double dZmn_drho = rmnZmnSet.Zmn_interp.getderiv(i, rho) ;

			R += rmn_val * cos_mtheta_nphi;
			Z += zmn_val * sin_mtheta_nphi;
			dR_drho   += dRmn_drho * cos_mtheta_nphi;
			dZ_drho   += dZmn_drho * sin_mtheta_nphi;
			dR_dtheta -= rmn_val * sin_mtheta_nphi * rmnZmnSet.m[i];
			dZ_dtheta += zmn_val * cos_mtheta_nphi * rmnZmnSet.m[i];
			dR_dphi   -= rmn_val * sin_mtheta_nphi * rmnZmnSet.n[i];
			dZ_dphi   += zmn_val * cos_mtheta_nphi * rmnZmnSet.n[i];
		}
		double determinant = dR_drho*dZ_dtheta - dZ_drho*dR_dtheta;
		
		// drho_dR
//		drho[0] = dZ_dtheta/determinant;
		drho[0] = 1.0/(dR_drho - dZ_drho*(dR_dtheta / dZ_dtheta));
		// drho_dZ
//		drho[1] = -dR_dtheta/determinant;
		drho[1] = 1.0/(dZ_drho - dR_drho*(dZ_dtheta/dR_dtheta));
		// drho_dphi
		drho[2] = (-dZ_dtheta*dR_dphi -dR_dtheta*dZ_dphi)/determinant;
	}

	inline void RhoTable::getRZax(const double phi, double &R, double &Z){
		rmnZmnSet.getRZ(0.0, 0.0, phi, R, Z);
		return;
	}
	inline void RhoTable::getRZfromRhoThetaPhi(const double rho, const double theta, const double phi, double& R, double&Z){
		rmnZmnSet.getRZ(rho, theta, phi, R, Z);
		return;
	}

	inline double RhoTable::getVpCheck(const double rho_, const size_t num_theta, const size_t num_phi){
		std::vector<double> R(rmnZmnSet.size), Z(rmnZmnSet.size);
		std::vector<int> m(rmnZmnSet.size), n(rmnZmnSet.size);
		for (size_t i = 0; i < rmnZmnSet.size; ++i){
			R[i] = rmnZmnSet.Rmn_interp.get(i, rho_);
			Z[i] = rmnZmnSet.Zmn_interp.get(i, rho_);
			m[i] = (int)floor(rmnZmnSet.m[i] + 0.5);
			n[i] = (int)floor(rmnZmnSet.n[i] + 0.5);
		}
		return getVpOne(R, Z, m, n);
	}

	inline double RmnZmnSet::getVp(const double rho_start, const double rho_end){
		size_t m_size = size;
		std::vector<double> Rmn(m_size), Zmn(m_size);
		std::vector<int>  m_tmp(m_size), n_tmp(m_size);
		for (size_t i = 0; i < m_size; ++i){
			m_tmp[i] = int(floor(m[i] + 0.5));
			n_tmp[i] = int(floor(n[i] + 0.5));
			Rmn[i] = Rmn_interp[i].get(rho_start);
			Zmn[i] = Zmn_interp[i].get(rho_start);
		}
		double Vp_start = RhoTable::getVpOne(Rmn, Zmn, m_tmp, n_tmp);
		for (size_t i = 0; i < m_size; ++i){
			Rmn[i] = Rmn_interp[i].get(rho_end);
			Zmn[i] = Zmn_interp[i].get(rho_end);
		}
		double Vp_end = RhoTable::getVpOne(Rmn, Zmn, m_tmp, n_tmp);
		return Vp_end - Vp_start;
	}
	inline double RhoTable::get_Vp(const double rho_start, const double rho_end){
		std::vector<double> Rmn(m_size), Zmn(m_size);
		std::vector<int>  m_tmp(m_size), n_tmp(m_size);
		for (size_t i = 0; i < m_size; ++i){
			m_tmp[i] = int(floor(rmnZmnSet.m[i] + 0.5));
			n_tmp[i] = int(floor(rmnZmnSet.n[i] + 0.5));
			Rmn[i] = rmnZmnSet.Rmn_interp[i].get(rho_start);
			Zmn[i] = rmnZmnSet.Zmn_interp[i].get(rho_start);
		}
		double Vp_start = RhoTable::getVpOne(Rmn, Zmn, m_tmp, n_tmp);
		for (size_t i = 0; i < m_size; ++i){
			Rmn[i] = rmnZmnSet.Rmn_interp[i].get(rho_end);
			Zmn[i] = rmnZmnSet.Zmn_interp[i].get(rho_end);
		}
		double Vp_end = RhoTable::getVpOne(Rmn, Zmn, m_tmp, n_tmp);
		return Vp_end - Vp_start;
	}

	inline double RhoTable::getVpOne(const std::vector<double>& R, const std::vector<double>&Z, const std::vector<int>& m_list, const std::vector<int>& n_list){
		double result = 0.0;
		size_t num = n_list.size();
		for (size_t i = 0; i < num; ++i){
			int m = m_list[i];
			int n = n_list[i];
			double rmn = R[i];
			for (size_t j = 0; j < num; ++j){
				int m_ = m_list[j];
				int n_ = n_list[j];
				double zmn_ = Z[j];
				for (size_t k = 0; k < num; ++k){
					int m__ = m_list[k];
					int n__ = n_list[k];
					double rmn__ = R[k];
					//	m" = m + m'  &&  n" = n + n'
					if (m__ == m + m_ && n__ == n + n_)	result += rmn * zmn_ * m__*rmn__;
					//	m" =-m + m'  &&  n" =-n + n'
					else if (m__ == m_ - m && n__ == n_ - n)	result += rmn * zmn_ * m__*rmn__;
					//	m" = m - m'  &&  n" = n - n'
					else if (m__ == m - m_ && n__ == n - n_)	result -= rmn * zmn_ * m__*rmn__;
				}
			}
		}
		return result * myconst::pi*myconst::pi * 2.0;
	}

};