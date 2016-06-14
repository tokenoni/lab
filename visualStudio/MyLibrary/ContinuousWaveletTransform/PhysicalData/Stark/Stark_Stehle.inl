namespace myPhysicalData{

	inline bool stark_Stehle::read(const std::string dir_, const enum Series series_, const size_t upperN_)
	{
		dir = dir_;
		series = series_;
		upperN = upperN_;
		std::stringstream ss;
		
		if(upperN > 30)
			return false;
		else if(upperN >= 10)
			ss << upperN;
		else ss << "0" << upperN;

		switch(series){
		case LymanSeries:
			dir += "ly\\ly" + ss.str() + "\\";
			lowerN = 1;
			break;
		case BalmerSeries:
			dir += "ba\\ba" + ss.str() + "\\";
			lowerN = 1;
			break;
		case PachenSeries:
			dir += "pa\\pa" + ss.str() + "\\";
			lowerN = 1;
			break;
		};

		data.resize(0);
		for(size_t i=1; 1; ++i){
			stark_Stehle_1Ne oneData;
			ss=std::stringstream("");
			ss << i;
			if(!oneData.read(dir + "index" + ss.str() + ".dat", dir + "profil" + ss.str() + ".dat"))
				break;
			else
				data.push_back(oneData);
		}

		for(size_t i=0; i<size_ne();++i)
			data[i].setInterp();

		return true;
	}
	inline size_t stark_Stehle::size_te(const size_t index)const{return data[index].size();}


	inline bool stark_Stehle_1Ne::read(const std::string index_filename_, const std::string profile_filename_){
		index_filename = index_filename_;
		profile_filename = profile_filename_;

		//---		index*.dat		---
		
		std::fstream ifs_index;
		ifs_index.open(index_filename, std::ifstream::in);
		if(!ifs_index.good()) return false;

		std::string sdata;
		getline(ifs_index, sdata);
		//	electronic density(1/cm**3) => converted to 1/m**3		
		getline(ifs_index, sdata);
		sdata = sdata.substr(sdata.find(")")+1);
		ne = atof(sdata.c_str())*1.0e6;
		//	conversion factor for wlengths
		getline(ifs_index, sdata);
		sdata = sdata.substr(sdata.find("=")+3);
		conversion_factor_for_wavelengths = atof(sdata.c_str());
		//	wings factor
		getline(ifs_index, sdata);
		sdata = sdata.substr(sdata.find("=")+1);
		wings_factor = atof(sdata.c_str());

		ifs_index.close();

		//---		profile*.dat		---
		std::fstream ifs_profile;
		ifs_profile.open(profile_filename, std::ifstream::in);
		if(!ifs_profile.good()) return false;
		data.resize(0);

		// header
		for(size_t i=0; i<5; ++i)
			getline(ifs_profile, sdata);
		
		do{
			double dalpha_tmp;
			std::vector<stark_Stehle_1Ne_1Te> data_tmp(3);
			std::vector<double> iStark_tmp(3), iDoppler_tmp(3);

			getline(ifs_profile, sdata);		//	T(K)
			if(sdata.size()>=20+9) data_tmp[0].te = atof(sdata.substr(20,10).c_str());
			if(sdata.size()>=45+9) data_tmp[1].te = atof(sdata.substr(45,10).c_str());
			if(sdata.size()>=70+9) data_tmp[2].te = atof(sdata.substr(70,10).c_str());

			getline(ifs_profile, sdata);		//	RO/Debye
			if(sdata.size()>=20+9) data_tmp[0].R0_Debye = atof(sdata.substr(20,10).c_str());
			if(sdata.size()>=45+9) data_tmp[1].R0_Debye = atof(sdata.substr(45,10).c_str());
			if(sdata.size()>=70+9) data_tmp[2].R0_Debye = atof(sdata.substr(70,10).c_str());

			getline(ifs_profile, sdata);		//	width Stark
			if(sdata.size()>=20+9) data_tmp[0].width_Stark = atof(sdata.substr(20,10).c_str());
			if(sdata.size()>=45+9) data_tmp[1].width_Stark = atof(sdata.substr(45,10).c_str());
			if(sdata.size()>=70+9) data_tmp[2].width_Stark = atof(sdata.substr(70,10).c_str());

			getline(ifs_profile, sdata);		//	width incl.Dop.
			if(sdata.size()>=20+9) data_tmp[0].width_incl_Dop = atof(sdata.substr(20,10).c_str());
			if(sdata.size()>=45+9) data_tmp[1].width_incl_Dop = atof(sdata.substr(45,10).c_str());
			if(sdata.size()>=70+9) data_tmp[2].width_incl_Dop = atof(sdata.substr(70,10).c_str());
			
			getline(ifs_profile, sdata);		//	blank

			getline(ifs_profile, sdata);		//	1st line
			dalpha_tmp = atof(sdata.substr(1,10).c_str());
			if(sdata.size()>=13+10) iDoppler_tmp[0] = atof(sdata.substr(13,10).c_str());
			if(sdata.size()>=38+10) iDoppler_tmp[1] = atof(sdata.substr(38,10).c_str());
			if(sdata.size()>=63+10) iDoppler_tmp[2] = atof(sdata.substr(63,10).c_str());
			
			if(sdata.size()>=25+10) iStark_tmp[0] = atof(sdata.substr(25,10).c_str());
			if(sdata.size()>=50+10) iStark_tmp[1] = atof(sdata.substr(50,10).c_str());
			if(sdata.size()>=75+10) iStark_tmp[2] = atof(sdata.substr(75,10).c_str());
			
			data_tmp[0].Dalpha.push_back(dalpha_tmp);
			data_tmp[1].Dalpha.push_back(dalpha_tmp);
			data_tmp[2].Dalpha.push_back(dalpha_tmp);
			
			data_tmp[0].iDoppler.push_back(iDoppler_tmp[0]);
			data_tmp[1].iDoppler.push_back(iDoppler_tmp[1]);
			data_tmp[2].iDoppler.push_back(iDoppler_tmp[2]);

			data_tmp[0].iStark.push_back(iStark_tmp[0]);
			data_tmp[1].iStark.push_back(iStark_tmp[1]);
			data_tmp[2].iStark.push_back(iStark_tmp[2]);

			getline(ifs_profile, sdata);		//	2nd line
			while(sdata.size() > 11){
				dalpha_tmp = atof(sdata.substr(1,10).c_str());
				if(sdata.size()>=13+10) iDoppler_tmp[0] = atof(sdata.substr(13,10).c_str());
				if(sdata.size()>=38+10) iDoppler_tmp[1] = atof(sdata.substr(38,10).c_str());
				if(sdata.size()>=63+10) iDoppler_tmp[2] = atof(sdata.substr(63,10).c_str());
				
				if(sdata.size()>=25+10) iStark_tmp[0] = atof(sdata.substr(25,10).c_str());
				if(sdata.size()>=50+10) iStark_tmp[1] = atof(sdata.substr(50,10).c_str());
				if(sdata.size()>=75+10) iStark_tmp[2] = atof(sdata.substr(75,10).c_str());

				data_tmp[0].Dalpha.push_back(dalpha_tmp);
				data_tmp[1].Dalpha.push_back(dalpha_tmp);
				data_tmp[2].Dalpha.push_back(dalpha_tmp);

				data_tmp[0].iDoppler.push_back((iDoppler_tmp[0]>0.0) ? iDoppler_tmp[0] : data_tmp[0].iDoppler[data_tmp[0].iDoppler.size()-1]);
				data_tmp[1].iDoppler.push_back((iDoppler_tmp[1]>0.0) ? iDoppler_tmp[1] : data_tmp[1].iDoppler[data_tmp[1].iDoppler.size()-1]);
				data_tmp[2].iDoppler.push_back((iDoppler_tmp[2]>0.0) ? iDoppler_tmp[2] : data_tmp[2].iDoppler[data_tmp[2].iDoppler.size()-1]);

				data_tmp[0].iStark.push_back((iStark_tmp[0]>0.0) ? iStark_tmp[0] : data_tmp[0].iStark[data_tmp[0].iStark.size()-1]);
				data_tmp[1].iStark.push_back((iStark_tmp[1]>0.0) ? iStark_tmp[1] : data_tmp[1].iStark[data_tmp[1].iStark.size()-1]);
				data_tmp[2].iStark.push_back((iStark_tmp[2]>0.0) ? iStark_tmp[2] : data_tmp[2].iStark[data_tmp[2].iStark.size()-1]);

				//	for next line
				if(ifs_profile.eof()) break;
				getline(ifs_profile, sdata);
			}

			if(data_tmp[0].te > 0.0)
				data.push_back(data_tmp[0]);
			if(data_tmp[1].te > 0.0)
				data.push_back(data_tmp[1]);
			if(data_tmp[2].te > 0.0)
				data.push_back(data_tmp[2]);

		}while(!ifs_profile.eof());
		return true;
	}

	inline void stark_Stehle_1Ne::setInterp(){
		for(size_t i=0; i<size();++i)
			data[i].setInterp();
	}
	
	inline void stark_Stehle_1Ne_1Te::setInterp(){
		iStark_interp.set(  Dalpha,   iStark, mygsl::interp_logx::Linear);
		iDoppler_interp.set(Dalpha, iDoppler, mygsl::interp_logx::Linear);
	}
};
