#ifndef __STARK_STEHLE_HPP__
#define __STARK_STEHLE_HPP__

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <MyGsl\log_interp.hpp>

namespace myPhysicalData{
	class stark_Stehle_1Ne;
	class stark_Stehle_1Ne_1Te;

	class stark_Stehle{
	public:
		enum Series{
			LymanSeries,
			BalmerSeries,
			PachenSeries
		};
		bool read(const std::string dir_, const enum Series series_, const size_t upperN_);
		
		size_t size_ne()const{return data.size();}
		size_t size_te(const size_t index)const;

		void setNeTe(const double Ne, const double Te);
		double getShape(const double dlambda_over_lambda);
	private:
		//	physical variables
		enum Series series;
		size_t upperN, lowerN;
		std::vector<stark_Stehle_1Ne> data;
		std::string dir;
		mygsl::interp_logxy iStark_tmp, iDoppler_tmp;
	};

	class stark_Stehle_1Ne{
	public:
		bool read(const std::string index_filename_, const std::string profile_filename_);
		size_t size()const {return data.size();}
		void setInterp();
	private:
		//	physical variables
		double ne;		//	in 1/m**3
		double conversion_factor_for_wavelengths;
		double wings_factor;
		//		
		std::string index_filename, profile_filename;
		std::vector<stark_Stehle_1Ne_1Te> data;
	};

	class stark_Stehle_1Ne_1Te{
	public:
	
		double te;
		double R0_Debye;
		double width_Stark, width_incl_Dop;

		//	stark data
		std::vector<double> Dalpha;		
			//    Detuning from central wavelength, expressed in units
            //    of {alpha} = lambda0 / E0 where
            //    lambda0 is the unperturbed wavelength (in {AA} = 0.1nm) and
            //    E0     is the normal Holstmark electric field in     uesVolts/cm

		std::vector<double> iDoppler;	//	Intensity with Doppler (convolved) in units of 1/{alpha}
		std::vector<double> iStark;		//	Intensity with Stark only (in brackets) in units of 1/{alpha}
		mygsl::interp_logxy iStark_interp, iDoppler_interp;
		void setInterp();
	};
};



#include "Stark_Stehle.inl"
#endif
