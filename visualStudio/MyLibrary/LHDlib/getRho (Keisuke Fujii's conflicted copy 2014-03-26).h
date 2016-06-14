#ifndef __LHDLIB_getRho_HPP__
#define __LHDLIB_getRho_HPP__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <MyGsl\interp.hpp>
#include <MyGsl\interp_vector2.hpp>
#include <MyGsl\constants.hpp>
#include <gsl\gsl_roots.h>
#include <gsl\gsl_multimin.h>
#include <IGOR\IGORitx.hpp>


namespace LHDlib{
	class RmnZmnSet{
	public:
		RmnZmnSet(void){};
		~RmnZmnSet(void){};
		mygsl::interp_vector2 Rmn_interp, Zmn_interp;
		mygsl::interp m, n;
		mygsl::interp dVdrho;
		mygsl::interp Vp;
		mygsl::interp rhoVsReff, reffVsRho;
		mygsl::interp e_helical, e_toroidal, q;
		void getRZ(const double rho, const double theta, const double phi, double& R, double& Z){
			R=0.0; Z=0.0;
			for(size_t i=0; i<size; ++i){
				R += Rmn_interp.get(i, rho) * cos(m[i]*theta + n[i]*phi);
				Z += Zmn_interp.get(i, rho) * sin(m[i]*theta + n[i]*phi);
			}
		};
		void getRZderive(const double rho, const double theta, const double phi, double&R, double&Z, double&dRdtheta, double&dZdtheta){
			R = 0.0, Z = 0.0, dRdtheta = 0.0, dZdtheta = 0.0;
			for (size_t i = 0; i<size; ++i){
				double rmn = Rmn_interp.get(i, rho);
				double zmn = Zmn_interp.get(i, rho);
				double cos_mtheta_nphi = cos(m[i] * theta + n[i] * phi);
				double sin_mtheta_nphi = sin(m[i] * theta + n[i] * phi);

				R += rmn * cos_mtheta_nphi;
				Z += zmn * sin_mtheta_nphi;
				dRdtheta += -rmn * m[i] * sin_mtheta_nphi;
				dZdtheta +=  zmn * m[i] * cos_mtheta_nphi;
			}
		};
		double getdRdPhi(const double rho, const double theta, const double phi){
			double dRdphi = 0.0;
			for (size_t i = 0; i < size; ++i){
				double rmn = Rmn_interp.get(i, rho);
				double sin_mtheta_nphi = sin(m[i] * theta + n[i] * phi);
				dRdphi += -rmn * n[i] * sin_mtheta_nphi;
			}
			return dRdphi;
		}

		double R0, Z0, phi;
		size_t size;
		double getVp(const double rho_start, const double rho_end);
		double getE_helical(const double rho){ return e_helical.get(rho); }
		double getE_toroidal(const double rho){ return e_toroidal.get(rho); }
	};

	class RhoTable
	{
	public:
		RhoTable(void);
		~RhoTable(void);
		bool read(const std::string filename);
		double getrho(const double *r);
		double getrho(const double r, const double z, const double phi);
		double getreff(const double *r){	return getReffFromRho(getrho(r));}
		double getreff(const double r, const double z, const double phi){	return getReffFromRho(getrho(r, z, phi));		};

		double getRhoFromReff(const double reff){ return rmnZmnSet.rhoVsReff.get(reff); }
		double getReffFromRho(const double rho_){ return rmnZmnSet.reffVsRho.get(rho_); }
		
		void getRZfromRhoThetaPhi(const double rho, const double theta, const double phi, double& R, double&Z);

		void get_rho_drho(const double *r, double* rho, double* theta, double *drho_dr);
		void get_rho_drho(const double r, const double z, const double phi, double* rho, double* theta, double* drho);
		void get_rho_drho_xyz(const double *x, double* rho, double* theta, double *drho_dxyz);
		void get_rho_drho_xyz(const double x, const double y, const double z, double* rho, double* theta, double* drho_xyz);
		void get_drho(const double rho, const double theta, const double phi, double* drho_dr);
		RmnZmnSet getRmnZmnSet()const{return rmnZmnSet;}

		double get_dVdrho( const double rho_){ return rmnZmnSet.dVdrho.get(rho_); }
		double get_dVpdrho(const double rho_){ return rmnZmnSet.Vp.get(rho_); }
		double get_Vp(const double rho_start, const double rho_end);
		double get_Vp_fromReff(const double reff_start, const double reff_end)
		{	return (reff_end*reff_end - reff_start*reff_start)*2.0*myconst::pi*myconst::pi*Rax_avg;		}
//		double get_reff(const double rho_){return sqrt(get_Vp(0.0, rho_)/(2.0*myconst::pi*myconst::pi*Rax_avg));}
		void getRZax(const double phi, double &R, double &Z);

		//	return the plasma volume area of the magnetic flux surface
		//	num_phi and num_theta are the number of the discretized point of the poloidal and toroidal directions
		double getVpCheck(const double rho_, const size_t num_theta = 200, const size_t num_phi = 200);
		static double getVpOne(const std::vector<double>& R, const std::vector<double>&Z, const std::vector<int>& m_list, const std::vector<int>& n_list);

		double getSafetyFactor(const double rho){ return rmnZmnSet.q.get(rho); }
		double getE_helical(const double rho){ return rmnZmnSet.getE_helical(rho); }
		double getE_toroidal(const double rho){ return rmnZmnSet.getE_helical(rho); }

	private:
		bool readFLXfile(const std::string filename);
		bool readErgofile(const std::string ergoRmnFilename, const std::string ergoZmnFilename);
		RmnZmnSet rmnZmnSet;
		bool prepare_interp();
		size_t m_size, rho_size;
		std::vector< std::vector<double> > Rmn, Zmn;
		std::vector< std::vector<int> > m_list, n_list;
		std::vector< double > rho_list;
		std::vector< double > dV_drho, Vp, e_helical, e_toroidal, q;

		mygsl::interp rho_vs_R2;
		double Rax_avg;

		//---	for multimin	---
	public:
		static double Psi(const gsl_vector* x, void* p);
		static void dPsi(const gsl_vector* x, void* p, gsl_vector *g);
		static void Psi_dPsi(const gsl_vector* x, void* p, double* f, gsl_vector *g);
	private:
		const gsl_multimin_fdfminimizer_type *T;
		gsl_multimin_fdfminimizer *s;
		gsl_vector* rho_theta;
		gsl_multimin_function_fdf my_func;
	};

	//---	constraint : rho > 0 should be taken into account   ---
	inline double RhoTable::Psi(const gsl_vector* x, void* p){
		RmnZmnSet* rzset = (RmnZmnSet*) p;
		double rho =  x->data[0]*x->data[0];//gsl_vector_get(x, 0);
		double theta =x->data[1];// gsl_vector_get(x, 1);
		size_t size = rzset->size;
		double R=0.0, Z=0.0;
		for(size_t i=0; i<size; ++i){
			R += rzset->Rmn_interp.get(i, rho) * cos(rzset->m[i]*theta + rzset->n[i]*rzset->phi);
			Z += rzset->Zmn_interp.get(i, rho) * sin(rzset->m[i]*theta + rzset->n[i]*rzset->phi);
		}
		return (rzset->R0 - R)*(rzset->R0 - R) + (rzset->Z0 - Z)*(rzset->Z0 - Z);
	}
	//---	constraint : rho > 0 should be taken into account   ---
	inline void RhoTable::dPsi(const gsl_vector* x, void* p, gsl_vector *g){
		RmnZmnSet* rzset = (RmnZmnSet*) p;
		double rho =  x->data[0]*x->data[0];//gsl_vector_get(x, 0);
		double sqrt_rho =  abs(x->data[0]);//gsl_vector_get(x, 0);
		double theta =x->data[1];// gsl_vector_get(x, 1);
		size_t size = rzset->size;
		double R=0.0, Z=0.0;
		double dPsi_drho_R=0.0, dPsi_drho_Z=0.0;
		double dPsi_dtheta_R=0.0, dPsi_dtheta_Z=0.0;
		for(size_t i=0; i<size; ++i){
			double cos_mtheta_nphi = cos(rzset->m[i]*theta + rzset->n[i]*rzset->phi);
			double sin_mtheta_nphi = sin(rzset->m[i]*theta + rzset->n[i]*rzset->phi);
			double rmn_val = rzset->Rmn_interp.get(i, rho);
			double zmn_val = rzset->Zmn_interp.get(i, rho);

			R += rmn_val * cos_mtheta_nphi;
			Z += zmn_val * sin_mtheta_nphi;
			dPsi_drho_R += rzset->Rmn_interp.getderiv(i, rho) * cos_mtheta_nphi;
			dPsi_drho_Z += rzset->Zmn_interp.getderiv(i, rho) * sin_mtheta_nphi;
			dPsi_dtheta_R -= rmn_val * sin_mtheta_nphi * rzset->m[i];
			dPsi_dtheta_Z += zmn_val * cos_mtheta_nphi * rzset->m[i];
		}
		dPsi_drho_R *= 2.0*sqrt_rho;
		dPsi_drho_Z *= 2.0*sqrt_rho;

		gsl_vector_set(g, 0, 2.0*(R - rzset->R0) * dPsi_drho_R   + 2.0*(Z - rzset->Z0) * dPsi_drho_Z);
		gsl_vector_set(g, 1, 2.0*(R - rzset->R0) * dPsi_dtheta_R + 2.0*(Z - rzset->Z0) * dPsi_dtheta_Z);
	}
	inline void RhoTable::Psi_dPsi(const gsl_vector* x, void* p, double* f, gsl_vector *g){
		RhoTable::dPsi(x, p, g);
		*f = RhoTable::Psi(x, p);
	}

};
#include "getRho.inl"
#endif