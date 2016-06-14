/*
 * $Id$
 * Copyright (C) 2010 John D Lamb
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef CCGSL_SF_COULOMB_HPP
#define CCGSL_SF_COULOMB_HPP

#include<gsl/gsl_sf_coulomb.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_hydrogenicR_1().
     * Normalized hydrogenic bound states, radial dependence. 
     * R_1 := 2Z sqrt(Z) exp(-Z r)
     * @param Z A real value
     * @param r A real value
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int hydrogenicR_1_e( double const Z, double const r, result& result ){
      return gsl_sf_hydrogenicR_1_e( Z, r, &result ); }
    /**
     * C++ version of gsl_sf_hydrogenicR_1().
     * Normalized hydrogenic bound states, radial dependence. 
     * R_1 := 2Z sqrt(Z) exp(-Z r)
     * @param Z A real value
     * @param r A real value
     * @return The function value
     */
    inline double hydrogenicR_1( double const Z, double const r ){
      return gsl_sf_hydrogenicR_1( Z, r ); } 
    /**
     * C++ version of gsl_sf_hydrogenicR_e().
     * R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
     * normalization such that psi(n,l,r) = R_n Y_{lm}
     * @param n An integer
     * @param l an integer
     * @param Z A real value
     * @param r A real value
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int hydrogenicR_e( int const n, int const l, double const Z, double const r, result& result ){
      return gsl_sf_hydrogenicR_e( n, l, Z, r, &result ); } 
    /**
     * C++ version of gsl_sf_hydrogenicR().
     * R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
     * normalization such that psi(n,l,r) = R_n Y_{lm}
     * @param n An integer
     * @param l An integer
     * @param Z A real value
     * @param r A real value
     * @return The function value
     */
    inline double hydrogenicR( int const n, int const l, double const Z, double const r ){
      return gsl_sf_hydrogenicR( n, l, Z, r ); } 
    /**
     * Namespace for @c gsl_sf_coulomb functions.
     */
    namespace coulomb {
      /**
       * C++ version of gsl_sf_coulomb_wave_FG_e().
       * Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
       * and their derivatives; lam_G := lam_F - k_lam_G
       *
       * lam_F, lam_G > -0.5
       * x > 0.0
       *
       * Conventions of Abramowitz+Stegun.
       *
       * Because there can be a large dynamic range of values,
       * overflows are handled gracefully. If an overflow occurs,
       * GSL_EOVRFLW is signalled and exponent(s) are returned
       * through exp_F, exp_G. These are such that
       *
       *   F_L(eta,x)  =  fc[k_L] * exp(exp_F)
       *   G_L(eta,x)  =  gc[k_L] * exp(exp_G)
       *   F_L'(eta,x) = fcp[k_L] * exp(exp_F)
       *   G_L'(eta,x) = gcp[k_L] * exp(exp_G)
       * @param eta A real value
       * @param x A real value
       * @param lam_F A real value
       * @param k_lam_G An integer
       * @param F A result as a @c gsl::sf::result object
       * @param Fp A result as a @c gsl::sf::result object
       * @param G A result as a @c gsl::sf::result object
       * @param Gp A result as a @c gsl::sf::result object
       * @param exp_F A double
       * @param exp_G A double
       * @return GSL_SUCCESS or GSLEOVERFLW
       */
      inline int wave_FG_e( double const eta, double const x, double const lam_F, int const k_lam_G,
			    result& F, result& Fp, result& G, result& Gp, double* exp_F, double* exp_G ){
	return gsl_sf_coulomb_wave_FG_e( eta, x, lam_F, k_lam_G, &F, &Fp, &G, &Gp, exp_F, exp_G ); } 
      /**
       * C++ version of gsl_sf_coulomb_wave_F_array().
       * @param lam_min A real value
       * @param kmax An integer
       * @param eta A real value
       * @param x A real value
       * @param fc_array An array
       * @param F_exponent A double
       * @return Error code on failure
       */
      inline int wave_F_array( double lam_min, int kmax, double eta, double x,
			       double* fc_array, double* F_exponent ){
	return gsl_sf_coulomb_wave_F_array( lam_min, kmax, eta, x, fc_array, F_exponent ); } 
      /**
       * C++ version of gsl_sf_coulomb_wave_FG_array().
       * @param lam_min A real value
       * @param kmax An integer
       * @param eta A real value
       * @param x A real value
       * @param fc_array An array
       * @param gc_array An array
       * @param F_exponent A double
       * @param G_exponent A double
       * @return Error code on failure
       */
      inline int wave_FG_array( double lam_min, int kmax, double eta, double x, double* fc_array,
				double* gc_array, double* F_exponent, double* G_exponent ){
	return gsl_sf_coulomb_wave_FG_array( lam_min, kmax, eta, x, fc_array, gc_array,
					     F_exponent, G_exponent ); } 
      /**
       * C++ version of gsl_sf_coulomb_wave_FGp_array().
       * @param lam_min A real value
       * @param kmax An integer
       * @param eta A real value
       * @param x A real value
       * @param fc_array An array
       * @param fcp_array An array
       * @param gc_array An array
       * @param gcp_array An array
       * @param F_exponent A double
       * @param G_exponent A double
       * @return Error code on failure
       */
      inline int wave_FGp_array( double lam_min, int kmax, double eta, double x, double* fc_array,
				 double* fcp_array, double* gc_array, double* gcp_array, double* F_exponent,
				 double* G_exponent ){
	return gsl_sf_coulomb_wave_FGp_array( lam_min, kmax, eta, x, fc_array, fcp_array,
					      gc_array, gcp_array, F_exponent, G_exponent ); } 
      /**
       * C++ version of gsl_sf_coulomb_wave_sphF_array().
       * @param lam_min A real value
       * @param kmax An integer
       * @param eta A real value
       * @param x A real value
       * @param fc_array An array
       * @param F_exponent A double
       * @return Error code on failure
       */
      inline int wave_sphF_array( double lam_min, int kmax, double eta, double x,
				  double* fc_array, double* F_exponent ){
	return gsl_sf_coulomb_wave_sphF_array( lam_min, kmax, eta, x, fc_array, F_exponent ); } 
      /**
       * C++ version of gsl_sf_coulomb_CL_e().
       * @param L A real value
       * @param eta A real value
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int CL_e( double L, double eta, result& result ){
	return gsl_sf_coulomb_CL_e( L, eta, &result ); } 
      /**
       * C++ version of gsl_sf_coulomb_CL_array().
       * @param Lmin A real value
       * @param kmax An integer
       * @param eta A real value
       * @param cl An array
       * @return Error code on failure
       */
      inline int CL_array( double Lmin, int kmax, double eta, double* cl ){
	return gsl_sf_coulomb_CL_array( Lmin, kmax, eta, cl ); }
    }
  }
}

#endif
