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

#ifndef CCGSL_SF_ELLINT_HPP
#define CCGSL_SF_ELLINT_HPP

#include<gsl/gsl_sf_ellint.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * Namespace for elliptic integrals.
     */
    namespace ellint {
      /**
       * C++ version of gsl_sf_ellint_Kcomp_e().
       * Legendre form of complete elliptic integrals
       *
       * K(k) = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
       *
       * @param k A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Kcomp_e( double k, mode_t mode, result& result ){
	return gsl_sf_ellint_Kcomp_e( k, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_Kcomp().
       * Legendre form of complete elliptic integrals
       *
       * K(k) = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
       *
       * @param k A real number
       * @param mode The mode
       * @return function value
       */
      inline double Kcomp( double k, mode_t mode ){
	return gsl_sf_ellint_Kcomp( k, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_Ecomp_e().
       * Legendre form of complete elliptic integrals
       *
       * E(k) = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
       *
       * @param k A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Ecomp_e( double k, mode_t mode, result& result ){
	return gsl_sf_ellint_Ecomp_e( k, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_Ecomp().
       * Legendre form of complete elliptic integrals
       *
       * E(k) = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, Pi/2}]
       *
       * @param k A real number
       * @param mode The mode
       * @return The function value
       */
      inline double Ecomp( double k, mode_t mode ){
	return gsl_sf_ellint_Ecomp( k, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_Pcomp_e().
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int Pcomp_e( double k, double n, mode_t mode, result& result ){
	return gsl_sf_ellint_Pcomp_e( k, n, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_Pcomp().
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @return The function value
       */
      inline double Pcomp( double k, double n, mode_t mode ){
	return gsl_sf_ellint_Pcomp( k, n, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_Dcomp_e().
       * @param k A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int Dcomp_e( double k, mode_t mode, result& result ){
	return gsl_sf_ellint_Dcomp_e( k, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_Dcomp().
       * @param k A real number
       * @param mode The mode
       * @return The function value
       */
      inline double Dcomp( double k, mode_t mode ){
	return gsl_sf_ellint_Dcomp( k, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_F_e().
       * F(phi,k)   = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int F_e( double phi, double k, mode_t mode, result& result ){
	return gsl_sf_ellint_F_e( phi, k, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_F().
       * F(phi,k)   = Integral[1/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param mode The mode
       * @return The function value
       */
      inline double F( double phi, double k, mode_t mode ){
	return gsl_sf_ellint_F( phi, k, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_E_e().
       * E(phi,k)   = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int E_e( double phi, double k, mode_t mode, result& result ){
	return gsl_sf_ellint_E_e( phi, k, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_E().
       * E(phi,k)   = Integral[  Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param mode The mode
       * @return The function value
       */
      inline double E( double phi, double k, mode_t mode ){
	return gsl_sf_ellint_E( phi, k, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_P_e().
       * P(phi,k,n) = Integral[(1 + n Sin[t]^2)^(-1)/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int P_e( double phi, double k, double n, mode_t mode, result& result ){
	return gsl_sf_ellint_P_e( phi, k, n, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_P().
       * P(phi,k,n) = Integral[(1 + n Sin[t]^2)^(-1)/Sqrt[1 - k^2 Sin[t]^2], {t, 0, phi}]
       * @param phi A real number
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @return The function value
       */
      inline double P( double phi, double k, double n, mode_t mode ){
	return gsl_sf_ellint_P( phi, k, n, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_D_e().
       * D(phi,k,n) = R_D(1-Sin[phi]^2, 1-k^2 Sin[phi]^2, 1.0)
       * @param phi A real number
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return Error code on failure
       */
      inline int D_e( double phi, double k, double n, mode_t mode, result& result ){
	return gsl_sf_ellint_D_e( phi, k, n, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_D().
       * D(phi,k,n) = R_D(1-Sin[phi]^2, 1-k^2 Sin[phi]^2, 1.0)
       * @param phi A real number
       * @param k A real number
       * @param n A real number
       * @param mode The mode
       * @return The function value
       */
      inline double D( double phi, double k, double n, mode_t mode ){
	return gsl_sf_ellint_D( phi, k, n, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_RC_e().
       * Carlson's symmetric basis of functions
       *
       * RC(x,y)   = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1)], {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int RC_e( double x, double y, mode_t mode, result& result ){
	return gsl_sf_ellint_RC_e( x, y, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_RC().
       * Carlson's symmetric basis of functions
       *
       * RC(x,y)   = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1)], {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param mode The mode
       * @return The function value
       */
      inline double RC( double x, double y, mode_t mode ){
	return gsl_sf_ellint_RC( x, y, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_RD_e().
       * Carlson's symmetric basis of functions
       *
       * RD(x,y,z) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int RD_e( double x, double y, double z, mode_t mode, result& result ){
	return gsl_sf_ellint_RD_e( x, y, z, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_RD().
       * Carlson's symmetric basis of functions
       *
       * RD(x,y,z) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param mode The mode
       * @return The function value
       */
      inline double RD( double x, double y, double z, mode_t mode ){
	return gsl_sf_ellint_RD( x, y, z, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_RF_e().
       * Carlson's symmetric basis of functions
       *
       * RF(x,y,z) = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int RF_e( double x, double y, double z, mode_t mode, result& result ){
	return gsl_sf_ellint_RF_e( x, y, z, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_RF().
       * Carlson's symmetric basis of functions
       *
       * RF(x,y,z) = 1/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param mode The mode
       * @return The function value
       */
      inline double RF( double x, double y, double z, mode_t mode ){
	return gsl_sf_ellint_RF( x, y, z, mode ); } 
      /**
       * C++ version of gsl_sf_ellint_RJ_e().
       * Carlson's symmetric basis of functions
       *
       * RJ(x,y,z,p) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param p A real number
       * @param mode The mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int RJ_e( double x, double y, double z, double p, mode_t mode, result& result ){
	return gsl_sf_ellint_RJ_e( x, y, z, p, mode, &result ); } 
      /**
       * C++ version of gsl_sf_ellint_RJ().
       * Carlson's symmetric basis of functions
       *
       * RJ(x,y,z,p) = 3/2 Integral[(t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1), {t,0,Inf}]
       * @param x A real number
       * @param y A real number
       * @param z A real number
       * @param p A real number
       * @param mode The mode
       * @return The function value
       */
      inline double RJ( double x, double y, double z, double p, mode_t mode ){
	return gsl_sf_ellint_RJ( x, y, z, p, mode ); } 
    }
  }
}

#endif
