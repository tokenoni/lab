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

#ifndef CCGSL_SF_FERMI_DIRAC_HPP
#define CCGSL_SF_FERMI_DIRAC_HPP

#include<gsl/gsl_sf_fermi_dirac.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_fermi_dirac_m1_e().
     * Complete integral F_{-1}(x) = e^x / (1 + e^x)
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW
     */
    inline int fermi_dirac_m1_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_m1_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_m1().
     * Complete integral F_{-1}(x) = e^x / (1 + e^x)
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_m1( double const x ){
      return gsl_sf_fermi_dirac_m1( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_0_e().
     * Complete integral F_0(x) = ln(1 + e^x)
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_0_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_0_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_0().
     * Complete integral F_0(x) = ln(1 + e^x)
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_0( double const x ){
      return gsl_sf_fermi_dirac_0( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_1_e().
     * F_1(x): F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_1_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_1_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_1().
     * F_1(x): F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_1( double const x ){
      return gsl_sf_fermi_dirac_1( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_2_e().
     * F_2(x): F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_2_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_2_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_2().
     * F_2(x): F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_2( double const x ){
      return gsl_sf_fermi_dirac_2( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_int_e().
     * F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param j An integer
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_int_e( int const j, double const x, result& result ){
      return gsl_sf_fermi_dirac_int_e( j, x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_int().
     * F_j(x)   := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param j An integer
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_int( int const j, double const x ){
      return gsl_sf_fermi_dirac_int( j, x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_mhalf_e().
     * F_{-1/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_mhalf_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_mhalf_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_mhalf().
     * F_{-1/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_mhalf( double const x ){
      return gsl_sf_fermi_dirac_mhalf( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_half_e().
     * F_{1/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_half_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_half_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_half().
     * F_{1/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_half( double const x ){
      return gsl_sf_fermi_dirac_half( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_3half_e().
     * F_{3/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int fermi_dirac_3half_e( double const x, result& result ){
      return gsl_sf_fermi_dirac_3half_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_3half().
     * F_{3/2}(x): F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double fermi_dirac_3half( double const x ){
      return gsl_sf_fermi_dirac_3half( x ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_inc_0_e().
     * Incomplete integral F_0(x,b) = ln(1 + e^(b-x)) - (b-x)
     * @param x A real number
     * @param b A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM
     */
    inline int fermi_dirac_inc_0_e( double const x, double const b, result& result ){
      return gsl_sf_fermi_dirac_inc_0_e( x, b, &result ); } 
    /**
     * C++ version of gsl_sf_fermi_dirac_inc_0().
     * Incomplete integral F_0(x,b) = ln(1 + e^(b-x)) - (b-x)
     * @param x A real number
     * @param b A real number
     * @return The function value
     */
    inline double fermi_dirac_inc_0( double const x, double const b ){
      return gsl_sf_fermi_dirac_inc_0( x, b );}
  }
}

#endif
