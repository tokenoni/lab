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

#ifndef CCGSL_SF_EXPINT_HPP
#define CCGSL_SF_EXPINT_HPP

#include<gsl/gsl_sf_expint.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_expint_E1_e().
     * E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int expint_E1_e( double x, result& result ){ return gsl_sf_expint_E1_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_E1().
     * E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param x A real number
     * @return The function value
     */
    inline double expint_E1( double const x ){ return gsl_sf_expint_E1( x ); } 
    /**
     * C++ version of gsl_sf_expint_E2_e().
     * E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int expint_E2_e( double const x, result& result ){ return gsl_sf_expint_E2_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_E2().
     * E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param x A real number
     * @return The function value
     */
    inline double expint_E2( double const x ){ return gsl_sf_expint_E2( x ); } 
    /**
     * C++ version of gsl_sf_expint_En_e().
     * E_n(x) := Re[ Integrate[ Exp[-xt]/t^n, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param n An integer
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int expint_En_e( int const n, double const x, result& result ){
      return gsl_sf_expint_En_e( n, x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_En().
     * E_n(x) := Re[ Integrate[ Exp[-xt]/t^n, {t,1,Infinity}] ]
     *
     * x != 0.0
     * @param n An integer
     * @param x A real number
     * @return The function value
     */
    inline double expint_En( int const n, double const x ){ return gsl_sf_expint_En( n, x ); } 
    /**
     * C++ version of gsl_sf_expint_E1_scaled_e().
     * E_1_scaled(x) := exp(x) E_1(x)
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int expint_E1_scaled_e( double const x, result& result ){
      return gsl_sf_expint_E1_scaled_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_E1_scaled().
     * E_1_scaled(x) := exp(x) E_1(x)
     *
     * x != 0.0
     * @param x A real number
     * @return The function value
     */
    inline double expint_E1_scaled( double const x ){ return gsl_sf_expint_E1_scaled( x ); } 
    /**
     * C++ version of gsl_sf_expint_E2_scaled_e().
     * E_2_scaled(x) := exp(x) E_2(x)
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return The function value
     */
    inline int expint_E2_scaled_e( double const x, result& result ){
      return gsl_sf_expint_E2_scaled_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_E2_scaled().
     * E_2_scaled(x) := exp(x) E_2(x)
     *
     * x != 0.0
     * @param x A real number
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline double expint_E2_scaled( double const x ){
      return gsl_sf_expint_E2_scaled( x ); } 
    /**
     * C++ version of gsl_sf_expint_En_scaled_e().
     * E_n_scaled(x) := exp(x) E_n(x)
     *
     * x != 0.0
     * @param n An integer
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return The function value
     */
    inline int expint_En_scaled_e( int const n, double const x, result& result ){
      return gsl_sf_expint_En_scaled_e( n, x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_En_scaled().
     * E_n_scaled(x) := exp(x) E_n(x)
     *
     * x != 0.0
     * @param n An integer
     * @param x A real number
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline double expint_En_scaled( int const n, double const x ){
      return gsl_sf_expint_En_scaled( n, x ); } 
    /**
     * C++ version of gsl_sf_expint_Ei_e().
     * Ei(x) := - PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
     *       :=   PV Integrate[ Exp[t]/t, {t,-Infinity,x}]
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return The function value
     */
    inline int expint_Ei_e( double const x, result& result ){
      return gsl_sf_expint_Ei_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_Ei().
     * Ei(x) := - PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
     *       :=   PV Integrate[ Exp[t]/t, {t,-Infinity,x}]
     *
     * x != 0.0
     * @param x A real number
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline double expint_Ei( double const x ){ return gsl_sf_expint_Ei( x ); } 
    /**
     * C++ version of gsl_sf_expint_Ei_scaled_e().
     * Ei_scaled(x) := exp(-x) Ei(x)
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int expint_Ei_scaled_e( double const x, result& result ){
      return gsl_sf_expint_Ei_scaled_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_Ei_scaled().
     * Ei_scaled(x) := exp(-x) Ei(x)
     *
     * x != 0.0
     * @param x A real number
     * @return The function value
     */
    inline double expint_Ei_scaled( double const x ){ return gsl_sf_expint_Ei_scaled( x ); } 
    /**
     * C++ version of gsl_sf_Shi_e().
     * Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int Shi_e( double const x, result& result ){ return gsl_sf_Shi_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_Shi().
     * Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double Shi( double const x ){ return gsl_sf_Shi( x ); } 
    /**
     * C++ version of gsl_sf_Chi_e().
     * Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
     *
     * x != 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int Chi_e( double const x, result& result ){ return gsl_sf_Chi_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_Chi().
     * Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
     *
     * x != 0.0
     * @param x A real number
     * @return The function value
     */
    inline double Chi( double const x ){ return gsl_sf_Chi( x ); } 
    /**
     * C++ version of gsl_sf_expint_3_e().
     * Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
     *
     * x >= 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int expint_3_e( double const x, result& result ){ return gsl_sf_expint_3_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expint_3().
     * Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
     *
     * x >= 0.0
     * @param x A real number
     * @return The function value
     */
    inline double expint_3( double x ){ return gsl_sf_expint_3( x ); } 
    /**
     * C++ version of gsl_sf_Si_e().
     * Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int Si_e( double const x, result& result ){ return gsl_sf_Si_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_Si().
     * Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double Si( double const x ){ return gsl_sf_Si( x ); } 
    /**
     * C++ version of gsl_sf_Ci_e().
     * Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
     *
     * x > 0.0
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int Ci_e( double const x, result& result ){ return gsl_sf_Ci_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_Ci().
     * Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
     *
     * x > 0.0
     * @param x A real number
     * @return The function value
     */
    inline double Ci( double const x ){ return gsl_sf_Ci( x ); } 
    /**
     * C++ version of gsl_sf_atanint_e().
     * AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int atanint_e( double const x, result& result ){ return gsl_sf_atanint_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_atanint().
     * AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double atanint( double const x ){ return gsl_sf_atanint( x ); }
  }
}

#endif
