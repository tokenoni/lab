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

#ifndef CCGSL_SF_DEBYE_HPP
#define CCGSL_SF_DEBYE_HPP

#include<gsl/gsl_sf_debye.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_debye_1_e().
     * D_1(x) := 1/x Integrate[t/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int debye_1_e( double const x, result& result ){ return gsl_sf_debye_1_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_1().
     * D_1(x) := n/x^n Integrate[t^1/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_1( double const x ){ return gsl_sf_debye_1( x ); } 
    /**
     * C++ version of gsl_sf_debye_2_e().
     * D_2(x) := 2/x^2 Integrate[t^2/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int debye_2_e( double const x, result& result ){ return gsl_sf_debye_2_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_2().
     * D_2(x) := 2/x^2 Integrate[t^2/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_2( double const x ){ return gsl_sf_debye_2( x ); } 
    /**
     * C++ version of gsl_sf_debye_3_e().
     * D_3(x) := 3/x^3 Integrate[t^3/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_UNDRFLW
     */
    inline int debye_3_e( double const x, result& result ){ return gsl_sf_debye_3_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_3().
     * D_3(x) := 3/x^3 Integrate[t^3/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_3( double const x ){ return gsl_sf_debye_3( x ); } 
    /**
     * C++ version of gsl_sf_debye_4_e().
     * D_4(x) := 4/x^4 Integrate[t^4/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_UNDRFLW
     */
    inline int debye_4_e( double const x, result& result ){ return gsl_sf_debye_4_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_4().
     * D_4(x) := 4/x^4 Integrate[t^4/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_4( double const x ){ return gsl_sf_debye_4( x ); } 
    /**
     * C++ version of gsl_sf_debye_5_e().
     * D_5(x) := 5/x^5 Integrate[t^5/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_UNDRFLW
     */
    inline int debye_5_e( double const x, result& result ){ return gsl_sf_debye_5_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_5().
     * D_5(x) := 5/x^5 Integrate[t^5/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_5( double const x ){ return gsl_sf_debye_5( x ); } 
    /**
     * C++ version of gsl_sf_debye_6_e().
     * D_6(x) := 6/x^6 Integrate[t^6/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_UNDRFLW
     */
    inline int debye_6_e( double const x, result& result ){ return gsl_sf_debye_6_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_debye_6().
     * D_6(x) := 6/x^6 Integrate[t^6/(e^t - 1), {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double debye_6( double const x ){ return gsl_sf_debye_6( x ); } 
  }
}

#endif
