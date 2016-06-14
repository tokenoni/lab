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

#ifndef CCGSL_SF_ERF_HPP
#define CCGSL_SF_ERF_HPP

#include<gsl/gsl_sf_erf.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_erfc_e().
     * Complementary Error Function
     * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int erfc_e( double x, result& result ){ return gsl_sf_erfc_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_erfc().
     * Complementary Error Function
     * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
     * @param x A real number
     * @return The function value
     */
    inline double erfc( double x ){ return gsl_sf_erfc( x ); } 
    /**
     * C++ version of gsl_sf_log_erfc_e().
     * Log Complementary Error Function
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int log_erfc_e( double x, result& result ){ return gsl_sf_log_erfc_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_log_erfc().
     * Log Complementary Error Function
     * @param x A real number
     * @return The function value
     */
    inline double log_erfc( double x ){ return gsl_sf_log_erfc( x ); } 
    /**
     * C++ version of gsl_sf_erf_e().
     * Error Function
     * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int erf_e( double x, result& result ){ return gsl_sf_erf_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_erf().
     * Error Function
     * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
     * @param x A real number
     * @return The function value
     */
    inline double erf( double x ){ return gsl_sf_erf( x ); } 
    /**
     * C++ version of gsl_sf_erf_Z_e().
     * Probability function
     * Z(x) :  Abramowitz+Stegun 26.2.1
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int erf_Z_e( double x, result& result ){ return gsl_sf_erf_Z_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_erf_Q_e().
     * Probability function
     * Q(x) :  Abramowitz+Stegun 26.2.3
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int erf_Q_e( double x, result& result ){ return gsl_sf_erf_Q_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_erf_Z().
     * Probability function
     * Z(x) :  Abramowitz+Stegun 26.2.1
     * @param x A real number
     * @return The function value
     */
    inline double erf_Z( double x ){ return gsl_sf_erf_Z( x ); } 
    /**
     * C++ version of gsl_sf_erf_Q().
     * Probability function
     * Q(x) :  Abramowitz+Stegun 26.2.3
     * @param x A real number
     * @return The function value
     */
    inline double erf_Q( double x ){ return gsl_sf_erf_Q( x ); } 
    /**
     * C++ version of gsl_sf_hazard_e().
     * Hazard function, also known as the inverse Mill's ratio.
     *
     *   H(x) := Z(x)/Q(x)
     *         = Sqrt[2/Pi] Exp[-x^2 / 2] / Erfc[x/Sqrt[2]]
     *
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW
     */
    inline int hazard_e( double x, result& result ){ return gsl_sf_hazard_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_hazard().
     * Hazard function, also known as the inverse Mill's ratio.
     *
     *   H(x) := Z(x)/Q(x)
     *         = Sqrt[2/Pi] Exp[-x^2 / 2] / Erfc[x/Sqrt[2]]
     *
     * @param x A real number
     * @return The function value
     */
    inline double hazard( double x ){ return gsl_sf_hazard( x ); } 
  }
}

#endif
