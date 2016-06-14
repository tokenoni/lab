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

#ifndef CCGSL_SF_LAGUERRE_HPP
#define CCGSL_SF_LAGUERRE_HPP

#include<gsl/gsl_sf_laguerre.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_laguerre_1_e().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int laguerre_1_e( double const a, double const x, result& result ){
      return gsl_sf_laguerre_1_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_laguerre_2_e().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int laguerre_2_e( double const a, double const x, result& result ){
      return gsl_sf_laguerre_2_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_laguerre_3_e().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int laguerre_3_e( double const a, double const x, result& result ){
      return gsl_sf_laguerre_3_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_laguerre_1().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @return The function value
     */
    inline double laguerre_1( double a, double x ){ return gsl_sf_laguerre_1( a, x ); } 
    /**
     * C++ version of gsl_sf_laguerre_2().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @return The function value
     */
    inline double laguerre_2( double a, double x ){ return gsl_sf_laguerre_2( a, x ); } 
    /**
     * C++ version of gsl_sf_laguerre_3().
     * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
     * @param a A real value
     * @param x A real value
     * @return The function value
     */
    inline double laguerre_3( double a, double x ){ return gsl_sf_laguerre_3( a, x ); } 
    /**
     * C++ version of gsl_sf_laguerre_n_e().
     * Evaluate generalized Laguerre polynomials.
     *
     * a > -1.0
     * n >= 0
     * @param n An integer
     * @param a A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int laguerre_n_e( int const n, double const a, double const x, result& result ){
      return gsl_sf_laguerre_n_e( n, a, x, &result ); } 
    /**
     * C++ version of gsl_sf_laguerre_n().
     * Evaluate generalized Laguerre polynomials.
     *
     * a > -1.0
     * n >= 0
     * @param n An integer
     * @param a A real value
     * @param x A real value
     * @return The function value
     */
    inline double laguerre_n( int n, double a, double x ){ return gsl_sf_laguerre_n( n, a, x ); }
  }
}

#endif
