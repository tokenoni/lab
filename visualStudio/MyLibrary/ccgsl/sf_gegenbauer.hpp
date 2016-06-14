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

#ifndef CCGSL_SF_GEGENBAUER_HPP
#define CCGSL_SF_GEGENBAUER_HPP

#include<gsl/gsl_sf_gegenbauer.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_gegenpoly_1_e().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int gegenpoly_1_e( double lambda, double x, result& result ){
      return gsl_sf_gegenpoly_1_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_2_e().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int gegenpoly_2_e( double lambda, double x, result& result ){
      return gsl_sf_gegenpoly_2_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_3_e().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int gegenpoly_3_e( double lambda, double x, result& result ){
      return gsl_sf_gegenpoly_3_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_1().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @return The function value
     */
    inline double gegenpoly_1( double lambda, double x ){
      return gsl_sf_gegenpoly_1( lambda, x ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_2().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @return The function value
     */
    inline double gegenpoly_2( double lambda, double x ){
      return gsl_sf_gegenpoly_2( lambda, x ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_3().
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @return The function value
     */
    inline double gegenpoly_3( double lambda, double x ){
      return gsl_sf_gegenpoly_3( lambda, x ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_n_e().
     * @param n A nonnegative integer
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gegenpoly_n_e( int n, double lambda, double x, result& result ){
      return gsl_sf_gegenpoly_n_e( n, lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_n().
     * @param n A nonnegative integer
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @return The function value
     */
    inline double gegenpoly_n( int n, double lambda, double x ){
      return gsl_sf_gegenpoly_n( n, lambda, x ); } 
    /**
     * C++ version of gsl_sf_gegenpoly_array().
     * @param nmax A nonnegative integer
     * @param lambda A real value greater than -0.5
     * @param x A real value
     * @param result_array An array of size @c nmax
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gegenpoly_array( int nmax, double lambda, double x, double* result_array ){
      return gsl_sf_gegenpoly_array( nmax, lambda, x, result_array );
    }
  }
}

#endif
