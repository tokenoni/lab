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

#ifndef CCGSL_SF_DILOG_HPP
#define CCGSL_SF_DILOG_HPP

#include<gsl/gsl_sf_dilog.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_dilog().
     * Real part of DiLogarithm(x), for real argument.
     * In Lewin's notation, this is Li_2(x).
     *
     * Li_2(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int dilog_e( double x, result& result ){ return gsl_sf_dilog_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_dilog().
     * Real part of DiLogarithm(x), for real argument.
     * In Lewin's notation, this is Li_2(x).
     *
     * Li_2(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
     * @param x A real number
     * @return Error code on failure
     */
    inline double dilog( double const x ){ return gsl_sf_dilog( x ); } 
    /**
     * C++ version of gsl_sf_complex_dilog_xy_e().
     * DiLogarithm(z), for complex argument z = x + i y.
     * Computes the principal branch.
     * @param x A real number
     * @param y A real number
     * @param result_re Real part of result as @c gsl::sf::result object
     * @param result_im Imaginary part of result as @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int complex_dilog_xy_e( double const x, double const y, result& result_re, result& result_im ){
      return gsl_sf_complex_dilog_xy_e( x, y, &result_re, &result_im ); } 
    /**
     * C++ version of gsl_sf_complex_dilog_e().
     * DiLogarithm(z), for complex argument z = r Exp[i theta].
     * Computes the principal branch, thereby assuming an
     * implicit reduction of theta to the range (-2 pi, 2 pi).
     * @param r A real number
     * @param theta A real number
     * @param result_re Real part of result as @c gsl::sf::result object
     * @param result_im Imaginary part of result as @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int complex_dilog_e( double const r, double const theta, result& result_re, result& result_im ){
      return gsl_sf_complex_dilog_e( r, theta, &result_re, &result_im ); } 
    /**
     * C++ version of gsl_sf_complex_spence_xy_e().
     * Spence integral; spence(s) := Li_2(1-s) 
     * @param x A real number
     * @param y A real number
     * @param real_sp Real part of result as @c gsl::sf::result object
     * @param imag_sp Imaginary part of result as @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int complex_spence_xy_e( double const x, double const y, result& real_sp, result& imag_sp ){
      return gsl_sf_complex_spence_xy_e( x, y, &real_sp, &imag_sp ); }
  }
}

#endif
