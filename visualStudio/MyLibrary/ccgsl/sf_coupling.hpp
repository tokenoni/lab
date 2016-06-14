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

#ifndef CCGSL_SF_COUPLING_HPP
#define CCGSL_SF_COUPLING_HPP

#include<gsl/gsl_sf_coupling.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_coupling_3j().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_ma Coupling coefficient in half-integer units
     * @param two_mb Coupling coefficient in half-integer units
     * @param two_mc Coupling coefficient in half-integer units
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVERFLW
     */
    inline int coupling_3j_e( int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc,
			      result& result ){
      return gsl_sf_coupling_3j_e( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, &result ); }
    /**
     * C++ version of gsl_sf_coupling_3j().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_ma Coupling coefficient in half-integer units
     * @param two_mb Coupling coefficient in half-integer units
     * @param two_mc Coupling coefficient in half-integer units
     * @return The Wigner 3-j coefficient
     */
    inline double coupling_3j( int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc ){
      return gsl_sf_coupling_3j( two_ja, two_jb, two_jc, two_ma, two_mb, two_mc ); } 
    /**
     * C++ version of gsl_sf_coupling_6j_e().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @param result Coupling coefficient in half-integer units
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVERFLW
     */
    inline int coupling_6j_e( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf,
			      result& result ){
      return gsl_sf_coupling_6j_e( two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &result ); } 
    /**
     * C++ version of gsl_sf_coupling_6j().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @return The Wigner 6-j coefficient
     */
    inline double coupling_6j( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf ){
      return gsl_sf_coupling_6j( two_ja, two_jb, two_jc, two_jd, two_je, two_jf ); } 
    /**
     * C++ version of gsl_sf_coupling_RacahW_e().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVERFLW
     */
    inline int coupling_RacahW_e( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf,
				  result& result ){
      return gsl_sf_coupling_RacahW_e( two_ja, two_jb, two_jc, two_jd, two_je, two_jf, &result ); } 
    /**
     * C++ version of gsl_sf_coupling_RacahW().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @return The function value
     */
    inline double coupling_RacahW( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf ){
      return gsl_sf_coupling_RacahW( two_ja, two_jb, two_jc, two_jd, two_je, two_jf ); } 
    /**
     * C++ version of gsl_sf_coupling_9j_e().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @param two_jg Coupling coefficient in half-integer units
     * @param two_jh Coupling coefficient in half-integer units
     * @param two_ji Coupling coefficient in half-integer units
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVERFLW
     */
    inline int coupling_9j_e( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg,
			      int two_jh, int two_ji, result& result ){
      return gsl_sf_coupling_9j_e( two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji,
				   &result ); } 
    /**
     * C++ version of gsl_sf_coupling_9j().
     * @param two_ja Coupling coefficient in half-integer units
     * @param two_jb Coupling coefficient in half-integer units
     * @param two_jc Coupling coefficient in half-integer units
     * @param two_jd Coupling coefficient in half-integer units
     * @param two_je Coupling coefficient in half-integer units
     * @param two_jf Coupling coefficient in half-integer units
     * @param two_jg Coupling coefficient in half-integer units
     * @param two_jh Coupling coefficient in half-integer units
     * @param two_ji Coupling coefficient in half-integer units
     * @return The Wigner 9-j coefficient
     */
    inline double coupling_9j( int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf,
			       int two_jg, int two_jh, int two_ji ){
      return gsl_sf_coupling_9j( two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji ); }
  }
}

#endif
