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

#ifndef CCGSL_SF_ELEMENTARY_HPP
#define CCGSL_SF_ELEMENTARY_HPP

#include<gsl/gsl_sf_elementary.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_multiply_e().
     * @param x A real number
     * @param y Another real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int multiply_e( double const x, double const y, result& result ){
      return gsl_sf_multiply_e( x, y, &result ); } 
    /**
     * C++ version of gsl_sf_multiply().
     * @param x A real number
     * @param y Another real number
     * @return GSL_SUCCESS or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline double multiply( double const x, double const y ){ return gsl_sf_multiply( x, y ); } 
    /**
     * C++ version of gsl_sf_multiply_err_e().
     * @param x A real number
     * @param dx Error in x
     * @param y Another real number
     * @param dy Error in y
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int multiply_err_e( double const x, double const dx, double const y, double const dy,
        result& result ){ return gsl_sf_multiply_err_e( x, dx, y, dy, &result ); } 
  }
}

#endif
