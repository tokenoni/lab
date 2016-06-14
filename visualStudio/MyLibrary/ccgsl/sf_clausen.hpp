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

#ifndef CCGSL_SF_CLAUSEN_HPP
#define CCGSL_SF_CLAUSEN_HPP

#include<gsl/gsl_sf_clausen.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_clausen().
     * Calculate the Clausen integral:
     *   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
     *
     * Relation to dilogarithm:
     *   Cl_2(theta) = Im[ Li_2(e^(i theta)) ] 
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int clausen_e( double const x, result& result ){
      return gsl_sf_clausen_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_clausen().
     * Calculate the Clausen integral:
     *   Cl_2(x) := Integrate[-Log[2 Sin[t/2]], {t,0,x}]
     *
     * Relation to dilogarithm:
     *   Cl_2(theta) = Im[ Li_2(e^(i theta)) ] 
     * @param x A real value
     * @return The function applied to @c x
     */
    inline double clausen( double const x ){ return gsl_sf_clausen( x ); } 
  }
}

#endif
