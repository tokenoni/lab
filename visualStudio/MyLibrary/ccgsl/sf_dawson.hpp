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

#ifndef CCGSL_SF_DAWSON_HPP
#define CCGSL_SF_DAWSON_HPP

#include<gsl/gsl_sf_dawson.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_dawson_e().
     * Dawson's integral:
     *
     *   Exp[-x^2] Integral[ Exp[t^2], {t,0,x}]
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int dawson_e( double x, result& result ){ return gsl_sf_dawson_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_dawson().
     * Dawson's integral:
     *
     *   Exp[-x^2] Integral[ Exp[t^2], {t,0,x}]
     * @param x A real number
     * @return The integral
     */
    inline double dawson( double x ){ return gsl_sf_dawson( x ); } 
  }
}

#endif
