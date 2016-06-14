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

#ifndef CCGSL_SF_RESULT_HPP
#define CCGSL_SF_RESULT_HPP

#include<gsl/gsl_sf_result.h>

namespace gsl {
  namespace sf {
    /**
     * Typedef for @c gsl_sf_result.
     */
    typedef gsl_sf_result result;
    /**
     * Typedef for @c gsl_sf_result_e10.
     */
    typedef gsl_sf_result_e10 result_e10;
    /**
     * C++ version of gsl_sf_result_smash_e().
     * @param re Object of type result_e10
     * @param r Object of type result
     * @return Error code
     */
    inline int result_smash_e( gsl_sf_result_e10 const& re, gsl_sf_result&  r){
      return gsl_sf_result_smash_e( &re, &r ); }
  }
}

#endif
