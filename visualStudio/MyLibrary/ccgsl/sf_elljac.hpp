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

#ifndef CCGSL_SF_ELLJAC_HPP
#define CCGSL_SF_ELLJAC_HPP

#include<gsl/gsl_sf_elljac.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_elljac_e().
     * Jacobian elliptic functions sn, dn, cn,
     * by descending Landen transformations
     * @param u A real number
     * @param m A real number
     * @param sn A real number
     * @param cn A real number
     * @param dn A real number
     * @return GSL_SUCCESS or GSL_EDOM
     */
    int gsl_sf_elljac_e( double const u, double const m, double* sn, double* cn, double* dn );
  }
}

#endif
