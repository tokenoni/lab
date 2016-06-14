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

#ifndef CCGSL_SF_LAMBERT_HPP
#define CCGSL_SF_LAMBERT_HPP

#include<gsl/gsl_sf_lambert.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * Lambert functions.
     */
    namespace lambert {
      /**
       * C++ version of gsl_sf_lambert_W0_e().
       * Lambert's Function W_0(x)
       *
       * W_0(x) is the principal branch of the
       * implicit function defined by W e^W = x.
       *
       * \f$-1/E < x < \infty\f$
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EMAXITER
       */
      inline int W0_e( double x, result& result ){
	return gsl_sf_lambert_W0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_lambert_W0().
       * Lambert's Function W_0(x)
       *
       * W_0(x) is the principal branch of the
       * implicit function defined by W e^W = x.
       *
       * \f$-1/E < x < \infty\f$
       * @param x A real value
       * @return The function value
       */
      inline double W0( double x ){ return gsl_sf_lambert_W0( x ); } 
      /**
       * C++ version of gsl_sf_lambert_Wm1_e().
       * Lambert's Function W_{-1}(x)
       *
       * W_{-1}(x) is the second real branch of the
       * implicit function defined by W e^W = x.
       * It agrees with W_0(x) when x >= 0.
       *
       * \f$-1/E < x < \infty\f$
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EMAXITER
       */
      inline int Wm1_e( double x, result& result ){
	return gsl_sf_lambert_Wm1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_lambert_Wm1().
       * Lambert's Function W_{-1}(x)
       *
       * W_{-1}(x) is the second real branch of the
       * implicit function defined by W e^W = x.
       * It agrees with W_0(x) when x >= 0.
       *
       * \f$-1/E < x < \infty\f$
       * @param x A real value
       * @return The function value
       */
      inline double Wm1( double x ){ return gsl_sf_lambert_Wm1( x ); }
    }
  }
}

#endif
