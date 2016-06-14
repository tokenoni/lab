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

#ifndef CCGSL_SF_EXP_HPP
#define CCGSL_SF_EXP_HPP

#include<gsl/gsl_sf_exp.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_exp_e().
     * Exponential function
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or  GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int exp_e( double x, result& result ){ return gsl_sf_exp_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_exp().
     * Exponential function
     * @param x A real number
     * @return The function value
     */
    inline double exp( double const x ){ return gsl_sf_exp( x ); } 
    /**
     * C++ version of gsl_sf_exp_e10_e().
     * Exp(x)
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int exp_e10_e( double const x, result_e10& result ){
      return gsl_sf_exp_e10_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_exp_mult_e().
     * Exponentiate and multiply by a given factor:  y * Exp(x)
     * @param x A real number
     * @param y A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int exp_mult_e( double const x, double const y, result& result ){
      return gsl_sf_exp_mult_e( x, y, &result ); } 
    /**
     * C++ version of gsl_sf_exp_mult().
     * Exponentiate and multiply by a given factor:  y * Exp(x)
     * @param x A real number
     * @param y A real number
     * @return The function value
     */
    inline double exp_mult( double const x, double const y ){
      return gsl_sf_exp_mult( x, y ); } 
    /**
     * C++ version of gsl_sf_exp_mult_e10_e().
     * Exponentiate and multiply by a given factor:  y * Exp(x)
     * @param x A real number
     * @param y A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int exp_mult_e10_e( double const x, double const y, result_e10& result ){
      return gsl_sf_exp_mult_e10_e( x, y, &result ); } 
    /**
     * C++ version of gsl_sf_expm1_e().
     * exp(x)-1
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW
     */
    inline int expm1_e( double const x, result& result ){
      return gsl_sf_expm1_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_expm1().
     * exp(x)-1
     * @param x A real number
     * @return The function value
     */
    inline double expm1( double const x ){ return gsl_sf_expm1( x ); } 
    /**
     * C++ version of gsl_sf_exprel_e().
     * (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW
     */
    inline int exprel_e( double const x, result& result ){
      return gsl_sf_exprel_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_exprel().
     * (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
     * @param x A real number
     * @return The function value
     */
    inline double exprel( double const x ){ return gsl_sf_exprel( x ); } 
    /**
     * C++ version of gsl_sf_exprel_2_e().
     * 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW
     */
    inline int exprel_2_e( double x, result& result ){
      return gsl_sf_exprel_2_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_exprel_2().
     * 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
     * @param x A real number
     * @return The function value
     */
    inline double exprel_2( double const x ){ return gsl_sf_exprel_2( x ); } 
    /**
     * C++ version of gsl_sf_exprel_n_e().
     * Similarly for the N-th generalization of
     * exprel_2. The so-called N-relative exponential
     *
     * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
     *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
     *             = 1F1(1,1+N,x)
     * @param n A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int exprel_n_e( int const n, double const x, result& result ){
      return gsl_sf_exprel_n_e( n, x, &result ); } 
    /**
     * C++ version of gsl_sf_exprel_n().
     * Similarly for the N-th generalization of
     * exprel_2. The so-called N-relative exponential
     *
     * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
     *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
     *             = 1F1(1,1+N,x)
     * @param n A real number
     * @param x A real number
     * @return The function value
     */
    inline double exprel_n( int const n, double const x ){
      return gsl_sf_exprel_n( n, x ); } 
    /**
     * C++ version of gsl_sf_exprel_n_CF_e().
     * @param n A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int exprel_n_CF_e( double const n, double const x, result& result ){
      return gsl_sf_exprel_n_CF_e( n, x, &result ); } 
    /**
     * C++ version of gsl_sf_exp_err_e().
     * Exponentiate a quantity with an associated error.
     * @param x A real number
     * @param dx A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int exp_err_e( double const x, double const dx, result& result ){
      return gsl_sf_exp_err_e( x, dx, &result ); } 
    /**
     * C++ version of gsl_sf_exp_err_e10_e().
     * Exponentiate a quantity with an associated error.
     * @param x A real number
     * @param dx A real number
     * @param result The result as a @c gsl::sf::result object
     * @return Error code on failure
     */
    inline int exp_err_e10_e( double const x, double const dx, result_e10& result ){
      return gsl_sf_exp_err_e10_e( x, dx, &result ); } 
    /**
     * C++ version of gsl_sf_exp_mult_err_e().
     * Exponentiate and multiply by a given factor:  y * Exp(x),
     * for quantities with associated errors.
     * @param x A real number
     * @param dx A real number
     * @param y A real number
     * @param dy A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int exp_mult_err_e( double const x, double const dx,
			       double const y, double const dy, result& result ){
      return gsl_sf_exp_mult_err_e( x, dx, y, dy, &result ); } 
    /**
     * C++ version of gsl_sf_exp_mult_err_e10_e().
     * Exponentiate and multiply by a given factor:  y * Exp(x),
     * for quantities with associated errors.
     * @param x A real number
     * @param dx A real number
     * @param y A real number
     * @param dy A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EOVRFLW
     */
    inline int exp_mult_err_e10_e( double const x, double const dx,
				   double const y, double const dy, result_e10& result ){
      return gsl_sf_exp_mult_err_e10_e( x, dx, y, dy, &result ); }
  }
}

#endif
