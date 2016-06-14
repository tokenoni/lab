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

#ifndef CCGSL_SF_HYPERG_HPP
#define CCGSL_SF_HYPERG_HPP

#include<gsl/gsl_sf_hyperg.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_hyperg_0F1_e().
     * Hypergeometric function related to Bessel functions
     * 0F1[c,x] =
     *            Gamma[c]    x^(1/2(1-c)) I_{c-1}(2 Sqrt[x])
     *            Gamma[c] (-x)^(1/2(1-c)) J_{c-1}(2 Sqrt[-x])
     * @param c A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int hyperg_0F1_e( double c, double x, result& result ){
      return gsl_sf_hyperg_0F1_e( c, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_0F1().
     * Hypergeometric function related to Bessel functions
     * 0F1[c,x] =
     *            Gamma[c]    x^(1/2(1-c)) I_{c-1}(2 Sqrt[x])
     *            Gamma[c] (-x)^(1/2(1-c)) J_{c-1}(2 Sqrt[-x])
     * @param c A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_0F1( double const c, double const x ){
      return gsl_sf_hyperg_0F1( c, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_1F1_int_e().
     * Confluent hypergeometric function  for integer parameters.
     * 1F1[m,n,x] = M(m,n,x)
     * @param m An integer
     * @param n An integer
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_1F1_int_e( int const m, int const n, double const x, result& result ){
      return gsl_sf_hyperg_1F1_int_e( m, n, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_1F1_int().
     * Confluent hypergeometric function  for integer parameters.
     * 1F1[m,n,x] = M(m,n,x)
     * @param m An integer
     * @param n An integer
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_1F1_int( int const m, int const n, double x ){
      return gsl_sf_hyperg_1F1_int( m, n, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_1F1_e().
     * Confluent hypergeometric function.
     * 1F1[a,b,x] = M(a,b,x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_1F1_e( double const a, double const b, double const x, result& result ){
      return gsl_sf_hyperg_1F1_e( a, b, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_1F1().
     * Confluent hypergeometric function.
     * 1F1[a,b,x] = M(a,b,x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_1F1( double a, double b, double x ){
      return gsl_sf_hyperg_1F1( a, b, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_U_int_e().
     * Confluent hypergeometric function for integer parameters.
     * U(m,n,x)
     * @param m An integer
     * @param n An integer
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_U_int_e( int const m, int const n, double const x, result& result ){
      return gsl_sf_hyperg_U_int_e( m, n, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_U_int().
     * Confluent hypergeometric function for integer parameters.
     * U(m,n,x)
     * @param m An integer
     * @param n An integer
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_U_int( int const m, int const n, double const x ){
      return gsl_sf_hyperg_U_int( m, n, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_U_int_e10_e().
     * Confluent hypergeometric function for integer parameters.
     * U(m,n,x)
     * @param m An integer
     * @param n An integer
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_U_int_e10_e( int const m, int const n, double const x, result_e10& result ){
      return gsl_sf_hyperg_U_int_e10_e( m, n, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_U_e().
     * Confluent hypergeometric function.
     * U(a,b,x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_U_e( double const a, double const b, double const x, result& result ){
      return gsl_sf_hyperg_U_e( a, b, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_U().
     * Confluent hypergeometric function.
     * U(a,b,x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_U( double const a, double const b, double const x ){
      return gsl_sf_hyperg_U( a, b, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_U_e10_e().
     * Confluent hypergeometric function.
     * U(a,b,x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_U_e10_e( double const a, double const b, double const x, result_e10& result ){
      return gsl_sf_hyperg_U_e10_e( a, b, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_e().
     * Gauss hypergeometric function 2F1[a,b,c,x]
     * |x| < 1
     * @param a A real value
     * @param b A real value
     * @param c A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_2F1_e( double a, double b, double const c, double const x, result& result ){
      return gsl_sf_hyperg_2F1_e( a, b, c, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1().
     * Gauss hypergeometric function 2F1[a,b,c,x]
     * |x| < 1
     * @param a A real value
     * @param b A real value
     * @param c A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_2F1( double a, double b, double c, double x ){
      return gsl_sf_hyperg_2F1( a, b, c, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_conj_e().
     * Gauss hypergeometric function
     * 2F1[aR + I aI, aR - I aI, c, x]
     * @param aR Real part of a complex number
     * @param aI Imaginary part of a complex number
     * @param c A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_2F1_conj_e( double const aR, double const aI, double const c, double const x,
				  result& result ){ return gsl_sf_hyperg_2F1_conj_e( aR, aI, c, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_conj().
     * Gauss hypergeometric function
     * 2F1[aR + I aI, aR - I aI, c, x]
     * @param aR Real part of a complex number
     * @param aI Imaginary part of a complex number
     * @param c A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_2F1_conj( double aR, double aI, double c, double x ){
      return gsl_sf_hyperg_2F1_conj( aR, aI, c, x ); }
    /**
     * C++ version of gsl_sf_hyperg_2F1_renorm_e().
     * Renormalized Gauss hypergeometric function
     * 2F1[a,b,c,x] / Gamma[c]
     * |x| < 1
     * @param a A real value
     * @param b A real value
     * @param c A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_2F1_renorm_e( double const a, double const b, double const c, double const x,
				    result& result ){ return gsl_sf_hyperg_2F1_renorm_e( a, b, c, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_renorm().
     * Renormalized Gauss hypergeometric function
     * 2F1[a,b,c,x] / Gamma[c]
     * |x| < 1
     * @param a A real value
     * @param b A real value
     * @param c A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_2F1_renorm( double a, double b, double c, double x ){
      return gsl_sf_hyperg_2F1_renorm( a, b, c, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_conj_renorm_e().
     * Renormalized Gauss hypergeometric function
     * 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
     * |x| < 1
     * @param aR Real part of a complex number
     * @param aI Imaginary part of a complex number
     * @param c A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int hyperg_2F1_conj_renorm_e( double const aR, double const aI, double const c,
					 double const x, result& result ){
      return gsl_sf_hyperg_2F1_conj_renorm_e( aR, aI, c, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F1_conj_renorm().
     * Renormalized Gauss hypergeometric function
     * 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
     * |x| < 1
     * @param aR Real part of a complex number
     * @param aI Imaginary part of a complex number
     * @param c A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_2F1_conj_renorm( double aR, double aI, double c, double x ){
      return gsl_sf_hyperg_2F1_conj_renorm( aR, aI, c, x ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F0_e().
     * Mysterious hypergeometric function. The series representation
     * is a divergent hypergeometric series. However, for x < 0 we
     * have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int hyperg_2F0_e( double const a, double const b, double const x, result& result ){
      return gsl_sf_hyperg_2F0_e( a, b, x, &result ); } 
    /**
     * C++ version of gsl_sf_hyperg_2F0().
     * Mysterious hypergeometric function. The series representation
     * is a divergent hypergeometric series. However, for x < 0 we
     * have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
     * @param a A real value
     * @param b A real value
     * @param x A real value
     * @return The function value
     */
    inline double hyperg_2F0( double const a, double const b, double const x ){
      return gsl_sf_hyperg_2F0( a, b, x ); }
  }
}

#endif
