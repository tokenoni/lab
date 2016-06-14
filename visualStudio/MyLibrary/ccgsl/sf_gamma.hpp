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

#ifndef CCGSL_SF_GAMMA_HPP
#define CCGSL_SF_GAMMA_HPP

#include<gsl/gsl_sf_gamma.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * C++ version of gsl_sf_lngamma_e().
     * Log[Gamma(x)], x not a negative integer
     * Uses real Lanczos method.
     * Returns the real part of Log[Gamma[x]] when x < 0,
     * i.e. Log[|Gamma[x]|].
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EROUND
     */
    inline int lngamma_e( double x, result& result ){
      return gsl_sf_lngamma_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_lngamma().
     * Log[Gamma(x)], x not a negative integer
     * Uses real Lanczos method.
     * Returns the real part of Log[Gamma[x]] when x < 0,
     * i.e. Log[|Gamma[x]|].
     * @param x A real number
     * @return The function value
     */
    inline double lngamma( double const x ){ return gsl_sf_lngamma( x ); } 
    /**
     * C++ version of gsl_sf_lngamma_sgn_e().
     * Log[Gamma(x)], x not a negative integer
     * Uses real Lanczos method.
     * Returns the real part of Log[Gamma[x]] when x < 0,
     * @param x A real number
     * @param result_lg Result as a @c gsl::sf::result object
     * @param sgn The sign as a return value
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EROUND
     */
    inline int lngamma_sgn_e( double x, result& result_lg, double* sgn ){
      return gsl_sf_lngamma_sgn_e( x, &result_lg, sgn ); } 
    /**
     * C++ version of gsl_sf_gamma_e().
     * Gamma(x), x not a negative integer
     * Uses real Lanczos method.
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EROUND
     */
    inline int gamma_e( double const x, result& result ){
      return gsl_sf_gamma_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_gamma().
     * Gamma(x), x not a negative integer
     * Uses real Lanczos method.
     * @param x A real number
     * @return The function value
     */
    inline double gamma( double const x ){ return gsl_sf_gamma( x ); } 
    /**
     * C++ version of gsl_sf_gammastar_e().
     * Regulated Gamma Function, x > 0
     * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
     *            = (1 + 1/(12x) + ...),  x->Inf
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gammastar_e( double const x, result& result ){
      return gsl_sf_gammastar_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_gammastar().
     * Regulated Gamma Function, x > 0
     * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
     *            = (1 + 1/(12x) + ...),  x->Inf
     * @param x A real number
     * @return The function value
     */
    inline double gammastar( double const x ){ return gsl_sf_gammastar( x ); } 
    /**
     * C++ version of gsl_sf_gammainv_e().
     * 1/Gamma(x)
     * Uses real Lanczos method.
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EROUND
     */
    inline int gammainv_e( double const x, result& result ){
      return gsl_sf_gammainv_e( x, &result ); } 
    /**
     * C++ version of gsl_sf_gammainv().
     * 1/Gamma(x)
     * Uses real Lanczos method.
     * @param x A real number
     * @return The function value
     */
    inline double gammainv( double const x ){ return gsl_sf_gammainv( x ); } 
    /**
     * C++ version of gsl_sf_lngamma_complex_e().
     * Log[Gamma(z)] for z complex, z not a negative integer
     * Uses complex Lanczos method. Note that the phase part (arg)
     * is not well-determined when |z| is very large, due
     * to inevitable roundoff in restricting to (-Pi,Pi].
     * This will raise the GSL_ELOSS exception when it occurs.
     * The absolute value part (lnr), however, never suffers.
     *
     * Calculates:
     *   lnr = log|Gamma(z)|
     *   arg = arg(Gamma(z))  in (-Pi, Pi]
     * @param zr The real part
     * @param zi The imaginary part
     * @param lnr Result as a @c gsl::sf::result object
     * @param arg Result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_ELOSS
     */
    inline int lngamma_complex_e( double zr, double zi, result& lnr, result& arg ){
      return gsl_sf_lngamma_complex_e( zr, zi, &lnr, &arg ); } 
    /**
     * C++ version of gsl_sf_taylorcoeff_e().
     * x^n / n!
     *
     * x >= 0.0, n >= 0
     * @param n An integer
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW or GSL_EUNDRFLW
     */
    inline int taylorcoeff_e( int const n, double const x, result& result ){
      return gsl_sf_taylorcoeff_e( n, x, &result ); } 
    /**
     * C++ version of gsl_sf_taylorcoeff().
     * x^n / n!
     *
     * x >= 0.0, n >= 0
     * @param n An integer
     * @param x A real number
     * @return The function value
     */
    inline double taylorcoeff( int const n, double const x ){
      return gsl_sf_taylorcoeff( n, x ); } 
    /**
     * C++ version of gsl_sf_fact_e().
     * n!
     * @param n An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
     */
    inline int fact_e( unsigned int const n, result& result ){
      return gsl_sf_fact_e( n, &result ); } 
    /**
     * C++ version of gsl_sf_fact().
     * n!
     * @param n An integer
     * @return The function value
     */
    inline double fact( unsigned int const n ){ return gsl_sf_fact( n ); } 
    /**
     * C++ version of gsl_sf_doublefact_e().
     * n!! = n(n-2)(n-4) ... 
     * @param n An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
     */
    inline int doublefact_e( unsigned int const n, result& result ){
      return gsl_sf_doublefact_e( n, &result ); } 
    /**
     * C++ version of gsl_sf_doublefact().
     * n!! = n(n-2)(n-4) ... 
     * @param n An integer
     * @return The function value
     */
    inline double doublefact( unsigned int const n ){ return gsl_sf_doublefact( n ); } 
    /**
     * C++ version of gsl_sf_lnfact_e().
     * log(n!) 
     * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
     * @param n An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int lnfact_e( unsigned int const n, result& result ){
      return gsl_sf_lnfact_e( n, &result ); } 
    /**
     * C++ version of gsl_sf_lnfact().
     * log(n!) 
     * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
     * @param n An integer
     * @return The function value
     */
    inline double lnfact( unsigned int const n ){ return gsl_sf_lnfact( n ); } 
    /**
     * C++ version of gsl_sf_lndoublefact_e().
     * log(n!!) 
     * @param n An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS
     */
    inline int lndoublefact_e( unsigned int const n, result& result ){
      return gsl_sf_lndoublefact_e( n, &result ); } 
    /**
     * C++ version of gsl_sf_lndoublefact().
     * log(n!!) 
     * @param n An integer
     * @return The function value
     */
    inline double lndoublefact( unsigned int const n ){
      return gsl_sf_lndoublefact( n ); } 
    /**
     * C++ version of gsl_sf_lnchoose_e().
     * log(n choose m)
     * @param n An integer
     * @param m An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int lnchoose_e( unsigned int n, unsigned int m, result& result ){
      return gsl_sf_lnchoose_e( n, m, &result ); } 
    /**
     * C++ version of gsl_sf_lnchoose().
     * log(n choose m)
     * @param n An integer
     * @param m An integer
     * @return The function value
     */
    inline double lnchoose( unsigned int n, unsigned int m ){
      return gsl_sf_lnchoose( n, m ); } 
    /**
     * C++ version of gsl_sf_choose_e().
     * n choose m
     * @param n An integer
     * @param m An integer
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
     */
    inline int choose_e( unsigned int n, unsigned int m, result& result ){
      return gsl_sf_choose_e( n, m, &result ); } 
    /**
     * C++ version of gsl_sf_choose().
     * n choose m
     * @param n An integer
     * @param m An integer
     * @return The function value
     */
    inline double choose( unsigned int n, unsigned int m ){ return gsl_sf_choose( n, m ); } 
    /**
     * C++ version of gsl_sf_lnpoch_e().
     * Logarithm of Pochhammer (Apell) symbol
     *   log( (a)_x )
     *   where (a)_x := Gamma[a + x]/Gamma[a]
     *
     * a > 0, a+x > 0
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int lnpoch_e( double const a, double const x, result& result ){
      return gsl_sf_lnpoch_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_lnpoch().
     * Logarithm of Pochhammer (Apell) symbol
     *   log( (a)_x )
     *   where (a)_x := Gamma[a + x]/Gamma[a]
     *
     * a > 0, a+x > 0
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double lnpoch( double const a, double const x ){
      return gsl_sf_lnpoch( a, x ); } 
    /**
     * C++ version of gsl_sf_lnpoch_sgn_e().
     * Logarithm of Pochhammer (Apell) symbol, with sign information.
     *   result = log( |(a)_x| )
     *   sgn    = sgn( (a)_x )
     *   where (a)_x := Gamma[a + x]/Gamma[a]
     *
     * a != neg integer, a+x != neg integer
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @param sgn Record the sign here
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int lnpoch_sgn_e( double const a, double const x, result& result, double* sgn ){
      return gsl_sf_lnpoch_sgn_e( a, x, &result, sgn ); } 
    /**
     * C++ version of gsl_sf_poch_e().
     * Pochhammer (Apell) symbol
     *   (a)_x := Gamma[a + x]/Gamma[x]
     *
     * a != neg integer, a+x != neg integer
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
     */
    inline int poch_e( double const a, double const x, result& result ){
      return gsl_sf_poch_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_poch().
     * Pochhammer (Apell) symbol
     *   (a)_x := Gamma[a + x]/Gamma[x]
     *
     * a != neg integer, a+x != neg integer
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double poch( double const a, double const x ){ return gsl_sf_poch( a, x ); } 
    /**
     * C++ version of gsl_sf_pochrel_e().
     * Relative Pochhammer (Apell) symbol
     *   ((a,x) - 1)/x
     *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int pochrel_e( double const a, double const x, result& result ){
      return gsl_sf_pochrel_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_pochrel().
     * Relative Pochhammer (Apell) symbol
     *   ((a,x) - 1)/x
     *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double pochrel( double const a, double const x ){ return gsl_sf_pochrel( a, x ); } 
    /**
     * C++ version of gsl_sf_gamma_inc_Q_e().
     * Normalized Incomplete Gamma Function
     *
     * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
     *
     * a >= 0, x >= 0
     *   Q(a,0) := 1
     *   Q(0,x) := 0, x != 0
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gamma_inc_Q_e( double const a, double const x, result& result ){
      return gsl_sf_gamma_inc_Q_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_gamma_inc_Q().
     * Normalized Incomplete Gamma Function
     *
     * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
     *
     * a >= 0, x >= 0
     *   Q(a,0) := 1
     *   Q(0,x) := 0, x != 0
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double gamma_inc_Q( double const a, double const x ){
      return gsl_sf_gamma_inc_Q( a, x ); } 
    /**
     * C++ version of gsl_sf_gamma_inc_P_e().
     * Complementary Normalized Incomplete Gamma Function
     *
     * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
     *
     * a > 0, x >= 0
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gamma_inc_P_e( double const a, double const x, result& result ){
      return gsl_sf_gamma_inc_P_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_gamma_inc_P().
     * Complementary Normalized Incomplete Gamma Function
     *
     * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
     *
     * a > 0, x >= 0
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double gamma_inc_P( double const a, double const x ){
      return gsl_sf_gamma_inc_P( a, x ); } 
    /**
     * C++ version of gsl_sf_gamma_inc_e().
     * Non-normalized Incomplete Gamma Function
     *
     * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
     *
     * x >= 0.0
     *   Gamma(a, 0) := Gamma(a)
     * @param a A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int gamma_inc_e( double const a, double const x, result& result ){
      return gsl_sf_gamma_inc_e( a, x, &result ); } 
    /**
     * C++ version of gsl_sf_gamma_inc().
     * Non-normalized Incomplete Gamma Function
     *
     * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
     *
     * x >= 0.0
     *   Gamma(a, 0) := Gamma(a)
     * @param a A real number
     * @param x A real number
     * @return The function value
     */
    inline double gamma_inc( double const a, double const x ){
      return gsl_sf_gamma_inc( a, x ); } 
    /**
     * C++ version of gsl_sf_lnbeta_e().
     * Logarithm of Beta Function
     * Log[B(a,b)]
     *
     * a > 0, b > 0
     * @param a A real number
     * @param b A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int lnbeta_e( double const a, double const b, result& result ){
      return gsl_sf_lnbeta_e( a, b, &result ); } 
    /**
     * C++ version of gsl_sf_lnbeta().
     * Logarithm of Beta Function
     * Log[B(a,b)]
     *
     * a > 0, b > 0
     * @param a A real number
     * @param b A real number
     * @return The function value
     */
    inline double lnbeta( double const a, double const b ){
      return gsl_sf_lnbeta( a, b ); } 
    /**
     * C++ version of gsl_sf_lnbeta_sgn_e().
     * Logarithm of Beta Function
     * Log[B(a,b)]
     *
     * a > 0, b > 0
     * @param x A real number
     * @param y A real number
     * @param result The result as a @c gsl::sf::result object
     * @param sgn Record the sign here
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int lnbeta_sgn_e( double const x, double const y, result& result, double* sgn ){
      return gsl_sf_lnbeta_sgn_e( x, y, &result, sgn ); } 
    /**
     * C++ version of gsl_sf_beta_e().
     * Beta Function
     * B(a,b)
     *
     * a > 0, b > 0
     * @param a A real number
     * @param b A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EOVRFLW or GSL_EUNDRFLW or GSL_EDOM
     */
    inline int beta_e( double const a, double const b, result& result ){
      return gsl_sf_beta_e( a, b, &result ); } 
    /**
     * C++ version of gsl_sf_beta().
     * Beta Function
     * B(a,b)
     *
     * a > 0, b > 0
     * @param a A real number
     * @param b A real number
     * @return The function value
     */
    inline double beta( double const a, double const b ){ return gsl_sf_beta( a, b ); } 
    /**
     * C++ version of gsl_sf_beta_inc_e().
     * Normalized Incomplete Beta Function
     * B_x(a,b)/B(a,b)
     *
     * a > 0, b > 0, 0 <= x <= 1
     * @param a A real number
     * @param b A real number
     * @param x A real number
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
     */
    inline int beta_inc_e( double const a, double const b, double const x, result& result ){
      return gsl_sf_beta_inc_e( a, b, x, &result ); } 
    /**
     * C++ version of gsl_sf_beta_inc().
     * Normalized Incomplete Beta Function
     * B_x(a,b)/B(a,b)
     *
     * a > 0, b > 0, 0 <= x <= 1
     * @param a A real number
     * @param b A real number
     * @param x A real number
     * @return The function value
     */
    inline double beta_inc( double const a, double const b, double const x ){
      return gsl_sf_beta_inc( a, b, x ); }
    /**
     * The maximum x such that gamma(x) is not
     * considered an overflow.
     */
    double const GAMMA_XMAX = GSL_SF_GAMMA_XMAX;
    /**
     * The maximum n such that fact(n) does not give an overflow.
     */
    double const FACT_NMAX = GSL_SF_FACT_NMAX;
    /**
     * The maximum n such that doublefact(n) does not give an overflow.
     */
    double const DOUBLEFACT_NMAX = GSL_SF_DOUBLEFACT_NMAX;
  }
}

#endif
