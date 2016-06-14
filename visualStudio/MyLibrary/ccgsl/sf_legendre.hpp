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

#ifndef CCGSL_SF_LEGENDRE_HPP
#define CCGSL_SF_LEGENDRE_HPP

#include<gsl/gsl_sf_legendre.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * Legendre functions.
     */
    namespace legendre {
      /**
       * C++ version of gsl_sf_legendre_Pl_e().
       * P_l(x)   l >= 0; |x| <= 1
       * @param l An integer
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Pl_e( int const l, double const x, result& result ){
	return gsl_sf_legendre_Pl_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_Pl().
       * P_l(x)   l >= 0; |x| <= 1
       * @param l An integer
       * @param x A real value
       * @return The function value
       */
      inline double Pl( int const l, double const x ){
	return gsl_sf_legendre_Pl( l, x ); } 
      /**
       * C++ version of gsl_sf_legendre_Pl_array().
       * P_l(x) for l=0,...,lmax; |x| <= 1
       * @param lmax An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Pl_array( int const lmax, double const x, double* result_array ){
	return gsl_sf_legendre_Pl_array( lmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_legendre_Pl_deriv_array().
       * P_l(x) and P_l'(x) for l=0,...,lmax; |x| <= 1
       * @param lmax An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @param result_deriv_array An array of size @c lmax
       * @return The function value
       */
      inline int Pl_deriv_array( int const lmax, double const x,
					  double* result_array, double* result_deriv_array ){
	return gsl_sf_legendre_Pl_deriv_array( lmax, x, result_array, result_deriv_array ); } 
      /**
       * C++ version of gsl_sf_legendre_P1_e().
       * P_l(x), l=1
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int P1_e( double x, result& result ){
	return gsl_sf_legendre_P1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_P2_e().
       * P_l(x), l=2
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int P2_e( double x, result& result ){
	return gsl_sf_legendre_P2_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_P3_e().
       * P_l(x), l=3
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int P3_e( double x, result& result ){
	return gsl_sf_legendre_P3_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_P1().
       * P_l(x), l=1
       * @param x A real value
       * @return The function value
       */
      inline double P1( double const x ){ return gsl_sf_legendre_P1( x ); } 
      /**
       * C++ version of gsl_sf_legendre_P2().
       * P_l(x), l=2
       * @param x A real value
       * @return The function value
       */
      inline double P2( double const x ){ return gsl_sf_legendre_P2( x ); } 
      /**
       * C++ version of gsl_sf_legendre_P3().
       * P_l(x), l=3
       * @param x A real value
       * @return The function value
       */
      inline double P3( double const x ){ return gsl_sf_legendre_P3( x ); } 
      /**
       * C++ version of gsl_sf_legendre_Q0_e().
       * Q_0(x), x > -1, x != 1
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Q0_e( double const x, result& result ){
	return gsl_sf_legendre_Q0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_Q0().
       * Q_0(x), x > -1, x != 1
       * @param x A real value
       * @return The function value
       */
      inline double Q0( double const x ){ return gsl_sf_legendre_Q0( x ); } 
      /**
       * C++ version of gsl_sf_legendre_Q1_e().
       * Q_1(x), x > -1, x != 1
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Q1_e( double const x, result& result ){
	return gsl_sf_legendre_Q1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_Q1().
       * Q_1(x), x > -1, x != 1
       * @param x A real value
       * @return The function value
       */
      inline double Q1( double const x ){ return gsl_sf_legendre_Q1( x ); } 
      /**
       * C++ version of gsl_sf_legendre_Ql_e().
       * Q_l(x), x > -1, x != 1, l >= 0
       * @param l An integer
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Ql_e( int const l, double const x, result& result ){
	return gsl_sf_legendre_Ql_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_Ql().
       * Q_l(x), x > -1, x != 1, l >= 0
       * @param l An integer
       * @param x A real value
       * @return The function value
       */
      inline double Ql( int const l, double const x ){ return gsl_sf_legendre_Ql( l, x ); }
      /**
       * C++ version of gsl_sf_legendre_Plm_e().
       * P_l^m(x)  m >= 0; l >= m; |x| <= 1.0
       *
       * Note that this function grows combinatorially with l.
       * Therefore we can easily generate an overflow for l larger
       * than about 150.
       *
       * There is no trouble for small m, but when m and l are both large,
       * then there will be trouble. Rather than allow overflows, these
       * functions refuse to calculate when they can sense that l and m are
       * too big.
       *
       * If you really want to calculate a spherical harmonic, then DO NOT
       * use this. Instead use legendre_sphPlm() below, which  uses a similar
       * recursion, but with the normalized functions.
       * @param l An integer
       * @param m An integer
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
       */
      inline int Plm_e( int const l, int const m, double const x, result& result ){
	return gsl_sf_legendre_Plm_e( l, m, x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_Plm().
       * P_l^m(x)  m >= 0; l >= m; |x| <= 1.0
       *
       * Note that this function grows combinatorially with l.
       * Therefore we can easily generate an overflow for l larger
       * than about 150.
       *
       * There is no trouble for small m, but when m and l are both large,
       * then there will be trouble. Rather than allow overflows, these
       * functions refuse to calculate when they can sense that l and m are
       * too big.
       *
       * If you really want to calculate a spherical harmonic, then DO NOT
       * use this. Instead use legendre_sphPlm() below, which  uses a similar
       * recursion, but with the normalized functions.
       * @param l An integer
       * @param m An integer
       * @param x A real value
       * @return The function value
       */
      inline double Plm( int const l, int const m, double const x ){
	return gsl_sf_legendre_Plm( l, m, x ); } 
      /**
       * C++ version of gsl_sf_legendre_Plm_array().
       * P_l^m(x)  m >= 0; l >= m; |x| <= 1.0
       * l=|m|,...,lmax
       * @param lmax An integer
       * @param m An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
       */
      inline int Plm_array( int const lmax, int const m, double const x, double* result_array ){
	return gsl_sf_legendre_Plm_array( lmax, m, x, result_array ); } 
      /**
       * C++ version of gsl_sf_legendre_Plm_deriv_array().
       * P_l^m(x)  and d(P_l^m(x))/dx;  m >= 0; lmax >= m; |x| <= 1.0
       * l=|m|,...,lmax
       * @param lmax An integer
       * @param m An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @param result_deriv_array An array of size @c lmax
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EOVRFLW
       */
      inline int Plm_deriv_array( int const lmax, int const m, double const x,
					   double* result_array, double* result_deriv_array ){
	return gsl_sf_legendre_Plm_deriv_array( lmax, m, x, result_array, result_deriv_array ); } 
      /**
       * C++ version of gsl_sf_legendre_sphPlm_e().
       * P_l^m(x), normalized properly for use in spherical harmonics
       * m >= 0; l >= m; |x| <= 1.0
       *
       * There is no overflow problem, as there is for the
       * standard normalization of P_l^m(x).
       *
       * Specifically, it returns:
       *
       *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
       * @param l An integer
       * @param m An integer
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int sphPlm_e( int const l, int m, double const x, result& result ){
	return gsl_sf_legendre_sphPlm_e( l, m, x, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_sphPlm().
       * P_l^m(x), normalized properly for use in spherical harmonics
       * m >= 0; l >= m; |x| <= 1.0
       *
       * There is no overflow problem, as there is for the
       * standard normalization of P_l^m(x).
       *
       * Specifically, it returns:
       *
       *        sqrt((2l+1)/(4pi)) sqrt((l-m)!/(l+m)!) P_l^m(x)
       * @param l An integer
       * @param m An integer
       * @param x A real value
       * @return The function value
       */
      inline double sphPlm( int const l, int const m, double const x ){
	return gsl_sf_legendre_sphPlm( l, m, x ); } 
      /**
       * C++ version of gsl_sf_legendre_sphPlm_array().
       * sphPlm(l,m,x) values
       * m >= 0; l >= m; |x| <= 1.0
       * l=|m|,...,lmax
       * @param lmax An integer
       * @param m An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int sphPlm_array( int const lmax, int m, double const x, double* result_array ){
	return gsl_sf_legendre_sphPlm_array( lmax, m, x, result_array ); } 
      /**
       * C++ version of gsl_sf_legendre_sphPlm_deriv_array().
       * sphPlm(l,m,x) and d(sphPlm(l,m,x))/dx values
       * m >= 0; l >= m; |x| <= 1.0
       * l=|m|,...,lmax
       * @param lmax An integer
       * @param m An integer
       * @param x A real value
       * @param result_array An array of size @c lmax
       * @param result_deriv_array An array of size @c lmax
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int sphPlm_deriv_array( int const lmax, int const m, double const x,
					      double* result_array, double* result_deriv_array ){
	return gsl_sf_legendre_sphPlm_deriv_array( lmax, m, x, result_array, result_deriv_array ); } 
      /**
       * C++ version of gsl_sf_legendre_array_size().
       * size of result_array[] needed for the array versions of Plm
       * (lmax - m + 1)
       * @param lmax An integer
       * @param m An integer
       * @return GSL_SUCCESS
       */
      inline int array_size( int const lmax, int const m ){
	return gsl_sf_legendre_array_size( lmax, m ); } 
      /**
       * C++ version of gsl_sf_conicalP_half_e().
       * Irregular Spherical Conical Function
       * P^{1/2}_{-1/2 + I lambda}(x)
       *
       * x > -1.0
       * @param lambda A real value
       * @param x A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
    }
    inline int conicalP_half_e( double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_half_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_half().
       * Irregular Spherical Conical Function
       * P^{1/2}_{-1/2 + I lambda}(x)
       *
       * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_half( double const lambda, double const x ){
      return gsl_sf_conicalP_half( lambda, x ); } 
    /**
     * C++ version of gsl_sf_conicalP_mhalf_e().
     * Regular Spherical Conical Function
     * P^{-1/2}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int conicalP_mhalf_e( double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_mhalf_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_mhalf().
     * Regular Spherical Conical Function
     * P^{-1/2}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_mhalf( double const lambda, double const x ){
      return gsl_sf_conicalP_mhalf( lambda, x ); } 
    /**
     * C++ version of gsl_sf_conicalP_0_e().
     * Conical Function
     * P^{0}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int conicalP_0_e( double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_0_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_0().
     * Conical Function
     * P^{0}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_0( double const lambda, double const x ){
      return gsl_sf_conicalP_0( lambda, x ); } 
    /**
     * C++ version of gsl_sf_conicalP_1_e().
     * Conical Function
     * P^{1}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int conicalP_1_e( double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_1_e( lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_1().
     * Conical Function
     * P^{1}_{-1/2 + I lambda}(x)
     *
     * x > -1.0
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_1( double const lambda, double const x ){
      return gsl_sf_conicalP_1( lambda, x ); } 
    /**
     * C++ version of gsl_sf_conicalP_sph_reg_e().
     * Regular Spherical Conical Function
     * P^{-1/2-l}_{-1/2 + I lambda}(x)
     *
     * x > -1.0, l >= -1
     * @param l An integer
     * @param lambda A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int conicalP_sph_reg_e( int const l, double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_sph_reg_e( l, lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_sph_reg().
     * Regular Spherical Conical Function
     * P^{-1/2-l}_{-1/2 + I lambda}(x)
     *
     * x > -1.0, l >= -1
     * @param l An integer
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_sph_reg( int const l, double const lambda, double const x ){
      return gsl_sf_conicalP_sph_reg( l, lambda, x ); } 
    /**
     * C++ version of gsl_sf_conicalP_cyl_reg_e().
     * Regular Cylindrical Conical Function
     * P^{-m}_{-1/2 + I lambda}(x)
     *
     * x > -1.0, m >= -1
     * @param m An integer
     * @param lambda A real value
     * @param x A real value
     * @param result The result as a @c gsl::sf::result object
     * @return GSL_SUCCESS or GSL_EDOM
     */
    inline int conicalP_cyl_reg_e( int const m, double const lambda, double const x, result& result ){
      return gsl_sf_conicalP_cyl_reg_e( m, lambda, x, &result ); } 
    /**
     * C++ version of gsl_sf_conicalP_cyl_reg().
     * Regular Cylindrical Conical Function
     * P^{-m}_{-1/2 + I lambda}(x)
     *
     * x > -1.0, m >= -1
     * @param m An integer
     * @param lambda A real value
     * @param x A real value
     * @return The function value
     */
    inline double conicalP_cyl_reg( int const m, double const lambda, double const x ){
      return gsl_sf_conicalP_cyl_reg( m, lambda, x ); } 
    namespace legendre {
      /**
       * C++ version of gsl_sf_legendre_H3d_0_e().
       * Zeroth radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * legendre_H3d_0(lambda,eta) := sin(lambda*eta)/(lambda*sinh(eta))
       * 
       * Normalization:
       * Flat-Lim legendre_H3d_0(lambda,eta) = j_0(lambda*eta)
       *
       * eta >= 0.0
       * @param lambda A real value
       * @param eta A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int H3d_0_e( double const lambda, double const eta, result& result ){
	return gsl_sf_legendre_H3d_0_e( lambda, eta, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d_0().
       * Zeroth radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * legendre_H3d_0(lambda,eta) := sin(lambda*eta)/(lambda*sinh(eta))
       * 
       * Normalization:
       * Flat-Lim legendre_H3d_0(lambda,eta) = j_0(lambda*eta)
       *
       * eta >= 0.0
       * @param lambda A real value
       * @param eta A real value
       * @return The function value
       */
      inline double H3d_0( double const lambda, double const eta ){
	return gsl_sf_legendre_H3d_0( lambda, eta ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d_1_e().
       * First radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * legendre_H3d_1(lambda,eta) :=
       *    1/sqrt(lambda^2 + 1) sin(lam eta)/(lam sinh(eta))
       *    (coth(eta) - lambda cot(lambda*eta))
       * 
       * Normalization:
       * Flat-Lim legendre_H3d_1(lambda,eta) = j_1(lambda*eta)
       *
       * eta >= 0.0
       * @param lambda A real value
       * @param eta A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int H3d_1_e( double const lambda, double const eta, result& result ){
	return gsl_sf_legendre_H3d_1_e( lambda, eta, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d_1().
       * First radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * legendre_H3d_1(lambda,eta) :=
       *    1/sqrt(lambda^2 + 1) sin(lam eta)/(lam sinh(eta))
       *    (coth(eta) - lambda cot(lambda*eta))
       * 
       * Normalization:
       * Flat-Lim legendre_H3d_1(lambda,eta) = j_1(lambda*eta)
       *
       * eta >= 0.0
       * @param lambda A real value
       * @param eta A real value
       * @return The function value
       */
      inline double H3d_1( double const lambda, double const eta ){
	return gsl_sf_legendre_H3d_1( lambda, eta ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d_e().
       * lth radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * Normalization:
       * Flat-Lim legendre_H3d_l(l,lambda,eta) = j_l(lambda*eta)
       *
       * eta >= 0.0, l >= 0
       * @param l An integer
       * @param lambda A real value
       * @param eta A real value
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int H3d_e( int const l, double const lambda, double const eta, result& result ){
	return gsl_sf_legendre_H3d_e( l, lambda, eta, &result ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d().
      * lth radial eigenfunction of the Laplacian on the
       * 3-dimensional hyperbolic space.
       *
       * Normalization:
       * Flat-Lim legendre_H3d_l(l,lambda,eta) = j_l(lambda*eta)
       *
       * eta >= 0.0, l >= 0
       * @param l An integer
       * @param lambda A real value
       * @param eta A real value
       * @return The function value
       */
      inline double H3d( int const l, double const lambda, double const eta ){
	return gsl_sf_legendre_H3d( l, lambda, eta ); } 
      /**
       * C++ version of gsl_sf_legendre_H3d_array().
       * Array of H3d(ell),  0 <= ell <= lmax
       * @param lmax An integer
       * @param lambda A real value
       * @param eta A real value
       * @param result_array An array of size @c lmax
       * @return GSL_SUCCESS
       */
      inline int H3d_array( int const lmax, double const lambda, double const eta,
				     double* result_array ){
	return gsl_sf_legendre_H3d_array( lmax, lambda, eta, result_array ); }
    }
  }
}

#endif
