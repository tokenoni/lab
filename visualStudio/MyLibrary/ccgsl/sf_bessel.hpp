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
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef CCGSL_SF_BESSEL_HPP
#define CCGSL_SF_BESSEL_HPP

#include<gsl/gsl_sf_bessel.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * Namespace for Bessel functions.
     */
    namespace bessel {
      /**
       * C++ version of gsl_sf_bessel_J0().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int J0_e( double const x, result& result ){ return gsl_sf_bessel_J0_e( x, &result ); }  
      /**
       * C++ version of gsl_sf_bessel_J0().
       * @param x A real number
       * @return Value of function at @c x
       */
      inline double J0( double const x ){ return gsl_sf_bessel_J0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_J1_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int J1_e( double const x, result& result ){ return gsl_sf_bessel_J1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_J1().
       * @param x A real number
       * @return Value of function at @c x
       */
      inline double J1( double const x ){ return gsl_sf_bessel_J1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_Jn_e().
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int Jn_e( int n, double x, result& result ){ return gsl_sf_bessel_Jn_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Jn().
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double Jn( int const n, double const x ){ return gsl_sf_bessel_Jn( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Jn_array().
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM
       */
      inline int Jn_array( int nmin, int nmax, double x, double* result_array ){
	return gsl_sf_bessel_Jn_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_Y0_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM 
       */
      inline int Y0_e( double const x, result& result ){ return gsl_sf_bessel_Y0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Y0().
       * @param x A real number
       * @return Value of function
       */
      inline double Y0( double const x ){ return gsl_sf_bessel_Y0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_Y1_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM or GSL_EOVERFLW
       */
      inline int Y1_e( double const x, result& result ){ return gsl_sf_bessel_Y1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Y1().
       * @param x A real number
       * @return Value of function
       */
      inline double Y1( double const x ){ return gsl_sf_bessel_Y1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_Yn_e().
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM or GSL_EOVERFLW
       */
      inline int Yn_e( int n, double const x, result& result ){ return gsl_sf_bessel_Yn_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Yn().
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double Yn( int const n, double const x ){ return gsl_sf_bessel_Yn( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Yn_array().
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM or GSL_EOVERFLW
       */
      inline int Yn_array( int const nmin, int const nmax, double const x, double* result_array ){
	return gsl_sf_bessel_Yn_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_I0_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EOVERFLW
       */
      inline int I0_e( double const x, result& result ){ return gsl_sf_bessel_I0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_I0().
       * @param x A real number
       * @return Value of function
       */
      inline double I0( double const x ){ return gsl_sf_bessel_I0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_I1_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EOVERFLW or GSL_EUNDRFLW
       */
      inline int I1_e( double const x, result& result ){ return gsl_sf_bessel_I1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_I1().
       * @param x A real number
       * @return Value of function
       */
      inline double I1( double const x ){ return gsl_sf_bessel_I1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_In_e().
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EOVERFLW or GSL_EUNDRFLW
       */
      inline int In_e( int const n, double const x, result& result ){ return gsl_sf_bessel_In_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_In().
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double In( int const n, double const x ){ return gsl_sf_bessel_In( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_In_array().
       * nmin >=0, nmax >= nmin
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EOVERFLW or GSL_EUNDRFLW
       */
      inline int In_array( int const nmin, int const nmax, double const x, double* result_array ){
	return gsl_sf_bessel_In_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_I0_scaled_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM or GSL_EOVERFLW
       */
      inline int I0_scaled_e( double const x, result& result ){ return gsl_sf_bessel_I0_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_I0_scaled().
       *  exp(-|x|) I_0(x)
       * @param x A real number
       * @return Value of function
       */
      inline double I0_scaled( double const x ){ return gsl_sf_bessel_I0_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_I1_scaled_e().
       *  exp(-|x|) I_0(x)
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int I1_scaled_e( double const x, result& result ){ return gsl_sf_bessel_I1_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_I1_scaled().
       *  exp(-|x|) I_1(x)
       * @param x A real number
       * @return Value of function
       */
      inline double I1_scaled( double const x ){ return gsl_sf_bessel_I1_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_In_scaled_e().
       *  exp(-|x|) I_1(x)
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int In_scaled_e( int n, double const x, result& result ){
	return gsl_sf_bessel_In_scaled_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_In_scaled().
       *  exp(-|x|) I_n(x)
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double In_scaled( int const n, double const x ){ return gsl_sf_bessel_In_scaled( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_In_scaled_array().
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or
       */
      inline int In_scaled_array( int const nmin, int const nmax, double const x, double* result_array ){
	return gsl_sf_bessel_In_scaled_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_K0_e().
       *  exp(-|x|) I_n(x)  for n=nmin,...,nmax
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int K0_e( double const x, result& result ){ return gsl_sf_bessel_K0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_K0().
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double K0( double const x ){ return gsl_sf_bessel_K0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_K1_e().
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int K1_e( double const x, result& result ){ return gsl_sf_bessel_K1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_K1().
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double K1( double const x ){ return gsl_sf_bessel_K1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn_e().
       * x > 0.0
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int Kn_e( int const n, double const x, result& result ){
	return gsl_sf_bessel_Kn_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn().
       * x > 0.0
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double Kn( int const n, double const x ){ return gsl_sf_bessel_Kn( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn_array().
       * x > 0.0, nmin >=0, nmax >= nmin
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int Kn_array( int const nmin, int const nmax, double const x, double* result_array ){
	return gsl_sf_bessel_Kn_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_K0_scaled_e().
       *  exp(x) K_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int K0_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_K0_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_K0_scaled().
       *  exp(x) K_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double K0_scaled( double const x ){ return gsl_sf_bessel_K0_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_K1_scaled_e().
       *  exp(x) K_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int K1_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_K1_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_K1_scaled().
       *  exp(x) K_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double K1_scaled( double const x ){ return gsl_sf_bessel_K1_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn_scaled_e().
       *  exp(x) K_n(x)
       *
       * x > 0.0
       * @param n The order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int Kn_scaled_e( int n, double const x, result& result ){
	return gsl_sf_bessel_Kn_scaled_e( n, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn_scaled().
       *  exp(x) K_n(x)
       *
       * x > 0.0
       * @param n The order
       * @param x A real number
       * @return Value of function
       */
      inline double Kn_scaled( int const n, double const x ){ return gsl_sf_bessel_Kn_scaled( n, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Kn_scaled_array().
       * Scaled irregular modified Bessel function  exp(x) K_n(x)  for n=nmin,...,nmax
       *
       * x > 0.0, nmin >=0, nmax >= nmin
       * @param nmin The lower value
       * @param nmax The upper value
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int Kn_scaled_array( int const nmin, int const nmax, double const x, double* result_array ){
	return gsl_sf_bessel_Kn_scaled_array( nmin, nmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_j0_e().
       * Regular spherical Bessel function j_0(x) = sin(x)/x
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int j0_e( double const x, result& result ){ return gsl_sf_bessel_j0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_j0().
       * Regular spherical Bessel function j_0(x) = sin(x)/x
       * @param x A real number
       * @return Value of function
       */
      inline double j0( double const x ){ return gsl_sf_bessel_j0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_j1_e().
       * Regular spherical Bessel function j_1(x) = (sin(x)/x - cos(x))/x
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int j1_e( double const x, result& result ){ return gsl_sf_bessel_j1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_j1().
       * Regular spherical Bessel function j_1(x) = (sin(x)/x - cos(x))/x
       * @param x A real number
       * @return Value of function
       */
      inline double j1( double const x ){ return gsl_sf_bessel_j1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_j2_e().
       * Regular spherical Bessel function j_2(x) = ((3/x^2 - 1)sin(x) - 3cos(x)/x)/x
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW or GSL_EDOM
       */
      inline int j2_e( double const x, result& result ){ return gsl_sf_bessel_j2_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_j2().
       * Regular spherical Bessel function j_2(x) = ((3/x^2 - 1)sin(x) - 3cos(x)/x)/x
       * @param x A real number
       * @return Value of function
       */
      inline double j2( double const x ){ return gsl_sf_bessel_j2( x ); } 
      /**
       * C++ version of gsl_sf_bessel_jl_e().
       * Regular spherical Bessel function j_l(x)
       *
       * l >= 0, x >= 0.0
       * @param l A nonnegative integer
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int jl_e( int const l, double const x, result& result ){
	return gsl_sf_bessel_jl_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_jl().
       * Regular spherical Bessel function j_l(x)
       *
       * l >= 0, x >= 0.0
       * @param l A nonnegative integer
       * @param x A real number
       * @return Value of function
       */
      inline double jl( int const l, double const x ){ return gsl_sf_bessel_jl( l, x ); } 
      /**
       * C++ version of gsl_sf_bessel_jl_array().
       * Regular spherical Bessel function j_l(x)
       *
       * l >= 0, x >= 0.0
       * @param lmax Maximum value of l
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int jl_array( int const lmax, double const x, double* result_array ){
	return gsl_sf_bessel_jl_array( lmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_jl_steed_array().
       * Regular spherical Bessel function j_l(x) for l=0,1,...,lmax
       * Uses Steed's method.
       * @param lmax Maximum value of l
       * @param x A real number
       * @param jl_x_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int jl_steed_array( int const lmax, double const x, double* jl_x_array ){
	return gsl_sf_bessel_jl_steed_array( lmax, x, jl_x_array ); } 
      /**
       * C++ version of gsl_sf_bessel_y0_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int y0_e( double const x, result& result ){
	return gsl_sf_bessel_y0_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_y0().
       * @param x A real number
       * @return Value of function
       */
      inline double y0( double const x ){ return gsl_sf_bessel_y0( x ); } 
      /**
       * C++ version of gsl_sf_bessel_y1_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int y1_e( double const x, result& result ){ return gsl_sf_bessel_y1_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_y1().
       * @param x A real number
       * @return Value of function
       */
      inline double y1( double const x ){ return gsl_sf_bessel_y1( x ); } 
      /**
       * C++ version of gsl_sf_bessel_y2_e().
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or  GSL_EUNDRFLW
       */
      inline int y2_e( double const x, result& result ){ return gsl_sf_bessel_y2_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_y2().
       * @param x A real number
       * @return Value of function
       */
      inline double y2( double const x ){ return gsl_sf_bessel_y2( x ); } 
      /**
       * C++ version of gsl_sf_bessel_yl_e().
       * @param l A nonnegative integer
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int yl_e( int l, double const x, result& result ){ return gsl_sf_bessel_yl_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_yl().
       * @param l A nonnegative integer
       * @param x A real number
       * @return Value of function
       */
      inline double yl( int const l, double const x ){ return gsl_sf_bessel_yl( l, x ); } 
      /**
       * C++ version of gsl_sf_bessel_yl_array().
       * @param lmax Maximum value of l
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int yl_array( int const lmax, double const x, double* result_array ){
	return gsl_sf_bessel_yl_array( lmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_i0_scaled_e().
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_0(x)
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int i0_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_i0_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_i0_scaled().
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_0(x)
       * @param x A real number
       * @return Value of function
       */
      inline double i0_scaled( double const x ){ return gsl_sf_bessel_i0_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_i1_scaled_e().
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_1(x)
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int i1_scaled_e( double const x, result& result ){ return gsl_sf_bessel_i1_scaled_e( x, &result ); } 
      /**
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_1(x)
       * C++ version of gsl_sf_bessel_i1_scaled().
       * @param x A real number
       * @return Value of function
       */
      inline double i1_scaled( double const x ){ return gsl_sf_bessel_i1_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_i2_scaled_e().
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_2(x)
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int i2_scaled_e( double const x, result& result ){ return gsl_sf_bessel_i2_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_i2_scaled().
       * Regular scaled modified spherical Bessel function
       *
       * Exp[-|x|] i_2(x)
       * @param x A real number
       * @return Value of function
       */
      inline double i2_scaled( double const x ){ return gsl_sf_bessel_i2_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_il_scaled_e().
       * Regular scaled modified spherical Bessel functions
       *
       * Exp[-|x|] i_l(x)
       *
       * i_l(x) = Sqrt[Pi/(2x)] BesselI[l+1/2,x]
       *
       * l >= 0
       * @param l A nonnegative integer
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW
       */
      inline int il_scaled_e( int const l, double x, result& result ){
	return gsl_sf_bessel_il_scaled_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_il_scaled().
       * Regular scaled modified spherical Bessel functions
       *
       * Exp[-|x|] i_l(x)
       *
       * i_l(x) = Sqrt[Pi/(2x)] BesselI[l+1/2,x]
       *
       * l >= 0
       * @param l A nonnegative integer
       * @param x A real number
       * @return Value of function
       */
      inline double il_scaled( int const l, double const x ){ return gsl_sf_bessel_il_scaled( l, x ); } 
      /**
       * C++ version of gsl_sf_bessel_il_scaled_array().
       * Regular scaled modified spherical Bessel functions
       *
       * Exp[-|x|] i_l(x)
       * for l=0,1,...,lmax
       * @param lmax Maximum value of l
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int il_scaled_array( int const lmax, double const x, double* result_array ){
	return gsl_sf_bessel_il_scaled_array( lmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_k0_scaled_e().
       * Irregular scaled modified spherical Bessel function
       * Exp[x] k_0(x)
       *
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int k0_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_k0_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_k0_scaled().
       * Irregular scaled modified spherical Bessel function
       * Exp[x] k_0(x)
       *
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double k0_scaled( double const x ){ return gsl_sf_bessel_k0_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_k1_scaled_e().
       * Irregular modified spherical Bessel function
       * Exp[x] k_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int k1_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_k1_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_k1_scaled().
       * Irregular modified spherical Bessel function
       * Exp[x] k_1(x)
       *
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double k1_scaled( double const x ){ return gsl_sf_bessel_k1_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_k2_scaled_e().
       * Irregular modified spherical Bessel function
       * Exp[x] k_2(x)
       *
       * x > 0.0
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EUNDRFLW or GSL_EOVRFLW
       */
      inline int k2_scaled_e( double const x, result& result ){
	return gsl_sf_bessel_k2_scaled_e( x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_k2_scaled().
       * Irregular modified spherical Bessel function
       * Exp[x] k_2(x)
       *
       * x > 0.0
       * @param x A real number
       * @return Value of function
       */
      inline double k2_scaled( double const x ){ return gsl_sf_bessel_k2_scaled( x ); } 
      /**
       * C++ version of gsl_sf_bessel_kl_scaled_e().
       * Irregular scaled modified spherical Bessel function
       * Exp[x] k_l(x)
       *
       * for l=0,1,...,lmax
       * @param l A nonnegative integer
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int kl_scaled_e( int l, double const x, result& result ){
	return gsl_sf_bessel_kl_scaled_e( l, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_kl_scaled().
       * Irregular scaled modified spherical Bessel function
       * Exp[x] k_l(x)
       *
       * for l=0,1,...,lmax
       * @param l A nonnegative integer
       * @param x A real number
       * @return Value of function
       */
      inline double kl_scaled( int const l, double const x ){ return gsl_sf_bessel_kl_scaled( l, x ); } 
      /**
       * C++ version of gsl_sf_bessel_kl_scaled_array().
       * Irregular scaled modified spherical Bessel function
       * Exp[x] k_l(x)
       *
       * for l=0,1,...,lmax
       * @param lmax Maximum value of l
       * @param x A real number
       * @param result_array Array of reals to store result in
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int kl_scaled_array( int const lmax, double const x, double* result_array ){
	return gsl_sf_bessel_kl_scaled_array( lmax, x, result_array ); } 
      /**
       * C++ version of gsl_sf_bessel_Jnu_e().
       * Regular cylindrical Bessel function J_nu(x)
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int Jnu_e( double const nu, double const x, result& result ){
	return gsl_sf_bessel_Jnu_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Jnu().
       * Regular cylindrical Bessel function J_nu(x)
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Jnu( double const nu, double const x ){ return gsl_sf_bessel_Jnu( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Ynu_e().
       * Irregular cylindrical Bessel function Y_nu(x)
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int Ynu_e( double nu, double x, result& result ){
	return gsl_sf_bessel_Ynu_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Ynu().
       * Irregular cylindrical Bessel function Y_nu(x)
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Ynu( double const nu, double const x ){ return gsl_sf_bessel_Ynu( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_sequence_Jnu_e().
       * Regular cylindrical Bessel function J_nu(x)
       * evaluated at a series of x values. The array
       * contains the x values. They are assumed to be
       * strictly ordered and positive. The array is
       * over-written with the values of J_nu(x_i).
       * @param nu The fractional order
       * @param mode A @c gsl::mod_t mode
       * @param size The size of the array
       * @param v The array
       * @return GSL_SUCCESS or GSL_EDOM or GSL_EINVAL
       */
      inline int sequence_Jnu_e( double nu, mode_t mode, size_t size, double* v ){
	return gsl_sf_bessel_sequence_Jnu_e( nu, mode, size, v ); } 
      /**
       * C++ version of gsl_sf_bessel_Inu_scaled_e().
       * Scaled modified cylindrical Bessel functions
       *
       * Exp[-|x|] BesselI[nu, x]
       * x >= 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int Inu_scaled_e( double nu, double x, result& result ){
	return gsl_sf_bessel_Inu_scaled_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Inu_scaled().
       * Scaled modified cylindrical Bessel functions
       *
       * Exp[-|x|] BesselI[nu, x]
       * x >= 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Inu_scaled( double nu, double x ){ return gsl_sf_bessel_Inu_scaled( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Inu_e().
       * Modified cylindrical Bessel functions
       *
       * BesselI[nu, x]
       * x >= 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or  GSL_EDOM, GSL_EOVERFLW
       */
      inline int Inu_e( double nu, double x, result& result ){
	return gsl_sf_bessel_Inu_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Inu().
       * Modified cylindrical Bessel functions
       *
       * BesselI[nu, x]
       * x >= 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Inu( double nu, double x ){ return gsl_sf_bessel_Inu( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Knu_scaled_e().
       * Scaled modified cylindrical Bessel functions
       *
       * Exp[+|x|] BesselK[nu, x]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int Knu_scaled_e( double const nu, double const x, result& result ){
	return gsl_sf_bessel_Knu_scaled_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Knu_scaled().
       * Scaled modified cylindrical Bessel functions
       *
       * Exp[+|x|] BesselK[nu, x]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Knu_scaled( double const nu, double const x ){
	return gsl_sf_bessel_Knu_scaled( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_Knu_e().
       * Modified cylindrical Bessel functions
       *
       * BesselK[nu, x]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM, GSL_EUNDRFLW
       */
      inline int Knu_e( double const nu, double const x, result& result ){
	return gsl_sf_bessel_Knu_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_Knu().
       * Modified cylindrical Bessel functions
       *
       * BesselK[nu, x]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double Knu( double const nu, double const x ){ return gsl_sf_bessel_Knu( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_lnKnu_e().
       * Logarithm of modified cylindrical Bessel functions.
       *
       * Log[BesselK[nu, x]]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EDOM
       */
      inline int lnKnu_e( double const nu, double const x, result& result ){
	return gsl_sf_bessel_lnKnu_e( nu, x, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_lnKnu().
       * Logarithm of modified cylindrical Bessel functions.
       *
       * Log[BesselK[nu, x]]
       * x > 0, nu >= 0
       * @param nu The fractional order
       * @param x A real number
       * @return Value of function
       */
      inline double lnKnu( double const nu, double const x ){
	return gsl_sf_bessel_lnKnu( nu, x ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_J0_e().
       * sth positive zero of the Bessel function J_0(x).
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int zero_J0_e( unsigned int s, result& result ){
	return gsl_sf_bessel_zero_J0_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_J0().
       * sth positive zero of the Bessel function J_0(x).
       * @param s An integer
       * @return Value of function
       */
      inline double zero_J0( unsigned int s ){ return gsl_sf_bessel_zero_J0( s ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_J1_e().
       * s'th positive zero of the Bessel function J_1(x).
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int zero_J1_e( unsigned int s, result& result ){
	return gsl_sf_bessel_zero_J1_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_J1().
       * sth positive zero of the Bessel function J_1(x).
       * @param s An integer
       * @return Value of function
       */
      inline double zero_J1( unsigned int s ){ return gsl_sf_bessel_zero_J1( s ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_Jnu_e().
       * s'th positive zero of the Bessel function J_nu(x).
       * @param nu The fractional order
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int zero_Jnu_e( double nu, unsigned int s, result& result ){
	return gsl_sf_bessel_zero_Jnu_e( nu, s, &result ); } 
      /**
       * C++ version of gsl_sf_bessel_zero_Jnu().
       * s'th positive zero of the Bessel function J_nu(x).
       * @param nu The fractional order
       * @param s An integer
       * @return Value of function
       */
      inline double zero_Jnu( double nu, unsigned int s ){
	return gsl_sf_bessel_zero_Jnu( nu, s ); }
    }
  }
}

#endif
