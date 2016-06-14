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

#ifndef CCGSL_SF_AIRY_HPP
#define CCGSL_SF_AIRY_HPP

#include<gsl/gsl_sf_airy.h>
#include"mode.hpp"
#include"sf_result.hpp"

namespace gsl {
  namespace sf {
    /**
     * Namespace for Airy functions.
     */
    namespace airy {
      /**
       * C++ version of gsl_sf_airy_Ai_e().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS or GSL_EUNDRFLW
       */
      inline int airy_Ai_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Ai_e( x, mode, &result ); }
      /**
       * C++ version of gsl_sf_airy_Ai().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Ai( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Ai( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_e().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS OR GSL_OVERFLW
       */
      inline int airy_Bi_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Bi_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Bi().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Bi( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Bi( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_scaled_e().
       * Scaled Ai(x): Ai(x) if  x < 0;
       * exp(+2/3 x^{3/2}) Ai(x) if  x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_Ai_scaled_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Ai_scaled_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_scaled().
       * Scaled Ai(x): Ai(x) if  x < 0;
       * exp(+2/3 x^{3/2}) Ai(x) if  x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Ai_scaled( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Ai_scaled( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_scaled_e().
       * scaled Bi(x): Bi(x) if x < 0;
       * exp(-2/3 x^{3/2}) Bi(x) if x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_Bi_scaled_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Bi_scaled_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_scaled().
       * scaled Bi(x): Bi(x) if x < 0;
       * exp(-2/3 x^{3/2}) Bi(x) if x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Bi_scaled( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Bi_scaled( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_deriv_e().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS OR GSL_EUNDRFLW
       */
      inline int airy_Ai_deriv_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Ai_deriv_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_deriv().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Ai_deriv( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Ai_deriv( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_deriv_e().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS OR GSL_EOVERFLW
       */
      inline int airy_Bi_deriv_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Bi_deriv_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_deriv().
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Bi_deriv( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Bi_deriv( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_deriv_scaled_e().
       * Scaled derivative Ai'(x): Ai'(x) if x < 0;
       * exp(+2/3 x^{3/2}) Ai'(x) if x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_Ai_deriv_scaled_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Ai_deriv_scaled_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Ai_deriv_scaled().
       * Scaled derivative Ai'(x): Ai'(x) if x < 0;
       * exp(+2/3 x^{3/2}) Ai'(x) if x > 0
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Ai_deriv_scaled( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Ai_deriv_scaled( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_deriv_scaled_e().
       * Scaled derivative: Bi'(x) if x < 0;
       * exp(-2/3 x^{3/2}) Bi'(x) if x > 0.
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_Bi_deriv_scaled_e( double const x, gsl::mode_t mode, result& result ){
	return gsl_sf_airy_Bi_deriv_scaled_e( x, mode, &result ); } 
      /**
       * C++ version of gsl_sf_airy_Bi_deriv_scaled().
       * Scaled derivative: Bi'(x) if x < 0;
       * exp(-2/3 x^{3/2}) Bi'(x) if x > 0.
       * @param x A real number
       * @param mode The @c gsl::mode_t mode
       * @return function result
       */
      inline double airy_Bi_deriv_scaled( double const x, gsl::mode_t mode ){
	return gsl_sf_airy_Bi_deriv_scaled( x, mode ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Ai_e().
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_zero_Ai_e( unsigned int s, result& result ){
	return gsl_sf_airy_zero_Ai_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Ai().
       * @param s An integer
       * @return function result
       */
      inline double airy_zero_Ai( unsigned int s ){ return gsl_sf_airy_zero_Ai( s ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Bi_e().
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_zero_Bi_e( unsigned int s, result& result ){
	return gsl_sf_airy_zero_Bi_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Bi().
       * @param s An integer
       * @return function result
       */
      inline double airy_zero_Bi( unsigned int s ){ return gsl_sf_airy_zero_Bi( s ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Ai_deriv_e().
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_zero_Ai_deriv_e( unsigned int s, result& result ){
	return gsl_sf_airy_zero_Ai_deriv_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Ai_deriv().
       * @param s An integer
       * @return function result
       */
      inline double airy_zero_Ai_deriv( unsigned int s ){ return gsl_sf_airy_zero_Ai_deriv( s ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Bi_deriv_e().
       * @param s An integer
       * @param result The result as a @c gsl::sf::result object
       * @return GSL_SUCCESS
       */
      inline int airy_zero_Bi_deriv_e( unsigned int s, result& result ){
	return gsl_sf_airy_zero_Bi_deriv_e( s, &result ); } 
      /**
       * C++ version of gsl_sf_airy_zero_Bi_deriv().
       * @param s An integer
       * @return function result
       */
      inline double airy_zero_Bi_deriv( unsigned int s ){ return gsl_sf_airy_zero_Bi_deriv( s ); } 
    }
  }
}

#endif
