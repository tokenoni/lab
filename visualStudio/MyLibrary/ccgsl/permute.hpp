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

#ifndef CCGSL_PERMUTE_HPP
#define CCGSL_PERMUTE_HPP

#include<cmath>
#include<gsl/gsl_permute.h>

namespace gsl {
  /**
   * This namespace handles the @c gsl_permute functions.
   */
  namespace permute {
    /**
     * C++ version of gsl_permute_complex_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int complex_inverse( size_t const* p, double* data, size_t const stride, size_t const n ){
      return gsl_permute_complex_inverse( p, data, stride, n ); }
    /**
     * C++ version of gsl_permute_complex_float_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int complex_float_inverse( size_t const* p, float* data, size_t const stride, size_t const n ){
      return gsl_permute_complex_float_inverse( p, data, stride, n ); }
    /**
     * C++ version of gsl_permute_complex_long_double_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int complex_long_double_inverse( size_t const* p, long double* data, size_t const stride,
        size_t const n ){
      return gsl_permute_complex_long_double_inverse( p, data, stride, n ); }
    /**
     * C++ version of gsl_permute_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int inverse( size_t const* p, double* data, size_t const stride, size_t const n ){
      return gsl_permute_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_float_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int float_inverse( size_t const* p, float* data, size_t const stride, size_t const n ){
      return gsl_permute_float_inverse( p, data, stride, n ); }
    /**
     * C++ version of gsl_permute_int_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int int_inverse( size_t const* p, int* data, size_t const stride, size_t const n ){
      return gsl_permute_int_inverse( p, data, stride, n ); }

    /**
     * C++ version of gsl_permute_long_double_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int long_double_inverse( size_t const* p, long double* data, size_t const stride,
        size_t const n ){
      return gsl_permute_long_double_inverse( p, data, stride, n ); }
    /**
     * C++ version of gsl_permute_long_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int long_inverse( size_t const* p, long* data, size_t const stride, size_t const n ){
      return gsl_permute_long_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_short_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int short_inverse( size_t const* p, short* data, size_t const stride, size_t const n ){
      return gsl_permute_short_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_uchar_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int uchar_inverse( size_t const* p, unsigned char* data, size_t const stride, size_t const n ){
      return gsl_permute_uchar_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_uint_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int uint_inverse( size_t const* p, unsigned int* data, size_t const stride, size_t const n ){
      return gsl_permute_uint_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_ulong_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int ulong_inverse( size_t const* p, unsigned long* data, size_t const stride, size_t const n ){
      return gsl_permute_ulong_inverse( p, data, stride, n ); } 
    /**
     * C++ version of gsl_permute_ushort_inverse().
     * @param p Integer representation of permutation
     * @param data An array to be permuted
     * @param stride The size of the stride to use within @c data
     * @param n The size of @c data
     * @return Error code on failure
     */
    inline int ushort_inverse( size_t const* p, unsigned short* data, size_t const stride, size_t const n ){
      return gsl_permute_ushort_inverse( p, data, stride, n ); }
  }
}

#endif
