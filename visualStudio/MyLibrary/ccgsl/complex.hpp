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

#ifndef CCGSL_COMPLEX_HPP
#define CCGSL_COMPLEX_HPP

#include<cmath>
#include<cstring>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include"exception.hpp"

namespace gsl {
  /**
   * This class handles complex numbers. There is an ambiguity in translating gsl_complex 
   * to C++: @c gsl_complex is a @c struct and so should be represented by a C++ @c class.
   * But @c gsl_complex is also used like a @c namespace and so it would be nice to have
   * a C++ @c namespace @c gsl:: complex. However, we cannot have both. The compromise is to
   * use @c static @c complex @c class functions. But, since we cannot use
   * @code
   * using namespace gsl::complex
   * @endcode
   * and similar namespace operations with @c gsl::complex, we also define a namespace
   * @c gsl::cpx with the same functions, again defined @c inline for efficiency.
   */
  class complex : private gsl_complex {
    friend class complex_ptr;
    friend class complex_ref;
    friend class vector_complex;
    friend class matrix_complex;
 public:
    // /**
    //  * Make sure this is convertible to @c gsl_complex. This is not ideal because it copies
    //  * twice. But it gets round the problem that a @c gsl_complex must contain constant data.
    //  * @return A @c gsl_complex
    //  */ 
    // operator gsl_complex() const {
    //   gsl_complex z; memcpy( z.dat, dat, 2 * sizeof( double ) ); return z; }
    // Default constructible: so that it can  be stored in a container
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex(){}
    /**
     * Allow construction from raw data. This is mainly for use by the container classes.
     * @param dat Raw data
     */
    explicit complex( double* dat ){
      memcpy( this->dat, dat, 2 * sizeof( double ) );
    }
    // copy constructor: default is OK
    // assignment operator: default is OK
    // destructor: default is OK
    /**
     * Constructor from base class.
     * @param z A gsl_complex
     */
    complex( gsl_complex const& z ){
      GSL_SET_REAL( this, GSL_REAL( z ) );
      GSL_SET_IMAG( this, GSL_IMAG( z ) );
    }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex& get(){ return dynamic_cast<gsl_complex&>( *this ); }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex const& get() const { return dynamic_cast<gsl_complex const&>( *this ); }
    /**
     * C++ version of gsl_complex_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex
     */
    static complex rect( double x, double y ){ complex z; GSL_SET_COMPLEX( &z, x, y ); return z; }
    /**
     * C++ version of gsl_complex_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex
     */
    static complex polar( double r, double theta ){ complex z;
      GSL_SET_COMPLEX( &z, r * std::cos( theta ), r * std::sin( theta ) ); return z; }
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex( double x, double y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    double real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    double imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( double x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( double y ){ GSL_SET_IMAG( this, y ); }
    /**
     * C++ version of gsl_complex_arg().
     * @param z A complex number
     * @return The argument of @c z
     */
    static double arg( complex const& z ){ return gsl_complex_arg( z ); }
    /**
     * C++ version of gsl_complex_abs().
     * @param z A complex number
     * @return The magnitude of @c z
     */
    static double abs( complex const& z ){ return gsl_complex_abs( z ); }
    /**
     * C++ version of gsl_complex_abs2().
     * @param z A complex number
     * @return The square of the magnitude of @c z
     */
    static double abs2( complex const& z ){ return gsl_complex_abs2( z ); }
    /**
     * C++ version of gsl_complex_logabs().
     * @param z A complex number
     * @return the natural logarithm of the magnitude of @c z
     */
    static double logabs( complex const& z ){ return gsl_complex_logabs( z ); }
    /**
     * C++ version of gsl_complex_add().
     * @param a A complex number
     * @param b A complex to add to @c a
     * @return The sum
     */
    static complex add( complex const& a, complex const& b )
    { return complex( gsl_complex_add( a, b ) ); }
    /**
     * C++ version of gsl_complex_sub().
     * @param a A complex number
     * @param b A complex to subtract from @c a
     * @return The difference
     */
    static complex sub( complex const& a, complex const& b )
    { return complex( gsl_complex_sub( a, b ) ); }
    /**
     * C++ version of gsl_complex_mul().
     * @param a A complex number
     * @param b A complex to multiply by @c a
     * @return The product
     */
    static complex mul( complex const& a, complex const& b )
    { return complex( gsl_complex_mul( a, b ) ); }
    /**
     * C++ version of gsl_complex_div().
     * @param a A complex number
     * @param b A complex to divide @c a by
     * @return The quotient
     */
    static complex div( complex const& a, complex const& b )
    { return complex( gsl_complex_div( a, b ) ); }
    /**
     * C++ version of gsl_complex_add_real().
     * @param a A complex number
     * @param x A double to add to @c a
     * @return The sum
     */
    static complex add_real( complex const& a, double x )
    { return complex( gsl_complex_add_real( a, x ) ); }
    /**
     * C++ version of gsl_complex_sub_real().
     * @param a A complex number
     * @param x A double to add to @c a
     * @return  The sum
     */
    static complex sub_real( complex const& a, double x )
    { return complex( gsl_complex_sub_real( a, x ) ); }
    /**
     * C++ version of gsl_complex_mul_real().
     * @param a A complex number
     * @param x A double to multiply @c a by
     * @return The product
     */
    static complex mul_real( complex const& a, double x )
    { return complex( gsl_complex_mul_real( a, x ) ); }
    /**
     * C++ version of gsl_complex_div_real().
     * @param a A complex number
     * @param x A double to divide @c a by
     * @return The quotient
     */
    static complex div_real( complex const& a, double x )
    { return complex( gsl_complex_div_real( a, x ) ); }
    /**
     * C++ version of gsl_complex_add_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to add to @c a
     * @return The sum
     */
    static complex add_imag( complex const& a, double y )
    { return complex( gsl_complex_add_imag( a, y ) ); }
    /**
     * C++ version of gsl_complex_sub_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to subtract from @c a
     * @return The difference
     */
    static complex sub_imag( complex const& a, double y )
    { return complex( gsl_complex_sub_imag( a, y ) ); }
    /**
     * C++ version of gsl_complex_mul_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to multiply @c a by
     * @return The product
     */
    static complex mul_imag( complex const& a, double y )
    { return complex( gsl_complex_mul_imag( a, y ) ); }
    /**
     * C++ version of gsl_complex_div_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to divide @c a by
     * @return The quotient
     */
    static complex div_imag( complex const& a, double y )
    { return complex( gsl_complex_div_imag( a, y ) ); }
    /**
     * C++ version of gsl_complex_conjugate().
     * @param z A complex number
     * @return The complex conjugate
     */
    static complex conjugate( complex const& z )
    { return complex( gsl_complex_conjugate( z ) ); } 
    /**
     * C++ version of gsl_complex_inverse().
     * @param a A complex number
     * @return The inverse
     */
    static complex inverse( complex const& a )
    { return complex( gsl_complex_inverse( a ) ); } 
    /**
     * C++ version of gsl_complex_negative().
     * @param a A complex number
     * @return The negative
     */
    static complex negative( complex const& a )
    { return complex( gsl_complex_negative( a ) ); } 
    /**
     * C++ version of gsl_complex_sqrt().
     * @param z A complex number
     * @return The square root
     */
    static complex sqrt( complex const& z )
    { return complex( gsl_complex_sqrt( z ) ); }
    /**
     * C++ version of gsl_complex_sqrt_real().
     * @param x A double
     * @return The complex square root
     */
    static complex sqrt_real( double x )
    { return complex( gsl_complex_sqrt_real( x ) ); }
    /**
     * C++ version of gsl_complex_pow().
     * @param a A complex number
     * @param b Another complex number
     * @return \f$a^b\f$
     */
    static complex pow( complex const& a, complex const& b )
    { return complex( gsl_complex_pow( a, b ) ); }
    /**
     * C++ version of gsl_complex_pow_real().
     * @param a A complex number
     * @param b A double
     * @return \f$a^b\f$
     */
    static complex pow_real( complex const& a, double b )
    { return complex( gsl_complex_pow_real( a, b ) ); }
    /**
     * C++ version of gsl_complex_exp().
     * @param a A complex number
     * @return The exponential function of @c z
     */
    static complex exp( complex const& a )
    { return complex( gsl_complex_exp( a ) ); } 
    /**
     * C++ version of gsl_complex_log().
     * @param a A complex number
     * @return The natural logarithm
     */
    static complex log( complex const& a )
    { return complex( gsl_complex_log( a ) ); } 
    /**
     * C++ version of gsl_complex_log10().
     * @param a A complex number
     * @return The base 10 logarithm
     */
    static complex log10( complex const& a )
    { return complex( gsl_complex_log10( a ) ); } 
    /**
     * C++ version of gsl_complex_log_b().
     * @param a A complex number
     * @param b A complex number
     * @return The complex base@c b logarithm
     */
    static complex log_b( complex const& a, complex const& b )
    { return complex( gsl_complex_log_b( a, b ) ); }
    /**
     * C++ version of gsl_complex_sin().
     * @param a A complex number
     * @return The sine
     */
    static complex sin( complex const& a )
    { return complex( gsl_complex_sin( a ) ); } 
    /**
     * C++ version of gsl_complex_cos().
     * @param a A complex number
     * @return The cosine
     */
    static complex cos( complex const& a )
    { return complex( gsl_complex_cos( a ) ); } 
    /**
     * C++ version of gsl_complex_sec().
     * @param a A complex number
     * @return The secant
     */
    static complex sec( complex const& a )
    { return complex( gsl_complex_sec( a ) ); } 
    /**
     * C++ version of gsl_complex_csc().
     * @param a A complex number
     * @return The cosecant
     */
    static complex csc( complex const& a )
    { return complex( gsl_complex_csc( a ) ); } 
    /**
     * C++ version of gsl_complex_tan().
     * @param a A complex number
     * @return The tangent
     */
    static complex tan( complex const& a )
    { return complex( gsl_complex_tan( a ) ); } 
    /**
     * C++ version of gsl_complex_cot().
     * @param a A complex number
     * @return The cotangent
     */
    static complex cot( complex const& a )
    { return complex( gsl_complex_cot( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsin().
     * @param a A complex number
     * @return The inverse sine
     */
    static complex arcsin( complex const& a )
    { return complex( gsl_complex_arcsin( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsin_real().
     * @param a A double
     * @return The inverse sine
     */
    static complex arcsin_real( double a )
    { return complex( gsl_complex_arcsin_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arccos().
     * @param a A complex number
     * @return The inverse cosine
     */
    static complex arccos( complex const& a )
    { return complex( gsl_complex_arccos( a ) ); } 
    /**
     * C++ version of gsl_complex_arccos_real().
     * @param a A double
     * @return The inverse cosine
     */
    static complex arccos_real( double a )
    { return complex( gsl_complex_arccos_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsec().
     * @param a A complex number
     * @return The inverse secant
     */
    static complex arcsec( complex const& a )
    { return complex( gsl_complex_arcsec( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsec_real().
     * @param a A double
     * @return The inverse secant
     */
    static complex arcsec_real( double a )
    { return complex( gsl_complex_arcsec_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arccsc().
     * @param a A complex number
     * @return The inverse cosecant
     */
    static complex arccsc( complex const& a )
    { return complex( gsl_complex_arccsc( a ) ); } 
    /**
     * C++ version of gsl_complex_arccsc_real().
     * @param a A double
     * @return The inverse cosecant
     */
    static complex arccsc_real( double a )
    { return complex( gsl_complex_arccsc_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arctan().
     * @param a A complex number
     * @return the inverse tangent
     */
    static complex arctan( complex const& a )
    { return complex( gsl_complex_arctan( a ) ); } 
    /**
     * C++ version of gsl_complex_arccot().
     * @param a A complex number
     * @return the inverse cotangent
     */
    static complex arccot( complex const& a )
    { return complex( gsl_complex_arccot( a ) ); } 
    /**
     * C++ version of gsl_complex_sinh().
     * @param a A complex number
     * @return The hyperbolic sine
     */
    static complex sinh( complex const& a )
    { return complex( gsl_complex_sinh( a ) ); } 
    /**
     * C++ version of gsl_complex_cosh().
     * @param a A complex number
     * @return The hyperbolic cosine
     */
    static complex cosh( complex const& a )
    { return complex( gsl_complex_cosh( a ) ); } 
    /**
     * C++ version of gsl_complex_sech().
     * @param a A complex number
     * @return The hyperbolic secant
     */
    static complex sech( complex const& a )
    { return complex( gsl_complex_sech( a ) ); } 
    /**
     * C++ version of gsl_complex_csch().
     * @param a A complex number
     * @return The hyperbolic cosecant
     */
    static complex csch( complex const& a )
    { return complex( gsl_complex_csch( a ) ); } 
    /**
     * C++ version of gsl_complex_tanh().
     * @param a A complex number
     * @return The hyperbolic tangent
     */
    static complex tanh( complex const& a )
    { return complex( gsl_complex_tanh( a ) ); } 
    /**
     * C++ version of gsl_complex_coth().
     * @param a A complex number
     * @return The hyperbolic cotangent
     */
    static complex coth( complex const& a )
    { return complex( gsl_complex_coth( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsinh().
     * @param a A complex number
     * @return The inverse hypebolic sine
     */
    static complex arcsinh( complex const& a )
    { return complex( gsl_complex_arcsinh( a ) ); } 
    /**
     * C++ version of gsl_complex_arccosh().
     * @param a A complex number
     * @return The inverse hypebolic cosine
     */
    static complex arccosh( complex const& a )
    { return complex( gsl_complex_arccosh( a ) ); } 
    /**
     * C++ version of gsl_complex_arccosh_real().
     * @param a A double
     * @return The inverse hypebolic cosine
     */
    static complex arccosh_real( double a )
    { return complex( gsl_complex_arccosh_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arcsech().
     * @param a A complex number
     * @return The inverse hyperbolic secant
     */
    static complex arcsech( complex const& a )
    { return complex( gsl_complex_arcsech( a ) ); } 
    /**
     * C++ version of gsl_complex_arccsch().
     * @param a A complex number
     * @return The inverse hyperbolic cosecant
     */
    static complex arccsch( complex const& a )
    { return complex( gsl_complex_arccsch( a ) ); } 
    /**
     * C++ version of gsl_complex_arctanh().
     * @param a A complex number
     * @return The inverse hyperbolic tangent
     */
    static complex arctanh( complex const& a )
    { return complex( gsl_complex_arctanh( a ) ); } 
    /**
     * C++ version of gsl_complex_arctanh_real().
     * @param a A double
     * @return The inverse hyperbolic tangent
     */
    static complex arctanh_real( double a )
    { return complex( gsl_complex_arctanh_real( a ) ); } 
    /**
     * C++ version of gsl_complex_arccoth().
     * @param a A complex number
     * @return The inverse hyperbolic cotangent
     */
    static complex arccoth( complex const& a )
    { return complex( gsl_complex_arccoth( a ) ); }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the equality operator.
     * @param z Another complex
     * @return @c true or @c false according as @c this equals @c z
     */
    bool operator==( complex const& z ) const { return dat[0] == z.dat[0] and dat[1] == z.dat[1]; }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the inequality operator.
     * @param z Another complex
     * @return @c false or @c true according as @c this equals @c z
     */
    bool operator!=( complex const& z ) const { return dat[0] != z.dat[0] or dat[1] != z.dat[1]; }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than operator.
     * @param z Another complex
     * @return @c true or @c false according as @c this is less than @c z
     */
    bool operator<( complex const& z ) const {
      return dat[0] < z.dat[0] or (dat[0] == z.dat[0] and dat[1] < z.dat[1]); }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than operator.
     * @param z Another complex
     * @return @c true or @c false according as @c this is greater than than @c z
     */
    bool operator>( complex const& z ) const {
      return dat[0] > z.dat[0] or (dat[0] == z.dat[0] and dat[1] > z.dat[1]); }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than or equal to operator.
     * @param z Another complex
     * @return @c true or @c false according as @c this is less than or equal to  @c z
     */
    bool operator<=( complex const& z ) const { return *this < z or *this == z; }
    /**
     * A complex object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than or equal to operator.
     * @param z Another complex
     * @return @c true or @c false according as @c this is greater than or equal to  @c z
     */
    bool operator>=( complex const& z ) const { return *this > z or *this == z; }
    
//--- operator for mathematical description ... by fujii  ---
	//-- summation --
	complex operator +(complex const& a){return complex::add(*this, a);}
	complex operator +( double const& a){return complex::add_real(*this, a);}
	void operator +=(complex const& a){	*this = complex::add(*this, a);}
	void operator +=(double const& a){	*this = complex::add_real(*this, a);}
	//-- subtraction --
	complex operator -(complex const& a){return complex::sub(*this, a);}
	complex operator -( double const& a){return complex::sub_real(*this, a);}
	void operator -=(complex const& a){	*this = complex::sub(*this, a);}
	void operator -=(double const& a){	*this = complex::sub_real(*this, a);}
	//-- multiplication-
	complex operator *(complex const& a){return complex::mul(*this, a);}
	complex operator *( double const& a){return complex::mul_real(*this, a);}
	void operator *=(complex const& a){	*this = complex::mul(*this, a);}
	void operator *=(double const& a){	*this = complex::mul_real(*this, a);}
	//-- multiplication-
	complex operator /(complex const& a){return complex::div(*this, a);}
	complex operator /( double const& a){return complex::div_real(*this, a);}
	void operator /=(complex const& a){	*this = complex::div(*this, a);}
	void operator /=(double const& a){	*this = complex::div_real(*this, a);}
	//-- substitution ---
	void operator = (double const &a){  *this = complex::rect(a, 0.0);}
  };

  /**
   * This class can be used like a reference for complex objects so that we can
   * iterate over a vector (for example) of them. 
   */
  class complex_ref {
  protected:
    /**
     * The data. Never initialised or destroyed within the class. Always precisely two elements.
     */
    double* dat;
  public:
    /**
     * We use this in constructing complex_ptr objects.
     * @param dat A pointer to the data
     */
    complex_ref( double* dat ) : dat( dat ){}
   /**
     * Make sure this is convertible to @c gsl_complex. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex must contain constant data.
     * @return A @c gsl_complex
     */ 
    operator gsl_complex() const {
      gsl_complex z; memcpy( z.dat, dat, 2 * sizeof( double ) ); return z; }
    /**
     * Make sure this is convertible to @c complex. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex must contain constant data.
     * @return A @c complex
     */ 
    operator complex() const {
      complex z; memcpy( z.dat, dat, 2 * sizeof( double ) ); return z; }
    // Default constructible
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex_ref() : dat( 0 ){}
    /**
     * Make sure we can construct from a complex.
     * @param z The complex to construct from
     */
    complex_ref( complex& z ) : dat( z.dat ){}
    // copy constructor: default is OK
    // assignment operator: default is OK
    // also need to copy from complex
    /**
     * Assignment from complex.
     * @param z The complex to assign from
     */
    complex_ref& operator=( complex const& z ){
      memcpy( dat, z.dat, 2 * sizeof( double ) ); return *this; }
    // destructor: default is OK
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex( double x, double y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    double real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    double imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( double x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( double y ){ GSL_SET_IMAG( this, y ); }
//--- operator for mathematical description ... by fujii  ---
	//-- summation --
	complex operator +(complex const& a){return complex::add(*this, a);}
	complex operator +( double const& a){return complex::add_real(*this, a);}
	void operator +=(complex const& a){	*this = complex::add(*this, a);}
	void operator +=(double const& a){	*this = complex::add_real(*this, a);}
	//-- subtraction --
	complex operator -(complex const& a){return complex::sub(*this, a);}
	complex operator -( double const& a){return complex::sub_real(*this, a);}
	void operator -=(complex const& a){	*this = complex::sub(*this, a);}
	void operator -=(double const& a){	*this = complex::sub_real(*this, a);}
	//-- multiplication-
	complex operator *(complex const& a){return complex::mul(*this, a);}
	complex operator *( double const& a){return complex::mul_real(*this, a);}
	void operator *=(complex const& a){	*this = complex::mul(*this, a);}
	void operator *=(double const& a){	*this = complex::mul_real(*this, a);}
	//-- multiplication-
	complex operator /(complex const& a){return complex::div(*this, a);}
	complex operator /( double const& a){return complex::div_real(*this, a);}
	void operator /=(complex const& a){	*this = complex::div(*this, a);}
	void operator /=(double const& a){	*this = complex::div_real(*this, a);}
	//-- substitution ---
	void operator = (double const &a){  *this = complex::rect(a, 0.0);}
   };

  /**
   * This class can be used like a pointer for complex objects so that we can
   * iterate over a vector (for example) of them.
   */
    class complex_ptr : private complex_ref {
  public:
    /**
     * Typically we are given a pointer to the data storing the complex and need
     * to construct a complex_ptr from it.
     * @param dat A pointer to the data
     */
    complex_ptr( double* dat ) : complex_ref( dat ){}
    /**
     * Dereference the pointer
     * @return the complex object
     */
    complex_ref operator*(){ return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex object
     */
    complex_ref* operator->(){ return this; }
    /**
     * Dereference the pointer
     * @return the complex object
     */
    complex_ref const operator*() const { return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex object
     */
    complex_ref const* operator->() const { return this; }
  };

  /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream, complex const& z ){
    double i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

    /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex_ref.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex_ref
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream,
					 complex_ref const& z ){
    double i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

  /**
   * The @c gsl::cpx namespace allows the @c gsl::complex @c static funtions to be called from
   * a @c namespace and so allows the usual C++ convenience of using code like the following.
   * @code
   * using namespace gsl;
   * using namespace gsl::cpx;
   * complex z = rect( 3, 1 );
   * complex w = sin( z );
   * @endcode
   */
  namespace cpx {
    /**
     * C++ version of gsl_complex_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex
     */
    inline complex rect( double x, double y ){ return complex::rect( x, y ); }
    /**
     * C++ version of gsl_complex_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex
     */
    inline complex polar( double r, double theta ){ return complex::polar( r, theta ); }
    /**
     * C++ version of gsl_complex_arg().
     * @param z A complex number
     * @return The argument of @c z
     */
    inline double arg( complex const& z ){ return complex::arg( z ); }
    /**
     * C++ version of gsl_complex_abs().
     * @param z A complex number
     * @return The magnitude of @c z
     */
    inline double abs( complex const& z ){ return complex::abs( z ); }
    /**
     * C++ version of gsl_complex_abs2().
     * @param z A complex number
     * @return The square of the magnitude of @c z
     */
    inline double abs2( complex const& z ){ return complex::abs2( z ); }
    /**
     * C++ version of gsl_complex_logabs().
     * @param z A complex number
     * @return the natural logarithm of the magnitude of @c z
     */
    inline double logabs( complex const& z ){ return complex::logabs( z ); }
    /**
     * C++ version of gsl_complex_add().
     * @param a A complex number
     * @param b A complex to add to @c a
     * @return The sum
     */
    inline complex add( complex const& a, complex const& b )
    { return complex::add( a, b ); }
    /**
     * C++ version of gsl_complex_sub().
     * @param a A complex number
     * @param b A complex to subtract from @c a
     * @return The difference
     */
    inline complex sub( complex const& a, complex const& b )
    { return complex::sub( a, b ); }
    /**
     * C++ version of gsl_complex_mul().
     * @param a A complex number
     * @param b A complex to multiply by @c a
     * @return The product
     */
    inline complex mul( complex const& a, complex const& b )
    { return complex::mul( a, b ); }
    /**
     * C++ version of gsl_complex_div().
     * @param a A complex number
     * @param b A complex to divide @c a by
     * @return The quotient
     */
    inline complex div( complex const& a, complex const& b )
    { return complex::div( a, b ); }
    /**
     * C++ version of gsl_complex_add_real().
     * @param a A complex number
     * @param x A double to add to @c a
     * @return The sum
     */
    inline complex add_real( complex const& a, double x )
    { return complex::add_real( a, x ); }
    /**
     * C++ version of gsl_complex_sub_real().
     * @param a A complex number
     * @param x A double to add to @c a
     * @return  The sum
     */
    inline complex sub_real( complex const& a, double x )
    { return complex::sub_real( a, x ); }
    /**
     * C++ version of gsl_complex_mul_real().
     * @param a A complex number
     * @param x A double to multiply @c a by
     * @return The product
     */
    inline complex mul_real( complex const& a, double x )
    { return complex::mul_real( a, x ); }
    /**
     * C++ version of gsl_complex_div_real().
     * @param a A complex number
     * @param x A double to divide @c a by
     * @return The quotient
     */
    inline complex div_real( complex const& a, double x )
    { return complex::div_real( a, x ); }
    /**
     * C++ version of gsl_complex_add_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to add to @c a
     * @return The sum
     */
    inline complex add_imag( complex const& a, double y )
    { return complex::add_imag( a, y ); }
    /**
     * C++ version of gsl_complex_sub_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to subtract from @c a
     * @return The difference
     */
    inline complex sub_imag( complex const& a, double y )
    { return complex::sub_imag( a, y ); }
    /**
     * C++ version of gsl_complex_mul_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to multiply @c a by
     * @return The product
     */
    inline complex mul_imag( complex const& a, double y )
    { return complex::mul_imag( a, y ); }
    /**
     * C++ version of gsl_complex_div_imag().
     * @param a A complex number
     * @param y The imaginary part of a double to divide @c a by
     * @return The quotient
     */
    inline complex div_imag( complex const& a, double y )
    { return complex::div_imag( a, y ); }
    /**
     * C++ version of gsl_complex_conjugate().
     * @param z A complex number
     * @return The complex conjugate
     */
    inline complex conjugate( complex const& z )
    { return complex::conjugate( z ); } 
    /**
     * C++ version of gsl_complex_inverse().
     * @param a A complex number
     * @return The inverse
     */
    inline complex inverse( complex const& a )
    { return complex::inverse( a ); } 
    /**
     * C++ version of gsl_complex_negative().
     * @param a A complex number
     * @return The negative
     */
    inline complex negative( complex const& a )
    { return complex::negative( a ); } 
    /**
     * C++ version of gsl_complex_sqrt().
     * @param z A complex number
     * @return The square root
     */
    inline complex sqrt( complex const& z )
    { return complex::sqrt( z ); }
    /**
     * C++ version of gsl_complex_sqrt_real().
     * @param x A double
     * @return The complex square root
     */
    inline complex sqrt_real( double x )
    { return complex::sqrt_real( x ); }
    /**
     * C++ version of gsl_complex_pow().
     * @param a A complex number
     * @param b Another complex number
     * @return \f$a^b\f$
     */
    inline complex pow( complex const& a, complex const& b )
    { return complex::pow( a, b ); }
    /**
     * C++ version of gsl_complex_pow_real().
     * @param a A complex number
     * @param b A double
     * @return \f$a^b\f$
     */
    inline complex pow_real( complex const& a, double b )
    { return complex::pow_real( a, b ); }
    /**
     * C++ version of gsl_complex_exp().
     * @param a A complex number
     * @return The exponential function of @c z
     */
    inline complex exp( complex const& a )
    { return complex::exp( a ); } 
    /**
     * C++ version of gsl_complex_log().
     * @param a A complex number
     * @return The natural logarithm
     */
    inline complex log( complex const& a )
    { return complex::log( a ); } 
    /**
     * C++ version of gsl_complex_log10().
     * @param a A complex number
     * @return The base 10 logarithm
     */
    inline complex log10( complex const& a )
    { return complex::log10( a ); } 
    /**
     * C++ version of gsl_complex_log_b().
     * @param a A complex number
     * @param b A complex number
     * @return The complex base@c b logarithm
     */
    inline complex log_b( complex const& a, complex const& b )
    { return complex::log_b( a, b ); }
    /**
     * C++ version of gsl_complex_sin().
     * @param a A complex number
     * @return The sine
     */
    inline complex sin( complex const& a )
    { return complex::sin( a ); } 
    /**
     * C++ version of gsl_complex_cos().
     * @param a A complex number
     * @return The cosine
     */
    inline complex cos( complex const& a )
    { return complex::cos( a ); } 
    /**
     * C++ version of gsl_complex_sec().
     * @param a A complex number
     * @return The secant
     */
    inline complex sec( complex const& a )
    { return complex::sec( a ); } 
    /**
     * C++ version of gsl_complex_csc().
     * @param a A complex number
     * @return The cosecant
     */
    inline complex csc( complex const& a )
    { return complex::csc( a ); } 
    /**
     * C++ version of gsl_complex_tan().
     * @param a A complex number
     * @return The tangent
     */
    inline complex tan( complex const& a )
    { return complex::tan( a ); } 
    /**
     * C++ version of gsl_complex_cot().
     * @param a A complex number
     * @return The cotangent
     */
    inline complex cot( complex const& a )
    { return complex::cot( a ); } 
    /**
     * C++ version of gsl_complex_arcsin().
     * @param a A complex number
     * @return The inverse sine
     */
    inline complex arcsin( complex const& a )
    { return complex::arcsin( a ); } 
    /**
     * C++ version of gsl_complex_arcsin_real().
     * @param a A double
     * @return The inverse sine
     */
    inline complex arcsin_real( double a )
    { return complex::arcsin_real( a ); } 
    /**
     * C++ version of gsl_complex_arccos().
     * @param a A complex number
     * @return The inverse cosine
     */
    inline complex arccos( complex const& a )
    { return complex::arccos( a ); } 
    /**
     * C++ version of gsl_complex_arccos_real().
     * @param a A double
     * @return The inverse cosine
     */
    inline complex arccos_real( double a )
    { return complex::arccos_real( a ); } 
    /**
     * C++ version of gsl_complex_arcsec().
     * @param a A complex number
     * @return The inverse secant
     */
    inline complex arcsec( complex const& a )
    { return complex::arcsec( a ); } 
    /**
     * C++ version of gsl_complex_arcsec_real().
     * @param a A double
     * @return The inverse secant
     */
    inline complex arcsec_real( double a )
    { return complex::arcsec_real( a ); } 
    /**
     * C++ version of gsl_complex_arccsc().
     * @param a A complex number
     * @return The inverse cosecant
     */
    inline complex arccsc( complex const& a )
    { return complex::arccsc( a ); } 
    /**
     * C++ version of gsl_complex_arccsc_real().
     * @param a A double
     * @return The inverse cosecant
     */
    inline complex arccsc_real( double a )
    { return complex::arccsc_real( a ); } 
    /**
     * C++ version of gsl_complex_arctan().
     * @param a A complex number
     * @return the inverse tangent
     */
    inline complex arctan( complex const& a )
    { return complex::arctan( a ); } 
    /**
     * C++ version of gsl_complex_arccot().
     * @param a A complex number
     * @return the inverse cotangent
     */
    inline complex arccot( complex const& a )
    { return complex::arccot( a ); } 
    /**
     * C++ version of gsl_complex_sinh().
     * @param a A complex number
     * @return The hyperbolic sine
     */
    inline complex sinh( complex const& a )
    { return complex::sinh( a ); } 
    /**
     * C++ version of gsl_complex_cosh().
     * @param a A complex number
     * @return The hyperbolic cosine
     */
    inline complex cosh( complex const& a )
    { return complex::cosh( a ); } 
    /**
     * C++ version of gsl_complex_sech().
     * @param a A complex number
     * @return The hyperbolic secant
     */
    inline complex sech( complex const& a )
    { return complex::sech( a ); } 
    /**
     * C++ version of gsl_complex_csch().
     * @param a A complex number
     * @return The hyperbolic cosecant
     */
    inline complex csch( complex const& a )
    { return complex::csch( a ); } 
    /**
     * C++ version of gsl_complex_tanh().
     * @param a A complex number
     * @return The hyperbolic tangent
     */
    inline complex tanh( complex const& a )
    { return complex::tanh( a ); } 
    /**
     * C++ version of gsl_complex_coth().
     * @param a A complex number
     * @return The hyperbolic cotangent
     */
    inline complex coth( complex const& a )
    { return complex::coth( a ); } 
    /**
     * C++ version of gsl_complex_arcsinh().
     * @param a A complex number
     * @return The inverse hypebolic sine
     */
    inline complex arcsinh( complex const& a )
    { return complex::arcsinh( a ); } 
    /**
     * C++ version of gsl_complex_arccosh().
     * @param a A complex number
     * @return The inverse hypebolic cosine
     */
    inline complex arccosh( complex const& a )
    { return complex::arccosh( a ); } 
    /**
     * C++ version of gsl_complex_arccosh_real().
     * @param a A double
     * @return The inverse hypebolic cosine
     */
    inline complex arccosh_real( double a )
    { return complex::arccosh_real( a ); } 
    /**
     * C++ version of gsl_complex_arcsech().
     * @param a A complex number
     * @return The inverse hyperbolic secant
     */
    inline complex arcsech( complex const& a )
    { return complex::arcsech( a ); } 
    /**
     * C++ version of gsl_complex_arccsch().
     * @param a A complex number
     * @return The inverse hyperbolic cosecant
     */
    inline complex arccsch( complex const& a )
    { return complex::arccsch( a ); } 
    /**
     * C++ version of gsl_complex_arctanh().
     * @param a A complex number
     * @return The inverse hyperbolic tangent
     */
    inline complex arctanh( complex const& a )
    { return complex::arctanh( a ); } 
    /**
     * C++ version of gsl_complex_arctanh_real().
     * @param a A double
     * @return The inverse hyperbolic tangent
     */
    inline complex arctanh_real( double a )
    { return complex::arctanh_real( a ); } 
    /**
     * C++ version of gsl_complex_arccoth().
     * @param a A complex number
     * @return The inverse hyperbolic cotangent
     */
    inline complex arccoth( complex const& a )
    { return complex::arccoth( a ); }
  }

}

#endif
