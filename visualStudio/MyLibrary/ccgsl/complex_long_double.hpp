/*
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

#ifndef CCGSL_COMPLEX_LONG_DOUBLE_HPP
#define CCGSL_COMPLEX_LONG_DOUBLE_HPP

#include<cmath>
#include<cstring>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

namespace gsl {
  /**
   * This class handles complex_long_double numbers. There is an ambiguity in translating gsl_complex_long_double 
   * to C++: @c gsl_complex_long_double is a @c struct and so should be represented by a C++ @c class.
   * But @c gsl_complex_long_double is also used like a @c namespace and so it would be nice to have
   * a C++ @c namespace @c gsl:: complex_long_double. However, we cannot have both. The compromise is to
   * use @c static @c complex_long_double @c class functions. But, since we cannot use
   * @code
   * using namespace gsl::complex_long_double
   * @endcode
   * and similar namespace operations with @c gsl::complex_long_double, we also define a namespace
   * @c gsl::cpx_long_double with the same functions, again defined @c inline for efficiency.
   */
  class complex_long_double : private gsl_complex_long_double {
    friend class complex_long_double_ptr;
    friend class complex_long_double_ref;
    friend class vector_complex_long_double;
    friend class matrix_complex_long_double;
 public:
    // /**
    //  * Make sure this is convertible to @c gsl_complex_long_double. This is not ideal because it copies
    //  * twice. But it gets round the problem that a @c gsl_comple must contain constant data.
    //  * @return A @c gsl_complex_long_double
    //  */ 
    // operator gsl_complex_long_double() const {
    //   gsl_complex_long_double z; memcpy( z.dat, dat, 2 * sizeof( long double ) ); return z; }
    // Default constructible: so that it can  be stored in a container
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex_long_double(){}
    /**
     * Allow construction from raw data. This is mainly for use by the container classes.
     * @param dat Raw data
     */
    explicit complex_long_double( long double* dat ){
      memcpy( this->dat, dat, 2 * sizeof( long double ) );
    }
    // copy constructor: default is OK
    // assignment operator: default is OK
    // destructor: default is OK
    /**
     * Constructor from base class.
     * @param z A gsl_complex_long_double
     */
    complex_long_double( gsl_complex_long_double const& z ){
      GSL_SET_REAL( this, GSL_REAL( z ) );
      GSL_SET_IMAG( this, GSL_IMAG( z ) );
    }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex_long_double& get(){ return dynamic_cast<gsl_complex_long_double&>( *this ); }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex_long_double const& get() const { return dynamic_cast<gsl_complex_long_double const&>( *this ); }
    /**
     * C++ version of gsl_complex_long_double_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex_long_double
     */
    static complex_long_double rect( long double x, long double y ){ complex_long_double z; GSL_SET_COMPLEX( &z, x, y ); return z; }
    /**
     * C++ version of gsl_complex_long_double_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex_long_double
     */
    static complex_long_double polar( long double r, long double theta ){ complex_long_double z;
      GSL_SET_COMPLEX( &z, r * std::cos( theta ), r * std::sin( theta ) ); return z; }
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex_long_double( long double x, long double y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    long double real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    long double imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( long double x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( long double y ){ GSL_SET_IMAG( this, y ); }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the equality operator.
     * @param z Another complex_long_double
     * @return @c true or @c false according as @c this equals @c z
     */
    bool operator==( complex_long_double const& z ) const { return dat[0] == z.dat[0] and dat[1] == z.dat[1]; }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the inequality operator.
     * @param z Another complex_long_double
     * @return @c false or @c true according as @c this equals @c z
     */
    bool operator!=( complex_long_double const& z ) const { return dat[0] != z.dat[0] or dat[1] != z.dat[1]; }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than operator.
     * @param z Another complex_long_double
     * @return @c true or @c false according as @c this is less than @c z
     */
    bool operator<( complex_long_double const& z ) const {
      return dat[0] < z.dat[0] or (dat[0] == z.dat[0] and dat[1] < z.dat[1]); }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than operator.
     * @param z Another complex_long_double
     * @return @c true or @c false according as @c this is greater than than @c z
     */
    bool operator>( complex_long_double const& z ) const {
      return dat[0] > z.dat[0] or (dat[0] == z.dat[0] and dat[1] > z.dat[1]); }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than or equal to operator.
     * @param z Another complex_long_double
     * @return @c true or @c false according as @c this is less than or equal to  @c z
     */
    bool operator<=( complex_long_double const& z ) const { return *this < z or *this == z; }
    /**
     * A complex_long_double object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than or equal to operator.
     * @param z Another complex_long_double
     * @return @c true or @c false according as @c this is greater than or equal to  @c z
     */
    bool operator>=( complex_long_double const& z ) const { return *this > z or *this == z; }
  };
    
  /**
   * This class can be used like a reference for complex_long_double objects so that we can
   * iterate over a vector (for example) of them. 
   */
  class complex_long_double_ref {
  protected:
    /**
     * The data. Never initialised or destroyed within the class. Always precisely two elements.
     */
    long double* dat;
  public:
    /**
     * We use this in constructing complex_long_double_ptr objects.
     * @param dat A pointer to the data
     */
    complex_long_double_ref( long double* dat ) : dat( dat ){}
   /**
     * Make sure this is convertible to @c gsl_complex_long_double. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex_long_double must contain constant data.
     * @return A @c gsl_complex_long_double
     */ 
    operator gsl_complex_long_double() const {
      gsl_complex_long_double z; memcpy( z.dat, dat, 2 * sizeof( long double ) ); return z; }
    /**
     * Make sure this is convertible to @c complex_long_double. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex_long_double must contain constant data.
     * @return A @c complex_long_double
     */ 
    operator complex_long_double() const {
      complex_long_double z; memcpy( z.dat, dat, 2 * sizeof( long double ) ); return z; }
    // Default constructible
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex_long_double_ref() : dat( 0 ){}
    /**
     * Make sure we can construct from a complex_long_double.
     * @param z The complex_long_double to construct from
     */
    complex_long_double_ref( complex_long_double& z ) : dat( z.dat ){}
    // copy constructor: default is OK
    // assignment operator: default is OK
    // also need to copy from complex_long_double
    /**
     * Assignment from complex_long_double.
     * @param z The complex_long_double to assign from
     */
    complex_long_double_ref& operator=( complex_long_double const& z ){
      memcpy( dat, z.dat, 2 * sizeof( long double ) ); return *this; }
    // destructor: default is OK
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex_long_double( long double x, long double y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    long double real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    long double imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( long double x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( long double y ){ GSL_SET_IMAG( this, y ); }
  };

  /**
   * This class can be used like a pointer for complex_long_double objects so that we can
   * iterate over a vector (for example) of them.
   */
    class complex_long_double_ptr : private complex_long_double_ref {
  public:
    /**
     * Typically we are given a pointer to the data storing the complex_long_double and need
     * to construct a complex_long_double_ptr from it.
     * @param dat A pointer to the data
     */
    complex_long_double_ptr( long double* dat ) : complex_long_double_ref( dat ){}
    /**
     * Dereference the pointer
     * @return the complex_long_double object
     */
    complex_long_double_ref operator*(){ return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex_long_double object
     */
    complex_long_double_ref* operator->(){ return this; }
    /**
     * Dereference the pointer
     * @return the complex_long_double object
     */
    complex_long_double_ref const operator*() const { return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex_long_double object
     */
    complex_long_double_ref const* operator->() const { return this; }
  };

  /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex_long_double.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex_long_double
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream, complex_long_double const& z ){
    long double i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

    /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex_long_double_ref.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex_long_double_ref
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream,
					 complex_long_double_ref const& z ){
    long double i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

  /**
   * The @c gsl::cpx_long_double namespace allows the @c gsl::complex_long_double @c static funtions to be called from
   * a @c namespace and so allows the usual C++ convenience of using code like the following.
   * @code
   * using namespace gsl;
   * using namespace gsl::cpx_long_double;
   * complex_long_double z = rect( 3, 1 );
   * complex_long_double w = sin( z );
   * @endcode
   */
  namespace cpx_long_double {
    /**
     * C++ version of gsl_complex_long_double_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex_long_double
     */
    inline complex_long_double rect( long double x, long double y ){ return complex_long_double::rect( x, y ); }
    /**
     * C++ version of gsl_complex_long_double_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex_long_double
     */
    inline complex_long_double polar( long double r, long double theta ){ return complex_long_double::polar( r, theta ); }  }

}

#endif
