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

#ifndef CCGSL_COMPLEX_FLOAT_HPP
#define CCGSL_COMPLEX_FLOAT_HPP

#include<cmath>
#include<cstring>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

namespace gsl {
  /**
   * This class handles complex_float numbers. There is an ambiguity in translating gsl_complex_float 
   * to C++: @c gsl_complex_float is a @c struct and so should be represented by a C++ @c class.
   * But @c gsl_complex_float is also used like a @c namespace and so it would be nice to have
   * a C++ @c namespace @c gsl:: complex_float. However, we cannot have both. The compromise is to
   * use @c static @c complex_float @c class functions. But, since we cannot use
   * @code
   * using namespace gsl::complex_float
   * @endcode
   * and similar namespace operations with @c gsl::complex_float, we also define a namespace
   * @c gsl::cpx_float with the same functions, again defined @c inline for efficiency.
   */
  class complex_float : private gsl_complex_float {
    friend class complex_float_ptr;
    friend class complex_float_ref;
    friend class vector_complex_float;
    friend class matrix_complex_float;
 public:
    // /**
    //  * Make sure this is convertible to @c gsl_complex_float. This is not ideal because it copies
    //  * twice. But it gets round the problem that a @c gsl_comple must contain constant data.
    //  * @return A @c gsl_complex_float
    //  */ 
    // operator gsl_complex_float() const {
    //   gsl_complex_float z; memcpy( z.dat, dat, 2 * sizeof( float ) ); return z; }
    // Default constructible: so that it can  be stored in a container
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex_float(){}
    /**
     * Allow construction from raw data. This is mainly for use by the container classes.
     * @param dat Raw data
     */
    explicit complex_float( float* dat ){
      memcpy( this->dat, dat, 2 * sizeof( float ) );
    }
    // copy constructor: default is OK
    // assignment operator: default is OK
    // destructor: default is OK
    /**
     * Constructor from base class.
     * @param z A gsl_complex_float
     */
    complex_float( gsl_complex_float const& z ){
      GSL_SET_REAL( this, GSL_REAL( z ) );
      GSL_SET_IMAG( this, GSL_IMAG( z ) );
    }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex_float& get(){ return dynamic_cast<gsl_complex_float&>( *this ); }
    /**
     * Get the base class object.
     * @return the base class object.
     */
    gsl_complex_float const& get() const { return dynamic_cast<gsl_complex_float const&>( *this ); }
    /**
     * C++ version of gsl_complex_float_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex_float
     */
    static complex_float rect( float x, float y ){ complex_float z; GSL_SET_COMPLEX( &z, x, y ); return z; }
    /**
     * C++ version of gsl_complex_float_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex_float
     */
    static complex_float polar( float r, float theta ){ complex_float z;
      GSL_SET_COMPLEX( &z, r * std::cos( theta ), r * std::sin( theta ) ); return z; }
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex_float( float x, float y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    float real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    float imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( float x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( float y ){ GSL_SET_IMAG( this, y ); }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the equality operator.
     * @param z Another complex_float
     * @return @c true or @c false according as @c this equals @c z
     */
    bool operator==( complex_float const& z ) const { return dat[0] == z.dat[0] and dat[1] == z.dat[1]; }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the inequality operator.
     * @param z Another complex_float
     * @return @c false or @c true according as @c this equals @c z
     */
    bool operator!=( complex_float const& z ) const { return dat[0] != z.dat[0] or dat[1] != z.dat[1]; }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than operator.
     * @param z Another complex_float
     * @return @c true or @c false according as @c this is less than @c z
     */
    bool operator<( complex_float const& z ) const {
      return dat[0] < z.dat[0] or (dat[0] == z.dat[0] and dat[1] < z.dat[1]); }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than operator.
     * @param z Another complex_float
     * @return @c true or @c false according as @c this is greater than than @c z
     */
    bool operator>( complex_float const& z ) const {
      return dat[0] > z.dat[0] or (dat[0] == z.dat[0] and dat[1] > z.dat[1]); }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) less than or equal to operator.
     * @param z Another complex_float
     * @return @c true or @c false according as @c this is less than or equal to  @c z
     */
    bool operator<=( complex_float const& z ) const { return *this < z or *this == z; }
    /**
     * A complex_float object must be less than comparable so that it can be  used as a container
     * value type. This is the (lexicographic) greater than or equal to operator.
     * @param z Another complex_float
     * @return @c true or @c false according as @c this is greater than or equal to  @c z
     */
    bool operator>=( complex_float const& z ) const { return *this > z or *this == z; }
  };
    
  /**
   * This class can be used like a reference for complex_float objects so that we can
   * iterate over a vector (for example) of them. 
   */
  class complex_float_ref {
  protected:
    /**
     * The data. Never initialised or destroyed within the class. Always precisely two elements.
     */
    float* dat;
  public:
    /**
     * We use this in constructing complex_float_ptr objects.
     * @param dat A pointer to the data
     */
    complex_float_ref( float* dat ) : dat( dat ){}
   /**
     * Make sure this is convertible to @c gsl_complex_float. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex_float must contain constant data.
     * @return A @c gsl_complex_float
     */ 
    operator gsl_complex_float() const {
      gsl_complex_float z; memcpy( z.dat, dat, 2 * sizeof( float ) ); return z; }
    /**
     * Make sure this is convertible to @c complex_float. This is not ideal because it copies
     * twice. But it gets round the problem that a @c gsl_complex_float must contain constant data.
     * @return A @c complex_float
     */ 
    operator complex_float() const {
      complex_float z; memcpy( z.dat, dat, 2 * sizeof( float ) ); return z; }
    // Default constructible
    /**
     * The default constructor is only really useful for assigning to.
     */
    complex_float_ref() : dat( 0 ){}
    /**
     * Make sure we can construct from a complex_float.
     * @param z The complex_float to construct from
     */
    complex_float_ref( complex_float& z ) : dat( z.dat ){}
    // copy constructor: default is OK
    // assignment operator: default is OK
    // also need to copy from complex_float
    /**
     * Assignment from complex_float.
     * @param z The complex_float to assign from
     */
    complex_float_ref& operator=( complex_float const& z ){
      memcpy( dat, z.dat, 2 * sizeof( float ) ); return *this; }
    // destructor: default is OK
    /**
     * C++ version of GSL_SET_COMPLEX().
     * @param x Real part 
     * @param y Imaginary part
     */
    void set_complex_float( float x, float y ){ GSL_SET_COMPLEX( this, x, y );}
    /**
     * C++ version of GSL_REAL().
     * @return The real part of @c this
     */
    float real() const { return GSL_REAL( *this ); }
    /**
     * C++ version of GSL_IMAG().
     * @return The real part of @c this
     */
    float imag() const { return GSL_IMAG( *this ); }
    /**
     * C++ version of GSL_SET_REAL().
     * @param x The new real part
     */
    void set_real( float x ){ GSL_SET_REAL( this, x ); }
    /**
     * C++ version of GSL_SET_IMAG().
     * @param y The new imaginary part
     */
    void set_imag( float y ){ GSL_SET_IMAG( this, y ); }
  };

  /**
   * This class can be used like a pointer for complex_float objects so that we can
   * iterate over a vector (for example) of them.
   */
    class complex_float_ptr : private complex_float_ref {
  public:
    /**
     * Typically we are given a pointer to the data storing the complex_float and need
     * to construct a complex_float_ptr from it.
     * @param dat A pointer to the data
     */
    complex_float_ptr( float* dat ) : complex_float_ref( dat ){}
    /**
     * Dereference the pointer
     * @return the complex_float object
     */
    complex_float_ref operator*(){ return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex_float object
     */
    complex_float_ref* operator->(){ return this; }
    /**
     * Dereference the pointer
     * @return the complex_float object
     */
    complex_float_ref const operator*() const { return *this; }
    /**
     * Dereference the pointer
     * @return a pointer to the complex_float object
     */
    complex_float_ref const* operator->() const { return this; }
  };

  /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex_float.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex_float
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream, complex_float const& z ){
    float i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

    /**
   * Define the @c << @c operator so that we can use expressions like
   * @code
   * std::cout << z << std::endl;
   * @endcode
   * where @c z is a complex_float_ref.
   * @tparam Ch The character class
   * @tparam Tr The character traits class
   * @param stream An output stream
   * @param z A complex_float_ref
   * @return A reference to the output stream
   */
  template<typename Ch,typename Tr>
  std::basic_ostream<Ch,Tr>& operator<<( std::basic_ostream<Ch,Tr>& stream,
					 complex_float_ref const& z ){
    float i = z.imag();
    if( i >= 0 ) stream << z.real() << "+" << i << "i";
    else stream << z.real() << "-" << -i << "i";
    return stream;
  }

  /**
   * The @c gsl::cpx_float namespace allows the @c gsl::complex_float @c static funtions to be called from
   * a @c namespace and so allows the usual C++ convenience of using code like the following.
   * @code
   * using namespace gsl;
   * using namespace gsl::cpx_float;
   * complex_float z = rect( 3, 1 );
   * complex_float w = sin( z );
   * @endcode
   */
  namespace cpx_float {
    /**
     * C++ version of gsl_complex_float_rect().
     * @param x Real part 
     * @param y Imaginary part
     * @return complex_float
     */
    inline complex_float rect( float x, float y ){ return complex_float::rect( x, y ); }
    /**
     * C++ version of gsl_complex_float_rect().
     * @param r The magnitude
     * @param theta The argument
     * @return complex_float
     */
    inline complex_float polar( float r, float theta ){ return complex_float::polar( r, theta ); }
  }

}

#endif
