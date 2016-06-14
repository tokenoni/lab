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

#ifndef CCGSL_MATRIX_COMPLEX_FLOAT_HPP
#define CCGSL_MATRIX_COMPLEX_FLOAT_HPP

#include<gsl/gsl_matrix.h>
#include<new>

#include"exception.hpp"
#include"vector_complex_float.hpp"

// This file is a template
#define CCGSL_MTY 2

namespace gsl {
  /**
   * This class handles matrix_complex_float objects as shared handles. It models a random access container
   * so that STL functions work with matrix_complex_float.
   *
   * Note that matrix_complex_float_views are implemented as matrix_complex_float objects here.
   */
  class matrix_complex_float {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    matrix_complex_float(){
      ccgsl_pointer = 0;
      // just plausibly we could fail to allocate count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // Refines random access container
    // Refines assignable
    /**
     * The default constructor creates a new matrix_complex_float with n elements
     * @param n1 The number of rows in the matrix_complex_float
     * @param n2 The number of columns in the matrix_complex_float
     */
    explicit matrix_complex_float( size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_complex_float_alloc( n1, n2 );
      // just plausibly we could allocate matrix_complex_float but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_matrix_complex_float_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_matrix_complex_float. This is not usually a good idea. In this case
     * we should not use gsl_matrix_complex_float_free() to deallocate the memory.
     * @param v The matrix_complex_float
     */
    explicit matrix_complex_float( gsl_matrix_complex_float* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the matrix_complex_float. Use clone() if you want a full copy.
     * @param v The matrix_complex_float to copy.
     */
    matrix_complex_float( matrix_complex_float const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // matrix_complex_float is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The matrix_complex_float to copy
     */
    matrix_complex_float& operator=( matrix_complex_float const& v ){
      // first, possibly delete anything pointed to by this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_matrix_complex_float_free( ccgsl_pointer );
	delete count;
      }
      // Then copy
      ccgsl_pointer = v.ccgsl_pointer;
      count = v.count;
      ++*count; // block_complex_float is now shared.
      return *this;
    }
    // clone()
    /**
     * The clone function. Use this if you want a copy of the block_complex_float that does
     * not share the underlying data.
     * @return a new copy of this.
     */
    matrix_complex_float clone() const {
      matrix_complex_float copy( get()->size1, get()->size2 );
      // Now copy
      gsl_matrix_complex_float_memcpy( copy.get(), get() );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~matrix_complex_float(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_matrix_complex_float_free( ccgsl_pointer );
	delete count;
      }
    }
    // Sizes
    /**
     * The number of rows of the matrix_complex_float
     * @return The number of rows of the matrix_complex_float
     */
    size_t size1() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size1; }
    /**
     * The number of columns of the matrix_complex_float
     * @return The number of columns of the matrix_complex_float
     */
    size_t size2() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size2; }
    /**
     * Swap two matrix_complex_float objects. This works even if the matrix_complex_float objects have different sizes
     * because it swaps pointers.
     * @param m The matrix_complex_float to swap with @c this.
     */
    void swap( matrix_complex_float& m ){
      gsl_matrix_complex_float* tmp = ccgsl_pointer; ccgsl_pointer = m.ccgsl_pointer; m.ccgsl_pointer = tmp;
      size_t* tmp2 = count; count = m.count; m.count = tmp2;
    }
    // view operations
    /**
     * C++ version of gsl_matrix_complex_float_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix_complex_float submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_row().
     * @param i A row index
     * @return A row as a vector_complex_float
     */
    vector_complex_float row( size_t const i ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_row( get(), i ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_column().
     * @param j A column index
     * @return A column as a vector_complex_float
     */
    vector_complex_float column( size_t const j ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_column( get(), j ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_diagonal().
     * @return The principal diagonal as a vector_complex_float
     */
    vector_complex_float
    diagonal(){ gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_diagonal( get() ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_complex_float
     */
    vector_complex_float subdiagonal( size_t const k ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_subdiagonal( get(), k ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_complex_float
     */
    vector_complex_float superdiagonal( size_t const k ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_superdiagonal( get(), k ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector_complex_float
     */
    vector_complex_float subrow( size_t const i, size_t const offset, size_t const n ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_subrow( get(), i, offset, n ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector_complex_float
     */
    vector_complex_float subcolumn( size_t const j, size_t const offset, size_t const n ){
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_subcolumn( get(), j, offset, n ).vector;
      return vector_complex_float( w );
    }
   /**
     * C++ version of gsl_matrix_complex_float_view_array().
     * @param base An array of type float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_complex_float
     */
    static matrix_complex_float view_array( float* base, size_t const n1, size_t const n2 ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_view_array( base, n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_view_array_with_tda().
     * @param base An array of type float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_complex_float
     */
    static matrix_complex_float view_array_with_tda( float* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m =  gsl_matrix_complex_float_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_view_vector().
     * @param v A vector_complex_float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_complex_float
     */
    static matrix_complex_float view_vector( vector_complex_float& v, size_t const n1, size_t const n2 ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m =  gsl_matrix_complex_float_view_vector( v.get(), n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_view_vector_with_tda().
     * @param v A vector_complex_float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_complex_float
     */
    static matrix_complex_float view_vector_with_tda( vector_complex_float& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix_complex_float( m );
    }
    // const versions ...
    /**
     * C++ version of gsl_matrix_complex_float_const_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix_complex_float const const_submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ) const {
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_const_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_row().
     * @param i A row index
     * @return A row as a vector_complex_float
     */
    vector_complex_float const const_row( size_t const i ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_row( get(), i ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_column().
     * @param j A column index
     * @return A column as a vector_complex_float
     */
    vector_complex_float const const_column( size_t const j ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_column( get(), j ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_diagonal().
     * @return The principal diagonal as a vector_complex_float
     */
    vector_complex_float const const_diagonal() const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_diagonal( get() ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_complex_float
     */
    vector_complex_float const const_subdiagonal( size_t const k ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_subdiagonal( get(), k ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_complex_float
     */
    vector_complex_float const const_superdiagonal( size_t const k ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_superdiagonal( get(), k ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector_complex_float
     */
    vector_complex_float const const_subrow( size_t const i, size_t const offset, size_t const n ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_subrow( get(), i, offset, n ).vector;
      return vector_complex_float( w );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector_complex_float
     */
    vector_complex_float const const_subcolumn( size_t const j, size_t const offset, size_t const n ) const {
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_const_subcolumn( get(), j, offset, n ).vector;
      return vector_complex_float( w );
    }
   /**
     * C++ version of gsl_matrix_complex_float_const_view_array().
     * @param base An array of type float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_complex_float
     */
    static matrix_complex_float const const_view_array( float* base, size_t const n1, size_t const n2 ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_const_view_array( base, n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_view_array_with_tda().
     * @param base An array of type float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_complex_float
     */
    static matrix_complex_float const
    const_view_array_with_tda( float* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m =  gsl_matrix_complex_float_const_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_view_vector().
     * @param v A vector_complex_float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_complex_float
     */
    static matrix_complex_float const const_view_vector( vector_complex_float& v, size_t const n1, size_t const n2 ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m =  gsl_matrix_complex_float_const_view_vector( v.get(), n1, n2 ).matrix;
      return matrix_complex_float( m );
    }
    /**
     * C++ version of gsl_matrix_complex_float_const_view_vector_with_tda().
     * @param v A vector_complex_float
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_complex_float
     */
    static matrix_complex_float const
    const_view_vector_with_tda( vector_complex_float& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_complex_float* m = static_cast<gsl_matrix_complex_float*>( malloc( sizeof( gsl_matrix_complex_float ) ) );
      *m = gsl_matrix_complex_float_const_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix_complex_float( m );
    }
  private:
    /**
     * The shared pointer
     */
    gsl_matrix_complex_float* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_matrix_complex_float.
     * @return the gsl_matrix_complex_float
     */
    gsl_matrix_complex_float* get(){
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Get the gsl_matrix_complex_float.
     * @return the gsl_matrix_complex_float
     */
    gsl_matrix_complex_float const* get() const {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Find if this is the only object sharing the gsl_matrix_complex_float.
     * @return @c true or @c falses according as 
     * this is the only matrix_complex_float object sharing the gsl_matrix_complex_float
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many matrix_complex_float objects share this pointer.
     * @return the number of matrix_complex_float objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as this contains a pointer
     * to a gsl_matrix_complex_float
     */
    operator bool() const { return ccgsl_pointer != 0; }
    // GSL functions
    /**
     * C++ version of gsl_matrix_complex_float_calloc(). This constructs a matrix_complex_float object with entries
     * initialised to zero.
     * @param n1 The number of rows in the matrix_complex_float
     * @param n2 The number of columns in the matrix_complex_float
     * @return A matrix_complex_float initialised to zero
     */
    static matrix_complex_float calloc( size_t const n1, size_t const n2 ){ return matrix_complex_float( gsl_matrix_complex_float_calloc( n1, n2 ) ); }
    /**
     * C++ version of gsl_matrix_complex_float_set_zero().
     */
    void set_zero(){ gsl_matrix_complex_float_set_zero( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_set_all().
     * @param x The value to which all elements are set
     */
    void set_all( complex_float x ){ gsl_matrix_complex_float_set_all( get(), x ); }
    /**
     * C++ version of gsl_matrix_complex_float_memcpy().
     * @param src source matrix_complex_float
     * @return error code on failure
     */
    int memcpy( matrix_complex_float const& src ){ return gsl_matrix_complex_float_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_add().
     * @param b matrix_complex_float to add to this
     * @return error code on failure
     */
    int add( matrix_complex_float const& b ){ return gsl_matrix_complex_float_add( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_sub().
     * @param b matrix_complex_float to subtract from this
     * @return error code on failure
     */
    int sub( matrix_complex_float const& b ){ return gsl_matrix_complex_float_sub( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_scale().
     * @param x constant to multiply this by
     * @return error code on failure
     */
    int scale( complex_float const x ){ return gsl_matrix_complex_float_scale( get(), x ); }
    /**
     * C++ version of gsl_matrix_complex_float_add_constant().
     * @param x constant to add to each element of this
     * @return error code on failure
     */
    int add_constant( complex_float const x ){ return gsl_matrix_complex_float_add_constant( get(), x ); }
    /**
     * C++ version of gsl_matrix_complex_float_isnull().
     * @return @c +1 or @c 0 according as elements are all zero or not
     */
    int isnull() const { return gsl_matrix_complex_float_isnull( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_ispos().
     * @return @c +1 or @c 0 according as elements are all positive or not
     */
    int ispos() const { return gsl_matrix_complex_float_ispos( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_isneg().
     * @return @c +1 or @c 0 according as elements are all negative or not
     */
    int isneg() const { return gsl_matrix_complex_float_isneg( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_isnonneg().
     * @return @c +1 or @c 0 according as elements are all nonnegative or not
     */
    int isnonneg() const { return gsl_matrix_complex_float_isnonneg( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_get().
     * @param i index of row
     * @param j index of column
     * @return value of element
     */
    complex_float get( size_t const i, size_t const j ) const { return gsl_matrix_complex_float_get( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_complex_float_set().
     * @param i index of row
     * @param j index of column
     * @param x new value for element
     */
    void set( size_t const i, size_t const j, complex_float x ){ gsl_matrix_complex_float_set( get(), i, j, x ); }
    /**
     * C++ version of gsl_matrix_complex_float_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    complex_float_ptr ptr( size_t const i, size_t const j ){
      if( i >= ccgsl_pointer->size1 )
	gsl_error( "Index out of range", __FILE__, __LINE__, exception::GSL_EINVAL );
      if( j >= ccgsl_pointer->size2 )
	gsl_error( "Index out of range", __FILE__, __LINE__, exception::GSL_EINVAL );
      return complex_float_ptr( get()->data + CCGSL_MTY * (i * get()->tda + j ) ); }
    /**
     * C++ version of gsl_matrix_complex_float_const_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    complex_float_ptr const const_ptr( size_t const i, size_t const j ) const {
      if( i >= ccgsl_pointer->size1 )
	gsl_error( "Index out of range", __FILE__, __LINE__, exception::GSL_EINVAL );
      if( j >= ccgsl_pointer->size2 )
	gsl_error( "Index out of range", __FILE__, __LINE__, exception::GSL_EINVAL );
      return complex_float_ptr( get()->data + CCGSL_MTY * (i * get()->tda + j ) ); }
    /**
     * C++ version of gsl_matrix_complex_float_fread().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fread( FILE* stream ){ return gsl_matrix_complex_float_fread( stream, get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_fwrite().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fwrite( FILE* stream ) const { return gsl_matrix_complex_float_fwrite( stream, get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_fscanf().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fscanf( FILE* stream ){ return gsl_matrix_complex_float_fscanf( stream, get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_fprintf().
     * @param stream A C file stream
     * @param format %d, %e, %f or %g
     * @return error code on failure
     */
    int fprintf( FILE* stream, char const* format ) const {
      return gsl_matrix_complex_float_fprintf( stream, get(), format ); }
    /**
     * C++ version of gsl_matrix_complex_float_alloc_from_block().
     * @param b The block_complex_float
     * @param offset The offset within the block_complex_float
     * @param n1 The number of rows in the matrix_complex_float
     * @param n2 The number of columns in the matrix_complex_float
     * @param d2 undocumented
     */
    matrix_complex_float( block_complex_float& b, size_t const offset, size_t const n1, size_t const n2, size_t const d2 ){
      ccgsl_pointer = gsl_matrix_complex_float_alloc_from_block( b.get(), offset, n1, n2, d2 );
      // just plausibly we could allocate vector_complex_float but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_complex_float_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * C++ version of gsl_matrix_complex_float_alloc_from_matrix().
     * @param m The matrix_complex_float
     * @param k1 the row of @c m to take as row zero
     * @param k2 the column of @c m to take as column zero
     * @param n1 The number of rows in the matrix_complex_float
     * @param n2 The number of columns in the matrix_complex_float
     */
    matrix_complex_float( matrix_complex_float& m, size_t const k1, size_t const k2, size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_complex_float_alloc_from_matrix( m.get(), k1, k2, n1, n2 );
      // just plausibly we could allocate matrix_complex_float but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_complex_float_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // More functions
    /**
     * C++ version of gsl_matrix_complex_float_set_identity().
     */
    void set_identity(){ gsl_matrix_complex_float_set_identity( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_swap_rows().
     * @param i Index of first row
     * @param j Index of second row
     * @return error code on failure
     */
    int swap_rows( size_t const i, size_t const j ){ return gsl_matrix_complex_float_swap_rows( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_complex_float_swap_columns().
     * @param i Index of first column
     * @param j Index of second column
     * @return error code on failure
     */
    int swap_columns( size_t const i, size_t const j ){
      return gsl_matrix_complex_float_swap_columns( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_complex_float_swap_rowcol(). Swap row and column in place.
     * Matrix must be square.
     * @param i index of row
     * @param j index of column
     * @return error code on failure
     */
    int swap_rowcol( size_t const i, size_t const j ){ return gsl_matrix_complex_float_swap_rowcol( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_complex_float_transpose().
     * @return error code on failure.
     */
    int transpose(){ return gsl_matrix_complex_float_transpose( get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_transpose_memcpy().
     * @param src matrix_complex_float whose transpose it to be copied to @c this
     * @return error code on failure
     */
    int transpose_memcpy( matrix_complex_float const& src ){
      return gsl_matrix_complex_float_transpose_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_mul_elements(). Multiply matrices elementwise.
     * @param b Another matrix_complex_float
     * @return error code on failure
     */
    int
    mul_elements( matrix_complex_float const& b ){
          return gsl_matrix_complex_float_mul_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_div_elements(). Divide each element of @c this by the
     * corrsponding element of @c b
     * @param b Another matrix_complex_float
     * @return error code on failure
     */
    int div_elements( matrix_complex_float const& b ){
      return gsl_matrix_complex_float_div_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_add_diagonal().
     * @param x A constant
     * @return error code on failure
     */
    int add_diagonal( complex_float const x ){
      return gsl_matrix_complex_float_add_diagonal( get(), x ); }
    /**
     * C++ version of gsl_matrix_complex_float_get_row().
     * @param v A vector_complex_float
     * @param i The index of the row
     * @return error code on failure
     */
    int get_row( vector_complex_float& v, size_t const i ) const {
      return gsl_matrix_complex_float_get_row( v.get(), get(), i ); }
    /**
     * C++ version of gsl_matrix_complex_float_get_col().
     * @param v A vector_complex_float
     * @param j The index of the column
     * @return error code on failure
     */
    int get_col( vector_complex_float& v, size_t const j ) const {
      return gsl_matrix_complex_float_get_col( v.get(), get(), j ); }
    /**
     * C++ version of gsl_matrix_complex_float_set_row().
     * @param i The index of the row
     * @param v A vector_complex_float
     * @return error code on failure
     */
    int set_row( size_t const i, vector_complex_float const& v ){
      return gsl_matrix_complex_float_set_row( get(), i, v.get() ); }
    /**
     * C++ version of gsl_matrix_complex_float_set_col().
     * @param j The index of the column
     * @param v A vector_complex_float
     * @return error code on failure
     */
    int set_col( size_t const j, vector_complex_float const& v ){
      return gsl_matrix_complex_float_set_col( get(), j, v.get() ); }
    // Extra functions for []
    /**
     * This function allows us to use a matrix_complex_float like an array. Use with caution.
     * Although @c matrix_complex_float[i][j] is possible, it is much less efficient than
     * matrix_complex_float::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector_complex_float representing a row
     */
    vector_complex_float operator[]( size_t const i ){
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix_complex_float is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector_complex_float();
      }
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_row( ccgsl_pointer, i ).vector;
      return vector_complex_float( w );
    }
    /**
     * This function allows us to use a matrix_complex_float like an array. Use with caution.
     * Although @c matrix_complex_float[i][j] is possible, it is much less efficient than
     * matrix_complex_float::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector_complex_float representing a row
     */
    vector_complex_float const operator[]( size_t const i ) const {
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix_complex_float is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector_complex_float();
      }
      gsl_vector_complex_float* w = static_cast<gsl_vector_complex_float*>( malloc( sizeof( gsl_vector_complex_float ) ) );
      *w = gsl_matrix_complex_float_row( ccgsl_pointer, i ).vector;
      return vector_complex_float( w );
    }
  };

  // Extra functions for vector_complex_float allocation from matrix_complex_float objects
  vector_complex_float vector_complex_float::alloc_row_from_matrix( matrix_complex_float& m, size_t const i ){
    return vector_complex_float ( gsl_vector_complex_float_alloc_row_from_matrix( m.get(), i ) ); }
   vector_complex_float vector_complex_float::alloc_col_from_matrix( matrix_complex_float& m, size_t const i ){
    return vector_complex_float ( gsl_vector_complex_float_alloc_col_from_matrix( m.get(), i ) ); }
}
#undef CCGSL_MTY
#endif
