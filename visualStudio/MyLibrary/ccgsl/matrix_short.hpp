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

#ifndef CCGSL_MATRIX_SHORT_HPP
#define CCGSL_MATRIX_SHORT_HPP

#include<gsl/gsl_matrix.h>
#include<new>

#include"exception.hpp"
#include"vector_short.hpp"

// This file is autogenerated

namespace gsl {
  /**
   * This class handles matrix_short objects as shared handles. It models a random access container
   * so that STL functions work with matrix_short.
   *
   * Note that matrix_short_views are implemented as matrix_short objects here.
   */
  class matrix_short {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    matrix_short(){
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
     * The default constructor creates a new matrix_short with n elements
     * @param n1 The number of rows in the matrix_short
     * @param n2 The number of columns in the matrix_short
     */
    explicit matrix_short( size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_short_alloc( n1, n2 );
      // just plausibly we could allocate matrix_short but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_matrix_short_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_matrix_short. This is not usually a good idea. In this case
     * we should not use gsl_matrix_short_free() to deallocate the memory.
     * @param v The matrix_short
     */
    explicit matrix_short( gsl_matrix_short* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the matrix_short. Use clone() if you want a full copy.
     * @param v The matrix_short to copy.
     */
    matrix_short( matrix_short const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // matrix_short is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The matrix_short to copy
     */
    matrix_short& operator=( matrix_short const& v ){
      // first, possibly delete anything pointed to by this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_matrix_short_free( ccgsl_pointer );
	delete count;
      }
      // Then copy
      ccgsl_pointer = v.ccgsl_pointer;
      count = v.count;
      ++*count; // block_short is now shared.
      return *this;
    }
    // clone()
    /**
     * The clone function. Use this if you want a copy of the block_short that does
     * not share the underlying data.
     * @return a new copy of this.
     */
    matrix_short clone() const {
      matrix_short copy( get()->size1, get()->size2 );
      // Now copy
      gsl_matrix_short_memcpy( copy.get(), get() );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~matrix_short(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_matrix_short_free( ccgsl_pointer );
	delete count;
      }
    }
    // Sizes
    /**
     * The number of rows of the matrix_short
     * @return The number of rows of the matrix_short
     */
    size_t size1() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size1; }
    /**
     * The number of columns of the matrix_short
     * @return The number of columns of the matrix_short
     */
    size_t size2() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size2; }
    /**
     * Swap two matrix_short objects. This works even if the matrix_short objects have different sizes
     * because it swaps pointers.
     * @param m The matrix_short to swap with @c this.
     */
    void swap( matrix_short& m ){
      gsl_matrix_short* tmp = ccgsl_pointer; ccgsl_pointer = m.ccgsl_pointer; m.ccgsl_pointer = tmp;
      size_t* tmp2 = count; count = m.count; m.count = tmp2;
    }
    // view operations
    /**
     * C++ version of gsl_matrix_short_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix_short submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_row().
     * @param i A row index
     * @return A row as a vector_short
     */
    vector_short row( size_t const i ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_row( get(), i ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_column().
     * @param j A column index
     * @return A column as a vector_short
     */
    vector_short column( size_t const j ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_column( get(), j ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_diagonal().
     * @return The principal diagonal as a vector_short
     */
    vector_short
    diagonal(){ gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_diagonal( get() ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_short
     */
    vector_short subdiagonal( size_t const k ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_subdiagonal( get(), k ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_short
     */
    vector_short superdiagonal( size_t const k ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_superdiagonal( get(), k ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector_short
     */
    vector_short subrow( size_t const i, size_t const offset, size_t const n ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_subrow( get(), i, offset, n ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector_short
     */
    vector_short subcolumn( size_t const j, size_t const offset, size_t const n ){
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_subcolumn( get(), j, offset, n ).vector;
      return vector_short( w );
    }
   /**
     * C++ version of gsl_matrix_short_view_array().
     * @param base An array of type short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_short
     */
    static matrix_short view_array( short* base, size_t const n1, size_t const n2 ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_view_array( base, n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_view_array_with_tda().
     * @param base An array of type short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_short
     */
    static matrix_short view_array_with_tda( short* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m =  gsl_matrix_short_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_view_vector().
     * @param v A vector_short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_short
     */
    static matrix_short view_vector( vector_short& v, size_t const n1, size_t const n2 ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m =  gsl_matrix_short_view_vector( v.get(), n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_view_vector_with_tda().
     * @param v A vector_short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_short
     */
    static matrix_short view_vector_with_tda( vector_short& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix_short( m );
    }
    // const versions ...
    /**
     * C++ version of gsl_matrix_short_const_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix_short const const_submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ) const {
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_const_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_const_row().
     * @param i A row index
     * @return A row as a vector_short
     */
    vector_short const const_row( size_t const i ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_row( get(), i ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_column().
     * @param j A column index
     * @return A column as a vector_short
     */
    vector_short const const_column( size_t const j ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_column( get(), j ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_diagonal().
     * @return The principal diagonal as a vector_short
     */
    vector_short const const_diagonal() const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_diagonal( get() ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_short
     */
    vector_short const const_subdiagonal( size_t const k ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_subdiagonal( get(), k ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector_short
     */
    vector_short const const_superdiagonal( size_t const k ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_superdiagonal( get(), k ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector_short
     */
    vector_short const const_subrow( size_t const i, size_t const offset, size_t const n ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_subrow( get(), i, offset, n ).vector;
      return vector_short( w );
    }
    /**
     * C++ version of gsl_matrix_short_const_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector_short
     */
    vector_short const const_subcolumn( size_t const j, size_t const offset, size_t const n ) const {
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_const_subcolumn( get(), j, offset, n ).vector;
      return vector_short( w );
    }
   /**
     * C++ version of gsl_matrix_short_const_view_array().
     * @param base An array of type short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_short
     */
    static matrix_short const const_view_array( short* base, size_t const n1, size_t const n2 ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_const_view_array( base, n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_const_view_array_with_tda().
     * @param base An array of type short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_short
     */
    static matrix_short const
    const_view_array_with_tda( short* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m =  gsl_matrix_short_const_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_const_view_vector().
     * @param v A vector_short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix_short
     */
    static matrix_short const const_view_vector( vector_short& v, size_t const n1, size_t const n2 ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m =  gsl_matrix_short_const_view_vector( v.get(), n1, n2 ).matrix;
      return matrix_short( m );
    }
    /**
     * C++ version of gsl_matrix_short_const_view_vector_with_tda().
     * @param v A vector_short
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix_short
     */
    static matrix_short const
    const_view_vector_with_tda( vector_short& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix_short* m = static_cast<gsl_matrix_short*>( malloc( sizeof( gsl_matrix_short ) ) );
      *m = gsl_matrix_short_const_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix_short( m );
    }
  private:
    /**
     * The shared pointer
     */
    gsl_matrix_short* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_matrix_short.
     * @return the gsl_matrix_short
     */
    gsl_matrix_short* get(){
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Get the gsl_matrix_short.
     * @return the gsl_matrix_short
     */
    gsl_matrix_short const* get() const {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Find if this is the only object sharing the gsl_matrix_short.
     * @return @c true or @c falses according as 
     * this is the only matrix_short object sharing the gsl_matrix_short
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many matrix_short objects share this pointer.
     * @return the number of matrix_short objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as this contains a pointer
     * to a gsl_matrix_short
     */
    operator bool() const { return ccgsl_pointer != 0; }
    // GSL functions
    /**
     * C++ version of gsl_matrix_short_calloc(). This constructs a matrix_short object with entries
     * initialised to zero.
     * @param n1 The number of rows in the matrix_short
     * @param n2 The number of columns in the matrix_short
     * @return A matrix_short initialised to zero
     */
    static matrix_short calloc( size_t const n1, size_t const n2 ){ return matrix_short( gsl_matrix_short_calloc( n1, n2 ) ); }
    /**
     * C++ version of gsl_matrix_short_set_zero().
     */
    void set_zero(){ gsl_matrix_short_set_zero( get() ); }
    /**
     * C++ version of gsl_matrix_short_set_all().
     * @param x The value to which all elements are set
     */
    void set_all( short x ){ gsl_matrix_short_set_all( get(), x ); }
    /**
     * C++ version of gsl_matrix_short_memcpy().
     * @param src source matrix_short
     * @return error code on failure
     */
    int memcpy( matrix_short const& src ){ return gsl_matrix_short_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_short_max().
     * @return maximum element of matrix_short
     */
    short max() const { return gsl_matrix_short_max( get() ); }
    /**
     * C++ version of gsl_matrix_short_min().
     * @return minimum element of matrix_short
     */
    short min() const { return gsl_matrix_short_min( get() ); }
    /**
     * C++ version of gsl_matrix_short_minmax().
     * @param min_out minimum element of matrix_short
     * @param max_out maximum element of matrix_short
     */
    void minmax( short* min_out, short* max_out ) const {
      gsl_matrix_short_minmax( get(), min_out, max_out ); }
    /**
     * C++ version of gsl_matrix_short_add().
     * @param b matrix_short to add to this
     * @return error code on failure
     */
    int add( matrix_short const& b ){ return gsl_matrix_short_add( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_short_sub().
     * @param b matrix_short to subtract from this
     * @return error code on failure
     */
    int sub( matrix_short const& b ){ return gsl_matrix_short_sub( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_short_scale().
     * @param x constant to multiply this by
     * @return error code on failure
     */
    int scale( short const x ){ return gsl_matrix_short_scale( get(), x ); }
    /**
     * C++ version of gsl_matrix_short_add_constant().
     * @param x constant to add to each element of this
     * @return error code on failure
     */
    int add_constant( short const x ){ return gsl_matrix_short_add_constant( get(), x ); }
    /**
     * C++ version of gsl_matrix_short_isnull().
     * @return @c +1 or @c 0 according as elements are all zero or not
     */
    int isnull() const { return gsl_matrix_short_isnull( get() ); }
    /**
     * C++ version of gsl_matrix_short_ispos().
     * @return @c +1 or @c 0 according as elements are all positive or not
     */
    int ispos() const { return gsl_matrix_short_ispos( get() ); }
    /**
     * C++ version of gsl_matrix_short_isneg().
     * @return @c +1 or @c 0 according as elements are all negative or not
     */
    int isneg() const { return gsl_matrix_short_isneg( get() ); }
    /**
     * C++ version of gsl_matrix_short_isnonneg().
     * @return @c +1 or @c 0 according as elements are all nonnegative or not
     */
    int isnonneg() const { return gsl_matrix_short_isnonneg( get() ); }
    /**
     * C++ version of gsl_matrix_short_get().
     * @param i index of row
     * @param j index of column
     * @return value of element
     */
    short get( size_t const i, size_t const j ) const { return gsl_matrix_short_get( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_set().
     * @param i index of row
     * @param j index of column
     * @param x new value for element
     */
    void set( size_t const i, size_t const j, short x ){ gsl_matrix_short_set( get(), i, j, x ); }
    /**
     * C++ version of gsl_matrix_short_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    short* ptr( size_t const i, size_t const j ){ return gsl_matrix_short_ptr( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_const_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    short const* const_ptr( size_t const i, size_t const j ) const {
      return gsl_matrix_short_const_ptr( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_fread().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fread( FILE* stream ){ return gsl_matrix_short_fread( stream, get() ); }
    /**
     * C++ version of gsl_matrix_short_fwrite().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fwrite( FILE* stream ) const { return gsl_matrix_short_fwrite( stream, get() ); }
    /**
     * C++ version of gsl_matrix_short_fscanf().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fscanf( FILE* stream ){ return gsl_matrix_short_fscanf( stream, get() ); }
    /**
     * C++ version of gsl_matrix_short_fprintf().
     * @param stream A C file stream
     * @param format %d, %e, %f or %g
     * @return error code on failure
     */
    int fprintf( FILE* stream, char const* format ) const {
      return gsl_matrix_short_fprintf( stream, get(), format ); }
    /**
     * C++ version of gsl_matrix_short_alloc_from_block().
     * @param b The block_short
     * @param offset The offset within the block_short
     * @param n1 The number of rows in the matrix_short
     * @param n2 The number of columns in the matrix_short
     * @param d2 undocumented
     */
    matrix_short( block_short& b, size_t const offset, size_t const n1, size_t const n2, size_t const d2 ){
      ccgsl_pointer = gsl_matrix_short_alloc_from_block( b.get(), offset, n1, n2, d2 );
      // just plausibly we could allocate vector_short but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_short_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * C++ version of gsl_matrix_short_alloc_from_matrix().
     * @param m The matrix_short
     * @param k1 the row of @c m to take as row zero
     * @param k2 the column of @c m to take as column zero
     * @param n1 The number of rows in the matrix_short
     * @param n2 The number of columns in the matrix_short
     */
    matrix_short( matrix_short& m, size_t const k1, size_t const k2, size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_short_alloc_from_matrix( m.get(), k1, k2, n1, n2 );
      // just plausibly we could allocate matrix_short but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_short_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // More functions
    /**
     * C++ version of gsl_matrix_short_set_identity().
     */
    void set_identity(){ gsl_matrix_short_set_identity( get() ); }
    /**
     * C++ version of gsl_matrix_short_swap_rows().
     * @param i Index of first row
     * @param j Index of second row
     * @return error code on failure
     */
    int swap_rows( size_t const i, size_t const j ){ return gsl_matrix_short_swap_rows( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_swap_columns().
     * @param i Index of first column
     * @param j Index of second column
     * @return error code on failure
     */
    int swap_columns( size_t const i, size_t const j ){
      return gsl_matrix_short_swap_columns( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_swap_rowcol(). Swap row and column in place. Matrix must be square.
     * @param i index of row
     * @param j index of column
     * @return error code on failure
     */
    int swap_rowcol( size_t const i, size_t const j ){ return gsl_matrix_short_swap_rowcol( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_short_transpose().
     * @return error code on failure.
     */
    int transpose(){ return gsl_matrix_short_transpose( get() ); }
    /**
     * C++ version of gsl_matrix_short_transpose_memcpy().
     * @param src matrix_short whose transpose it to be copied to @c this
     * @return error code on failure
     */
    int transpose_memcpy( matrix_short const& src ){
      return gsl_matrix_short_transpose_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_short_max_index().
     * @param imax row index of the first maximum element in the matrix_short
     * @param jmax column index of the first maximum element in the matrix_short
     */
    void max_index( size_t* imax, size_t* jmax ) const {
      gsl_matrix_short_max_index( get(), imax, jmax ); }
    /**
     * C++ version of gsl_matrix_short_min_index().
     * @param imin row index of the first minimum element in the matrix_short
     * @param jmin column index of the first minimum element in the matrix_short
     */
    void min_index( size_t* imin, size_t* jmin ) const {
      gsl_matrix_short_min_index( get(), imin, jmin ); }
    /**
     * C++ version of gsl_matrix_short_minmax_index().
     * @param imin row index of the first minimum element in the matrix_short
     * @param jmin column index of the first minimum element in the matrix_short
     * @param imax row index of the first maximum element in the matrix_short
     * @param jmax column index of the first maximum element in the matrix_short
     */
    void minmax_index( size_t* imin, size_t* jmin, size_t* imax, size_t* jmax ) const {
      gsl_matrix_short_minmax_index( get(), imin, jmin, imax, jmax ); }
    /**
     * C++ version of gsl_matrix_short_mul_elements(). Multiply matrices elementwise.
     * @param b Another matrix_short
     * @return error code on failure
     */
    int
    mul_elements( matrix_short const& b ){
          return gsl_matrix_short_mul_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_short_div_elements(). Divide each element of @c this by the
     * corrsponding element of @c b
     * @param b Another matrix_short
     * @return error code on failure
     */
    int div_elements( matrix_short const& b ){
      return gsl_matrix_short_div_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_short_add_diagonal().
     * @param x A constant
     * @return error code on failure
     */
    int add_diagonal( short const x ){
      return gsl_matrix_short_add_diagonal( get(), x ); }
    /**
     * C++ version of gsl_matrix_short_get_row().
     * @param v A vector_short
     * @param i The index of the row
     * @return error code on failure
     */
    int get_row( vector_short& v, size_t const i ) const {
      return gsl_matrix_short_get_row( v.get(), get(), i ); }
    /**
     * C++ version of gsl_matrix_short_get_col().
     * @param v A vector_short
     * @param j The index of the column
     * @return error code on failure
     */
    int get_col( vector_short& v, size_t const j ) const {
      return gsl_matrix_short_get_col( v.get(), get(), j ); }
    /**
     * C++ version of gsl_matrix_short_set_row().
     * @param i The index of the row
     * @param v A vector_short
     * @return error code on failure
     */
    int set_row( size_t const i, vector_short const& v ){
      return gsl_matrix_short_set_row( get(), i, v.get() ); }
    /**
     * C++ version of gsl_matrix_short_set_col().
     * @param j The index of the column
     * @param v A vector_short
     * @return error code on failure
     */
    int set_col( size_t const j, vector_short const& v ){
      return gsl_matrix_short_set_col( get(), j, v.get() ); }
    // Extra functions for []
    /**
     * This function allows us to use a matrix_short like an array. Use with caution.
     * Although @c matrix_short[i][j] is possible, it is much less efficient than
     * matrix_short::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector_short representing a row
     */
    vector_short operator[]( size_t const i ){
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix_short is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector_short();
      }
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_row( ccgsl_pointer, i ).vector;
      return vector_short( w );
    }
    /**
     * This function allows us to use a matrix_short like an array. Use with caution.
     * Although @c matrix_short[i][j] is possible, it is much less efficient than
     * matrix_short::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector_short representing a row
     */
    vector_short const operator[]( size_t const i ) const {
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix_short is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector_short();
      }
      gsl_vector_short* w = static_cast<gsl_vector_short*>( malloc( sizeof( gsl_vector_short ) ) );
      *w = gsl_matrix_short_row( ccgsl_pointer, i ).vector;
      return vector_short( w );
    }
  };

  // Extra functions for vector_short allocation from matrix_short objects
  vector_short vector_short::alloc_row_from_matrix( matrix_short& m, size_t const i ){
    return vector_short ( gsl_vector_short_alloc_row_from_matrix( m.get(), i ) ); }
   vector_short vector_short::alloc_col_from_matrix( matrix_short& m, size_t const i ){
    return vector_short ( gsl_vector_short_alloc_col_from_matrix( m.get(), i ) ); }
}
#endif