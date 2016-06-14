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

#ifndef CCGSL_MATRIX_HPP
#define CCGSL_MATRIX_HPP

#include<gsl/gsl_matrix.h>
#include<new>

#include"exception.hpp"
#include"vector.hpp"

// This file is used as a template

namespace gsl {
  /**
   * This class handles matrix objects as shared handles. It models a random access container
   * so that STL functions work with matrix.
   *
   * Note that matrix_views are implemented as matrix objects here.
   */
  class matrix {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    matrix(){
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
     * This constructor creates a new matrix with n1 rows and n2 columns
     * @param n1 The number of rows in the matrix
     * @param n2 The number of columns in the matrix
     */
    explicit matrix( size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_alloc( n1, n2 );
      // just plausibly we could allocate matrix but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_matrix_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_matrix. This is not usually a good idea. In this case
     * we should not use gsl_matrix_free() to deallocate the memory.
     * @param v The matrix
     */
    explicit matrix( gsl_matrix* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the matrix. Use clone() if you want a full copy.
     * @param v The matrix to copy.
     */
    matrix( matrix const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // matrix is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The matrix to copy
     */
    matrix& operator=( matrix const& v ){
      // first, possibly delete anything pointed to by this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_matrix_free( ccgsl_pointer );
	delete count;
      }
      // Then copy
      ccgsl_pointer = v.ccgsl_pointer;
      count = v.count;
      ++*count; // block is now shared.
      return *this;
    }
    // clone()
    /**
     * The clone function. Use this if you want a copy of the block that does
     * not share the underlying data.
     * @return a new copy of this.
     */
    matrix clone() const {
      matrix copy( get()->size1, get()->size2 );
      // Now copy
      gsl_matrix_memcpy( copy.get(), get() );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~matrix(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_matrix_free( ccgsl_pointer );
	delete count;
      }
    }
    // Sizes
    /**
     * The number of rows of the matrix
     * @return The number of rows of the matrix
     */
    size_t size1() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size1; }
    /**
     * The number of columns of the matrix
     * @return The number of columns of the matrix
     */
    size_t size2() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size2; }
    /**
     * Swap two matrix objects. This works even if the matrix objects have different sizes
     * because it swaps pointers.
     * @param m The matrix to swap with @c this.
     */
    void swap( matrix& m ){
      gsl_matrix* tmp = ccgsl_pointer; ccgsl_pointer = m.ccgsl_pointer; m.ccgsl_pointer = tmp;
      size_t* tmp2 = count; count = m.count; m.count = tmp2;
    }
    // view operations
    /**
     * C++ version of gsl_matrix_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_row().
     * @param i A row index
     * @return A row as a vector
     */
    vector row( size_t const i ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_row( get(), i ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_column().
     * @param j A column index
     * @return A column as a vector
     */
    vector column( size_t const j ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_column( get(), j ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_diagonal().
     * @return The principal diagonal as a vector
     */
    vector
    diagonal(){ gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_diagonal( get() ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector
     */
    vector subdiagonal( size_t const k ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_subdiagonal( get(), k ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector
     */
    vector superdiagonal( size_t const k ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_superdiagonal( get(), k ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector
     */
    vector subrow( size_t const i, size_t const offset, size_t const n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_subrow( get(), i, offset, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector
     */
    vector subcolumn( size_t const j, size_t const offset, size_t const n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_subcolumn( get(), j, offset, n ).vector;
      return vector( w );
    }
   /**
     * C++ version of gsl_matrix_view_array().
     * @param base An array of type double
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix
     */
    static matrix view_array( double* base, size_t const n1, size_t const n2 ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_view_array( base, n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_view_array_with_tda().
     * @param base An array of type double
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix
     */
    static matrix view_array_with_tda( double* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m =  gsl_matrix_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_view_vector().
     * @param v A vector
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix
     */
    static matrix view_vector( vector& v, size_t const n1, size_t const n2 ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m =  gsl_matrix_view_vector( v.get(), n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_view_vector_with_tda().
     * @param v A vector
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix
     */
    static matrix view_vector_with_tda( vector& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix( m );
    }
    // const versions ...
    /**
     * C++ version of gsl_matrix_const_submatrix().
     * @param i Index in @c this of first row of submatrix
     * @param j Index in @c this of first column of submatrix
     * @param n1 Number of rows of submatrix
     * @param n2 Number of columns of submatrix
     * @return The submatrix
     */
    matrix const const_submatrix( size_t const i, size_t const j, size_t const n1, size_t const n2 ) const {
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_const_submatrix( get(), i, j, n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_const_row().
     * @param i A row index
     * @return A row as a vector
     */
    vector const const_row( size_t const i ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_row( get(), i ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_column().
     * @param j A column index
     * @return A column as a vector
     */
    vector const const_column( size_t const j ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_column( get(), j ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_diagonal().
     * @return The principal diagonal as a vector
     */
    vector const const_diagonal() const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_diagonal( get() ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_subdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector
     */
    vector const const_subdiagonal( size_t const k ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_subdiagonal( get(), k ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_superdiagonal().
     * @param k An index
     * @return Subdiagonal @c k as a vector
     */
    vector const const_superdiagonal( size_t const k ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_superdiagonal( get(), k ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_subrow().
     * @param i A row index
     * @param offset A column offset
     * @param n The number of elements
     * @return A subrow as a vector
     */
    vector const const_subrow( size_t const i, size_t const offset, size_t const n ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_subrow( get(), i, offset, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_matrix_const_subcolumn().
     * @param j A column index
     * @param offset A row offset
     * @param n The number of elements
     * @return A subcolumn as a vector
     */
    vector const const_subcolumn( size_t const j, size_t const offset, size_t const n ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_const_subcolumn( get(), j, offset, n ).vector;
      return vector( w );
    }
   /**
     * C++ version of gsl_matrix_const_view_array().
     * @param base An array of type double
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix
     */
    static matrix const const_view_array( double* base, size_t const n1, size_t const n2 ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_const_view_array( base, n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_const_view_array_with_tda().
     * @param base An array of type double
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix
     */
    static matrix const
    const_view_array_with_tda( double* base, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m =  gsl_matrix_const_view_array_with_tda( base, n1, n2, tda ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_const_view_vector().
     * @param v A vector
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @return A matrix
     */
    static matrix const const_view_vector( vector& v, size_t const n1, size_t const n2 ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m =  gsl_matrix_const_view_vector( v.get(), n1, n2 ).matrix;
      return matrix( m );
    }
    /**
     * C++ version of gsl_matrix_const_view_vector_with_tda().
     * @param v A vector
     * @param n1 The number of rows
     * @param n2 The number of columns
     * @param tda The number of columns in memory
     * @return A matrix
     */
    static matrix const
    const_view_vector_with_tda( vector& v, size_t const n1, size_t const n2, size_t const tda ){
      gsl_matrix* m = static_cast<gsl_matrix*>( malloc( sizeof( gsl_matrix ) ) );
      *m = gsl_matrix_const_view_vector_with_tda( v.get(), n1, n2, tda ).matrix;
      return matrix( m );
    }
  private:
    /**
     * The shared pointer
     */
    gsl_matrix* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_matrix.
     * @return the gsl_matrix
     */
    gsl_matrix* get(){
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Get the gsl_matrix.
     * @return the gsl_matrix
     */
    gsl_matrix const* get() const {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Find if this is the only object sharing the gsl_matrix.
     * @return @c true or @c falses according as 
     * this is the only matrix object sharing the gsl_matrix
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many matrix objects share this pointer.
     * @return the number of matrix objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as this contains a pointer
     * to a gsl_matrix
     */
    operator bool() const { return ccgsl_pointer != 0; }
    // GSL functions
    /**
     * C++ version of gsl_matrix_calloc(). This constructs a matrix object with entries
     * initialised to zero.
     * @param n1 The number of rows in the matrix
     * @param n2 The number of columns in the matrix
     * @return A matrix initialised to zero
     */
    static matrix calloc( size_t const n1, size_t const n2 ){ return matrix( gsl_matrix_calloc( n1, n2 ) ); }
    /**
     * C++ version of gsl_matrix_set_zero().
     */
    void set_zero(){ gsl_matrix_set_zero( get() ); }
    /**
     * C++ version of gsl_matrix_set_all().
     * @param x The value to which all elements are set
     */
    void set_all( double x ){ gsl_matrix_set_all( get(), x ); }
    /**
     * C++ version of gsl_matrix_memcpy().
     * @param src source matrix
     * @return error code on failure
     */
    int memcpy( matrix const& src ){ return gsl_matrix_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_max().
     * @return maximum element of matrix
     */
    double max_val() const { return gsl_matrix_max( get() ); }
    /**
     * C++ version of gsl_matrix_min().
     * @return minimum element of matrix
     */
    double min_val() const { return gsl_matrix_min( get() ); }
    /**
     * C++ version of gsl_matrix_minmax().
     * @param min_out minimum element of matrix
     * @param max_out maximum element of matrix
     */
    void minmax( double* min_out, double* max_out ) const {
      gsl_matrix_minmax( get(), min_out, max_out ); }
    /**
     * C++ version of gsl_matrix_add().
     * @param b matrix to add to this
     * @return error code on failure
     */
    int add( matrix const& b ){ return gsl_matrix_add( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_sub().
     * @param b matrix to subtract from this
     * @return error code on failure
     */
    int sub( matrix const& b ){ return gsl_matrix_sub( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_scale().
     * @param x constant to multiply this by
     * @return error code on failure
     */
    int scale( double const x ){ return gsl_matrix_scale( get(), x ); }
    /**
     * C++ version of gsl_matrix_add_constant().
     * @param x constant to add to each element of this
     * @return error code on failure
     */
    int add_constant( double const x ){ return gsl_matrix_add_constant( get(), x ); }
    /**
     * C++ version of gsl_matrix_isnull().
     * @return @c +1 or @c 0 according as elements are all zero or not
     */
    int isnull() const { return gsl_matrix_isnull( get() ); }
    /**
     * C++ version of gsl_matrix_ispos().
     * @return @c +1 or @c 0 according as elements are all positive or not
     */
    int ispos() const { return gsl_matrix_ispos( get() ); }
    /**
     * C++ version of gsl_matrix_isneg().
     * @return @c +1 or @c 0 according as elements are all negative or not
     */
    int isneg() const { return gsl_matrix_isneg( get() ); }
    /**
     * C++ version of gsl_matrix_isnonneg().
     * @return @c +1 or @c 0 according as elements are all nonnegative or not
     */
    int isnonneg() const { return gsl_matrix_isnonneg( get() ); }
    /**
     * C++ version of gsl_matrix_get().
     * @param i index of row
     * @param j index of column
     * @return value of element
     */
    double get( size_t const i, size_t const j ) const { return gsl_matrix_get( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_set().
     * @param i index of row
     * @param j index of column
     * @param x new value for element
     */
    void set( size_t const i, size_t const j, double x ){ gsl_matrix_set( get(), i, j, x ); }
    /**
     * C++ version of gsl_matrix_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    double* ptr( size_t const i, size_t const j ){ return gsl_matrix_ptr( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_const_ptr().
     * @param i index of row
     * @param j index of column
     * @return pointer to element
     */
    double const* const_ptr( size_t const i, size_t const j ) const {
      return gsl_matrix_const_ptr( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_fread().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fread( FILE* stream ){ return gsl_matrix_fread( stream, get() ); }
    /**
     * C++ version of gsl_matrix_fwrite().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fwrite( FILE* stream ) const { return gsl_matrix_fwrite( stream, get() ); }
    /**
     * C++ version of gsl_matrix_fscanf().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fscanf( FILE* stream ){ return gsl_matrix_fscanf( stream, get() ); }
    /**
     * C++ version of gsl_matrix_fprintf().
     * @param stream A C file stream
     * @param format %d, %e, %f or %g
     * @return error code on failure
     */
    int fprintf( FILE* stream, char const* format ) const {
      return gsl_matrix_fprintf( stream, get(), format ); }
    /**
     * C++ version of gsl_matrix_alloc_from_block().
     * @param b The block
     * @param offset The offset within the block
     * @param n1 The number of rows in the matrix
     * @param n2 The number of columns in the matrix
     * @param d2 undocumented
     */
    matrix( block& b, size_t const offset, size_t const n1, size_t const n2, size_t const d2 ){
      ccgsl_pointer = gsl_matrix_alloc_from_block( b.get(), offset, n1, n2, d2 );
      // just plausibly we could allocate vector but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * C++ version of gsl_matrix_alloc_from_matrix().
     * @param m The matrix
     * @param k1 the row of @c m to take as row zero
     * @param k2 the column of @c m to take as column zero
     * @param n1 The number of rows in the matrix
     * @param n2 The number of columns in the matrix
     */
    matrix( matrix& m, size_t const k1, size_t const k2, size_t const n1, size_t const n2 ){
      ccgsl_pointer = gsl_matrix_alloc_from_matrix( m.get(), k1, k2, n1, n2 );
      // just plausibly we could allocate matrix but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
    	// try to tidy up before rethrowing
    	gsl_matrix_free( ccgsl_pointer );
    	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // More functions
    /**
     * C++ version of gsl_matrix_set_identity().
     */
    void set_identity(){ gsl_matrix_set_identity( get() ); }
    /**
     * C++ version of gsl_matrix_swap_rows().
     * @param i Index of first row
     * @param j Index of second row
     * @return error code on failure
     */
    int swap_rows( size_t const i, size_t const j ){ return gsl_matrix_swap_rows( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_swap_columns().
     * @param i Index of first column
     * @param j Index of second column
     * @return error code on failure
     */
    int swap_columns( size_t const i, size_t const j ){
      return gsl_matrix_swap_columns( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_swap_rowcol(). Swap row and column in place. Matrix must be square.
     * @param i index of row
     * @param j index of column
     * @return error code on failure
     */
    int swap_rowcol( size_t const i, size_t const j ){ return gsl_matrix_swap_rowcol( get(), i, j ); }
    /**
     * C++ version of gsl_matrix_transpose().
     * @return error code on failure.
     */
    int transpose(){ return gsl_matrix_transpose( get() ); }
    /**
     * C++ version of gsl_matrix_transpose_memcpy().
     * @param src matrix whose transpose it to be copied to @c this
     * @return error code on failure
     */
    int transpose_memcpy( matrix const& src ){
      return gsl_matrix_transpose_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_matrix_max_index().
     * @param imax row index of the first maximum element in the matrix
     * @param jmax column index of the first maximum element in the matrix
     */
    void max_index( size_t* imax, size_t* jmax ) const {
      gsl_matrix_max_index( get(), imax, jmax ); }
    /**
     * C++ version of gsl_matrix_min_index().
     * @param imin row index of the first minimum element in the matrix
     * @param jmin column index of the first minimum element in the matrix
     */
    void min_index( size_t* imin, size_t* jmin ) const {
      gsl_matrix_min_index( get(), imin, jmin ); }
    /**
     * C++ version of gsl_matrix_minmax_index().
     * @param imin row index of the first minimum element in the matrix
     * @param jmin column index of the first minimum element in the matrix
     * @param imax row index of the first maximum element in the matrix
     * @param jmax column index of the first maximum element in the matrix
     */
    void minmax_index( size_t* imin, size_t* jmin, size_t* imax, size_t* jmax ) const {
      gsl_matrix_minmax_index( get(), imin, jmin, imax, jmax ); }
    /**
     * C++ version of gsl_matrix_mul_elements(). Multiply matrices elementwise.
     * @param b Another matrix
     * @return error code on failure
     */
    int
    mul_elements( matrix const& b ){
          return gsl_matrix_mul_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_div_elements(). Divide each element of @c this by the
     * corrsponding element of @c b
     * @param b Another matrix
     * @return error code on failure
     */
    int div_elements( matrix const& b ){
      return gsl_matrix_div_elements( get(), b.get() ); }
    /**
     * C++ version of gsl_matrix_add_diagonal().
     * @param x A constant
     * @return error code on failure
     */
    int add_diagonal( double const x ){
      return gsl_matrix_add_diagonal( get(), x ); }
    /**
     * C++ version of gsl_matrix_get_row().
     * @param v A vector
     * @param i The index of the row
     * @return error code on failure
     */
    int get_row( vector& v, size_t const i ) const {
      return gsl_matrix_get_row( v.get(), get(), i ); }
    /**
     * C++ version of gsl_matrix_get_col().
     * @param v A vector
     * @param j The index of the column
     * @return error code on failure
     */
    int get_col( vector& v, size_t const j ) const {
      return gsl_matrix_get_col( v.get(), get(), j ); }
    /**
     * C++ version of gsl_matrix_set_row().
     * @param i The index of the row
     * @param v A vector
     * @return error code on failure
     */
    int set_row( size_t const i, vector const& v ){
      return gsl_matrix_set_row( get(), i, v.get() ); }
    /**
     * C++ version of gsl_matrix_set_col().
     * @param j The index of the column
     * @param v A vector
     * @return error code on failure
     */
    int set_col( size_t const j, vector const& v ){
      return gsl_matrix_set_col( get(), j, v.get() ); }
    // Extra functions for []
    /**
     * This function allows us to use a matrix like an array. Use with caution.
     * Although @c matrix[i][j] is possible, it is much less efficient than
     * matrix::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector representing a row
     */
    vector operator[]( size_t const i ){
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector();
      }
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_row( ccgsl_pointer, i ).vector;
      return vector( w );
    }
    /**
     * This function allows us to use a matrix like an array. Use with caution.
     * Although @c matrix[i][j] is possible, it is much less efficient than
     * matrix::set(). The effect is the same as row()
     * @param i The index of the row
     * @return A vector representing a row
     */
    vector const operator[]( size_t const i ) const {
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "matrix is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return vector();
      }
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_matrix_row( ccgsl_pointer, i ).vector;
      return vector( w );
    }

	//---	tip functions added by fujii	---
	void mul(const matrix &a, const matrix& b){
		size_t size1 = a.size1();
		size_t size2 = a.size2();
		size_t size3 = b.size2();
		if(size2 != b.size1())
  		gsl_error( "size doesn't match between two matrices ", __FILE__, __LINE__, exception::GSL_EFAILED );
		else{
			for(size_t i=0; i<size1; i++){
				for(size_t k=0; k<size3; k++){
					set(i,k,0.0);
					for(size_t j=0; j<size2; j++)
						set(i,k,get(i,k) + a.get(i,j)*b.get(j,k));
					}}
		}return;
	}

	matrix operator * (const matrix& a)const{
		matrix product(size1(), a.size2());
		product.mul(*this, a);
		return product;
	}

	vector operator * (const vector& b)const{
		vector product(size1());
		if(size2() != b.size())
  		gsl_error( "size doesn't match between two matrices ", __FILE__, __LINE__, exception::GSL_EFAILED );
		else{
			product.set_zero();
			for(size_t i=0; i<size1();++i){
				for(size_t j=0; j<size2();++j){
					product.set(i, product.get(i) + get(i, j)*b.get(j));
				}
			}
		}
		return product;
	}
	};

  #ifndef CCGSL_VECTOR_HPP
  // Extra functions for vector allocation from matrix objects
   vector vector::alloc_row_from_matrix( matrix& m, size_t const i ){
    return vector ( gsl_vector_alloc_row_from_matrix( m.get(), i ) ); }
   vector vector::alloc_col_from_matrix( matrix& m, size_t const i ){
    return vector ( gsl_vector_alloc_col_from_matrix( m.get(), i ) ); }
  #endif
   };
#endif
