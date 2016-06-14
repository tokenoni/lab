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

#ifndef CCGSL_PERMUTATION_HPP
#define CCGSL_PERMUTATION_HPP

#include<gsl/gsl_permutation.h>

#include"exception.hpp"


namespace gsl {
  /**
   * This class handles GSL permutation objects.
   */
  class permutation {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    permutation(){
      ccgsl_pointer = 0;
      // just plausibly we could fail to allocate count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * This constructor creates a new permutation with n elements
     * @param n The number of elements in the permutation
     */
    explicit permutation( size_t const n ){
      ccgsl_pointer = gsl_permutation_alloc( n );
      // just plausibly we could allocate permutation but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_permutation_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_permutation. This is not usually a good idea. In this case
     * we should not use gsl_permutation_free() to deallocate the memory.
     * @param v The permutation
     */
    explicit permutation( gsl_permutation* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the permutation. Use clone() if you want a full copy.
     * @param v The permutation to copy.
     */
    permutation( permutation const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // permutation is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The permutation to copy
     */
    permutation& operator=( permutation const& v ){
      // first, possibly delete anything pointed to by this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_permutation_free( ccgsl_pointer );
	delete count;
      }
      // Then copy
      ccgsl_pointer = v.ccgsl_pointer;
      count = v.count;
      ++*count; // block is now shared.
      return *this;
    }

//------------	modified by fujii	-------------------
	size_t operator[](size_t i){return gsl_permutation_get(ccgsl_pointer, i);}
//-----------------------------------------------------

	// clone()
    /**
     * The clone function. Use this if you want a copy of the block that does
     * not share the underlying data.
     * @return a new copy of this.
     */
    permutation clone() const {
      permutation copy( get()->size );
      // Now copy
      gsl_permutation_memcpy( copy.get(), get() );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~permutation(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_permutation_free( ccgsl_pointer );
	delete count;
      }
    }
  private:
    /**
     * The shared pointer
     */
    gsl_permutation* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_permutation.
     * @return the gsl_permutation
     */
    gsl_permutation* get(){
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Get the gsl_permutation.
     * @return the gsl_permutation
     */
    gsl_permutation const* get() const {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Find if this is the only object sharing the gsl_permutation.
     * @return @c true or @c falses according as 
     * this is the only permutation object sharing the gsl_permutation
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many permutation objects share this pointer.
     * @return the number of permutation objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as this contains a pointer
     * to a gsl_permutation
     */
    operator bool() const { return ccgsl_pointer != 0; }
  public:
    // GSL functions
    /**
     * C++ version of gsl_permutation_calloc().
     * @param n The size of the permutation
     * @return A new permutation, initialised to the identity
     */
    static permutation calloc( size_t const n ){ return permutation( gsl_permutation_calloc( n ) ); }
    /**
     * C++ version of gsl_permutation_init().
     */
    void init(){ gsl_permutation_init( get() ); }
    /**
     * C++ version of gsl_permutation_memcpy().
     * @param src The permutation to copy
     * @return Error code on failure
     */
    int memcpy( permutation const& src ){ return gsl_permutation_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_permutation_fread().
     * @param stream An output stream
     * @return Error code on failure
     */
    int fread( FILE* stream ){ return gsl_permutation_fread( stream, get() ); }
    /**
     * C++ version of gsl_permutation_fwrite().
     * @param stream An output stream
     * @return Error code on failure
     */
    int fwrite( FILE* stream ) const { return gsl_permutation_fwrite( stream, get() ); }
    /**
     * C++ version of gsl_permutation_fscanf().
     * @param stream An output stream
     * @return Error code on failure
     */
    int fscanf( FILE* stream ){ return gsl_permutation_fscanf( stream, get() ); }
    /**
     * C++ version of gsl_permutation_fprintf().
     * @param stream An output stream
     * @param format A print format: typically "%zu\n"
     * @return Error code on failure
     */
    int fprintf( FILE* stream, char const* format ) const {
      return gsl_permutation_fprintf( stream, get(), format ); }
    /**
     * C++ version of gsl_permutation_size().
     * @return The size
     */
    size_t size() const { return gsl_permutation_size( get() ); }
    /**
     * C++ version of gsl_permutation_data().
     * @return Pointer to the array of elements
     */
    size_t* data() const { return gsl_permutation_data( get() ); }
    /**
     * C++ version of gsl_permutation_swap().
     * @param i An element
     * @param j An element to swap with @c i
     * @return Error code on failure
     */
    int swap( size_t const i, size_t const j ){ return gsl_permutation_swap( get(), i, j ); }
    /**
     * C++ version of gsl_permutation_valid().
     * @return Error code for valid
     */
    int valid() const { return gsl_permutation_valid( get() ); }
    /**
     * C++ version of gsl_permutation_reverse().
     */
    void reverse(){ gsl_permutation_reverse( get() ); }
    /**
     * C++ version of gsl_permutation_inverse().
     * @param p Another permutation
     * @return Error code on failure
     */
    int inverse( permutation const& p ){ return gsl_permutation_inverse( get(), p.get() ); }
    /**
     * C++ version of gsl_permutation_next().
     * @return Error code on failure
     */
    int next(){ return gsl_permutation_next( get() ); }
    /**
     * C++ version of gsl_permutation_prev().
     * @return Error code on failure
     */
    int prev(){ return gsl_permutation_prev( get() ); }
    /**
     * C++ version of gsl_permutation_mul().
     * @param pa The first permutation
     * @param pb The second permutation (applied after @c pb)
     * @return Error code on failure
     */
    int mul( permutation const& pa, permutation const& pb ){
      return gsl_permutation_mul( get(), pa.get(), pb.get() ); }

    /**
     * C++ version of gsl_permutation_linear_to_canonical().
     * @param p A permutation
     * @return Error code on failure
     */
    int linear_to_canonical( permutation const& p ){
      return gsl_permutation_linear_to_canonical( get(), p.get() ); }
    /**
     * C++ version of gsl_permutation_canonical_to_linear().
     * @param q Apermutation
     * @return Error code on failure
     */
    int canonical_to_linear( permutation const& q ){
      return gsl_permutation_canonical_to_linear( get(), q.get() ); }
    /**
     * C++ version of gsl_permutation_inversions().
     * @return The number of inversions (pairs not in order)
     */
    size_t inversions() const { return gsl_permutation_inversions( get() ); }
    /**
     * C++ version of gsl_permutation_linear_cycles().
     * @return The number of cycles (linear form)
     */
    size_t linear_cycles() const { return gsl_permutation_linear_cycles( get() ); }
    /**
     * C++ version of gsl_permutation_canonical_cycles().
     * @return The number of cycles (canonical form)
     */
    size_t canonical_cycles(  ) const { return gsl_permutation_canonical_cycles( get() ); }
    /**
     * C++ version of gsl_permutation_get().
     * @param i index of element
     * @return Element at index @c i
     */
    size_t get( size_t const i ) const { return gsl_permutation_get( get(), i ); }
  };

}
#endif
