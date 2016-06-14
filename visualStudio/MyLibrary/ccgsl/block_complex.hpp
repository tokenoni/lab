/*
 * Copyright (C) 2010 John D Lamb
 * Enum copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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

#ifndef CCGSL_BLOCK_COMPLEX_HPP
#define CCGSL_BLOCK_COMPLEX_HPP

#include<gsl/gsl_block.h>
#include<new>
#include<iterator>
#include<cstring>

#include"exception.hpp"

// This file is used as a template
#include"complex.hpp"
#define CCGSL_MTY 2

namespace gsl {
  /**
   * This class handles vector_complexs as shared handles. It models a random access container
   * so that STL functions work with block_complex.
   */
  class block_complex {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    block_complex(){
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
     * The default constructor creates a new block_complex with n elements
     * @param n The number of elements in the block_complex
     */
    explicit block_complex( size_t const n ){
      ccgsl_pointer = gsl_block_complex_calloc( n );
      // just plausibly we could allocate block_complex but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_block_complex_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_block_complex. This is not usually a good idea. In this case
     * we should not use gsl_block_complex_free() to deallocate the memory.
     * @param v The block_complex
     */
    explicit block_complex( gsl_block_complex* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the block_complex. Use clone() if you want a full copy.
     * @param v The block_complex to copy.
     */
    block_complex( block_complex const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // block_complex is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The block_complex to copy
     */
    block_complex& operator=( block_complex const& v ){
      // first, possibly delete anything pointed to by this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_block_complex_free( ccgsl_pointer );
	delete count;
      }
      // Then copy
      ccgsl_pointer = v.ccgsl_pointer;
      count = v.count;
      ++*count; // block_complex is now shared.
      return *this;
    }
    // clone()
    /**
     * The clone function. Use this if you want a copy of the block_complex that does
     * not share the underlying data.
     * @return a new copy of this.
     */
    block_complex clone() const {
      block_complex copy( ccgsl_pointer->size );
      // Now copy
      memcpy( copy.ccgsl_pointer->data, ccgsl_pointer->data,
	      CCGSL_MTY * ccgsl_pointer->size * sizeof( double ) );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~block_complex(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_block_complex_free( ccgsl_pointer );
	delete count;
      }
    }
    // Refines equality comparable
    // == operator
    /**
     * Two block_complex objects are identically equal if their elements are identical.
     * @param v The block_complex to be compared with @c this
     * @return @c true or @c false according as @c this and @c v have
     * identical elements or not
     */
    bool operator==( block_complex const& v ) const {
      // trivially equal if gsl_block_complex*s are identical
      if( ccgsl_pointer == v.ccgsl_pointer ) return true;
      // trivially not equal if one is zero: != should be same as xor here
      if( (ccgsl_pointer == 0) != (v.ccgsl_pointer == 0 ) ) return false;
      // trivially not equal if sizes are different
      if( ccgsl_pointer->size != v.ccgsl_pointer->size ) return false;
      // check elementwise for equality
      for( size_t i = 0; i < CCGSL_MTY * ccgsl_pointer->size; ++i )
  	if( ccgsl_pointer->data[i] != v.ccgsl_pointer->data[i] ) return false; 
      return true;
    }
    // != operator
    /**
     * Two block_complex objects are different equal if their elements are not identical.
     * @param v The block_complex to be compared with @c this
     * @return @c false or @c true according as @c this and @c v have
     * identical elements or not
     */
    bool operator!=( block_complex const& v ) const { return not operator==( v ); }
    // Refines forward container
    // Refines less than comparable
    // operator<
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a block_complex is nonnegative.
     * @param v The block_complex to be compared with @c this
     * @return @c false or @c true according as @c this is less than @c v
     * lexicographically
     */
    bool operator<( block_complex const& v ) const {
      // null block_complex comes first
      if( ccgsl_pointer == 0 ) return v.ccgsl_pointer != 0;
      // if v is null then this > v
      if( v.ccgsl_pointer == 0 ) return false;
      // Now compare elementwise
      size_t const  size = CCGSL_MTY * ccgsl_pointer->size;
      size_t const  v_size = CCGSL_MTY * v.ccgsl_pointer->size;
      size_t const min = size > v_size ? size : v_size;
      for( size_t i = 0; i < min; ++i ){
	double const t = ccgsl_pointer->data[i];
	double const u =v.ccgsl_pointer->data[i];
	if( t < u ) return true;
	if( u < t ) return false;
      }
      // elements match.
      return size < v_size;
    }
    // operator>
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a block_complex is nonnegative.
     * @param v The block_complex to be compared with @c this
     * @return @c false or @c true according as @c this is greater than @c v
     * lexicographically
     */
    bool operator>( block_complex const& v ) const {
      // null block_complex comes first
      if( ccgsl_pointer == 0 ) return false;
      // if v is null then this > v
      if( v.ccgsl_pointer == 0 ) return true;
      // Now compare elementwise
      size_t const  size = CCGSL_MTY * ccgsl_pointer->size;
      size_t const  v_size = CCGSL_MTY * v.ccgsl_pointer->size;
      size_t const min = size > v_size ? size : v_size;
      for( size_t i = 0; i < min; ++i ){
	double const t = ccgsl_pointer->data[i];
	double const u =v.ccgsl_pointer->data[i];
	if( t > u ) return true;
	if( u > t ) return false;
      }
      // elements match.
      return size > v_size;
    }
    // operator<=
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a block_complex is nonnegative.
     * @param v The block_complex to be compared with @c this
     * @return @c false or @c true according as @c this is less than
     * or equal to @c v lexicographically
     */
    bool operator<=( block_complex const& v ) const {
      return operator<( v ) or operator==( v );
    }
    // operator>=
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a block_complex is nonnegative.
     * @param v The block_complex to be compared with @c this
     * @return @c false or @c true according as @c this is no 
     * less than @c v lexicographically
     */ 
    bool operator>=( block_complex const& v ) const {
      return operator>( v ) or operator==( v );
    }
    // Refines container
    // type value_type
    /**
     * A container must have a value_type.
     */
    typedef complex value_type;
    // type reference
    /**
     * A container must have a reference type.
     */
    typedef  complex_ref reference;
    // type const_reference
    /**
     * A container must have a constant reference type.
     */
    typedef  reference const const_reference;
    // type pointer
    /**
     * A container must have a pointer type.
     */
    typedef  complex_ptr pointer;
    // type const_pointer
    /**
     * A container must have a constant pointer type.
     */
    typedef  complex_ptr const const_pointer;
    // type iterator
  private:
    /**
     * The container must have iterator types. We create a suitable class
     * here.
     */
    template<typename container, typename content,bool reverse>
    class iterator_base {
      friend class block_complex;
    public:
      /**
       * An iterator must have a value_type.
       */
      typedef complex value_type;
      /**
       * An iterator must have a reference type.
       */
      typedef complex_ref reference;
      /**
       * An iterator must have a pointer type.
       */
      typedef complex_ptr pointer;
      /**
       * An iterator must have a pointer type.
       */
      typedef std::random_access_iterator_tag iterator_category;
      // // type iterator_traits<block_complex>::difference_type
      /**
       * An iterator must have a difference_type.
       */
      typedef ptrdiff_t difference_type;
    public:
      // // operator*
      /**
       * Dereference the pointer.
       * @return a reference to content
       */
      reference operator*() const {
  	// First check that iterator is initialised.
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
	} else if( v->ccgsl_pointer == 0 ){
	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	// Check that position make sense
  	if( position >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference beyond rbegin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	if( position <= -1 ){
  	  gsl_error( "trying to dereference beyond begin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	// position and v are valid: return data
	return complex_ref( v->ccgsl_pointer->data + CCGSL_MTY * position );
      }
      // // operator->
      /**
       * Dereference the pointer.
       * @return a pointer to content
       */
      pointer operator->() const {
  	// First check that iterator is initialised.
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return pointer( 0 );
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return pointer( 0 );
  	}
  	// Check that position make sense
  	if( position >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference end()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return pointer( 0 );
  	}
  	if( position <= -1 ){
  	  gsl_error( "trying to dereference rend()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return pointer( 0 );
  	}
  	// position and v are valid: return data
  	return pointer( v->ccgsl_pointer->data + CCGSL_MTY * position );
      }
      // // operator[]
      /**
       * Get element at i + n by reference ([] operator).
       * @param n The offset from i
       * @return a reference to content
       */
      reference operator[]( difference_type const n ) const {
  	// First check that iterator is initialised.
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	// Check that position make sense
  	difference_type const p = reverse ? position - n : position + n;
  	if( p >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference beyond rbegin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	if( p <= -1 ){
  	  gsl_error( "trying to dereference beyond begin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return complex_ref( 0 );
  	}
  	// p is a valid position
  	return complex_ref( v->ccgsl_pointer->data + CCGSL_MTY * p );
      }
      // // operator-: find distance between two iterators
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_base<container,content,reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse ? i.position - position : position - i.position;
      }
      // // operator!=
      // // operator<
      /**
       * The == operator.
       * @param i The iterator being compared
       * @return @c true or @c false according as i != *this
       */
      bool operator==( iterator_base<container,content,reverse> const& i ) const {
	return this->v == i.v and this->position == i.position;
      }
      /**
       * The != operator.
       * @param i The iterator being compared
       * @return @c true or @c false according as i != *this
       */
      bool operator!=( iterator_base<container,content,reverse> const& i ) const {
	return not this->operator==( i );
      }
      /**
       * The < operator is used to compare iterators. This only makes sense
       * if the iterators iterate over the same block_complex and the function calls
       * a GSL error handler and returns @c false if they do not.
       * @param i The iterator being compared
       * @return @c true or @c false according as i < j
       */
      bool operator<( iterator_base<container,content,reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse ? i.position < position : position < i.position;
      }
    protected:
      /**
       * Increment the iterator.
       * @return 0 for success, anything else for failure
       */
      void increment(){
  	// Only makes sense if v points to a block_complex
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	// increment position and check against size
  	if( reverse ){ if( position >= 0 ) --position; }
  	else { if( position < static_cast<difference_type>( v->size() ) ) ++position; }
      }
      /**
       * Derement the iterator.
       * @return 0 for success, anything else for failure
       */
      void decrement(){
  	// Only makes sense if v points to a block_complex
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	// decrement position and check against size
  	if( reverse ){ if( position < static_cast<difference_type>( v->size() ) ) ++position; }
  	else { if( position >= 0 ) --position; }
      }
      /**
       * Shift iterator n places
       * @param n A difference_type value to be added to position of iterator
       */
      void shift( difference_type const n ){
  	// Warn if iterator undefined
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	position += reverse ? -n : n;
      }
      /**
       * The iterator is default constructible.
       */
      iterator_base<container,content,reverse>(){ v = 0; }
      /**
       * This constructor allows block_complex to create non-default iterators.
       * @param v The block_complex that creates @c this
       * @param position The initial postion of the iterator
       */
      iterator_base<container,content,reverse>( container* v, difference_type position )
  	: v( v ), position( position ){}
      /**
       * Store a pointer to a block_complex we can iterate over: 0 if no block_complex.
       */
      container* v;
      /**
       * Mark position of iterator within block_complex.
       */
      difference_type position;
    };
    // Need to know about const_iterator_t
    template<bool reverse> class const_iterator_t;
    /**
     * A class template for the two non-const iterators.
     */
    template<bool reverse> class iterator_t : public iterator_base<block_complex,double,reverse>{
    public:
      // // Refines output iterator
      // // operator=
      /**
       * We can assign one output iterator from another.
       * @param i The iterator to copy
       */
      iterator_t<reverse>& operator=( iterator_t<reverse> const& i ){
  	iterator_base<block_complex,double,reverse>::v = i.v;
  	iterator_base<block_complex,double,reverse>::position = i.position;
  	return *this;
      }
      // // Refines forward iterator
      // // operator++ (both)
      /**
       * The prefix ++ operator.
       * @return @c this
       */
      iterator_t<reverse>& operator++(){
  	iterator_base<block_complex,double,reverse>::increment();
  	return *this;
      }
      /**
       * The postfix ++ operator.
       * @return @c this
       */
      iterator_t<reverse> operator++( int ){
  	// return value
  	iterator_t<reverse> result( *this );
  	iterator_base<block_complex,complex,reverse>::increment();
  	return result;
      }
      // // Refines bidirectional iterator
      // // operator-- (both)
      /**
       * The prefix -- operator.
       * @return @c this
       */
      iterator_t<reverse>& operator--(){
  	iterator_base<block_complex,double,reverse>::decrement();
  	return *this;
      }
      /**
       * The postfix -- operator.
       * @return @c this;
       */
      iterator_t<reverse> operator--( int ){
  	// return value
  	iterator_t<reverse> result( *this );
  	iterator_base<block_complex,double,reverse>::decrement();
  	return result;
      }
      /**
       * Difference type.
       */
      typedef typename iterator_base<block_complex,double,reverse>::difference_type
  	difference_type;
      // // operator+=
      /**
       * The += operator.
       * @param n A difference_type value to be added to position of iterator
       * @return @c this
       */
      iterator_t<reverse>& operator+=( difference_type const n ){
  	shift( n );
  	return *this;
      }
      // // operator-=
      /**
       * The -= operator.
       * @param n A difference_type value to be subtracted from position of iterator
       * @return @c this
       */
      iterator_t<reverse>& operator-=( difference_type const n ){
  	shift( -n );
  	return *this;
      }
      // // operator+ (n+i)(i+n)
      /**
       * The + operator.
       * @param n A difference_type value to be added
       * @return A new iterator
       * @see ccgsl::operator+();
       */
      iterator_t<reverse> operator+( difference_type const n ) const {
  	iterator_t<reverse> result( *this );
  	result.shift( n );
  	return result;
      }
      // // operator- (n-i)(i-n)(i-j)
      /**
       * The - operator: subtract distance from iterator
       * @param n A difference_type value to be subtracted
       * @return A new iterator
       * @see ccgsl::operator-();
       */
      iterator_t<reverse> operator-( difference_type const n ) const {
  	iterator_t<reverse> result( *this );
  	result.shift( -n );
  	return result;
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_t<reverse> const& i ) const {
  	return iterator_base<block_complex,double,reverse>::operator-( i );
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A const iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( const_iterator_t<reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( this->v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse ? i.position - this->position : this->position - i.position;
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c true or @c false according as this points to same element as @c i
       */
      bool operator==( const_iterator_t<reverse> const& i ) const {
       	return this->v == i.v and this->position == i.position;
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c false or @c true according as this points to same element as @c i
       */
      bool operator!=( const_iterator_t<reverse> const& i ) const {
       	return not this->operator==( i );
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c true or @c false according as this points to earlier element than @c i
       */
      bool operator<( const_iterator_t<reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse ? i.position < this->position : this->position < i.position;
      }
      /**
       * The default constructor.
       */
      iterator_t<reverse>() : iterator_base<block_complex,complex,reverse>(){}
    protected:
      friend class block_complex;
      // We need a constructor for block_complex
      /**
       * This constructor allows block_complex to create non-default iterators.
       * @param v The block_complex that creates @c this
       * @param position The initial postion of the iterator
       */
      iterator_t<reverse>( block_complex* v, difference_type position )
      : iterator_base<block_complex,double,reverse>( v, position ){}
    };
    /**
     * A class template for the const iterators.
     */
    template<bool reverse> class const_iterator_t
      : public iterator_base<block_complex const,double,reverse>{
    public:
      // // Refines output iterator
      // // operator=
      /**
       * We can assign one output iterator from another.
       * @param i The iterator to copy
       */
      const_iterator_t<reverse>& operator=( const_iterator_t<reverse> const& i ){
  	iterator_base<block_complex const,double,reverse>::v = i.v;
  	iterator_base<block_complex const,double,reverse>::position = i.position;
  	return *this;
      }
      // // Refines forward iterator
      // // operator++ (both)
      /**
       * The prefix ++ operator.
       * @return @c this
       */
      const_iterator_t<reverse>& operator++(){
  	iterator_base<block_complex const,double,reverse>::increment();
  	return *this;
      }
      /**
       * The postfix ++ operator.
       * @return @c this
       */
      const_iterator_t<reverse> operator++( int ){
  	// return value
  	const_iterator_t<reverse> result( *this );
  	iterator_base<block_complex const,double,reverse>::increment();
  	return result;
      }
      // // Refines bidirectional iterator
      // // operator-- (both)
      /**
       * The prefix -- operator.
       * @return @c this
       */
      const_iterator_t<reverse>& operator--(){
  	iterator_base<block_complex const,double,reverse>::decrement();
  	return *this;
      }
      /**
       * The postfix -- operator.
       * @return @c this
       */
      const_iterator_t<reverse> operator--( int ){
  	// return value
  	const_iterator_t<reverse> result( *this );
  	iterator_base<block_complex const,double,reverse>::decrement();
  	return result;
      }
      /**
       * Difference type.
       */
      typedef typename iterator_base<block_complex const,double,reverse>::difference_type
  	difference_type;
      // // operator+=
      /**
       * The += operator.
       * @param n A difference_type value to be added to position of iterator
       * @return @c this
       */
      const_iterator_t<reverse>& operator+=( difference_type const n ){
  	shift( n );
  	return *this;
      }
      // // operator-=
      /**
       * The -= operator.
       * @param n A difference_type value to be subtracted from position of iterator
       * @return @c this
       */
      const_iterator_t<reverse>& operator-=( difference_type const n ){
  	shift( -n );
  	return *this;
      }
      // // operator+ (n+i)(i+n)
      /**
       * The + operator.
       * @param n A difference_type value to be added
       * @return A new iterator
       * @see ccgsl::operator+();
       */
      const_iterator_t<reverse> operator+( difference_type const n ) const {
  	const_iterator_t<reverse> result( *this );
  	result += n;
  	return result;
      }
      // // operator- (n-i)(i-n)(i-j)
      /**
       * The - operator: subtract distance from iterator
       * @param n A difference_type value to be subtracted
       * @return A new iterator
       * @see ccgsl::operator-();
       */
      const_iterator_t<reverse> operator-( difference_type const n ) const {
  	const_iterator_t<reverse> result( *this );
  	result -= n;
  	return result;
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( const_iterator_t<reverse> const& i ) const {
  	return iterator_base<block_complex const,double,reverse>::operator-( i );
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_t<reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( this->v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "block_complex not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse ? i.position - this->position : this->position - i.position;
      }
      /**
       * The default constructor.
       */
      const_iterator_t<reverse>() : iterator_base<block_complex,double,reverse>(){}
      /**
       * A copy constructor
       * @param i The non-const iterator to copy
       */
      const_iterator_t<reverse>( iterator_t<reverse> const&i ){
  	const_iterator_t<reverse>::v = i.v;
  	const_iterator_t<reverse>::position = i.position;
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c true or @c false according as this points to same element as @c i
       */
      bool operator==( iterator_t<reverse> const& i ) const {
       	return this->v == i.v and this->position == i.position;
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c false or @c true according as this points to same element as @c i
       */
      bool operator!=( iterator_t<reverse> const& i ) const {
       	return not this->operator==( i );
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c true or @c false according as this points to earlier element than @c i
       */
      bool operator<( iterator_t<reverse> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same block_complex
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different block_complexs", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse ? i.position < this->position : this->position < i.position;
      }
    protected:
      // We need a constructor for block_complex
      friend class block_complex;
      /**
       * This constructor allows block_complex to create non-default iterators.
       * @param v The block_complex that creates @c this
       * @param position The initial postion of the iterator
       */
      const_iterator_t<reverse>( block_complex const* v, difference_type position )
      : iterator_base<block_complex const,double,reverse>( v, position ){}
    };
  public:
    /**
     * The const_iterator type.
     */
    typedef const_iterator_t<false> const_iterator;
    /**
     * The iterator type.
     */
    typedef iterator_t<false> iterator;
    /**
     * The const_reverse_iterator type.
     */
    typedef const_iterator_t<true> const_reverse_iterator;
    /**
     * The reverse_iterator type.
     */
    typedef iterator_t<true> reverse_iterator;
    // type difference_type
    /**
     * A container must have a difference_type.
     */
    typedef  const_iterator::difference_type difference_type;
    // type size_type
    /**
     * A container must have a size_type.
     */
    typedef  size_t size_type;
    // begin()
    /**
     * Get iterator pointing to first block_complex element.
     * @return iterator pointing to first block_complex element
     */
    iterator begin(){
      return iterator( this, 0 );
    }
    /**
     * Get iterator pointing to first block_complex element.
     * @return iterator pointing to first block_complex element
     */
    const_iterator begin() const {
      return const_iterator( this, 0 );
    }
    // end()
    /**
     * Get iterator pointing beyond last block_complex element.
     * @return iterator pointing beyond last block_complex element
     */
    iterator end(){
      if( ccgsl_pointer == 0 ) return iterator( this, 0 );
      return iterator( this, size() );
    }
    /**
     * Get iterator pointing beyond last block_complex element.
     * @return iterator pointing beyond last block_complex element
     */
    const_iterator end() const {
      if( ccgsl_pointer == 0 ) return const_iterator( this, 0 );
      return const_iterator( this, size() );
    }
    // size()
    /**
     * The size (number of elements) of the block_complex.
     * @return The size of the block_complex
     */
    size_type size() const { return ccgsl_pointer->size; }
    // max_size()
    /**
     * The max size (number of elements) of the block_complex. Identical to
     * size but required for a container.
     * @return The size of the block_complex
     */
    size_type max_size() const { return ccgsl_pointer->size; }
    // empty()
    /**
     * Find if the block_complex is empty.
     * @return @c true if has size zero; otherwise @c false
     */
    bool empty() const { return ccgsl_pointer == 0 or ccgsl_pointer->size == 0; }
    // swap() --- should work even if sizes don't match
    /**
     * Swap two block_complex objects. This works even if the block_complex objects have
     * different sizes because it swaps pointers.
     * @param v The block_complex to swap with @c this.
     */
    void swap( block_complex& v ){
      gsl_block_complex* tmp = ccgsl_pointer; ccgsl_pointer = v.ccgsl_pointer; v.ccgsl_pointer = tmp;
      size_t* tmp2 = count; count = v.count; v.count = tmp2;
    }
    // Refines reversible container
    // rbegin()
    /**
     * Get iterator pointing to first block_complex element.
     * @return iterator pointing to first block_complex element
     */
    reverse_iterator rbegin(){
      if( ccgsl_pointer ==0 ) return reverse_iterator( this, 0 );
      return reverse_iterator( this, size() - 1 );
    }
    /**
     * Get iterator pointing to first block_complex element.
     * @return iterator pointing to first block_complex element
     */
    const_reverse_iterator rbegin() const {
      if( ccgsl_pointer ==0 ) return const_reverse_iterator( this, 0 );
      return const_reverse_iterator( this, size() - 1 );
    }
    // rend()
    /**
     * Get iterator pointing beyon last block_complex element.
     * @return iterator pointing beyond last block_complex element
     */
    reverse_iterator rend(){
      return reverse_iterator( this, -1 );
    }
    /**
     * Get iterator pointing beyon last block_complex element.
     * @return iterator pointing beyond last block_complex element
     */
    const_reverse_iterator rend() const {
      return const_reverse_iterator( this, -1 );
    }
    // operator[]
    /**
     * Get element at position @c n by reference ([] operator).
     * @param n The position of the element
     * @return a reference to a complex
     */
    reference operator[]( size_t const n ){
      // Always try to return something
      static complex something;
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "block_complex is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // Check that position make sense
      if( n >= size() ){
  	gsl_error( "trying to read beyond end of block_complex", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // n is a valid position
      return complex_ref( ccgsl_pointer->data + CCGSL_MTY * n );
    }
    /**
     * Get element at position @c n by reference ([] operator).
     * @param n The position of the element
     * @return a reference to a complex
     */
    const_reference operator[]( size_t const n ) const {
      // Always try to return something
      static complex something;
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "block_complex is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // Check that position make sense
      if( n >= size() ){
  	gsl_error( "trying to read beyond end of block_complex", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // n is a valid position
      return complex_ref( ccgsl_pointer->data + CCGSL_MTY * n );
    }
  private:
    /**
     * The shared pointer
     */
    gsl_block_complex* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_block_complex.
     * @return the gsl_block_complex
     */
    gsl_block_complex* get() const { return ccgsl_pointer; }
    /**
     * Find if this is the only object sharing the gsl_block_complex.
     * @return @c true or @c falses according as 
     * this is the only block_complex object sharing the gsl_block_complex
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many block_complex objects share this pointer.
     * @return the number of block_complex objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as this contains a pointer
     * to a gsl_block_complex
     */
    operator bool() const { return ccgsl_pointer != 0; }
    // gsl functions
  };
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see block_complex::iterator::operator+()
   */
  inline block_complex::iterator operator+
  ( block_complex::iterator::difference_type const n,
  			      block_complex::iterator const& i ){
    return i + n;
  }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see block_complex::iterator::operator+()
   */
  inline block_complex::const_iterator operator+
  ( block_complex::const_iterator::difference_type const n,
  				    block_complex::const_iterator const& i ){
    return i + n;
  }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see block_complex::iterator::operator+()
   */
  inline block_complex::reverse_iterator operator+
  ( block_complex::reverse_iterator::difference_type const n,
  			      block_complex::reverse_iterator const& i ){
    return i + n;
  }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see block_complex::iterator::operator+()
   */
  inline block_complex::const_reverse_iterator operator+
  ( block_complex::const_reverse_iterator::difference_type const n,
  			      block_complex::const_reverse_iterator const& i ){
    return i + n;
  }
  
}

#undef CCGSL_MTY
#endif