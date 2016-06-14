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

#ifndef CCGSL_VECTOR_HPP
#define CCGSL_VECTOR_HPP

#include<gsl/gsl_vector.h>
#include<new>
#include<iterator>

#include"exception.hpp"
#include"block.hpp"

// This file is used as a template

namespace gsl {
  // declare matrix class
  class matrix;
  /**
   * This class handles vector objects as shared handles. It models a random access container
   * so that STL functions work with vector.
   *
   * Note that vector_views are implemented as vector objects here.
   */
  class vector {
  public:
    /**
     * The default constructor is only really useful for assigning to.
     */
    vector(){
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
     * The default constructor creates a new vector with n elements
     * @param n The number of elements in the vector
     */
    explicit vector( size_t const n ){
      ccgsl_pointer = gsl_vector_alloc( n );
      // just plausibly we could allocate vector but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_vector_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * Could construct from a gsl_vector. This is not usually a good idea. In this case
     * we should not use gsl_vector_free() to deallocate the memory.
     * @param v The vector
     */
    explicit vector( gsl_vector* v ){
      ccgsl_pointer = v;
      // just plausibly we could fail to allocate count: no further action needed.
      count = new size_t;
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    // copy constructor
    /**
     * The copy constructor. This shares the vector. Use clone() if you want a full copy.
     * @param v The vector to copy.
     */
    vector( vector const& v ) : ccgsl_pointer( v.ccgsl_pointer ), count( v.count ){
      ++*count; // vector is now shared.
    }
    // assignment operator
    /**
     * The assignment operator. This makes a shared copy.
     * @param v The vector to copy
     */
    vector& operator=( vector const& v ){
      // first, possibly delete anything pointed to by @c this
      if( --*count == 0 ){
	if( ccgsl_pointer != 0 ) gsl_vector_free( ccgsl_pointer );
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
     * @return a new copy of @c this.
     */
    vector clone() const {
      vector copy( get()->size );
      // Now copy
      gsl_vector_memcpy( copy.get(), get() );
      // return new object
      return copy;
    }
    // destructor
    /**
     * The destructor only deletes the pointers if count reaches zero.
     */
    ~vector(){
      if( --*count == 0 ){
	// could have allocated null pointer
	if( ccgsl_pointer != 0 ) gsl_vector_free( ccgsl_pointer );
	delete count;
      }
    }
    // Refines equality comparable
    // == operator
    /**
     * Two vector objects are identically equal if their elements are identical.
     * @param v The vector to be compared with @c this
     * @return @c true or @c false according as @c this and @c v have
     * identical elements or not
     */
    bool operator==( vector const& v ) const {
      // trivially equal if gsl_vector*s are identical
      if( ccgsl_pointer == v.ccgsl_pointer ) return true;
      // trivially not equal if one is zero: != should be same as xor here
      if( (ccgsl_pointer == 0) != (v.ccgsl_pointer == 0 ) ) return false;
      // trivially not equal if sizes are different
      if( ccgsl_pointer->size != v.ccgsl_pointer->size ) return false;
      // check elementwise for equality
      for( size_t i = 0; i < ccgsl_pointer->size; ++i )
  	if( gsl_vector_get( ccgsl_pointer, i ) != gsl_vector_get( v.ccgsl_pointer, i ) ) return false; 
      return true;
    }
    // != operator
    /**
     * Two vector objects are different equal if their elements are not identical.
     * @param v The vector to be compared with @c this
     * @return @c false or @c true according as @c this and @c v have
     * identical elements or not
     */
    bool operator!=( vector const& v ) const { return not operator==( v ); }
    // Refines forward container
    // Refines less than comparable
    // operator<
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a vector is nonnegative.
     * @param v The vector to be compared with @c this
     * @return @c false or @c true according as @c this is less than @c v
     * lexicographically
     */
    bool operator<( vector const& v ) const {
      // null vector comes first
      if( ccgsl_pointer == 0 ) return v.ccgsl_pointer != 0;
      // if v is null then this > v
      if( v.ccgsl_pointer == 0 ) return false;
      // Now compare elementwise
      size_t const  size = ccgsl_pointer->size;
      size_t const  v_size = v.ccgsl_pointer->size;
      size_t const min = size > v_size ? size : v_size;
      for( size_t i = 0; i < min; ++i ){
	double const t = gsl_vector_get( ccgsl_pointer, i );
	double const u =gsl_vector_get( v.ccgsl_pointer, i );
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
     * for example, for checking, that a vector is nonnegative.
     * @param v The vector to be compared with @c this
     * @return @c false or @c true according as @c this is greater than @c v
     * lexicographically
     */
    bool operator>( vector const& v ) const {
      // null vector comes first
      if( ccgsl_pointer == 0 ) return false;
      // if v is null then this > v
      if( v.ccgsl_pointer == 0 ) return true;
      // Now compare elementwise
      size_t const  size = ccgsl_pointer->size;
      size_t const  v_size = v.ccgsl_pointer->size;
      size_t const min = size > v_size ? size : v_size;
      for( size_t i = 0; i < min; ++i ){
	double const t = gsl_vector_get( ccgsl_pointer, i );
	double const u =gsl_vector_get( v.ccgsl_pointer, i );
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
     * for example, for checking, that a vector is nonnegative.
     * @param v The vector to be compared with @c this
     * @return @c false or @c true according as @c this is less than
     * or equal to @c v lexicographically
     */
    bool operator<=( vector const& v ) const {
      return operator<( v ) or operator==( v );
    }
    // operator>=
    /**
     * A container needs to define an ordering for sorting. This uses
     * standard lexicographical ordering and so is not useful,
     * for example, for checking, that a vector is nonnegative.
     * @param v The vector to be compared with @c this
     * @return @c false or @c true according as @c this is no 
     * less than @c v lexicographically
     */ 
    bool operator>=( vector const& v ) const {
      return operator>( v ) or operator==( v );
    }

    // Refines container
    // type value_type
    /**
     * A container must have a value_type.
     */
    typedef double value_type;
    // type reference
    /**
     * A container must have a reference type.
     */
    typedef  value_type& reference;
    // type const_reference
    /**
     * A container must have a constant reference type.
     */
    typedef  value_type const& const_reference;
    // type pointer
    /**
     * A container must have a pointer type.
     */
    typedef  value_type* pointer;
    // type const_pointer
    /**
     * A container must have a constant pointer type.
     */
    typedef  value_type const* const_pointer;
    // type iterator
  private:
    /**
     * The container must have iterator types. We create a suitable class
     * here.
     */
    template<typename container, typename content,bool reverse_t>
    class iterator_base
      : public std::iterator<std::random_access_iterator_tag,ptrdiff_t,container*,container&> {
      friend class vector;
    public:
      // // type iterator_traits<vector>::difference_type
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
  	// Always try to return something
  	static content something = 0;
  	// First check that iterator is initialised.
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	// Check that position make sense
  	if( position >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference beyond rbegin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	if( position <= -1 ){
  	  gsl_error( "trying to dereference beyond begin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	// position and v are valid: return data
  	return *(v->ccgsl_pointer->data + position * v->ccgsl_pointer->stride);
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
  	  return 0;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Check that position make sense
  	if( position >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference end()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	if( position <= -1 ){
  	  gsl_error( "trying to dereference rend()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// position and v are valid: return data
  	return v->ccgsl_pointer->data + position * v->ccgsl_pointer->stride;
      }
      // // operator[]
      /**
       * Get element at i + n by reference ([] operator).
       * @param n The offset from i
       * @return a reference to content
       */
      reference operator[]( difference_type const n ) const {
  	// Always try to return something
  	static content something = 0;
  	// First check that iterator is initialised.
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	// Check that position make sense
  	difference_type const p = reverse_t ? position - n : position + n;
  	if( p >= static_cast<difference_type>( v->size() ) ){
  	  gsl_error( "trying to dereference beyond rbegin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	if( p <= -1 ){
  	  gsl_error( "trying to dereference beyond begin()", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return something;
  	}
  	// p is a valid position
  	return *(v->ccgsl_pointer->data + p * v->ccgsl_pointer->stride);
      }
      // // operator-: find distance between two iterators
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_base<container,content,reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same vector
  	if( v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse_t ? i.position - position : position - i.position;
      }
      // // operator!=
      // // operator<
      /**
       * The == operator.
       * @param i The iterator being compared
       * @return @c true or @c false according as i != *this
       */
      bool operator==( iterator_base<container,content,reverse_t> const& i ) const {
	return this->v == i.v and this->position == i.position;
      }
      /**
       * The != operator.
       * @param i The iterator being compared
       * @return @c true or @c false according as i != *this
       */
      bool operator!=( iterator_base<container,content,reverse_t> const& i ) const {
	return not this->operator==( i );
      }
      /**
       * The < operator is used to compare iterators. This only makes sense
       * if the iterators iterate over the same vector and the function calls
       * a GSL error handler and returns @c false if they do not.
       * @param i The iterator being compared
       * @return @c true or @c false according as i < j
       */
      bool operator<( iterator_base<container,content,reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same vector
  	if( v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse_t ? i.position < position : position < i.position;
      }
    protected:
      /**
       * Increment the iterator.
       * @return 0 for success, anything else for failure
       */
      void increment(){
  	// Only makes sense if v points to a vector
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	// increment position and check against size
  	if( reverse_t ){ if( position >= 0 ) --position; }
  	else { if( position < static_cast<difference_type>( v->size() ) ) ++position; }
      }
      /**
       * Derement the iterator.
       * @return 0 for success, anything else for failure
       */
      void decrement(){
  	// Only makes sense if v points to a vector
  	if( v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	} else if( v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	// decrement position and check against size
  	if( reverse_t ){ if( position < static_cast<difference_type>( v->size() ) ) ++position; }
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
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return;
  	}
  	position += reverse_t ? -n : n;
      }
      /**
       * The iterator is default constructible.
       */
      iterator_base<container,content,reverse_t>(){ v = 0; }
      /**
       * This constructor allows vector to create non-default iterators.
       * @param v The vector that creates @c this
       * @param position The initial postion of the iterator
       */
      iterator_base<container,content,reverse_t>( container* v, difference_type position )
  	: v( v ), position( position ){}
      /**
       * Store a pointer to a vector we can iterate over: 0 if no vector.
       */
      container* v;
      /**
       * Mark position of iterator within vector.
       */
      difference_type position;
    };
    // Need to know about const_iterator_t
    template<bool reverse_t> class const_iterator_t;
    /**
     * A class template for the two non-const iterators.
     */
    template<bool reverse_t> class iterator_t : public iterator_base<vector,double,reverse_t>{
    public:
      // // Refines output iterator
      // // operator=
      /**
       * We can assign one output iterator from another.
       * @param i The iterator to copy
       */
      iterator_t<reverse_t>& operator=( iterator_t<reverse_t> const& i ){
  	iterator_base<vector,double,reverse_t>::v = i.v;
  	iterator_base<vector,double,reverse_t>::position = i.position;
  	return *this;
      }
      // // Refines forward iterator
      // // operator++ (both)
      /**
       * The prefix ++ operator.
       * @return @c this
       */
      iterator_t<reverse_t>& operator++(){
  	iterator_base<vector,double,reverse_t>::increment();
  	return *this;
      }
      /**
       * The postfix ++ operator.
       * @return @c this
       */
      iterator_t<reverse_t> operator++( int ){
  	// return value
  	iterator_t<reverse_t> result( *this );
  	iterator_base<vector,double,reverse_t>::increment();
  	return result;
      }
      // // Refines bidirectional iterator
      // // operator-- (both)
      /**
       * The prefix -- operator.
       * @return @c this
       */
      iterator_t<reverse_t>& operator--(){
  	iterator_base<vector,double,reverse_t>::decrement();
  	return *this;
      }
      /**
       * The postfix -- operator.
       * @return @c this;
       */
      iterator_t<reverse_t> operator--( int ){
  	// return value
  	iterator_t<reverse_t> result( *this );
  	iterator_base<vector,double,reverse_t>::decrement();
  	return result;
      }
      /**
       * Difference type.
       */
      typedef typename iterator_base<vector,double,reverse_t>::difference_type
  	difference_type;
      // // operator+=
      /**
       * The += operator.
       * @param n A difference_type value to be added to position of iterator
       * @return @c this
       */
      iterator_t<reverse_t>& operator+=( difference_type const n ){
  	shift( n );
  	return *this;
      }
      // // operator-=
      /**
       * The -= operator.
       * @param n A difference_type value to be subtracted from position of iterator
       * @return @c this
       */
      iterator_t<reverse_t>& operator-=( difference_type const n ){
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
      iterator_t<reverse_t> operator+( difference_type const n ) const {
  	iterator_t<reverse_t> result( *this );
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
      iterator_t<reverse_t> operator-( difference_type const n ) const {
  	iterator_t<reverse_t> result( *this );
  	result.shift( -n );
  	return result;
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_t<reverse_t> const& i ) const {
  	return iterator_base<vector,double,reverse_t>::operator-( i );
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A const iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( const_iterator_t<reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( this->v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same vector
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse_t ? i.position - this->position : this->position - i.position;
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c true or @c false according as @c this points to same element as @c i
       */
      bool operator==( const_iterator_t<reverse_t> const& i ) const {
       	return this->v == i.v and this->position == i.position;
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c false or @c true according as @c this points to same element as @c i
       */
      bool operator!=( const_iterator_t<reverse_t> const& i ) const {
       	return not this->operator==( i );
      }
      /**
       * Comparison with const iterator.
       * @param i another iterator
       * @return @c true or @c false according as @c this points to earlier element than @c i
       */
      bool operator<( const_iterator_t<reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same vector
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse_t ? i.position < this->position : this->position < i.position;
      }
      /**
       * The default constructor.
       */
      iterator_t<reverse_t>() : iterator_base<vector,double,reverse_t>(){}
    protected:
      friend class vector;
      // We need a constructor for vector
      /**
       * This constructor allows vector to create non-default iterators.
       * @param v The vector that creates @c this
       * @param position The initial postion of the iterator
       */
      iterator_t<reverse_t>( vector* v, difference_type position )
      : iterator_base<vector,double,reverse_t>( v, position ){}
    };
    /**
     * A class template for the const iterators.
     */
    template<bool reverse_t> class const_iterator_t
      : public iterator_base<vector const,double,reverse_t>{
    public:
      // // Refines output iterator
      // // operator=
      /**
       * We can assign one output iterator from another.
       * @param i The iterator to copy
       */
      const_iterator_t<reverse_t>& operator=( const_iterator_t<reverse_t> const& i ){
  	iterator_base<vector const,double,reverse_t>::v = i.v;
  	iterator_base<vector const,double,reverse_t>::position = i.position;
  	return *this;
      }
      // // Refines forward iterator
      // // operator++ (both)
      /**
       * The prefix ++ operator.
       * @return @c this
       */
      const_iterator_t<reverse_t>& operator++(){
  	iterator_base<vector const,double,reverse_t>::increment();
  	return *this;
      }
      /**
       * The postfix ++ operator.
       * @return @c this
       */
      const_iterator_t<reverse_t> operator++( int ){
  	// return value
  	const_iterator_t<reverse_t> result( *this );
  	iterator_base<vector const,double,reverse_t>::increment();
  	return result;
      }
      // // Refines bidirectional iterator
      // // operator-- (both)
      /**
       * The prefix -- operator.
       * @return @c this
       */
      const_iterator_t<reverse_t>& operator--(){
  	iterator_base<vector const,double,reverse_t>::decrement();
  	return *this;
      }
      /**
       * The postfix -- operator.
       * @return @c this
       */
      const_iterator_t<reverse_t> operator--( int ){
  	// return value
  	const_iterator_t<reverse_t> result( *this );
  	iterator_base<vector const,double,reverse_t>::decrement();
  	return result;
      }
      /**
       * Difference type.
       */
      typedef typename iterator_base<vector const,double,reverse_t>::difference_type
  	difference_type;
      // // operator+=
      /**
       * The += operator.
       * @param n A difference_type value to be added to position of iterator
       * @return @c this
       */
      const_iterator_t<reverse_t>& operator+=( difference_type const n ){
  	shift( n );
  	return *this;
      }
      // // operator-=
      /**
       * The -= operator.
       * @param n A difference_type value to be subtracted from position of iterator
       * @return @c this
       */
      const_iterator_t<reverse_t>& operator-=( difference_type const n ){
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
      const_iterator_t<reverse_t> operator+( difference_type const n ) const {
  	const_iterator_t<reverse_t> result( *this );
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
      const_iterator_t<reverse_t> operator-( difference_type const n ) const {
  	const_iterator_t<reverse_t> result( *this );
  	result -= n;
  	return result;
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( const_iterator_t<reverse_t> const& i ) const {
  	return iterator_base<vector const,double,reverse_t>::operator-( i );
      }
      /**
       * The - operator: find distance between two iterators
       * @param i A second iterator
       * @return (signed) distance between @c this and @c i
       */
      difference_type operator-( iterator_t<reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	} else if( this->v->ccgsl_pointer == 0 or i.v->ccgsl_pointer == 0 ){
  	  gsl_error( "vector not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return 0;
  	}
  	// Warn if iterators do not point to same vector
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return 0;
  	}
  	return reverse_t ? i.position - this->position : this->position - i.position;
      }
      /**
       * The default constructor.
       */
      const_iterator_t<reverse_t>() : iterator_base<vector,double,reverse_t>(){}
      /**
       * A copy constructor
       * @param i The non-const iterator to copy
       */
      const_iterator_t<reverse_t>( iterator_t<reverse_t> const&i ){
  	const_iterator_t<reverse_t>::v = i.v;
  	const_iterator_t<reverse_t>::position = i.position;
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c true or @c false according as @c this points to same element as @c i
       */
      bool operator==( iterator_t<reverse_t> const& i ) const {
       	return this->v == i.v and this->position == i.position;
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c false or @c true according as @c this points to same element as @c i
       */
      bool operator!=( iterator_t<reverse_t> const& i ) const {
       	return not this->operator==( i );
      }
      /**
       * Comparison with non-const iterator.
       * @param i another iterator
       * @return @c true or @c false according as @c this points to earlier element than @c i
       */
      bool operator<( iterator_t<reverse_t> const& i ) const {
  	// Warn if either iterator undefined
  	if( this->v == 0 or i.v == 0 ){
  	  gsl_error( "iterator not initialised", __FILE__, __LINE__, exception::GSL_EFAILED );
  	  return false;
  	}
  	// Warn if iterators do not point to same vector
  	if( this->v->ccgsl_pointer != i.v->ccgsl_pointer ){
  	  gsl_error( "trying to take difference of iterators for different vector objects", __FILE__, __LINE__,
  		     exception::GSL_EFAILED );
  	  return false;
  	}
  	return reverse_t ? i.position < this->position : this->position < i.position;
      }
    protected:
      // We need a constructor for vector
      friend class vector;
      /**
       * This constructor allows vector to create non-default iterators.
       * @param v The vector that creates @c this
       * @param position The initial postion of the iterator
       */
      const_iterator_t<reverse_t>( vector const* v, difference_type position )
      : iterator_base<vector const,double,reverse_t>( v, position ){}
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
     * The const_reverse_t type.
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
     * Get iterator pointing to first vector element.
     * @return iterator pointing to first vector element
     */
    iterator begin(){
      return iterator( this, 0 );
    }
    /**
     * Get iterator pointing to first vector element.
     * @return iterator pointing to first vector element
     */
    const_iterator begin() const {
      return const_iterator( this, 0 );
    }
    // end()
    /**
     * Get iterator pointing beyond last vector element.
     * @return iterator pointing beyond last vector element
     */
    iterator end(){
      if( ccgsl_pointer == 0 ) return iterator( this, 0 );
      return iterator( this, size() );
    }
    /**
     * Get iterator pointing beyond last vector element.
     * @return iterator pointing beyond last vector element
     */
    const_iterator end() const {
      if( ccgsl_pointer == 0 ) return const_iterator( this, 0 );
      return const_iterator( this, size() );
    }
    // size()
    /**
     * The size (number of elements) of the vector.
     * @return The size of the vector
     */
    size_type size() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size; }
    // max_size()
    /**
     * The max size (number of elements) of the vector. Identical to
     * size but required for a container.
     * @return The size of the vector
     */
    size_type max_size() const { return ccgsl_pointer == 0 ? 0 : ccgsl_pointer->size; }
    // empty()
    /**
     * Find if the vector is empty.
     * @return @c true if has size zero; otherwise @c false
     */
    bool empty() const { return ccgsl_pointer == 0 or ccgsl_pointer->size == 0; }
    // swap() --- should work even if sizes don't match
    /**
     * Swap two vector objects. This works even if the vector objects have different sizes
     * because it swaps pointers.
     * @param v The vector to swap with @c this.
     */
    void swap( vector& v ){
      gsl_vector* tmp = ccgsl_pointer; ccgsl_pointer = v.ccgsl_pointer; v.ccgsl_pointer = tmp;
      size_t* tmp2 = count; count = v.count; v.count = tmp2;
    }
    // Refines reversible container
    // rbegin()
    /**
     * Get iterator pointing to first vector element.
     * @return iterator pointing to first vector element
     */
    reverse_iterator rbegin(){
      if( ccgsl_pointer ==0 ) return reverse_iterator( this, 0 );
      return reverse_iterator( this, size() - 1 );
    }
    /**
     * Get iterator pointing to first vector element.
     * @return iterator pointing to first vector element
     */
    const_reverse_iterator rbegin() const {
      if( ccgsl_pointer ==0 ) return const_reverse_iterator( this, 0 );
      return const_reverse_iterator( this, size() - 1 );
    }
    // rend()
    /**
     * Get iterator pointing beyon last vector element.
     * @return iterator pointing beyond last vector element
     */
    reverse_iterator rend(){
      return reverse_iterator( this, -1 );
    }
    /**
     * Get iterator pointing beyon last vector element.
     * @return iterator pointing beyond last vector element
     */
    const_reverse_iterator rend() const {
      return const_reverse_iterator( this, -1 );
    }
    // operator[]
    /**
     * Get element at position @c n by reference ([] operator).
     * @param n The position of the element
     * @return a reference to a double
     */
    double& operator[]( size_t const n ){
      // Always try to return something
      static double something = 0;
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "vector is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // Check that position make sense
      if( n >= size() ){
  	gsl_error( "trying to read beyond end of vector", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // n is a valid position
      return *(ccgsl_pointer->data + n * ccgsl_pointer->stride);
    }
    /**
     * Get element at position @c n by reference ([] operator).
     * @param n The position of the element
     * @return a reference to a double
     */
    double const& operator[]( size_t const n ) const {
      // Always try to return something
      static double something = 0;
      // First check that iterator is initialised.
      if( ccgsl_pointer == 0 ){
  	gsl_error( "vector is null", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // Check that position make sense
      if( n >= size() ){
  	gsl_error( "trying to read beyond end of vector", __FILE__, __LINE__, exception::GSL_EFAILED );
  	return something;
      }
      // n is a valid position
      return *(ccgsl_pointer->data + n * ccgsl_pointer->stride);
    }
  private:
    /**
     * The shared pointer
     */
    gsl_vector* ccgsl_pointer;
    /**
     * The shared reference count
     */
    size_t* count;
  public:
    // shared reference functions
    /**
     * Get the gsl_vector.
     * @return the gsl_vector
     */
    gsl_vector* get() {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Get the gsl_vector.
     * @return the gsl_vector
     */
    gsl_vector const* get() const {
#ifndef GSL_RANGE_CHECK_OFF
      if( ccgsl_pointer == 0 )
  	gsl_error( "Object not initialised", __FILE__, __LINE__, exception::GSL_EFAULT );
#endif
      return ccgsl_pointer; }
    /**
     * Find if @c this is the only object sharing the gsl_vector.
     * @return @c true or @c falses according as 
     * this is the only vector object sharing the gsl_vector
     */
    bool unique() const { return *count == 1; }
    /**
     * Find how many vector objects share this pointer.
     * @return the number of vector objects that share this pointer
     */
    size_t use_count() const { return *count; }
    /**
     * Allow conversion to bool.
     * @return @c true or @c false according as @c this contains a pointer
     * to a gsl_vector
     */
    operator bool() const { return ccgsl_pointer != 0; }

    // GSL functions
    /**
     * C++ version of gsl_vector_calloc(). This constructs a vector object with entries
     * initialised to zero.
     * @param n The size of the vector
     * @return A vector initialised to zero
     */
    static vector calloc( size_t const n ){ return vector( gsl_vector_calloc( n ) ); }
    /**
     * C++ version of gsl_vector_set_zero().
     */
    void set_zero(){ gsl_vector_set_zero( get() ); }
    /**
     * C++ version of gsl_vector_set_all().
     * @param x The value to which all elements are set
     */
    void set_all( double x ){ gsl_vector_set_all( get(), x ); }
    /**
     * C++ version of gsl_vector_set_basis(). Creates a basis vector with one nonzero element.
     * @param i The element to be set to 1.
     * @return error code on failure
     */
    int set_basis( size_t i ){ return gsl_vector_set_basis( get(), i ); }
    /**
     * C++ version of gsl_vector_memcpy().
     * @param src source vector
     * @return error code on failure
     */
    int memcpy( vector const& src ){ return gsl_vector_memcpy( get(), src.get() ); }
    /**
     * C++ version of gsl_vector_reverse().
     * @return error code on failure
     */
    int reverse(){ return gsl_vector_reverse( get() ); }
    /**
     * C++ version of gsl_vector_swap_elements().
     * @param i first element
     * @param j second element
     * @return error code on failure
     */
    int swap_elements( size_t const i, size_t const j ){
      return gsl_vector_swap_elements( get(), i, j ); }
    /**
     * C++ version of gsl_vector_max().
     * @return maximum element of vector
     */
    double max_val() const { return gsl_vector_max( get() ); }
    /**
     * C++ version of gsl_vector_min().
     * @return minimum element of vector
     */
    double min_val() const { return gsl_vector_min( get() ); }
    /**
     * C++ version of gsl_vector_minmax().
     * @param min_out minimum element of vector
     * @param max_out maximum element of vector
     */
    void minmax( double* min_out, double* max_out ) const {
      gsl_vector_minmax( get(), min_out, max_out ); }
    /**
     * C++ version of gsl_vector_max_index().
     * @return index of maximum value of vector
     */
    size_t max_index() const { return gsl_vector_max_index( get() ); }
    /**
     * C++ version of gsl_vector_min_index().
     * @return index of minimum value of vector
     */
    size_t min_index() const { return gsl_vector_min_index( get() ); }
    /**
     * C++ version of gsl_vector_minmax_index().
     * @param imin index of minimum value of vector
     * @param imax index of maximum value of vector
     */
    void minmax_index( size_t* imin, size_t* imax ) const {
      gsl_vector_minmax_index( get(), imin, imax ); }
    /**
     * C++ version of gsl_vector_add().
     * @param b vector to add to @c this
     * @return error code on failure
     */
    int add( vector const& b ){ return gsl_vector_add( get(), b.get() ); }
    /**
     * C++ version of gsl_vector_sub().
     * @param b vector to subtract from @c this
     * @return error code on failure
     */
    int sub( vector const& b ){ return gsl_vector_sub( get(), b.get() ); }
    /**
     * C++ version of gsl_vector_mul().
     * @param b vector to multiply elementwise with @c this
     * @return error code on failure
     */
    int mul( vector const& b ){ return gsl_vector_mul( get(), b.get() ); }
    /**
     * C++ version of gsl_vector_div().
     * @param b vector to dividev @c this by elementwise
     * @return error code on failure
     */
    int div( vector const& b ){ return gsl_vector_div( get(), b.get() ); }
    /**
     * C++ version of gsl_vector_scale().
     * @param x constant to multiply @c this by
     * @return error code on failure
     */
    int scale( double const x ){ return gsl_vector_scale( get(), x ); }
    /**
     * C++ version of gsl_vector_add_constant().
     * @param x constant to add to each element of @c this
     * @return error code on failure
     */
    int add_constant( double const x ){ return gsl_vector_add_constant( get(), x ); }
    /**
     * C++ version of gsl_vector_isnull().
     * @return @c +1 or @c 0 according as elements are all zero or not
     */
    int isnull() const { return gsl_vector_isnull( get() ); }
    /**
     * C++ version of gsl_vector_ispos().
     * @return @c +1 or @c 0 according as elements are all positive or not
     */
    int ispos() const { return gsl_vector_ispos( get() ); }
    /**
     * C++ version of gsl_vector_isneg().
     * @return @c +1 or @c 0 according as elements are all negative or not
     */
    int isneg() const { return gsl_vector_isneg( get() ); }
    /**
     * C++ version of gsl_vector_isnonneg().
     * @return @c +1 or @c 0 according as elements are all nonnegative or not
     */
    int isnonneg() const { return gsl_vector_isnonneg( get() ); }
    /**
     * C++ version of gsl_vector_get().
     * @param i index of element to get
     * @return value of element
     */
    double get( size_t const i ) const { return gsl_vector_get( get(), i ); }
    /**
     * C++ version of gsl_vector_set().
     * @param i index to set
     * @param x new value for element
     */
    void set( size_t const i, double x ){ gsl_vector_set( get(), i, x ); }
    /**
     * C++ version of gsl_vector_ptr().
     * @param i index of element to get
     * @return pointer to element
     */
    double* ptr( size_t const i ){ return gsl_vector_ptr( get(), i ); }
    double* ptr(){ return get()->data ; }
    /**
     * C++ version of gsl_vector_const_ptr().
     * @param i index of element to get
     * @return pointer to element
     */
    double const* const_ptr( size_t const i ) const { return gsl_vector_const_ptr( get(), i ); }
    /**
     * C++ version of gsl_vector_fread().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fread( FILE* stream ){ return gsl_vector_fread( stream, get() ); }
    /**
     * C++ version of gsl_vector_fwrite().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fwrite( FILE* stream ) const { return gsl_vector_fwrite( stream, get() ); }
    /**
     * C++ version of gsl_vector_fscanf().
     * @param stream A C file stream
     * @return error code on failure
     */
    int fscanf( FILE* stream ){ return gsl_vector_fscanf( stream, get() ); }
    /**
     * C++ version of gsl_vector_fprintf().
     * @param stream A C file stream
     * @param format %d, %e, %f or %g
     * @return error code on failure
     */
    int fprintf( FILE* stream, char const* format ) const {
      return gsl_vector_fprintf( stream, get(), format ); }
    /**
     * C++ version of gsl_vector_alloc_from_block().
     * @param b The block
     * @param offset The offset within the block
     * @param n The number of elements
     * @param stride The stride
     */
    vector( block& b, size_t const offset, size_t const n, size_t const stride ){
      ccgsl_pointer = gsl_vector_alloc_from_block( b.get(), offset, n, stride );
      // just plausibly we could allocate vector but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_vector_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * C++ version of gsl_vector_alloc_from_vector().
     * @param v The vector
     * @param offset The offset
     * @param n The number of elements
     * @param stride The stride
     */
    vector( vector& v, size_t const offset, size_t const n, size_t const stride ){
      ccgsl_pointer = gsl_vector_alloc_from_vector( v.get(), offset, n, stride );
      // just plausibly we could allocate vector but not count
      try { count = new size_t; } catch( std::bad_alloc& e ){
	// try to tidy up before rethrowing
	gsl_vector_free( ccgsl_pointer );
	throw e;
      }
      *count = 1; // initially there is just one reference to ccgsl_pointer
    }
    /**
     * C++ version of gsl_vector_view_array().
     * @param v An array of type double
     * @param n The size of the vector
     * @return A vector
     */
    static vector view_array( double* v, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_view_array( v, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_view_array_with_stride().
     * @param base An array of type double
     * @param stride The stride
     * @param n The size of the vector
     * @return A vector
     */
    static vector view_array_with_stride( double* base, size_t stride, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_view_array_with_stride( base, stride, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_const_view_array().
     * @param v An array of type double
     * @param n The size of the vector
     * @return A vector
     */
    static vector const const_view_array( double const* v, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_const_view_array( v, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_const_view_array_with_stride().
     * @param base An array of type double
     * @param stride The stride
     * @param n The size of the vector
     * @return A vector
     */
    static vector const const_view_array_with_stride( double const* base, size_t stride, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_const_view_array_with_stride( base, stride, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_subvector().
     * @param i The offset
     * @param n The size
     * @return A subvector
     */
    vector subvector( size_t i, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_subvector( get(), i, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_subvector_with_stride().
     * @param i The offset
     * @param stride The stride
     * @param n The size
     * @return A subvector
     */
    vector subvector_with_stride( size_t i, size_t stride, size_t n ){
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_subvector_with_stride( get(), i, stride, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_subvector_with_stride().
     * @param i The offset
     * @param n The size
     * @return A subvector
     */
    vector const const_subvector( size_t i, size_t n ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_const_subvector( get(), i, n ).vector;
      return vector( w );
    }
    /**
     * C++ version of gsl_vector_subvector_with_stride().
     * @param i The offset
     * @param stride The stride
     * @param n The size
     * @return A subvector
     */
    vector const const_subvector_with_stride( size_t i, size_t stride, size_t n ) const {
      gsl_vector* w = static_cast<gsl_vector*>( malloc( sizeof( gsl_vector ) ) );
      *w = gsl_vector_const_subvector_with_stride( get(), i, stride, n ).vector;
      return vector( w );
    }
    // Extra allocators from matrix objects. The definition must come after matrix.
    /**
     * C++ version of gsl_vector_alloc_row_from_matrix().
     * @param m A matrix
     * @param i A row
     * @return A vector
     */
    static vector alloc_row_from_matrix( matrix& m, size_t const i );
    /**
     * C++ version of gsl_vector_alloc_col_from_matrix().
     * @param m A amtrix
     * @param j A column
     * @return A vector
     */
    static vector alloc_col_from_matrix( matrix& m, size_t const j );

/*//-------	operator for mathematical expression added by fujii	-------
	double operator * (vector a);
	vector operator * (matrix a);
	//	* operator returning a inner product defined be (this * a)
	//	* operator returning a inner product defined be (this * a)
	double operator * (vector a){
		size_t size = size();
		double inner_product = 0.0;
		if(size != a.size())
	  	  gsl_error( "inner product cannot be defined between different size vectors ", __FILE__, __LINE__, exception::GSL_EFAILED );
		else{
			for(size_t i=0; i<size; i++)
				inner_product += get()->data[i]*a[i];
		}
		return inner_product;
	}
   //	* operator returning a product defined be (this * a)
	vector operator * (matrix a){
		size_t size1 a.size1();
		size_t size2 a.size2();
		vector product(size2);
		if(get()->size != size1)
	  	  gsl_error( "size doesn't match between vector and matrix ", __FILE__, __LINE__, exception::GSL_EFAILED );
		else{
			for(size_t j=0; j<size2; j++){
				product[j] = 0.0;
				for(size_t i=0; i<size1; i++)
					product[j] += get->data[i]*a[i][j];
			}
		}
		return product;
	}
*/
  };

  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see vector::iterator::operator+()
   */
  inline vector::iterator operator+
  ( vector::iterator::difference_type const n, vector::iterator const& i ){ return i + n; }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see vector::iterator::operator+()
   */
  inline vector::const_iterator operator+
  ( vector::const_iterator::difference_type const n, vector::const_iterator const& i ){ return i + n; }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see vector::iterator::operator+()
   */
  inline vector::reverse_iterator operator+
  ( vector::reverse_iterator::difference_type const n, vector::reverse_iterator const& i ){
    return i + n; }
  
  /**
   * Allows constant to be added to iterator.
   * @param n The constant
   * @param i The iterator
   * @see vector::iterator::operator+()
   */
  inline vector::const_reverse_iterator operator+
  ( vector::const_reverse_iterator::difference_type const n, vector::const_reverse_iterator const& i ){
    return i + n; }
  
};
#endif
