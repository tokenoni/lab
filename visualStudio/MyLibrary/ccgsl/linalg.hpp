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

#ifndef CCGSL_LINALG_HPP
#define CCGSL_LINALG_HPP

#include<cmath>
#include<gsl/gsl_linalg.h>
#include"vector.hpp"
#include"matrix.hpp"
#include"complex.hpp"
#include"permutation.hpp"
#include"vector_complex.hpp"
#include"matrix_complex.hpp"

namespace gsl {
  /**
   *  This namespace handles the GSL linear algebra functions.
   */
  namespace linalg {
    /**
     * C++ version of gsl_linalg_householder_transform().
     * @param v A vector
     * @return The Householder transform
     */
    inline double householder_transform( vector& v ){
      return gsl_linalg_householder_transform( v.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_householder_transform().
     * @param v A vector
     * @return The Householder transform
     */
    inline complex complex_householder_transform( vector_complex& v ){
      return gsl_linalg_complex_householder_transform( v.get() ); } 
    /**
     * C++ version of gsl_linalg_householder_hm().
     * @param tau A scalar
     * @param v A vector
     * @param A A matrix
     * @return Error code on failure
     */
    inline int householder_hm( double tau, vector const& v, matrix& A ){
      return gsl_linalg_householder_hm( tau, v.get(), A.get() ); } 
    /**
     * C++ version of gsl_linalg_householder_mh().
     * @param tau A scalar
     * @param v A vector
     * @param A A matrix
     * @return Error code on failure
     */
    inline int householder_mh( double tau, vector const& v, matrix& A ){
      return gsl_linalg_householder_mh( tau, v.get(), A.get() ); } 
    /**
     * C++ version of gsl_linalg_householder_hv().
     * @param tau A scalar
     * @param v A vector
     * @param w Another vector
     * @return Error code on failure
     */
    inline int householder_hv( double tau, vector const& v, vector& w ){
      return gsl_linalg_householder_hv( tau, v.get(), w.get() ); } 
    /**
     * C++ version of gsl_linalg_householder_hm1().
     * @param tau A scalar
     * @param A A matrix
     * @return Error code on failure
     */
    inline int householder_hm1( double tau, matrix& A ){
      return gsl_linalg_householder_hm1( tau, A.get() ); }
    /**
     * C++ version of gsl_linalg_complex_householder_hm().
     * @param tau A scalar
     * @param v A vector
     * @param A A matrix
     * @return Error code on failure
     */
    inline int complex_householder_hm( complex& tau, vector_complex const& v, matrix_complex& A ){
      return gsl_linalg_complex_householder_hm( tau.get(), v.get(), A.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_householder_mh().
     * @param tau A scalar
     * @param v A vector
     * @param A A matrix
     * @return Error code on failure
     */
    inline int complex_householder_mh( complex& tau, vector_complex const& v, matrix_complex& A ){
      return gsl_linalg_complex_householder_mh( tau.get(), v.get(), A.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_householder_hv().
     * @param tau A scalar
     * @param v A vector
     * @param w Another vector
     * @return Error code on failure
     */
    inline int complex_householder_hv( complex& tau, vector_complex const& v, vector_complex& w ){
      return gsl_linalg_complex_householder_hv( tau.get(), v.get(), w.get() ); } 
    /**
     * C++ version of gsl_linalg_hessenberg_decomp().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int hessenberg_decomp( matrix& A, vector& tau ){
          return gsl_linalg_hessenberg_decomp( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_hessenberg_unpack().
     * @param H An upper Hessenberg matrix
     * @param tau A vector
     * @param U An orthogonal matrix
     * @return Error code on failure
     */
    inline int hessenberg_unpack( matrix& H, vector& tau, matrix& U ){
      return gsl_linalg_hessenberg_unpack( H.get(), tau.get(), U.get() ); } 
    /**
     * C++ version of gsl_linalg_hessenberg_unpack_accum().
     * @param H An upper Hessenberg matrix
     * @param tau A vector
     * @param U An orthogonal matrix
     * @return Error code on failure
     */
    inline int hessenberg_unpack_accum( matrix& H, vector& tau, matrix& U ){
      return gsl_linalg_hessenberg_unpack_accum( H.get(), tau.get(), U.get() ); } 
    /**
     * C++ version of gsl_linalg_hessenberg_set_zero().
     * @param H An upper Hessenberg matrix
     * @return Error code on failure
     */
    inline int hessenberg_set_zero( matrix& H ){ return gsl_linalg_hessenberg_set_zero( H.get() ); } 
    /**
     * C++ version of gsl_linalg_hessenberg_submatrix().
     * @param M A matrix
     * @param A A matrix
     * @param top An integer
     * @param tau A vector
     * @return Error code on failure
     */
    inline int hessenberg_submatrix( matrix& M, matrix& A, size_t top, vector& tau ){
      return gsl_linalg_hessenberg_submatrix( M.get(), A.get(), top, tau.get() ); }
    /**
     * C++ version of gsl_linalg_hessenberg().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int hessenberg( matrix& A, vector& tau ){ return gsl_linalg_hessenberg( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_hesstri_decomp().
     * @param A A matrix
     * @param B A matrix
     * @param U An orthogonal matrix
     * @param V An orthogonal matrix
     * @param work A vector
     * @return Error code on failure
     */
    inline int hesstri_decomp( matrix& A, matrix& B, matrix& U, matrix& V, vector& work ){
      return gsl_linalg_hesstri_decomp( A.get(), B.get(), U.get(), V.get(), work.get() ); } 
    /**
     * C++ version of gsl_linalg_SV_decomp().
     * @param A A matrix
     * @param V A matrix (part of SVD)
     * @param S A vector
     * @param work A vector
     * @return Error code on failure
     */
    inline int SV_decomp( matrix& A, matrix& V, vector& S, vector& work ){
      return gsl_linalg_SV_decomp( A.get(), V.get(), S.get(), work.get() ); } 
    /**
     * C++ version of gsl_linalg_SV_decomp_mod().
     * @param A A matrix
     * @param X A matrix
     * @param V A matrix (part of SVD)
     * @param S A vector
     * @param work A vector
     * @return Error code on failure
     */
    inline int SV_decomp_mod( matrix& A, matrix& X, matrix& V, vector& S, vector& work ){
      return gsl_linalg_SV_decomp_mod( A.get(), X.get(), V.get(), S.get(), work.get() ); } 
    /**
     * C++ version of gsl_linalg_SV_decomp_jacobi().
     * @param A A matrix
     * @param Q A matrix
     * @param S A vector
     * @return Error code on failure
     */
    inline int SV_decomp_jacobi( matrix& A, matrix& Q, vector& S ){
      return gsl_linalg_SV_decomp_jacobi( A.get(), Q.get(), S.get() ); } 
    /**
     * C++ version of gsl_linalg_SV_solve().
     * @param U An orthogonal matrix
     * @param Q A matrix
     * @param S A vector
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int SV_solve( matrix const& U, matrix const& Q, vector const& S, vector const& b, vector& x ){
      return gsl_linalg_SV_solve( U.get(), Q.get(), S.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_decomp().
     * @param A A matrix
     * @param p A permutation
     * @param signum The sign of the permutation
     * @return Error code on failure
     */
    inline int LU_decomp( matrix& A, permutation& p, int* signum ){
      return gsl_linalg_LU_decomp( A.get(), p.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_LU_solve().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int LU_solve( matrix const& LU, permutation const& p, vector const& b, vector& x ){
      return gsl_linalg_LU_solve( LU.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_svx().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int LU_svx( matrix const& LU, permutation const& p, vector& x ){
      return gsl_linalg_LU_svx( LU.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_refine().
     * @param A A matrix
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @param residual A residual vector
     * @return Error code on failure
     */
    inline int LU_refine( matrix const& A, matrix const& LU, permutation const& p, vector const& b,
		   vector& x, vector& residual ){
      return gsl_linalg_LU_refine( A.get(), LU.get(), p.get(), b.get(), x.get(), residual.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_invert().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param inverse The inverse matrix
     * @return Error code on failure
     */
    inline int LU_invert( matrix const& LU, permutation const& p, matrix& inverse ){
      return gsl_linalg_LU_invert( LU.get(), p.get(), inverse.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_det().
     * @param LU An LU decomposition matrix
     * @param signum The sign of the permutation
     * @return The determinant of the matrix with LU decomposition @c LU
     */
    inline double LU_det( matrix& LU, int signum ){ return gsl_linalg_LU_det( LU.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_LU_lndet().
     * @param LU An LU decomposition matrix
     * @return The logarithm of the absolute value of
     * the determinant of the matrix with LU decomposition @c LU
     */
    inline double LU_lndet( matrix& LU ){ return gsl_linalg_LU_lndet( LU.get() ); } 
    /**
     * C++ version of gsl_linalg_LU_sgndet().
     * @param lu An LU decomposition marix
     * @param signum The sign of the permutation
     * @return The sign of the determinant of the matrix with LU decomposition @c LU
     */
    inline int LU_sgndet( matrix& lu, int signum ){ return gsl_linalg_LU_sgndet( lu.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_decomp().
     * @param A A matrix
     * @param p A permutation
     * @param signum The sign of the permutation
     * @return Error code on failure
     */
    inline int complex_LU_decomp( matrix_complex& A, permutation& p, int* signum ){
      return gsl_linalg_complex_LU_decomp( A.get(), p.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_solve().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int complex_LU_solve( matrix_complex const& LU, permutation const& p,
				 vector_complex const& b, vector_complex& x ){
      return gsl_linalg_complex_LU_solve( LU.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_svx().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int complex_LU_svx( matrix_complex const& LU, permutation const& p, vector_complex& x ){
      return gsl_linalg_complex_LU_svx( LU.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_refine().
     * @param A A matrix
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @param residual A residual vector
     * @return Error code on failure
     */
    inline int complex_LU_refine( matrix_complex const& A, matrix_complex const& LU,
				  permutation const& p, vector_complex const& b,
				  vector_complex& x, vector_complex& residual ){
      return gsl_linalg_complex_LU_refine( A.get(), LU.get(), p.get(), b.get(), x.get(), residual.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_invert().
     * @param LU An LU decomposition matrix
     * @param p A permutation
     * @param inverse The inverse matrix
     * @return Error code on failure
     */
    inline int complex_LU_invert( matrix_complex const& LU, permutation const& p,
				  matrix_complex& inverse ){
      return gsl_linalg_complex_LU_invert( LU.get(), p.get(), inverse.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_det().
     * @param LU An LU decomposition matrix
     * @param signum The sign of the permutation
     * @return Error code on failure
     */
    inline complex complex_LU_det( matrix_complex& LU, int signum ){
      return gsl_linalg_complex_LU_det( LU.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_complex_LU_lndet().
     * @param LU An LU decomposition matrix
     * @return The determinant of the matrix with LU decomposition @c LU 
     */
    inline double complex_LU_lndet( matrix_complex& LU ){
      return gsl_linalg_complex_LU_lndet( LU.get() ); }
    /**
     * C++ version of gsl_linalg_complex_LU_sgndet().
     * @param LU An LU decomposition matrix
     * @param signum The sign of the permutation
     * @return The logarithm of the determinant of the matrix with LU decomposition @c LU
     */
    inline complex complex_LU_sgndet( matrix_complex& LU, int signum ){
      return gsl_linalg_complex_LU_sgndet( LU.get(), signum ); } 
    /**
     * C++ version of gsl_linalg_QR_decomp().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int QR_decomp( matrix& A, vector& tau ){ return gsl_linalg_QR_decomp( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_solve().
     * @param QR 
     * @param tau A vector
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QR_solve( matrix const& QR, vector const& tau, vector const& b, vector& x ){
      return gsl_linalg_QR_solve( QR.get(), tau.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_svx().
     * @param QR A QR decomposition
     * @param tau A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QR_svx( matrix const& QR, vector const& tau, vector& x ){
      return gsl_linalg_QR_svx( QR.get(), tau.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_lssolve().
     * @param QR A QR decomposition
     * @param tau A vector
     * @param b A vector
     * @param x A vector
     * @param residual A residual vector
     * @return Error code on failure
     */
    inline int QR_lssolve( matrix const& QR, vector const& tau, vector const& b, vector& x,
			   vector& residual ){
      return gsl_linalg_QR_lssolve( QR.get(), tau.get(), b.get(), x.get(), residual.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_QRsolve().
     * @param Q A Matrix
     * @param R A Matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QR_QRsolve( matrix& Q, matrix& R, vector const& b, vector& x ){
      return gsl_linalg_QR_QRsolve( Q.get(), R.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_Rsolve().
     * @param QR A QR decomposition matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QR_Rsolve( matrix const& QR, vector const& b, vector& x ){
      return gsl_linalg_QR_Rsolve( QR.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_Rsvx().
     * @param QR A QR decomposition matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int QR_Rsvx( matrix const& QR, vector& x ){ return gsl_linalg_QR_Rsvx( QR.get(), x.get() ); }
    /**
     * C++ version of gsl_linalg_QR_update().
     * @param Q A matrix
     * @param R A matrix
     * @param w A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int QR_update( matrix& Q, matrix& R, vector& w, vector const& v ){
      return gsl_linalg_QR_update( Q.get(), R.get(), w.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_QTvec().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int QR_QTvec( matrix const& QR, vector const& tau, vector& v ){
      return gsl_linalg_QR_QTvec( QR.get(), tau.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_Qvec().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int QR_Qvec( matrix const& QR, vector const& tau, vector& v ){
      return gsl_linalg_QR_Qvec( QR.get(), tau.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_QTmat().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param A A matrix
     * @return Error code on failure
     */
    inline int QR_QTmat( matrix const& QR, vector const& tau, matrix& A ){
      return gsl_linalg_QR_QTmat( QR.get(), tau.get(), A.get() ); } 
    /**
     * C++ version of gsl_linalg_QR_unpack().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param Q A matrix
     * @param R A matrix
     * @return Error code on failure
     */
    inline int QR_unpack( matrix const& QR, vector const& tau, matrix& Q, matrix& R ){
      return gsl_linalg_QR_unpack( QR.get(), tau.get(), Q.get(), R.get() ); } 
    /**
     * C++ version of gsl_linalg_R_solve().
     * @param R A matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int R_solve( matrix const& R, vector const& b, vector& x ){
      return gsl_linalg_R_solve( R.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_R_svx().
     * @param R A matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int R_svx( matrix const& R, vector& x ){ return gsl_linalg_R_svx( R.get(), x.get() ); }
    /**
     * C++ version of gsl_linalg_QRPT_decomp().
     * @param A A matrix
     * @param tau A vector
     * @param p A permutation
     * @param signum The sign of the permutation
     * @param norm A vector used for column pivoting
     * @return Error code on failure
     */
    inline int QRPT_decomp( matrix& A, vector& tau, permutation& p, int* signum, vector& norm ){
      return gsl_linalg_QRPT_decomp( A.get(), tau.get(), p.get(), signum, norm.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_decomp2().
     * @param A A matrix
     * @param q A matrix
     * @param r A matrix
     * @param tau A vector
     * @param p A permutation
     * @param signum The sign of the permutation
     * @param norm A vector used for column pivoting
     * @return Error code on failure
     */
    inline int QRPT_decomp2( matrix const& A, matrix& q, matrix& r, vector& tau,
			     permutation& p, int* signum, vector& norm ){
      return gsl_linalg_QRPT_decomp2( A.get(), q.get(), r.get(), tau.get(), p.get(), signum, norm.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_solve().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QRPT_solve( matrix const& QR, vector const& tau, permutation const& p,
			   vector const& b, vector& x ){
      return gsl_linalg_QRPT_solve( QR.get(), tau.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_svx().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int QRPT_svx( matrix const& QR, vector const& tau, permutation const& p, vector& x ){
      return gsl_linalg_QRPT_svx( QR.get(), tau.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_QRsolve().
     * @param Q A matrix
     * @param R A matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QRPT_QRsolve( matrix const& Q, matrix const& R, permutation const& p,
			     vector const& b, vector& x ){
      return gsl_linalg_QRPT_QRsolve( Q.get(), R.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_Rsolve().
     * @param QR A QR decomposition matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int QRPT_Rsolve( matrix const& QR, permutation const& p, vector const& b, vector& x ){
      return gsl_linalg_QRPT_Rsolve( QR.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_Rsvx().
     * @param QR A QR decomposition matrix
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int QRPT_Rsvx( matrix const& QR, permutation const& p, vector& x ){
      return gsl_linalg_QRPT_Rsvx( QR.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_QRPT_update().
     * @param Q A matrix
     * @param R A matrix
     * @param p A permutation
     * @param u A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int QRPT_update( matrix& Q, matrix& R, permutation const& p, vector& u, vector const& v ){
      return gsl_linalg_QRPT_update( Q.get(), R.get(), p.get(), u.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_decomp().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int LQ_decomp( matrix& A, vector& tau ){
      return gsl_linalg_LQ_decomp( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_solve_T().
     * @param LQ A matrix
     * @param tau A vector
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int LQ_solve_T( matrix const& LQ, vector const& tau, vector const& b, vector& x ){
      return gsl_linalg_LQ_solve_T( LQ.get(), tau.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_svx_T().
     * @param LQ A matrix
     * @param tau A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int LQ_svx_T( matrix const& LQ, vector const& tau, vector& x ){
      return gsl_linalg_LQ_svx_T( LQ.get(), tau.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_lssolve_T().
     * @param LQ A matrix
     * @param tau A vector
     * @param b A vector
     * @param x A vector
     * @param residual A residual vector
     * @return Error code on failure
     */
    inline int LQ_lssolve_T( matrix const& LQ, vector const& tau, vector const& b,
			     vector& x, vector& residual ){
      return gsl_linalg_LQ_lssolve_T( LQ.get(), tau.get(), b.get(), x.get(), residual.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_Lsolve_T().
     * @param LQ A matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int LQ_Lsolve_T( matrix const& LQ, vector const& b, vector& x ){
      return gsl_linalg_LQ_Lsolve_T( LQ.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_Lsvx_T().
     * @param LQ A matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int LQ_Lsvx_T( matrix const& LQ, vector& x ){
      return gsl_linalg_LQ_Lsvx_T( LQ.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_L_solve_T().
     * @param L A matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int L_solve_T( matrix const& L, vector const& b, vector& x ){
      return gsl_linalg_L_solve_T( L.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_vecQ().
     * @param LQ A matrix
     * @param tau A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int LQ_vecQ( matrix const& LQ, vector const& tau, vector& v ){
      return gsl_linalg_LQ_vecQ( LQ.get(), tau.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_vecQT().
     * @param LQ A matrix
     * @param tau A vector
     * @param v A vector
     * @return Error code on failure
     */
    inline int LQ_vecQT( matrix const& LQ, vector const& tau, vector& v ){
      return gsl_linalg_LQ_vecQT( LQ.get(), tau.get(), v.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_unpack().
     * @param LQ A matrix
     * @param tau A vector
     * @param Q A matrix
     * @param L A matrix
     * @return Error code on failure
     */
    inline int LQ_unpack( matrix const& LQ, vector const& tau, matrix& Q, matrix& L ){
      return gsl_linalg_LQ_unpack( LQ.get(), tau.get(), Q.get(), L.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_update().
     * @param Q A matrix
     * @param R A matrix
     * @param v A vector
     * @param w A vector
     * @return Error code on failure
     */
    inline int LQ_update( matrix& Q, matrix& R, vector const& v, vector& w ){
      return gsl_linalg_LQ_update( Q.get(), R.get(), v.get(), w.get() ); } 
    /**
     * C++ version of gsl_linalg_LQ_LQsolve().
     * @param Q A matrix
     * @param L A matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int LQ_LQsolve( matrix& Q, matrix& L, vector const& b, vector& x ){
      return gsl_linalg_LQ_LQsolve( Q.get(), L.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_decomp().
     * @param A A matrix
     * @param tau A vector
     * @param p A permutation
     * @param signum The sign of the permutation
     * @param norm A norm vector
     * @return Error code on failure
     */
    inline int PTLQ_decomp( matrix& A, vector& tau, permutation& p, int* signum, vector& norm ){
      return gsl_linalg_PTLQ_decomp( A.get(), tau.get(), p.get(), signum, norm.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_decomp2().
     * @param A A matrix
     * @param q A matrix
     * @param r A matrix
     * @param tau A vector
     * @param p A permutation
     * @param signum The sign of the permutation
     * @param norm A vector used as workspace
     * @return Error code on failure
     */
    inline int PTLQ_decomp2( matrix const& A, matrix& q, matrix& r, vector& tau,
			     permutation& p, int* signum, vector& norm ){
      return gsl_linalg_PTLQ_decomp2( A.get(), q.get(), r.get(), tau.get(), p.get(), signum, norm.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_solve_T().
     * @param QR A QR decomposition matrix
     * @param tau A vector
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int PTLQ_solve_T( matrix const& QR, vector const& tau, permutation const& p,
			     vector const& b, vector& x ){
      return gsl_linalg_PTLQ_solve_T( QR.get(), tau.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_svx_T().
     * @param LQ A matrix
     * @param tau A vector
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int PTLQ_svx_T( matrix const& LQ, vector const& tau, permutation const& p, vector& x ){
      return gsl_linalg_PTLQ_svx_T( LQ.get(), tau.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_LQsolve_T().
     * @param Q A matrix
     * @param L A matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int PTLQ_LQsolve_T( matrix const& Q, matrix const& L, permutation const& p,
			       vector const& b, vector& x ){
      return gsl_linalg_PTLQ_LQsolve_T( Q.get(), L.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_Lsolve_T().
     * @param LQ A matrix
     * @param p A permutation
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int PTLQ_Lsolve_T( matrix const& LQ, permutation const& p, vector const& b, vector& x ){
      return gsl_linalg_PTLQ_Lsolve_T( LQ.get(), p.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_Lsvx_T().
     * @param LQ A matrix
     * @param p A permutation
     * @param x A vector
     * @return Error code on failure
     */
    inline int PTLQ_Lsvx_T( matrix const& LQ, permutation const& p, vector& x ){
      return gsl_linalg_PTLQ_Lsvx_T( LQ.get(), p.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_PTLQ_update().
     * @param Q A matrix
     * @param L A matrix
     * @param p A permutation
     * @param v A vector
     * @param w A vector
     * @return Error code on failure
     */
    inline int PTLQ_update( matrix& Q, matrix& L, permutation const& p, vector const& v, vector& w ){
      return gsl_linalg_PTLQ_update( Q.get(), L.get(), p.get(), v.get(), w.get() ); } 
    /**
     * C++ version of gsl_linalg_cholesky_decomp().
     * @param A A matrix
     * @return Error code on failure
     */
    inline int cholesky_decomp( matrix& A ){ return gsl_linalg_cholesky_decomp( A.get() ); } 
    /**
     * C++ version of gsl_linalg_cholesky_solve().
     * @param cholesky A Cholesky decomposition matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int cholesky_solve( matrix const& cholesky, vector const& b, vector& x ){
      return gsl_linalg_cholesky_solve( cholesky.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_cholesky_svx().
     * @param cholesky A Cholesky decomposition matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int cholesky_svx( matrix const& cholesky, vector& x ){
      return gsl_linalg_cholesky_svx( cholesky.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_cholesky_invert().
     * @param cholesky A Cholesky decomposition matrix
     * @return Error code on failure
     */
    inline int cholesky_invert( matrix& cholesky ){ return gsl_linalg_cholesky_invert( cholesky.get() ); } 
    /**
     * C++ version of gsl_linalg_cholesky_decomp_unit().
     * @param A A matrix
     * @param D A vector
     * @return Error code on failure
     */
    inline int cholesky_decomp_unit( matrix& A, vector& D ){
      return gsl_linalg_cholesky_decomp_unit( A.get(), D.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_cholesky_decomp().
     * @param A A matrix
     * @return Error code on failure
     */
    inline int complex_cholesky_decomp( matrix_complex& A ){
      return gsl_linalg_complex_cholesky_decomp( A.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_cholesky_solve().
     * @param cholesky A Cholesky decomposition matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int complex_cholesky_solve( matrix_complex const& cholesky, vector_complex const& b,
				       vector_complex& x ){
      return gsl_linalg_complex_cholesky_solve( cholesky.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_complex_cholesky_svx().
     * @param cholesky A Cholesky decomposition matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int complex_cholesky_svx( matrix_complex const& cholesky, vector_complex& x ){
      return gsl_linalg_complex_cholesky_svx( cholesky.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_symmtd_decomp().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int symmtd_decomp( matrix& A, vector& tau ){
      return gsl_linalg_symmtd_decomp( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_symmtd_unpack().
     * @param A A matrix
     * @param tau A vector
     * @param Q A matrix
     * @param diag A vector of diagonal elements
     * @param subdiag The vector of subdiagonal elements
     * @return Error code on failure
     */
    inline int symmtd_unpack( matrix const& A, vector const& tau, matrix& Q,
			      vector& diag, vector& subdiag ){
      return gsl_linalg_symmtd_unpack( A.get(), tau.get(), Q.get(), diag.get(), subdiag.get() ); } 
    /**
     * C++ version of gsl_linalg_symmtd_unpack_T().
     * @param A A matrix
     * @param diag A vector of diagonal elements
     * @param subdiag The vector of subdiagonal elements
     * @return Error code on failure
     */
    inline int symmtd_unpack_T( matrix const& A, vector& diag, vector& subdiag ){
      return gsl_linalg_symmtd_unpack_T( A.get(), diag.get(), subdiag.get() ); } 
    /**
     * C++ version of gsl_linalg_hermtd_decomp().
     * @param A A matrix
     * @param tau A vector
     * @return Error code on failure
     */
    inline int hermtd_decomp( matrix_complex& A, vector_complex& tau ){
      return gsl_linalg_hermtd_decomp( A.get(), tau.get() ); } 
    /**
     * C++ version of gsl_linalg_hermtd_unpack().
     * @param A A matrix
     * @param tau A vector
     * @param U A unitary matrix
     * @param diag A vector of diagonal elements
     * @param sudiag The vector of subdiagonal elements
     * @return Error code on failure
     */
    inline int hermtd_unpack( matrix_complex const& A, vector_complex const& tau,
			      matrix_complex& U, vector& diag, vector& sudiag ){
      return gsl_linalg_hermtd_unpack( A.get(), tau.get(), U.get(), diag.get(), sudiag.get() ); } 
    /**
     * C++ version of gsl_linalg_hermtd_unpack_T().
     * @param A A matrix
     * @param diag A vector of diagonal elements
     * @param subdiag The vector of subdiagonal elements
     * @return Error code on failure
     */
    inline int hermtd_unpack_T( matrix_complex const& A, vector& diag, vector& subdiag ){
      return gsl_linalg_hermtd_unpack_T( A.get(), diag.get(), subdiag.get() ); } 
    /**
     * C++ version of gsl_linalg_HH_solve().
     * @param A A matrix
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int HH_solve( matrix& A, vector const& b, vector& x ){
      return gsl_linalg_HH_solve( A.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_HH_svx().
     * @param A A matrix
     * @param x A vector
     * @return Error code on failure
     */
    inline int HH_svx( matrix& A, vector& x ){ return gsl_linalg_HH_svx( A.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_solve_symm_tridiag().
     * @param diag A vector of diagonal elements
     * @param offdiag Off-diagonal vector (one element shorte than @c diag)
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int solve_symm_tridiag( vector const& diag, vector const& offdiag, vector const& b,
				   vector& x ){
      return gsl_linalg_solve_symm_tridiag( diag.get(), offdiag.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_solve_tridiag().
     * @param diag A vector of diagonal elements
     * @param abovediag Off-diagonal vector (one element shorte than @c diag)
     * @param belowdiag Off-diagonal vector (one element shorte than @c diag)
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int solve_tridiag( vector const& diag, vector const& abovediag, vector const& belowdiag,
			      vector const& b, vector& x ){
      return gsl_linalg_solve_tridiag( diag.get(), abovediag.get(), belowdiag.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_solve_symm_cyc_tridiag().
     * @param diag A vector of diagonal elements
     * @param offdiag Off-diagonal vector (one element shorte than @c diag)
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int solve_symm_cyc_tridiag( vector const& diag, vector const& offdiag,
				       vector const& b, vector& x ){
      return gsl_linalg_solve_symm_cyc_tridiag( diag.get(), offdiag.get(), b.get(), x.get() ); } 
    /**
     * C++ version of gsl_linalg_solve_cyc_tridiag().
     * @param diag A vector of diagonal elements
     * @param abovediag Off-diagonal vector (one element shorte than @c diag)
     * @param belowdiag Off-diagonal vector (one element shorte than @c diag)
     * @param b A vector
     * @param x A vector
     * @return Error code on failure
     */
    inline int solve_cyc_tridiag( vector const& diag, vector const& abovediag, vector const& belowdiag,
				  vector const& b, vector& x ){
      return gsl_linalg_solve_cyc_tridiag( diag.get(), abovediag.get(), belowdiag.get(), b.get(), x.get() ); }
    /**
     * C++ version of gsl_linalg_bidiag_decomp().
     * @param A A matrix
     * @param tau_U A vector
     * @param tau_V A vector
     * @return Error code on failure
     */
    inline int bidiag_decomp( matrix& A, vector& tau_U, vector& tau_V ){
      return gsl_linalg_bidiag_decomp( A.get(), tau_U.get(), tau_V.get() ); } 
    /**
     * C++ version of gsl_linalg_bidiag_unpack().
     * @param A A matrix
     * @param tau_U A vector
     * @param U An orthogonal matrix
     * @param tau_V A vector
     * @param V AN orthogonal matrix
     * @param diag A vector of diagonal elements
     * @param superdiag Off-diagonal vector (one element shorte than @c diag)
     * @return Error code on failure
     */
    inline int bidiag_unpack( matrix const& A, vector const& tau_U, matrix& U,
			      vector const& tau_V, matrix& V, vector& diag, vector& superdiag ){
      return gsl_linalg_bidiag_unpack( A.get(), tau_U.get(), U.get(), tau_V.get(),
				       V.get(), diag.get(), superdiag.get() ); }
    /**
     * C++ version of gsl_linalg_bidiag_unpack2().
     * @param A A matrix
     * @param tau_U A vector
     * @param tau_V A vector
     * @param V AN orthogonal matrix
     * @return Error code on failure
     */
    inline int bidiag_unpack2( matrix& A, vector& tau_U, vector& tau_V, matrix& V ){
      return gsl_linalg_bidiag_unpack2( A.get(), tau_U.get(), tau_V.get(), V.get() ); } 
    /**
     * C++ version of gsl_linalg_bidiag_unpack_B().
     * @param A A matrix
     * @param diag A vector of diagonal elements
     * @param superdiag Off-diagonal vector (one element shorte than @c diag)
     * @return Error code on failure
     */
    inline int bidiag_unpack_B( matrix const& A, vector& diag, vector& superdiag ){
      return gsl_linalg_bidiag_unpack_B( A.get(), diag.get(), superdiag.get() ); } 
    /**
     * C++ version of gsl_linalg_balance_matrix().
     * @param A A matrix
     * @param D A vector
     * @return Error code on failure
     */
    inline int balance_matrix( matrix& A, vector& D ){
      return gsl_linalg_balance_matrix( A.get(), D.get() ); } 
    /**
     * C++ version of gsl_linalg_balance_accum().
     * @param A A matrix
     * @param D A vector
     * @return Error code on failure
     */
    inline int balance_accum( matrix& A, vector& D ){
      return gsl_linalg_balance_accum( A.get(), D.get() ); } 
    /**
     * C++ version of gsl_linalg_balance_columns().
     * @param A A matrix
     * @param D A vector
     * @return Error code on failure
     */
    inline int balance_columns( matrix& A, vector& D ){
      return gsl_linalg_balance_columns( A.get(), D.get() ); } 
  }
}

#endif
