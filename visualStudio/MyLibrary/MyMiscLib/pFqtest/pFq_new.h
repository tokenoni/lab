#pragma once

//#include "f2c.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <MyMiscLib\f2c\f2c\f2c.h>


doublereal z_abs(const doublecomplex* c_){
	return sqrt(c_->r*c_->r + c_->i*c_->i );
};

//	z_div
void z_div(doublecomplex* z_1, doublecomplex*z_2, doublecomplex* z__)
{
	gsl_complex* z_1_ = (gsl_complex* )z_1;
	gsl_complex* z_2_ = (gsl_complex* )z_2;
	gsl_complex* z___ = (gsl_complex* )z__;
	*z_1_ = gsl_complex_div(*z_2_, *z___);
};
// pow_zz
void pow_zz(doublecomplex* z_6, doublecomplex* z_7, doublecomplex* z_8)
{
	gsl_complex* z_6_ = (gsl_complex*)z_6;
	gsl_complex* z_7_ = (gsl_complex*)z_7;
	gsl_complex* z_8_ = (gsl_complex*)z_8;
	*z_6_ = gsl_complex_pow(*z_7_, *z_8_);
};
//	d_imag
doublereal d_imag(doublecomplex* z)
{ return z->i;};
	
doublereal i_dnnt(doublereal* r)
{return floor(*r + 0.5);};

/* Common Block Declarations */
struct consts_1_ {
    doublereal zero, half, one, two, ten, eps;
};

#define consts_1 (*(struct consts_1_ *) &consts_)

struct {
    integer nout;
} io_;

#define io_1 io_

/* Initialized data */

struct {
    doublereal e_1[6];
    } consts_ = { 0., .5, 1., 2., 10., 1e-10 };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__25000 = 25000;
static doublereal c_b64 = 10.;

/* Double Complex */ 
VOID pfq_(doublecomplex * ret_val, doublecomplex *a, 
	doublecomplex *b, integer *ip, integer *iq, doublecomplex *z__, 
	integer *lnpfq, integer *ix, integer *nsigfig);

doublereal bits_(void);

/* Double Complex */ 
VOID hyper_(doublecomplex * ret_val, doublecomplex *a, 
	doublecomplex *b, integer *ip, integer *iq, doublecomplex *z__, 
	integer *lnpfq, integer *ix, integer *nsigfig);

/* Subroutine */ 
int aradd_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *z__, integer *l, doublereal *rmax);

/* Subroutine */ 
int arsub_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *wk1, doublereal *wk2, integer *l, doublereal *rmax);

/* Subroutine */ 
int armult_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *z__, integer *l, doublereal *rmax);

/* Subroutine */ 
int cmpadd_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, doublereal *wk1, 
	integer *l, doublereal *rmax);

/* Subroutine */ 
int cmpsub_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, doublereal *wk1, 
	doublereal *wk2, integer *l, doublereal *rmax);

/* Subroutine */ 
int cmpmul_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, doublereal *wk1, 
	doublereal *wk2, doublereal *cr2, doublereal *d1, doublereal *d2, 
	doublereal *wk6, integer *l, doublereal *rmax);

/* Subroutine */ 
int arydiv_(doublereal *ar, doublereal *ai, doublereal *br, 
	doublereal *bi, doublecomplex *c__, integer *l, integer *lnpfq, 
	doublereal *rmax, integer *ibit);

/* Subroutine */ 
int emult_(doublereal *n1, doublereal *e1, doublereal *n2, 
	doublereal *e2, doublereal *nf, doublereal *ef);

/* Subroutine */ 
int ediv_(doublereal *n1, doublereal *e1, doublereal *n2, 
	doublereal *e2, doublereal *nf, doublereal *ef);

/* Subroutine */ 
int eadd_(doublereal *n1, doublereal *e1, doublereal *n2, 
	doublereal *e2, doublereal *nf, doublereal *ef);

/* Subroutine */ 
int esub_(doublereal *n1, doublereal *e1, doublereal *n2, 
	doublereal *e2, doublereal *nf, doublereal *ef);

/* Subroutine */ 
int conv12_(doublecomplex *cn, doublereal *cae);

/* Subroutine */ 
int conv21_(doublereal *cae, doublecomplex *cn);

/* Subroutine */ 
int ecpmul_(doublereal *a, doublereal *b, doublereal *c__);

/* Subroutine */ int ecpdiv_(doublereal *a, doublereal *b, doublereal *c__);

integer ipremax_(doublecomplex *a, doublecomplex *b, integer *ip, integer *iq,
	 doublecomplex *z__);

/* Double Complex */ 
VOID factor_(doublecomplex * ret_val, doublecomplex *z__);

/* Double Complex */ 
VOID cgamma_(doublecomplex * ret_val, doublecomplex *arg,
	 integer *lnpfq);

/* Subroutine */ 
int bldat1_(void);