/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * hiplar_matprod.c is based on src/main/array.c from R
 * Modifications for the HiPLAR project copyright (C) 2012-2013  Vendel Szeremi
 * HiPLAR is distributed under the same licence as R.
 *  
 * src/main/array.c 
 *      R : A Computer Language for Statistical Data Analysis
 *      Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *      Copyright (C) 1998-2012   The R Core Team
 *      Copyright (C) 2002-2008   The R Foundation
 *  
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */     


#define USE_RINTERNALS

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <R_ext/Applic.h> /* for dgemm */

#include "hiplar_at.h"
#include "hiplar_dbg.h"
#include "plasma_matprod.h"
#include "magma_matprod.h"


#define _(String) (String)
int PRIMVAL(SEXP);
Rboolean R_has_methods(SEXP);
SEXP R_possible_dispatch(SEXP, SEXP, SEXP, SEXP, Rboolean);
typedef SEXP (*CCODE)();
CCODE (PRIMFUN)(SEXP);
void (SET_PRIMFUN)(SEXP, CCODE);


/*******************************************************************************
 * This requires src/main/objects.c to be patched:
 * attribute_hidden needs to be removed from R_possible_dispatch.
 * 
 * SEXP attribute_hidden
 * R_possible_dispatch(SEXP call, SEXP op, SEXP args, SEXP rho,
 *            Rboolean promisedArgs)
 ******************************************************************************/


SEXP hiplar_matprod(SEXP call, SEXP op, SEXP args, SEXP rho) {

    int ldx, ldy, nrx, ncx, nry, ncy, mode;
    SEXP x = CAR(args), y = CADR(args), xdims, ydims, ans;
    Rboolean sym;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_matprod");
#endif

    if(PRIMVAL(op) == 0 && /* %*% is primitive, the others are .Internal() */
       (IS_S4_OBJECT(x) || IS_S4_OBJECT(y))
       && R_has_methods(op)) {
	SEXP s, value;
	/* Remove argument names to ensure positional matching */
	for(s = args; s != R_NilValue; s = CDR(s)) SET_TAG(s, R_NilValue);
	value = R_possible_dispatch(call, op, args, rho, FALSE);
	if(value) return value;
    }

    sym = isNull(y);
    if (sym && (PRIMVAL(op) > 0)) y = x;
    if ( !(isNumeric(x) || isComplex(x)) || !(isNumeric(y) || isComplex(y)) )
	errorcall(call, _("requires numeric/complex matrix/vector arguments"));

    xdims = getAttrib(x, R_DimSymbol);
    ydims = getAttrib(y, R_DimSymbol);
    ldx = length(xdims);
    ldy = length(ydims);

    if (ldx != 2 && ldy != 2) {		/* x and y non-matrices */
	if (PRIMVAL(op) == 0) {
	    nrx = 1;
	    ncx = LENGTH(x);
	}
	else {
	    nrx = LENGTH(x);
	    ncx = 1;
	}
	nry = LENGTH(y);
	ncy = 1;
    }
    else if (ldx != 2) {		/* x not a matrix */
	nry = INTEGER(ydims)[0];
	ncy = INTEGER(ydims)[1];
	nrx = 0;
	ncx = 0;
	if (PRIMVAL(op) == 0) {
	    if (LENGTH(x) == nry) {	/* x as row vector */
		nrx = 1;
		ncx = nry; /* == LENGTH(x) */
	    }
	    else if (nry == 1) {	/* x as col vector */
		nrx = LENGTH(x);
		ncx = 1;
	    }
	}
	else if (PRIMVAL(op) == 1) { /* crossprod() */
	    if (LENGTH(x) == nry) {	/* x is a col vector */
		nrx = nry; /* == LENGTH(x) */
		ncx = 1;
	    }
	    /* else if (nry == 1) ... not being too tolerant
	       to treat x as row vector, as t(x) *is* row vector */
	}
	else { /* tcrossprod */
	    if (LENGTH(x) == ncy) {	/* x as row vector */
		nrx = 1;
		ncx = ncy; /* == LENGTH(x) */
	    }
	    else if (ncy == 1) {	/* x as col vector */
		nrx = LENGTH(x);
		ncx = 1;
	    }
	}
    }
    else if (ldy != 2) {		/* y not a matrix */
	nrx = INTEGER(xdims)[0];
	ncx = INTEGER(xdims)[1];
	nry = 0;
	ncy = 0;
	if (PRIMVAL(op) == 0) {
	    if (LENGTH(y) == ncx) {	/* y as col vector */
		nry = ncx;
		ncy = 1;
	    }
	    else if (ncx == 1) {	/* y as row vector */
		nry = 1;
		ncy = LENGTH(y);
	    }
	}
	else if (PRIMVAL(op) == 1) { /* crossprod() */
	    if (LENGTH(y) == nrx) {	/* y is a col vector */
		nry = nrx;
		ncy = 1;
	    }
	}
	else { /* tcrossprod --		y is a col vector */
	    nry = LENGTH(y);
	    ncy = 1;
	}
    }
    else {				/* x and y matrices */
	nrx = INTEGER(xdims)[0];
	ncx = INTEGER(xdims)[1];
	nry = INTEGER(ydims)[0];
	ncy = INTEGER(ydims)[1];
    }
    /* nr[ow](.) and nc[ol](.) are now defined for x and y */

    if (PRIMVAL(op) == 0) {
	/* primitive, so use call */
	if (ncx != nry)
	    errorcall(call, _("non-conformable arguments"));
    }
    else if (PRIMVAL(op) == 1) {
	if (nrx != nry)
	    error(_("non-conformable arguments"));
    }
    else {
	if (ncx != ncy)
	    error(_("non-conformable arguments"));
    }


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((nrx < xover_matprod) || (nrx > maxmagma_matprod))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_matprod(call, op, args, rho);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_matprod(call, op, args, rho);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif


}


SEXP hiplar_replace(void) {

    SEXP p_symbol = install("%*%");
    SEXP c_symbol = install("crossprod");
    SEXP t_symbol = install("tcrossprod");
    SEXP p = SYMVALUE(p_symbol);
    SEXP c = INTERNAL(c_symbol);
    SEXP t = INTERNAL(t_symbol);

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_replace");
#endif

    SET_PRIMFUN(p, hiplar_matprod);
    SET_PRIMFUN(c, hiplar_matprod);
    SET_PRIMFUN(t, hiplar_matprod);

    return R_NilValue;
}

