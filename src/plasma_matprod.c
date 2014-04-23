/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * plasma_matprod.c is based on src/main/array.c from R
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

#include "hiplar_dbg.h"
#include "P.h"

#include <R_ext/Applic.h> /* for dgemm */

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


static void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0;
    long double sum;
    Rboolean have_na = FALSE;
    int info;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
	 * The test is only O(n) here
	 */
	for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
	if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
		if (ISNAN(y[i])) {have_na = TRUE; break;}
	if (have_na) {
	    for (i = 0; i < nrx; i++)
		for (k = 0; k < ncy; k++) {
		    sum = 0.0;
		    for (j = 0; j < ncx; j++)
			sum += x[i + j * nrx] * y[j + k * nry];
		    z[i + k * nrx] = sum;
		}
	} else {
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: matprod");
#endif
/*
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
*/
        info = P_dgemm(transa, transb, nrx, ncy, ncx, one,
			x, nrx, y, nry, zero, z, nrx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
        }
    }

    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
#endif
}

static void cmatprod(Rcomplex *x, int nrx, int ncx,
		     Rcomplex *y, int nry, int ncy, Rcomplex *z) {
#ifdef HIPLAR_WITH_PLASMA
/* #ifdef HAVE_FORTRAN_DOUBLE_COMPLEX */
#if(1)
    char *transa = "N", *transb = "N";
    int i, info;
    Rcomplex one, zero;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: cmatprod");
#endif

    one.r = 1.0; one.i = zero.r = zero.i = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
/*
	F77_CALL(zgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
*/
        info = P_zgemm(transa, transb, nrx, ncy, ncx, one,
			x, nrx, y, nry, zero, z, nrx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_zgemm", info);
        }
    } else { /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i].r = z[i].i = 0;
    }
#else
    int i, j, k;
    double xij_r, xij_i, yjk_r, yjk_i;
    long double sum_i, sum_r;

    for (i = 0; i < nrx; i++)
	for (k = 0; k < ncy; k++) {
	    z[i + k * nrx].r = NA_REAL;
	    z[i + k * nrx].i = NA_REAL;
	    sum_r = 0.0;
	    sum_i = 0.0;
	    for (j = 0; j < ncx; j++) {
		xij_r = x[i + j * nrx].r;
		xij_i = x[i + j * nrx].i;
		yjk_r = y[j + k * nry].r;
		yjk_i = y[j + k * nry].i;
		if (ISNAN(xij_r) || ISNAN(xij_i)
		    || ISNAN(yjk_r) || ISNAN(yjk_i))
		    goto next_ik;
		sum_r += (xij_r * yjk_r - xij_i * yjk_i);
		sum_i += (xij_r * yjk_i + xij_i * yjk_r);
	    }
	    z[i + k * nrx].r = sum_r;
	    z[i + k * nrx].i = sum_i;
	next_ik:
	    ;
	}
#endif
#endif
}

static void symcrossprod(double *x, int nr, int nc, double *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *trans = "T", *uplo = "U";
    double one = 1.0, zero = 0.0;
    int i, j;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: symcrossprod");
#endif

    if (nr > 0 && nc > 0) {
/*
	F77_CALL(dsyrk)(uplo, trans, &nc, &nr, &one, x, &nr, &zero, z, &nc);
*/
    info = P_dsyrk(uplo, trans, nc, nr, one, x, nr, zero, z, nc);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dsyrk", info);
    }

	for (i = 1; i < nc; i++)
	    for (j = 0; j < i; j++) z[i + nc *j] = z[j + nc * i];
    } else { /* zero-extent operations should return zeroes */
	for(i = 0; i < nc*nc; i++) z[i] = 0;
    }

#endif
}

static void crossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: crossprod");
#endif

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
/*
	F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
*/
        info = P_dgemm(transa, transb, ncx, ncy, nrx, one,
			x, nrx, y, nry, zero, z, ncx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
        }

    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i] = 0;
    }
#endif
}

static void ccrossprod(Rcomplex *x, int nrx, int ncx,
		       Rcomplex *y, int nry, int ncy, Rcomplex *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *transa = "T", *transb = "N";
    Rcomplex one, zero;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: ccrossprod");
#endif

    one.r = 1.0; one.i = zero.r = zero.i = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
/*
	F77_CALL(zgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
*/
        info = P_zgemm(transa, transb, ncx, ncy, nrx, one,
			x, nrx, y, nry, zero, z, ncx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_zgemm", info);
        }
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i].r = z[i].i = 0;
    }
#endif
}

static void symtcrossprod(double *x, int nr, int nc, double *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *trans = "N", *uplo = "U";
    double one = 1.0, zero = 0.0;
    int i, j;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: symtcrossprod");
#endif

    if (nr > 0 && nc > 0) {
/*
	F77_CALL(dsyrk)(uplo, trans, &nr, &nc, &one, x, &nr, &zero, z, &nr);
*/
    info = P_dsyrk(uplo, trans, nr, nc, one, x, nr, zero, z, nr);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dsyrk", info);
    }

	for (i = 1; i < nr; i++)
	    for (j = 0; j < i; j++) z[i + nr *j] = z[j + nr * i];
    } else { /* zero-extent operations should return zeroes */
	for(i = 0; i < nr*nr; i++) z[i] = 0;
    }

#endif
}

static void tcrossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: tcrossprod");
#endif

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
/*
	F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
*/
        info = P_dgemm(transa, transb, nrx, nry, ncx, one,
			x, nrx, y, nry, zero, z, nrx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
        }

    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
    }
#endif
}

static void tccrossprod(Rcomplex *x, int nrx, int ncx,
			Rcomplex *y, int nry, int ncy, Rcomplex *z) {
#ifdef HIPLAR_WITH_PLASMA
    char *transa = "N", *transb = "T";
    Rcomplex one, zero;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: tccrossprod");
#endif

    one.r = 1.0; one.i = zero.r = zero.i = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
/*
	F77_CALL(zgemm)(transa, transb, &nrx, &nry, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
*/
        info = P_zgemm(transa, transb, nrx, nry, ncx, one,
			x, nrx, y, nry, zero, z, nrx);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_zgemm", info);
        }

    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < nrx*nry; i++) z[i].r = z[i].i = 0;
    }
#endif
}


SEXP plasma_matprod(SEXP call, SEXP op, SEXP args, SEXP rho) {
#ifdef HIPLAR_WITH_PLASMA

    int ldx, ldy, nrx, ncx, nry, ncy, mode;
    SEXP x = CAR(args), y = CADR(args), xdims, ydims, ans;
    Rboolean sym;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_matprod");
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

    if (isComplex(CAR(args)) || isComplex(CADR(args)))
	mode = CPLXSXP;
    else
	mode = REALSXP;
    SETCAR(args, coerceVector(CAR(args), mode));
    SETCADR(args, coerceVector(CADR(args), mode));

    if (PRIMVAL(op) == 0) {			/* op == 0 : matprod() */

	PROTECT(ans = allocMatrix(mode, nrx, ncy));
	if (mode == CPLXSXP)
	    cmatprod(COMPLEX(CAR(args)), nrx, ncx,
		     COMPLEX(CADR(args)), nry, ncy, COMPLEX(ans));
	else
	    matprod(REAL(CAR(args)), nrx, ncx,
		    REAL(CADR(args)), nry, ncy, REAL(ans));

	PROTECT(xdims = getAttrib(CAR(args), R_DimNamesSymbol));
	PROTECT(ydims = getAttrib(CADR(args), R_DimNamesSymbol));

	if (xdims != R_NilValue || ydims != R_NilValue) {
	    SEXP dimnames, dimnamesnames, dnx=R_NilValue, dny=R_NilValue;

	    /* allocate dimnames and dimnamesnames */

	    PROTECT(dimnames = allocVector(VECSXP, 2));
	    PROTECT(dimnamesnames = allocVector(STRSXP, 2));
	    if (xdims != R_NilValue) {
		if (ldx == 2 || ncx == 1) {
		    SET_VECTOR_ELT(dimnames, 0, VECTOR_ELT(xdims, 0));
		    dnx = getAttrib(xdims, R_NamesSymbol);
		    if(!isNull(dnx))
			SET_STRING_ELT(dimnamesnames, 0, STRING_ELT(dnx, 0));
		}
	    }

#define YDIMS_ET_CETERA							\
	    if (ydims != R_NilValue) {					\
		if (ldy == 2) {						\
		    SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(ydims, 1));	\
		    dny = getAttrib(ydims, R_NamesSymbol);		\
		    if(!isNull(dny))					\
			SET_STRING_ELT(dimnamesnames, 1, STRING_ELT(dny, 1)); \
		} else if (nry == 1) {					\
		    SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(ydims, 0));	\
		    dny = getAttrib(ydims, R_NamesSymbol);		\
		    if(!isNull(dny))					\
			SET_STRING_ELT(dimnamesnames, 1, STRING_ELT(dny, 0)); \
		}							\
	    }								\
									\
	    /* We sometimes attach a dimnames attribute			\
	     * whose elements are all NULL ...				\
	     * This is ugly but causes no real damage.			\
	     * Now (2.1.0 ff), we don't anymore: */			\
	    if (VECTOR_ELT(dimnames,0) != R_NilValue ||			\
		VECTOR_ELT(dimnames,1) != R_NilValue) {			\
		if (dnx != R_NilValue || dny != R_NilValue)		\
		    setAttrib(dimnames, R_NamesSymbol, dimnamesnames);	\
		setAttrib(ans, R_DimNamesSymbol, dimnames);		\
	    }								\
	    UNPROTECT(2)

	    YDIMS_ET_CETERA;
	}
    }

    else if (PRIMVAL(op) == 1) {	/* op == 1: crossprod() */

	PROTECT(ans = allocMatrix(mode, ncx, ncy));
	if (mode == CPLXSXP)
	    if(sym)
		ccrossprod(COMPLEX(CAR(args)), nrx, ncx,
			   COMPLEX(CAR(args)), nry, ncy, COMPLEX(ans));
	    else
		ccrossprod(COMPLEX(CAR(args)), nrx, ncx,
			   COMPLEX(CADR(args)), nry, ncy, COMPLEX(ans));
	else {
	    if(sym)
		symcrossprod(REAL(CAR(args)), nrx, ncx, REAL(ans));
	    else
		crossprod(REAL(CAR(args)), nrx, ncx,
			  REAL(CADR(args)), nry, ncy, REAL(ans));
	}

	PROTECT(xdims = getAttrib(CAR(args), R_DimNamesSymbol));
	if (sym)
	    PROTECT(ydims = xdims);
	else
	    PROTECT(ydims = getAttrib(CADR(args), R_DimNamesSymbol));

	if (xdims != R_NilValue || ydims != R_NilValue) {
	    SEXP dimnames, dimnamesnames, dnx=R_NilValue, dny=R_NilValue;

	    /* allocate dimnames and dimnamesnames */

	    PROTECT(dimnames = allocVector(VECSXP, 2));
	    PROTECT(dimnamesnames = allocVector(STRSXP, 2));

	    if (xdims != R_NilValue) {
		if (ldx == 2) {/* not nrx==1 : .. fixed, ihaka 2003-09-30 */
		    SET_VECTOR_ELT(dimnames, 0, VECTOR_ELT(xdims, 1));
		    dnx = getAttrib(xdims, R_NamesSymbol);
		    if(!isNull(dnx))
			SET_STRING_ELT(dimnamesnames, 0, STRING_ELT(dnx, 1));
		}
	    }

	    YDIMS_ET_CETERA;
	}

    }
    else {					/* op == 2: tcrossprod() */

	PROTECT(ans = allocMatrix(mode, nrx, nry));
	if (mode == CPLXSXP)
	    if(sym)
		tccrossprod(COMPLEX(CAR(args)), nrx, ncx,
			    COMPLEX(CAR(args)), nry, ncy, COMPLEX(ans));
	    else
		tccrossprod(COMPLEX(CAR(args)), nrx, ncx,
			    COMPLEX(CADR(args)), nry, ncy, COMPLEX(ans));
	else {
	    if(sym)
		symtcrossprod(REAL(CAR(args)), nrx, ncx, REAL(ans));
	    else
		tcrossprod(REAL(CAR(args)), nrx, ncx,
			   REAL(CADR(args)), nry, ncy, REAL(ans));
	}

	PROTECT(xdims = getAttrib(CAR(args), R_DimNamesSymbol));
	if (sym)
	    PROTECT(ydims = xdims);
	else
	    PROTECT(ydims = getAttrib(CADR(args), R_DimNamesSymbol));

	if (xdims != R_NilValue || ydims != R_NilValue) {
	    SEXP dimnames, dimnamesnames, dnx=R_NilValue, dny=R_NilValue;

	    /* allocate dimnames and dimnamesnames */

	    PROTECT(dimnames = allocVector(VECSXP, 2));
	    PROTECT(dimnamesnames = allocVector(STRSXP, 2));

	    if (xdims != R_NilValue) {
		if (ldx == 2) {
		    SET_VECTOR_ELT(dimnames, 0, VECTOR_ELT(xdims, 0));
		    dnx = getAttrib(xdims, R_NamesSymbol);
		    if(!isNull(dnx))
			SET_STRING_ELT(dimnamesnames, 0, STRING_ELT(dnx, 0));
		}
	    }
	    if (ydims != R_NilValue) {
		if (ldy == 2) {
		    SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(ydims, 0));
		    dny = getAttrib(ydims, R_NamesSymbol);
		    if(!isNull(dny))
			SET_STRING_ELT(dimnamesnames, 1, STRING_ELT(dny, 0));
		}
	    }
	    if (VECTOR_ELT(dimnames,0) != R_NilValue ||
		VECTOR_ELT(dimnames,1) != R_NilValue) {
		if (dnx != R_NilValue || dny != R_NilValue)
		    setAttrib(dimnames, R_NamesSymbol, dimnamesnames);
		setAttrib(ans, R_DimNamesSymbol, dimnames);
	    }

	    UNPROTECT(2);
	}
    }
    UNPROTECT(3);
    return ans;
#endif
}
#undef YDIMS_ET_CETERA


