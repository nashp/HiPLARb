/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * hiplar.c is based on src/modules/lapack/Lapack.c, src/appl/bakslv.c
 * Modifications for the HiPLAR project copyright
 *    (C) 2012-2013  Vendel Szeremi, HiPLAR Team
 * HiPLAR is distributed under the same license as R.
 *
 * src/modules/lapack/Lapack.c
 *     R : A Computer Language for Statistical Data Analysis
 *     Copyright (C) 2001--2012  The R Core Team.
 *     Copyright (C) 2003--2010  The R Foundation
 *
 * src/appl/bakslv.c
 *     R : A Computer Language for Statistical Data Analysis
 *     Copyright (C) 1997-1998   Robert Gentleman, Ross Ihaka and the
 *                            R Core Team
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


#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>

#include <ctype.h> /* for toupper */

#include "hiplar_at.h"
#include "hiplar_dbg.h"

#include "plasma_wrapper.h"
#include "magma_wrapper.h"




SEXP hiplar_zgeqrf(SEXP Ain) {
    int m, n, *Adims;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_zgeqrf");
#endif


    if (!(isMatrix(Ain) && isComplex(Ain))) {
		error("'a' must be a complex matrix");
	}

    Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_zgeqrf) || (m > maxmagma_zgeqrf))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_zgeqrf(Ain);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_zgeqrf(Ain);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_coef_cmplx(SEXP Q, SEXP Bin) {
    int n, nrhs, *Bdims, *Qdims;
    SEXP qr=VECTOR_ELT(Q, 0);
    int useP = asLogical(getAttrib(Q, ScalarString(mkChar("usePLASMA"))));
    int useM = asLogical(getAttrib(Q, ScalarString(mkChar("useMAGMA"))));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_coef_cmplx");
#endif


    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];


	if (useP == 1) {
		return plasma_wrapper_coef_cmplx(Q, Bin);
	} else if (useM == 1) {
		return magma_wrapper_coef_cmplx(Q, Bin);
	} else {
        R_ShowMessage("hiplar_coef_cmplx: either PLASMA or MAGMA qr object required");
    }

    return R_NilValue;
}


SEXP hiplar_qr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans) {
    int n, nrhs, *Bdims, *Qdims, tr;
    SEXP qr=VECTOR_ELT(Q, 0);
    int useP = asLogical(getAttrib(Q, ScalarString(mkChar("usePLASMA"))));
    int useM = asLogical(getAttrib(Q, ScalarString(mkChar("useMAGMA"))));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_qr_qy_cmplx");
#endif


    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    tr = asLogical(trans);
    if (tr == NA_LOGICAL) {
		error("invalid '%s' argument", "trans");
	}

    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));
    if(Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];


	if (useP == 1) {
		return plasma_wrapper_qr_qy_cmplx(Q, Bin, trans);
	} else if (useM == 1) {
		return magma_wrapper_qr_qy_cmplx(Q, Bin, trans);
	} else {
        R_ShowMessage("hiplar_qr_qy_cmplx: either PLASMA or MAGMA qr object required");
    }

    return R_NilValue;
}


SEXP hiplar_dgeqrf(SEXP Ain) {
    int m, n, *Adims;
    
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgeqrf");
#endif


    if (!(isMatrix(Ain) && isReal(Ain))) {
		error("'a' must be a numeric matrix");
	}

    Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_dgeqrf) || (m > maxmagma_dgeqrf))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_dgeqrf(Ain);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_dgeqrf(Ain);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_coef_real(SEXP Q, SEXP Bin) {
    int n, nrhs, *Bdims, *Qdims;
    SEXP qr=VECTOR_ELT(Q, 0);
    int useP = asLogical(getAttrib(Q, ScalarString(mkChar("usePLASMA"))));
    int useM = asLogical(getAttrib(Q, ScalarString(mkChar("useMAGMA"))));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_coef_real");
#endif


    if (!(isMatrix(Bin) && isReal(Bin))) {
		error("'b' must be a numeric matrix");
	}

    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];


	if (useP == 1) {
		return plasma_wrapper_coef_real(Q, Bin);
	} else if (useM == 1) {
		return magma_wrapper_coef_real(Q, Bin);
	} else {
        R_ShowMessage("hiplar_coef_real: either PLASMA or MAGMA qr object required");
	}

    return R_NilValue;
}


SEXP hiplar_qr_qy_real(SEXP Q, SEXP Bin, SEXP trans) {
    int n, nrhs, *Bdims, *Qdims, tr;
    SEXP qr=VECTOR_ELT(Q, 0);
    int useP = asLogical(getAttrib(Q, ScalarString(mkChar("usePLASMA"))));
    int useM = asLogical(getAttrib(Q, ScalarString(mkChar("useMAGMA"))));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_qr_qy_real");
#endif


    if (!(isMatrix(Bin) && isReal(Bin))) {
		error("'b' must be a numeric matrix");
	}

    tr = asLogical(trans);
    if (tr == NA_LOGICAL) {
		error("invalid '%s' argument", "trans");
	}

    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];


	if (useP == 1) {
		return plasma_wrapper_qr_qy_real(Q, Bin, trans);
	} else if (useM == 1) {
		return magma_wrapper_qr_qy_real(Q, Bin, trans);
	} else {
        R_ShowMessage("hiplar_qr_qy_real: either PLASMA or MAGMA qr object required");
	}

    return R_NilValue;
}


SEXP hiplar_chol(SEXP A, SEXP pivot) {

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_chol");
#endif


    if (isMatrix(A)) {
		SEXP adims = getAttrib(A, R_DimSymbol);
		int m = INTEGER(adims)[0];
		int n = INTEGER(adims)[1];

		if (m != n) {
			error("'a' must be a square matrix");
		}
		if (m <= 0) {
			error("'a' must have dims > 0");
		}

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_chol) || (m > maxmagma_chol))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
			return plasma_wrapper_chol(A, pivot);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
			return magma_wrapper_chol(A, pivot);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		}
#endif

    } else {
		error("'a' must be a numeric matrix");
	}

    return R_NilValue; /* -Wall */
}


SEXP hiplar_chol2inv(SEXP A, SEXP size) {
    int sz = asInteger(size);

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_chol2inv");
#endif


    if (sz == NA_INTEGER || sz < 1) {
		error("'size' argument must be a positive integer");
		return R_NilValue; /* -Wall */
    } else {
		int m = 1, n = 1;

		if (sz == 1 && !isMatrix(A) && isReal(A)) {
		    /* nothing to do; m = n = 1; ... */
		} else if (isMatrix(A)) {
		    SEXP adims = getAttrib(A, R_DimSymbol);
		    m = INTEGER(adims)[0];
		    n = INTEGER(adims)[1];
		} else {
			error("'a' must be a numeric matrix");
		}

		if (sz > n) {
			error("'size' cannot exceed ncol(x) = %d", n);
		}
		if (sz > m) {
			error("'size' cannot exceed nrow(x) = %d", m);
		}

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_chol2inv) || (m > maxmagma_chol2inv))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
			return plasma_wrapper_chol2inv(A, size);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
			return magma_wrapper_chol2inv(A, size);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
		}
#endif

    }
}


SEXP hiplar_zgesv(SEXP A, SEXP Bin) {
    int n, p, *Adims, *Bdims;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_zgesv");
#endif


    if (!(isMatrix(A) && isComplex(A))) {
		error("'a' must be a complex matrix");
	}
    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));

    n = Adims[0];
    if (n == 0) {
		error("'a' is 0-diml");
	}

    p = Bdims[1];
    if (p == 0) {
		error("no right-hand side in 'b'");
	}

    if(Adims[1] != n) {
		error("'a' (%d x %d) must be square", n, Adims[1]);
	}

    if(Bdims[0] != n) {
		error("'b' (%d x %d) must be compatible with 'a' (%d x %d)", Bdims[0], p, n, n);
	}

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_zgesv) || (n > maxmagma_zgesv))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_zgesv(A, Bin);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_zgesv(A, Bin);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_dgesv(SEXP A, SEXP Bin, SEXP tolin) {
    int n, p, *Adims, *Bdims;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgesv");
#endif


    if (!(isMatrix(A) && isReal(A)))
		error("'a' must be a numeric matrix");
    if (!(isMatrix(Bin) && isReal(Bin)))
		error("'b' must be a numeric matrix");

    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));

    n = Adims[0];
    if (n == 0) {
		error("'a' is 0-diml");
	}

    p = Bdims[1];
    if (p == 0) {
		error("no right-hand side in 'b'");
	}

    if (Adims[1] != n) {
		error("'a' (%d x %d) must be square", n, Adims[1]);
	}

    if (Bdims[0] != n) {
		error("'b' (%d x %d) must be compatible with 'a' (%d x %d)", Bdims[0], p, n, n);
	}


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_dgesv) || (n > maxmagma_dgesv))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_dgesv(A, Bin, tolin);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_dgesv(A, Bin, tolin);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_dlange(SEXP A, SEXP type) {
	char tp;
	int m, n;
	int *xdims;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dlange");
#endif


    if (!isString(type)) {
		error("'type' must be a character string");
	}

	tp = *(CHAR(asChar(type)));
	tp = toupper(tp);
	if (tp == '1') {
		tp = 'O';
	}
	if (tp == 'E') {
		tp = 'F';
	}

    if (tp != 'M' && tp != 'O' && tp != 'I' && tp != 'F') {
		error("argument type must be one of 'M','1','O','I','F' or 'E'");
	}

    if (!(isMatrix(A) && (isReal(A) || isNumeric(A)) )) {
		error("'A' must be a numeric matrix");
    }

    xdims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = xdims[0];
    n = xdims[1];


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_dlange) || (m > maxmagma_dlange))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_dlange(A, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_dlange(A, type);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


void hiplar_bakslv(
double *t, int *ldt, int *n,
double *b, int *ldb, int *nb,
double *x, int *job, int *info)
{

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_bakslv");
#endif

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((*n < xover_bakslv) || (*n > maxmagma_bakslv))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		plasma_wrapper_bakslv(t, ldt, n, b, ldb, nb, x, job, info);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		magma_wrapper_bakslv(t, ldt, n, b, ldb, nb, x, job, info);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_bakslv2(SEXP r, SEXP x, SEXP k, SEXP nb, SEXP utri, SEXP trans) {
	int N, NRHS, LDA, LDB;
    int *rdims, *xdims;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_bakslv2");
#endif


    N = asInteger(k);
	NRHS = asInteger(nb);

    rdims = INTEGER(coerceVector(getAttrib(r, R_DimSymbol), INTSXP));
    LDA = rdims[0];
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    LDB = rdims[0];


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((N < xover_bakslv) || (N > maxmagma_bakslv))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_bakslv2(r, x, k, nb, utri, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_bakslv2(r, x, k, nb, utri, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v) {
    int *xdims, n, p;
	int ldu, ldvt;
	char ju, jvt;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_svd");
#endif


    if (!(isString(jobu) && isString(jobv))) {
		error("'jobu' and 'jobv' must be character strings");
	}

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
	p = xdims[1];

	ldu = INTEGER(getAttrib(u, R_DimSymbol))[0];
	ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];

	ju = *(CHAR(STRING_ELT(jobu, 0)));
	jvt = *(CHAR(STRING_ELT(jobv, 0)));


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_svd) || (n > maxmagma_svd))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_svd(jobu, jobv, x, s, u, v);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_svd(jobu, jobv, x, s, u, v);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v) {
    int *xdims, n, p;
	char ju, jvt;
	int ldu, ldvt;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_svd_cmplx");
#endif


    if (!(isString(jobu) && isString(jobv))) {
		error("'jobu' and 'jobv' must be character strings");
	}

    xdims = INTEGER(coerceVector(getAttrib(xin, R_DimSymbol), INTSXP));
    n = xdims[0];
	p = xdims[1];

	ldu = INTEGER(getAttrib(u, R_DimSymbol))[0];
	ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];

	ju = *(CHAR(STRING_ELT(jobu, 0)));
	jvt = *(CHAR(STRING_ELT(jobv, 0)));


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_svd_cmplx) || (n > maxmagma_svd_cmplx))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_svd_cmplx(jobu, jobv, xin, s, u, v);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_svd_cmplx(jobu, jobv, xin, s, u, v);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_rs(SEXP xin, SEXP only_values) {
    int *xdims, n, ov;
    char jobv[1], uplo[1];

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_rs");
#endif


    xdims = INTEGER(coerceVector(getAttrib(xin, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
		error("'x' must be a square numeric matrix");
	}

    ov = asLogical(only_values);
    if (ov == NA_LOGICAL) {
		error("invalid '%s' argument", "only.values");
	}

    if (ov) {
		jobv[0] = 'N';
	} else {
		jobv[0] = 'V';
	}


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_rs) || (n > maxmagma_rs))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_rs(xin, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_rs(xin, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_rs_cmplx(SEXP xin, SEXP only_values) {
    int *xdims, n, ov;
    char jobv[1], uplo[1];

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_rs_cmplx");
#endif


    xdims = INTEGER(coerceVector(getAttrib(xin, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
		error("'x' must be a square numeric matrix");
	}

    ov = asLogical(only_values);
    if (ov == NA_LOGICAL) {
		error("invalid '%s' argument", "only.values");
	}
    if (ov) {
		jobv[0] = 'N';
	} else {
		jobv[0] = 'V';
	}

    uplo[0] = 'L';


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_rs_cmplx) || (n > maxmagma_rs_cmplx))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_rs_cmplx(xin, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_rs_cmplx(xin, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_rg(SEXP x, SEXP only_values) {
    int n, *xdims, ov;
    char jobVL[1], jobVR[1];

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_rg");
#endif

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
    	error("'x' must be a square numeric matrix");
    }


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_rg) || (n > maxmagma_rg))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_rg(x, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_rg(x, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_rg_cmplx(SEXP x, SEXP only_values) {
    int  n, *xdims, ov;
    char jobVL[1], jobVR[1];

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_rg_cmplx");
#endif

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
	    error("'x' must be a square numeric matrix");
    }


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_rg_cmplx) || (n > maxmagma_rg_cmplx))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_rg_cmplx(x, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_rg_cmplx(x, only_values);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_det_ge_real(SEXP Ain, SEXP logarithm) {
    int n, *Adims, useLog;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_det_ge_real");
#endif

    if (!(isMatrix(Ain) && isReal(Ain))) {
    	error("'a' must be a numeric matrix");
    }

    useLog = asLogical(logarithm);
    if (useLog == NA_LOGICAL) {
        error("argument 'logarithm' must be logical");
    }

    Adims = INTEGER(coerceVector(getAttrib(Ain, R_DimSymbol), INTSXP));
    n = Adims[0];
    if (Adims[1] != n) {
	    error("'a' must be a square matrix");
    }


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_det_ge_real) || (n > maxmagma_det_ge_real))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_det_ge_real(Ain, logarithm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_det_ge_real(Ain, logarithm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_dgecon(SEXP A, SEXP norm) {
    int *xdims, m, n;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dgecon");
#endif

    if (!isString(norm)) {
    	error("'norm' must be a character string");
    }

    if (!(isMatrix(A) && (isReal(A) || isNumeric(A)) )) {
    	error("'A' must be a numeric matrix");
    }

    xdims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = xdims[0];
    n = xdims[1];


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((m < xover_dgecon) || (m > maxmagma_dgecon))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_dgecon(A, norm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_dgecon(A, norm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


SEXP hiplar_zgecon(SEXP A, SEXP norm) {
    int *dims, n;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_zgecon");
#endif

    if (!isString(norm)) {
	    error("'norm' must be a character string");
    }

    if (!(isMatrix(A) && isComplex(A))) {
	    error("'A' must be a complex matrix");
    }

    dims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    n = dims[0];
    if(n != dims[1]) {
	    error("'A' must be a *square* matrix");
    }


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	if ((hiplar_library == HIPLAR_USE_PLASMA) || ((hiplar_library == HIPLAR_USE_AUTO) && ((n < xover_zgecon) || (n > maxmagma_zgecon))) ) {
#endif
#ifdef HIPLAR_WITH_PLASMA
		return plasma_wrapper_zgecon(A, norm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	} else {
#endif
#ifdef HIPLAR_WITH_MAGMA
		return magma_wrapper_zgecon(A, norm);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
	}
#endif

}


