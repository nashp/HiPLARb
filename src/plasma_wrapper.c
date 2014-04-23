/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * plasma_wrapper.c is based on src/modules/lapack/Lapack.c, src/appl/bakslv.c
 * Modifications for the HiPLAR project copyright
 *    (C) 2012-2013  Vendel Szeremi
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

#include "hiplar_dbg.h"
#include "plasma_wrapper_init.h"
#include "plasma_wrapper.h"
#include "P.h"


#include <R_ext/libextern.h>
typedef struct {
    int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
    double eps, epsneg, xmin, xmax;
} AccuracyInfo;
LibExtern AccuracyInfo R_AccuracyInfo;


SEXP plasma_wrapper_zgeqrf(SEXP Ain) {
#ifdef HIPLAR_WITH_PLASMA
    int i, m, n, *Adims, info, lwork;
    Rcomplex *work, tmp;
    double *rwork;
    SEXP val, nm, jpvt, tau, rank, A;
	SEXP T;
	int T_size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_zgeqrf");
#endif


    if (!(isMatrix(Ain) && isComplex(Ain))) {
		error("'a' must be a complex matrix");
	}

    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];

	T_size = factorQR_T_size(m, n);
	T = PROTECT(allocVector(CPLXSXP, T_size));
	tmp.r = 0.0;
	tmp.i = 0.0;
    for (i = 0; i < T_size; i++) {
		COMPLEX(T)[i] = tmp;
	}

	jpvt = PROTECT(allocVector(INTSXP, n));
	for (i=0; i<n; i++) {
		INTEGER(jpvt)[i] = i+1;
	}

	info = P_zgeqrf(m, n, COMPLEX(A), COMPLEX(T));
    if (info < 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_zgeqrf");
	}

    val = PROTECT(allocVector(VECSXP, 4));
    nm = PROTECT(allocVector(STRSXP, 4));
    rank = PROTECT(ScalarInteger(m < n ? m : n));
    SET_STRING_ELT(nm, 0, mkChar("qr"));
    SET_STRING_ELT(nm, 1, mkChar("rank"));
    SET_STRING_ELT(nm, 2, mkChar("qraux"));
    SET_STRING_ELT(nm, 3, mkChar("pivot"));

    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, A);
    SET_VECTOR_ELT(val, 1, rank);
    SET_VECTOR_ELT(val, 2, T);
    SET_VECTOR_ELT(val, 3, jpvt);

    setAttrib(val, ScalarString(mkChar("usePLASMA")), ScalarLogical(1));

    UNPROTECT(6);
    return val;
#endif
}


SEXP plasma_wrapper_coef_cmplx(SEXP Q, SEXP Bin) {
#ifdef HIPLAR_WITH_PLASMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims;
    SEXP B, qr=VECTOR_ELT(Q, 0), T=VECTOR_ELT(Q, 2);
	SEXP rank=VECTOR_ELT(Q, 1);
    Rcomplex *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_coef_cmplx");
#endif


    k = INTEGER(rank)[0];

    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(Bin, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];


	info = P_zunmqr("L", "C", n, nrhs, k,
		COMPLEX(qr), n, COMPLEX(T), COMPLEX(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_zunmqr");
	}


    /* ztrtrs checks for singularity */
    for (info=0; info<k; info++) {
        if ((COMPLEX(qr)[info + info*n].r == 0.0) && (COMPLEX(qr)[info + info*n].i == 0.0)) {
		    error("error code %d from singularity check", (info+1));
        }
    }

	tmp.r = 1.0;
	tmp.i = 0.0;
	info = P_ztrsm("L", "U", "N", "N", k, nrhs, tmp,
		COMPLEX(qr), n, COMPLEX(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_ztrsm");
	}

    UNPROTECT(1);
    return B;
#endif
}


SEXP plasma_wrapper_qr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans) {
#ifdef HIPLAR_WITH_PLASMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims, tr;
    SEXP B, qr=VECTOR_ELT(Q, 0), T=VECTOR_ELT(Q, 2);
	SEXP rank=VECTOR_ELT(Q, 1);
    Rcomplex *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_qr_qy_cmplx");
#endif


/* MIN(M,N) of A */
    k = INTEGER(rank)[0];

    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    tr = asLogical(trans);
    if (tr == NA_LOGICAL) {
		error("invalid '%s' argument", "trans");
	}

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if(Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];

	info = P_zunmqr("L", tr ? "C" : "N", n, nrhs, k,
		COMPLEX(qr), n, COMPLEX(T), COMPLEX(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_zunmqr");
	}

    UNPROTECT(1);
    return B;
#endif
}


SEXP plasma_wrapper_dgeqrf(SEXP Ain) {
#ifdef HIPLAR_WITH_PLASMA
    int i, m, n, *Adims, info, lwork;
    double *work, tmp;
    SEXP val, nm, jpvt, tau, rank, A;
	SEXP T;
	int T_size;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_dgeqrf");
#endif


    if (!(isMatrix(Ain) && isReal(Ain))) {
		error("'a' must be a numeric matrix");
	}

    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];

	T_size = factorQR_T_size(m, n);
	T = PROTECT(allocVector(REALSXP, T_size));
    for (i = 0; i < T_size; i++) {
		REAL(T)[i] = 0.0;
	}

	jpvt = PROTECT(allocVector(INTSXP, n));
	for (i=0; i<n; i++) {
		INTEGER(jpvt)[i] = i+1;
	}

	info = P_dgeqrf(m, n, REAL(A), REAL(T));
    if (info < 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dgeqrf");
	}

    val = PROTECT(allocVector(VECSXP, 4));
    nm = PROTECT(allocVector(STRSXP, 4));
    rank = PROTECT(ScalarInteger(m < n ? m : n));
    SET_STRING_ELT(nm, 0, mkChar("qr"));
    SET_STRING_ELT(nm, 1, mkChar("rank"));
    SET_STRING_ELT(nm, 2, mkChar("qraux"));
    SET_STRING_ELT(nm, 3, mkChar("pivot"));

    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, A);
    SET_VECTOR_ELT(val, 1, rank);
    SET_VECTOR_ELT(val, 2, T);
    SET_VECTOR_ELT(val, 3, jpvt);

    setAttrib(val, ScalarString(mkChar("usePLASMA")), ScalarLogical(1));

    UNPROTECT(6);
    return val;
#endif
}


SEXP plasma_wrapper_coef_real(SEXP Q, SEXP Bin) {
#ifdef HIPLAR_WITH_PLASMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims;
    SEXP B, qr=VECTOR_ELT(Q, 0), T=VECTOR_ELT(Q, 2);
	SEXP rank=VECTOR_ELT(Q, 1);
    double *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_coef_real");
#endif


/* MIN(M,N) of A */
    k = INTEGER(rank)[0];

    if (!(isMatrix(Bin) && isReal(Bin))) {
		error("'b' must be a numeric matrix");
	}

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];

	info = P_dormqr("L", "T" , n, nrhs, k,
		REAL(qr), n, REAL(T), REAL(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dormqr");
	}

    /* dtrtrs checks for singularity */
    for (info=0; info<k; info++) {
        if (REAL(qr)[info + info*n] == 0.0) {
		    error("error code %d from singularity check", (info+1));
        }
    }

	info = P_dtrsm("L", "U", "N", "N", k, nrhs, 1.0,
		REAL(qr), n, REAL(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dtrsm");
	}

    UNPROTECT(1);
    return B;
#endif
}

SEXP plasma_wrapper_qr_qy_real(SEXP Q, SEXP Bin, SEXP trans) {
#ifdef HIPLAR_WITH_PLASMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims, tr;
    SEXP B, qr=VECTOR_ELT(Q, 0), T=VECTOR_ELT(Q, 2);
	SEXP rank=VECTOR_ELT(Q, 1);
    double *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_qr_qy_real");
#endif


/* MIN(M,N) of A */
    k = INTEGER(rank)[0];

    if (!(isMatrix(Bin) && isReal(Bin))) {
		error("'b' must be a numeric matrix");
	}

    tr = asLogical(trans);
    if (tr == NA_LOGICAL) {
		error("invalid '%s' argument", "trans");
	}

    PROTECT(B = duplicate(Bin));
    Qdims = INTEGER(coerceVector(getAttrib(qr, R_DimSymbol), INTSXP));
    n = Qdims[0];
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
    if (Bdims[0] != n) {
		error("right-hand side should have %d not %d rows", n, Bdims[0]);
	}
    nrhs = Bdims[1];

	info = P_dormqr("L", tr ? "T" : "N", n, nrhs, k,
		REAL(qr), n, REAL(T), REAL(B), n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dormqr");
	}

    UNPROTECT(1);
    return B;
#endif
}


SEXP plasma_wrapper_chol(SEXP A, SEXP pivot) {
#ifdef HIPLAR_WITH_PLASMA
	int piv;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_chol");
#endif


    if (isMatrix(A)) {
		SEXP ans = PROTECT((TYPEOF(A) == REALSXP)?duplicate(A):
				   coerceVector(A, REALSXP));
		SEXP adims = getAttrib(A, R_DimSymbol);
		int m = INTEGER(adims)[0];
		int n = INTEGER(adims)[1];
		int i, j;

		if (m != n) {
			error("'a' must be a square matrix");
		}
		if (m <= 0) {
			error("'a' must have dims > 0");
		}

		for (j = 0; j < n; j++) {	/* zero the lower triangle */
		    for (i = j+1; i < n; i++) {
				REAL(ans)[i + j * n] = 0.;
		    }
		}

		piv = asInteger(pivot);
		if (piv != 0 && piv != 1) error("invalid '%s' value", "pivot");
		if(!piv) {
			i = P_dpotrf("U", m, REAL(ans), m);

			if (i != 0) {
			    if (i > 0) {
					error("the leading minor of order %d is not positive definite", i);
				}
				error("error code %d from PLASMA routine '%s'", i, "PLASMA_dpotrf");
			}
		} else {
		    double tol = -1;
		    SEXP piv = PROTECT(allocVector(INTSXP, m));
		    int *ip = INTEGER(piv);
		    double *work = (double *) R_alloc(2 * (size_t)m, sizeof(double));
		    int rank, info;
		    F77_CALL(dpstrf)("U", &m, REAL(ans), &m, ip, &rank, &tol, work, &info);
		    if (info != 0) {
			if (info > 0)
			    warning("the matrix is either rank-deficient or indefinite");
			else
			    error("argument %d of Lapack routine %s had invalid value",
				  -info, "dpstrf");
		    }
		    setAttrib(ans, install("pivot"), piv);
		    setAttrib(ans, install("rank"), ScalarInteger(rank));
		    SEXP cn, dn = getAttrib(ans, R_DimNamesSymbol);
		    if (!isNull(dn) && !isNull(cn = VECTOR_ELT(dn, 1))) {
			// need to pivot the colnames
			SEXP dn2 = PROTECT(duplicate(dn));
			SEXP cn2 = VECTOR_ELT(dn2, 1);
			for(int i = 0; i < m; i++) 
			    SET_STRING_ELT(cn2, i, STRING_ELT(cn, ip[i] - 1)); // base 1
			setAttrib(ans, R_DimNamesSymbol, dn2);
			UNPROTECT(1);
		    }
		    UNPROTECT(1);
		}

		unprotect(1);
		return ans;
    } else {
		error("'a' must be a numeric matrix");
	}

    return R_NilValue; /* -Wall */
#endif
}


SEXP plasma_wrapper_chol2inv(SEXP A, SEXP size) {
#ifdef HIPLAR_WITH_PLASMA
    int sz = asInteger(size);

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_chol2inv");
#endif


    if (sz == NA_INTEGER || sz < 1) {
		error("'size' argument must be a positive integer");
		return R_NilValue; /* -Wall */
    } else {
		SEXP ans, Amat = A; /* -Wall: we initialize here as for the 1x1 case */
		int m = 1, n = 1, i, j, nprot = 0;

		if (sz == 1 && !isMatrix(A) && isReal(A)) {
		    /* nothing to do; m = n = 1; ... */
		} else if (isMatrix(A)) {
		    SEXP adims = getAttrib(A, R_DimSymbol);
		    Amat = PROTECT(coerceVector(A, REALSXP));
			nprot++;
		    m = INTEGER(adims)[0];
		    n = INTEGER(adims)[1];
		} else {
			error("'a' must be a numeric matrix");
		}

		if (sz > n) {
			UNPROTECT(nprot);
			error("'size' cannot exceed ncol(x) = %d", n);
		}
		if (sz > m) {
			UNPROTECT(nprot);
			error("'size' cannot exceed nrow(x) = %d", m);
		}

		ans = PROTECT(allocMatrix(REALSXP, sz, sz));
		nprot++;
		for (j = 0; j < sz; j++) {
		    for (i = 0; i <= j; i++) {
				REAL(ans)[i + j * sz] = REAL(Amat)[i + j * m];
			}
		}

		i = P_dpotri("U", sz, REAL(ans), sz);

		if (i != 0) {
		    UNPROTECT(nprot);
		    if (i > 0) {
				error("element (%d, %d) is zero, so the inverse cannot be computed", i, i);
			}
			error("error code %d from PLASMA routine '%s'", i, "PLASMA_dpotri");
		}

		for (j = 0; j < sz; j++) {
		    for (i = j+1; i < sz; i++) {
				REAL(ans)[i + j * sz] = REAL(ans)[j + i * sz];
			}
		}

		UNPROTECT(nprot);
		return ans;
    }
#endif
}


SEXP plasma_wrapper_zgesv(SEXP A, SEXP Bin) {
#ifdef HIPLAR_WITH_PLASMA
    int n, p, info, *ipiv, *Adims, *Bdims;
    Rcomplex *avals;
    SEXP B, Ad, Bd;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_zgesv");
#endif


    if (!(isMatrix(A) && isComplex(A))) {
		error("'a' must be a complex matrix");
	}
    if (!(isMatrix(Bin) && isComplex(Bin))) {
		error("'b' must be a complex matrix");
	}

    PROTECT(B = duplicate(Bin));
	/* This could allocate so needs protection
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
	*/
	PROTECT(Ad = coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
	PROTECT(Bd = coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
	Adims = INTEGER(Ad);
	Bdims = INTEGER(Bd);

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

    ipiv = (int *) R_alloc(n, sizeof(int));

    avals = (Rcomplex *) R_alloc(n * n, sizeof(Rcomplex));
    /* work on a copy of x */
    Memcpy(avals, COMPLEX(A), n * n);

	info = P_zgesv(n, p, avals, n, ipiv, COMPLEX(B), n);

    if (info < 0) {
		error("argument %d of PLASMA routine %s had invalid value", -info, "PLASMA_zgesv");
	}

    if (info > 0) {
		error("PLASMA routine PLASMA_zgesv: system is exactly singular");
	}

    UNPROTECT(3);
    return B;
#endif
}


SEXP plasma_wrapper_dgesv(SEXP A, SEXP Bin, SEXP tolin) {
#ifdef HIPLAR_WITH_PLASMA
    int n, p, info, *ipiv, *Adims, *Bdims;
    double *avals, anorm, rcond, tol = asReal(tolin), *work, *work2;
    SEXP B, Ad, Bd;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_dgesv");
#endif


    if (!(isMatrix(A) && isReal(A)))
		error("'a' must be a numeric matrix");
    if (!(isMatrix(Bin) && isReal(Bin)))
		error("'b' must be a numeric matrix");

    PROTECT(B = duplicate(Bin));
	/* This could allocate so needs protection
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    Bdims = INTEGER(coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
	*/
	PROTECT(Ad = coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
	PROTECT(Bd = coerceVector(getAttrib(B, R_DimSymbol), INTSXP));
	Adims = INTEGER(Ad);
	Bdims = INTEGER(Bd);

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

    ipiv = (int *) R_alloc(n, sizeof(int));

    avals = (double *) R_alloc(n * n, sizeof(double));
				/* work on a copy of A */
    Memcpy(avals, REAL(A), n * n);

	info = P_dgesv(n, p, avals, n, ipiv, REAL(B), n);

    if (info < 0) {
		error("argument %d of PLASMA routine %s had invalid value", -info, "PLASMA_dgesv");
	}
    if (info > 0) {
		error("PLASMA routine dgesv: system is exactly singular");
	}

	if (tol > 0) {
	    work2 = (double *) R_alloc(R_PLASMA_NUM_THREADS, sizeof(double));
   		anorm = P_dlange("O", n, n, REAL(A), n, work2);
    	work = (double *) R_alloc(4*n, sizeof(double));
    	F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
    	if (rcond < tol) {
            error("system is computationally singular: reciprocal condition number = %g", rcond);
		}
	}

    UNPROTECT(3);
    return B;
#endif
}


SEXP plasma_wrapper_dlange(SEXP A, SEXP type) {
#ifdef HIPLAR_WITH_PLASMA
	SEXP x, val;
	char tp;
	int m, n;
	int *xdims, nprot = 0;
	double *work;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_dlange");
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

    if ((tp == 'F') && (CHECK_VVERSION_BEQ(2,4,5))) {
		error("not implemented");
    }


    if (!isReal(A) && isNumeric(A)) {
		x = PROTECT(coerceVector(A, REALSXP));
		nprot++;
    } else {
		x = A;
    }

    if (!(isMatrix(x) && isReal(x))) {
		UNPROTECT(nprot);
		error("'A' must be a numeric matrix");
    }

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    m = xdims[0];
    n = xdims[1];

/*
	if (tp == 'I') {
		work = (double *) R_alloc(m, sizeof(double));
	} else {
		work = NULL;
	}
*/
	if (tp == 'F') {
		work = (double *) R_alloc(2*R_PLASMA_NUM_THREADS, sizeof(double));
	} else {
		work = (double *) R_alloc(R_PLASMA_NUM_THREADS, sizeof(double));
    }

    val = PROTECT(allocVector(REALSXP, 1));
	nprot++;
	REAL(val)[0] = P_dlange(&tp, m, n, REAL(x), m, work);

	UNPROTECT(nprot);
	return val;
#endif
}


void plasma_wrapper_bakslv(
double *t, int *ldt, int *n,
double *b, int *ldb, int *nb,
double *x, int *job, int *info) {
#ifdef HIPLAR_WITH_PLASMA
	char *side = "L", *uplo, *transa, *diag = "N";
	int i, ione = 1, j, nn = *n;
	double one = 1.0;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_bakslv");
#endif

	*info = 0;
	for(i = 0; i < nn; i++) {   /* check for zeros on diagonal */
		if (t[i * (*ldt + 1)] == 0.0) {
			*info = i + 1;
			return;
		}
	}

	for(j = 0; j < *nb; j++) {  /* copy b to x */
		F77_CALL(dcopy)(n, &b[j * *ldb], &ione, &x[j * *ldb], &ione);
	}

	transa = ((*job) / 10) ? "T" : "N";
	uplo = ((*job) % 10) ? "U" : "L";
	if (*n > 0 && *nb > 0 && *ldt > 0 && *ldb > 0) {
		*info = P_dtrsm(side, uplo, transa, diag, *n, *nb, one, t, *ldt, x, *ldb);
	}
#endif
}


SEXP plasma_wrapper_bakslv2(SEXP r, SEXP x, SEXP k, SEXP nb, SEXP utri, SEXP trans) {
#ifdef HIPLAR_WITH_PLASMA
	int info;
	int N, NRHS, LDA, LDB;
    int *rdims, *xdims;
	SEXP b;
	char u, t;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_bakslv2");
#endif


    N = asInteger(k);
	NRHS = asInteger(nb);

    rdims = INTEGER(coerceVector(getAttrib(r, R_DimSymbol), INTSXP));
    LDA = rdims[0];
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    LDB = rdims[0];

    PROTECT(b = duplicate(x));

	if (asInteger(utri)) {
		u = 'U';
	} else {
		u = 'L';
	}
	
	if (asInteger(trans)) {
		t = 'T';
	} else {
		t = 'N';
	}


    for (info=0; info<N; info++) {
        if (REAL(r)[info + info*LDA] == 0.0) {
		    error("error code %d from singularity check", (info+1));
        }
    }

	info = P_dtrsm("L", &u, &t, "N", N, NRHS, 1.0, REAL(r), LDA, REAL(b), LDB);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dtrsm");
	}

	UNPROTECT(1);
	return b;
#endif
}


SEXP plasma_wrapper_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v) {
#ifdef HIPLAR_WITH_PLASMA
    int *xdims, n, p, info = 0;
    double *xvals;
    SEXP val, nm;
	int ldu, ldvt;
	char ju, jvt;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_svd");
#endif


    if (!(isString(jobu) && isString(jobv))) {
		error("'jobu' and 'jobv' must be character strings");
	}

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
	p = xdims[1];

    xvals = (double *) R_alloc(n * p, sizeof(double));

    /* work on a copy of x */
    Memcpy(xvals, REAL(x), n * p);

	ldu = INTEGER(getAttrib(u, R_DimSymbol))[0];
	ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];

	ju = *(CHAR(STRING_ELT(jobu, 0)));
	jvt = *(CHAR(STRING_ELT(jobv, 0)));

    /* PLASMA 2.4.6 only supports singular values */
    if ((ju == 'N') && (jvt == 'N')) {
    	info = P_dgesvd(&ju, &jvt,
    		n, p, xvals, n, REAL(s), REAL(u), ldu, REAL(v), ldvt);
    	if (info != 0) {
    	    error("error code %d from PLASMA routine '%s'", info, "PLASMA_dgesvd");
        }
    } else {
        /* fallback to Lapack */ 
        int lwork;
        double tmp, *work;
        int *iwork= (int *) R_alloc(8*(n<p ? n : p), sizeof(int));

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_wrapper_svd fallback to LAPACK");
#endif

        /* ask for optimal size of work array */
        lwork = -1;
        F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
             &n, &p, xvals, &n, REAL(s),
             REAL(u), &ldu,
             REAL(v), &ldvt,
             &tmp, &lwork, iwork, &info); 
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "dgesdd");
        }
        lwork = (int) tmp;
        work = (double *) R_alloc(lwork, sizeof(double));
        F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
             &n, &p, xvals, &n, REAL(s),
             REAL(u), &ldu,
             REAL(v), &ldvt,
             work, &lwork, iwork, &info); 
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "dgesdd");
        }
    }


    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(2);
    return val;
#endif
}


SEXP plasma_wrapper_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v) {
#ifdef HIPLAR_WITH_PLASMA
    int *xdims, n, p, info;
    Rcomplex *work, tmp;
    SEXP x, val, nm;
	char ju, jvt;
	int ldu, ldvt;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_svd_cmplx");
#endif


    if (!(isString(jobu) && isString(jobv))) {
		error("'jobu' and 'jobv' must be character strings");
	}

    PROTECT(x = duplicate(xin));
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
	p = xdims[1];

	ldu = INTEGER(getAttrib(u, R_DimSymbol))[0];
	ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];

	ju = *(CHAR(STRING_ELT(jobu, 0)));
	jvt = *(CHAR(STRING_ELT(jobv, 0)));

    /* PLASMA 2.4.6 only supports singular values */
    if ((ju == 'N') && (jvt == 'N')) {
    	info = P_zgesvd(&ju, &jvt, n, p, COMPLEX(x), n, REAL(s),
		     COMPLEX(u), ldu, COMPLEX(v), ldvt);
        if (info != 0) {
		    error("error code %d from PLASMA routine '%s'", info, "PLASMA_zgesvd");
	    }
    } else {
        /* fallback to LAPACK */
        int lwork;
        double *rwork;
        Rcomplex *work;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_wrapper_svd_cmplx fallback to LAPACK");
#endif

        rwork = (double *) R_alloc(5*(n < p ? n:p), sizeof(double));
        /* ask for optimal size of work array */
        lwork = -1;
        F77_CALL(zgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
             &n, &p, COMPLEX(x), &n, REAL(s),
             COMPLEX(u), INTEGER(getAttrib(u, R_DimSymbol)),
             COMPLEX(v), INTEGER(getAttrib(v, R_DimSymbol)),
             &tmp, &lwork, rwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "zgesvd");
        }
        lwork = (int) tmp.r;
        work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
        F77_CALL(zgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
             &n, &p, COMPLEX(x), &n, REAL(s),
             COMPLEX(u), INTEGER(getAttrib(u, R_DimSymbol)),
             COMPLEX(v), INTEGER(getAttrib(v, R_DimSymbol)),
             work, &lwork, rwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "zgesvd");
        }
    }

    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(3);
    return val;
#endif
}


SEXP plasma_wrapper_rs(SEXP xin, SEXP only_values) {
#ifdef HIPLAR_WITH_PLASMA
    int *xdims, n, info = 0, ov;
    char jobv[1], uplo[1];
    SEXP values, ret, nm, x, z = R_NilValue;
    double *rz = NULL;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_rs");
#endif


    PROTECT(x = duplicate(xin));
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
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

    PROTECT(values = allocVector(REALSXP, n));
    if (!ov) {
		PROTECT(z = allocMatrix(REALSXP, n, n));
		rz = REAL(z);
    }

	info = P_dsyev(jobv, uplo, n, REAL(x), n, REAL(values), rz, n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_dsyev");
	}

    if (!ov) {
		ret = PROTECT(allocVector(VECSXP, 2));
		nm = PROTECT(allocVector(STRSXP, 2));
		SET_STRING_ELT(nm, 1, mkChar("vectors"));
		SET_VECTOR_ELT(ret, 1, z);
		UNPROTECT_PTR(z);
    } else {
		ret = PROTECT(allocVector(VECSXP, 1));
		nm = PROTECT(allocVector(STRSXP, 1));
    }

    SET_STRING_ELT(nm, 0, mkChar("values"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 0, values);
    UNPROTECT(4);
    return ret;
#endif
}


SEXP plasma_wrapper_rs_cmplx(SEXP xin, SEXP only_values) {
#ifdef HIPLAR_WITH_PLASMA
    int *xdims, n, lwork, info, ov;
    char jobv[1], uplo[1];
    SEXP values, ret, nm, x, z = R_NilValue;
    Rcomplex *rz = NULL;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_rs_cmplx");
#endif


    PROTECT(x = duplicate(xin));
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
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

    PROTECT(values = allocVector(REALSXP, n));
    if (!ov) {
		PROTECT(z = allocMatrix(CPLXSXP, n, n));
		rz = COMPLEX(z);
    }

	info = P_zheev(jobv, uplo, n, COMPLEX(x), n, REAL(values), rz, n);
    if (info != 0) {
		error("error code %d from PLASMA routine '%s'", info, "PLASMA_zheev");
	}

    if (!ov) {
		ret = PROTECT(allocVector(VECSXP, 2));
		nm = PROTECT(allocVector(STRSXP, 2));
		SET_STRING_ELT(nm, 1, mkChar("vectors"));
		SET_VECTOR_ELT(ret, 1, z);
		UNPROTECT_PTR(z);
    }
    else {
		ret = PROTECT(allocVector(VECSXP, 1));
		nm = PROTECT(allocVector(STRSXP, 1));
    }

    SET_STRING_ELT(nm, 0, mkChar("values"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 0, values);
    UNPROTECT(4);
    return ret;
#endif
}


static SEXP unscramble(const double* imaginary, int n,
               const double* vecs) {
    int i, j;
    SEXP s = allocMatrix(CPLXSXP, n, n);

    for (j = 0; j < n; j++) {
    if (imaginary[j] != 0) { 
        int j1 = j + 1;
        for (i = 0; i < n; i++) {
        COMPLEX(s)[i+n*j].r = COMPLEX(s)[i+n*j1].r = vecs[i + j * n]; 
        COMPLEX(s)[i+n*j1].i = -(COMPLEX(s)[i+n*j].i = vecs[i + j1 * n]);
        }       
        j = j1; 
    } else {
        for (i = 0; i < n; i++) {
        COMPLEX(s)[i+n*j].r = vecs[i + j * n]; 
        COMPLEX(s)[i+n*j].i = 0.0;
        }       
    }
    }
    return s;
}
SEXP plasma_wrapper_rg(SEXP x, SEXP only_values) {
#ifdef HIPLAR_WITH_PLASMA
    Rboolean vectors, complexValues;
    int i, n, lwork, info, *xdims, ov;
    double *work, *wR, *wI, *left, *right, *xvals, tmp;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, val;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_rg");
#endif

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_wrapper_rg fallback to Lapack");
#endif

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
    	error("'x' must be a square numeric matrix");
    }

    xvals = (double *) R_alloc(n * n, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), n * n);

    ov = asLogical(only_values);
    if (ov == NA_LOGICAL) {
        error("invalid '%s' argument", "only.values");
    }
    vectors = !ov;
    jobVL[0] = jobVR[0] = 'N';
    left = right = (double *) 0;
    if (vectors) {
    	jobVR[0] = 'V';
    	right = (double *) R_alloc(n * n, sizeof(double));
    }
    wR = (double *) R_alloc(n, sizeof(double));
    wI = (double *) R_alloc(n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, &tmp, &lwork, &info);
    if (info != 0) {
    	error("error code %d from Lapack routine '%s'", info, "dgeev");
    }

    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		    left, &n, right, &n, work, &lwork, &info);
    if (info != 0) {
    	error("error code %d from Lapack routine '%s'", info, "dgeev");
    }

    complexValues = FALSE;
    for (i = 0; i < n; i++)
	/* This test used to be !=0 for R < 2.3.0.  This is OK for 0+0i */
	if (fabs(wI[i]) >  10 * R_AccuracyInfo.eps * fabs(wR[i])) {
	    complexValues = TRUE;
	    break;
	}
    ret = PROTECT(allocVector(VECSXP, 2));
    nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_STRING_ELT(nm, 1, mkChar("vectors"));
    setAttrib(ret, R_NamesSymbol, nm);
    SET_VECTOR_ELT(ret, 1, R_NilValue);
    if (complexValues) {
    	val = allocVector(CPLXSXP, n);
    	for (i = 0; i < n; i++) {
    	    COMPLEX(val)[i].r = wR[i];
    	    COMPLEX(val)[i].i = wI[i];
    	}
    	SET_VECTOR_ELT(ret, 0, val);

    	if (vectors) {
	        SET_VECTOR_ELT(ret, 1, unscramble(wI, n, right));
        }
    } else {
    	val = allocVector(REALSXP, n);
       	for (i = 0; i < n; i++) {
       	    REAL(val)[i] = wR[i];
        }
    	SET_VECTOR_ELT(ret, 0, val);
    	if(vectors) {
    	    val = allocMatrix(REALSXP, n, n);
    	    for (i = 0; i < (n * n); i++) {
    		    REAL(val)[i] = right[i];
            }
    	    SET_VECTOR_ELT(ret, 1, val);
    	}
    }
    UNPROTECT(2);
    return ret;
#endif
}


SEXP plasma_wrapper_rg_cmplx(SEXP x, SEXP only_values) {
#ifdef HIPLAR_WITH_PLASMA
    int  n, lwork, info, *xdims, ov;
    Rcomplex *work, *left, *right, *xvals, tmp;
    double *rwork;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, values, val = R_NilValue;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_rg_cmplx");
#endif

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_wrapper_rg_cmplx fallback to Lapack");
#endif

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0];
    if (n != xdims[1]) {
	    error("'x' must be a square numeric matrix");
    }

    xvals = (Rcomplex *) R_alloc(n * n, sizeof(Rcomplex));
    /* work on a copy of x */
    Memcpy(xvals, COMPLEX(x), n * n);
    ov = asLogical(only_values);
    if (ov == NA_LOGICAL) {
        error("invalid '%s' argument", "only.values");
    }

    jobVL[0] = jobVR[0] = 'N';
    left = right = (Rcomplex *) 0;
    if (!ov) {
	jobVR[0] = 'V';
	PROTECT(val = allocMatrix(CPLXSXP, n, n));
	right = COMPLEX(val);
    }
    PROTECT(values = allocVector(CPLXSXP, n));
    rwork = (double *) R_alloc(2*n, sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(zgeev)(jobVL, jobVR, &n, xvals, &n, COMPLEX(values),
		    left, &n, right, &n, &tmp, &lwork, rwork, &info);
    if (info != 0) {
    	error("error code %d from Lapack routine '%s'", info, "zgeev");
    }

    lwork = (int) tmp.r;
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
    F77_CALL(zgeev)(jobVL, jobVR, &n, xvals, &n, COMPLEX(values),
		    left, &n, right, &n, work, &lwork, rwork, &info);
    if (info != 0) {
	    error("error code %d from Lapack routine '%s'", info, "zgeev");
    }

    if(!ov){
    	ret = PROTECT(allocVector(VECSXP, 2));
    	nm = PROTECT(allocVector(STRSXP, 2));
    	SET_STRING_ELT(nm, 1, mkChar("vectors"));
    	SET_VECTOR_ELT(ret, 1, val);
    } else {
    	ret = PROTECT(allocVector(VECSXP, 1));
    	nm = PROTECT(allocVector(STRSXP, 1));
    }
    SET_STRING_ELT(nm, 0, mkChar("values"));
    SET_VECTOR_ELT(ret, 0, values);
    setAttrib(ret, R_NamesSymbol, nm);
    UNPROTECT(ov ? 3 : 4);
    return ret;
#endif
}


SEXP plasma_wrapper_det_ge_real(SEXP Ain, SEXP logarithm) {
#ifdef HIPLAR_WITH_PLASMA
    int i, n, *Adims, info, *jpvt, sign, useLog;
    double modulus = 0.0; /* -Wall */
    SEXP val, nm, A;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_det_ge_real");
#endif

    if (!(isMatrix(Ain) && isReal(Ain))) {
    	error("'a' must be a numeric matrix");
    }

    useLog = asLogical(logarithm);
    if (useLog == NA_LOGICAL) {
        error("argument 'logarithm' must be logical");
    }

    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    n = Adims[0];
    if (Adims[1] != n) {
	    error("'a' must be a square matrix");
    }

    jpvt = (int *) R_alloc(n, sizeof(int));
    info = P_dgetrf(n, n, REAL(A), n, jpvt);

    sign = 1;
    if (info < 0) {
    	error("error code %d from PLASMA routine '%s'", info, "PLASMA_dgetrf");
    } else if (info > 0) {
        /* Singular matrix:  U[i,i] (i := info) is 0 */
    	/*warning("Lapack dgetrf(): singular matrix: U[%d,%d]=0", info,info);*/
	    modulus = (useLog ? R_NegInf : 0.);
    } else {
    	for (i = 0; i < n; i++) if (jpvt[i] != (i + 1)) {
	        sign = -sign;
        }
    	if (useLog) {
    	    modulus = 0.0;
    	    for (i = 0; i < n; i++) {
        		double dii = REAL(A)[i*(n + 1)]; /* ith diagonal element */
        		modulus += log(dii < 0 ? -dii : dii);
        		if (dii < 0) {
                    sign = -sign;
                }
    	    }
       	} else {
    	    modulus = 1.0;
    	    for (i = 0; i < n; i++)
    		modulus *= REAL(A)[i*(n + 1)];
    	    if (modulus < 0) {
    		    modulus = -modulus;
    		    sign = -sign;
	        }
	    }
    }

    val = PROTECT(allocVector(VECSXP, 2));
    nm = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(nm, 0, mkChar("modulus"));
    SET_STRING_ELT(nm, 1, mkChar("sign"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, ScalarReal(modulus));
    setAttrib(VECTOR_ELT(val, 0), install("logarithm"), ScalarLogical(useLog));
    SET_VECTOR_ELT(val, 1, ScalarInteger(sign));
    setAttrib(val, R_ClassSymbol, ScalarString(mkChar("det")));
    UNPROTECT(3);
    return val;
#endif
}


SEXP plasma_wrapper_dgecon(SEXP A, SEXP norm) {
#ifdef HIPLAR_WITH_PLASMA
    SEXP x, val;
    int *xdims, m, n, info, *iwork;
    double anorm, *work, *work2;
    char typNorm[] = {'\0', '\0'};

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_dgecon");
#endif

    if (!isString(norm)) {
    	error("'norm' must be a character string");
    }

    if (!isReal(A) && isNumeric(A)) {
    	x = PROTECT(coerceVector(A, REALSXP));
    } else {
	    x = PROTECT(duplicate(A)); /* will be overwritten by LU */
    }
    if (!(isMatrix(x) && isReal(x))) {
    	UNPROTECT(1);
    	error("'A' must be a numeric matrix");
    }

    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    m = xdims[0];
    n = xdims[1];

    typNorm[0] = La_rcond_type(CHAR(asChar(norm)));

    val = PROTECT(allocVector(REALSXP, 1));

//    work = (double *) R_alloc((*typNorm == 'I' && m > 4*n) ? m : 4*n,
//			      sizeof(double));
    work = (double *) R_alloc(4*n, sizeof(double));

	if (typNorm[0] == 'F') {
		work2 = (double *) R_alloc(2*R_PLASMA_NUM_THREADS, sizeof(double));
	} else {
		work2 = (double *) R_alloc(R_PLASMA_NUM_THREADS, sizeof(double));
    }
    /* La_rcond_type only does "O" or "I"; otherwise check for PLASMA version */
    anorm = P_dlange(typNorm, m, n, REAL(x), m, work2);

    iwork = (int *) R_alloc(m, sizeof(int));

    /* Compute the LU-decomposition and overwrite 'x' with result :*/
    info = P_dgetrf(m, n, REAL(x), m, iwork);
    if (info) {
    	if (info < 0) {
    	    UNPROTECT(2);
    	    error("error [%d] from PLASMA_dgetrf()", info);
    	}
    	else { /* i := info > 0:  LU decomp. is completed, but  U[i,i] = 0
    		* <==> singularity */
    	    REAL(val)[0] = 0.; /* rcond = 0 <==> singularity */
    	    UNPROTECT(2);
    	    return val;
    	}
    }

    F77_CALL(dgecon)(typNorm, &n, REAL(x), &n, &anorm,
		     /* rcond = */ REAL(val),
		     work, iwork, &info);

    UNPROTECT(2);
    if (info) {
        error("error [%d] from Lapack 'dgecon()'", info);
    }
    return val;
#endif
}


SEXP plasma_wrapper_zgecon(SEXP A, SEXP norm) {
#ifdef HIPLAR_WITH_PLASMA
    SEXP val;
    Rcomplex *avals; /* copy of A, to be modified */
    int *dims, n, info;
    double anorm, *rwork;
    char typNorm[] = {'\0', '\0'};

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_wrapper_zgecon");
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

    typNorm[0] = La_rcond_type(CHAR(asChar(norm)));

    val = PROTECT(allocVector(REALSXP, 1));

    rwork = (double *) R_alloc(2*n, sizeof(Rcomplex));
    anorm = F77_CALL(zlange)(typNorm, &n, &n, COMPLEX(A), &n, rwork);

    /* Compute the LU-decomposition and overwrite 'x' with result;
     * working on a copy of A : */
    avals = (Rcomplex *) R_alloc(n * n, sizeof(Rcomplex));
    Memcpy(avals, COMPLEX(A), n * n);

    info = P_zgetrf(n, n, avals, n, (int *) R_alloc(n, sizeof(int)) );

    if (info) {
    	UNPROTECT(1);
    	error("error [%d] from PLASMA_zgetrf()'", info);
    }

    F77_CALL(zgecon)(typNorm, &n, avals, &n, &anorm,
		     /* rcond = */ REAL(val),
		     /* work : */ (Rcomplex *) R_alloc(4*n, sizeof(Rcomplex)),
		     rwork, &info);
    UNPROTECT(1);
    if (info) {
        error("error [%d] from Lapack 'zgecon()'", info);
    }
    return val;
#endif
}


SEXP plasma_wrapper_dplrnt(SEXP A, SEXP seed) {
#ifdef HIPLAR_WITH_PLASMA
    int *dims;
    int info;
    SEXP val;

    dims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));

    info = P_dplrnt(dims[0], dims[1], REAL(A), dims[0], INTEGER(coerceVector(seed, INTSXP))[0]);
    if (info) {
        error("error [%d] from PLASMA_dplrnt", info);
    }

    val = PROTECT(allocVector(INTSXP, 1));
    INTEGER(val)[0] = info;
   	UNPROTECT(1);
    return val;

#endif
}


SEXP plasma_wrapper_zplrnt(SEXP A, SEXP seed) {
#ifdef HIPLAR_WITH_PLASMA
    int *dims;
    int info;
    SEXP val;

    dims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));

    info = P_zplrnt(dims[0], dims[1], COMPLEX(A), dims[0], INTEGER(coerceVector(seed, INTSXP))[0]);
    if (info) {
        error("error [%d] from PLASMA_zplrnt", info);
    }

    val = PROTECT(allocVector(INTSXP, 1));
    INTEGER(val)[0] = info;
   	UNPROTECT(1);
    return val;

#endif
}


SEXP plasma_wrapper_dplgsy(SEXP A, SEXP seed) {
#ifdef HIPLAR_WITH_PLASMA
    int *dims;
    int info;
    SEXP val;

    dims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));

    info = P_dplgsy(dims[0], dims[0], REAL(A), dims[0], INTEGER(coerceVector(seed, INTSXP))[0]);
    if (info) {
        error("error [%d] from PLASMA_dplrnt", info);
    }

    val = PROTECT(allocVector(INTSXP, 1));
    INTEGER(val)[0] = info;
   	UNPROTECT(1);
    return val;

#endif
}


