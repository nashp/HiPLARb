/*
 * HiPLAR - High Performance Linear Algebra in R
 *
 * magma_wrapper.c is based on src/modules/lapack/Lapack.c, src/appl/bakslv.c
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

#ifdef HIPLAR_WITH_MAGMA
#include <magma.h>
#endif
#include "hiplar_at.h"
#include "hiplar_dbg.h"


#include <R_ext/libextern.h>
typedef struct {
    int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
    double eps, epsneg, xmin, xmax;
} AccuracyInfo;
LibExtern AccuracyInfo R_AccuracyInfo;


SEXP magma_wrapper_zgeqrf(SEXP Ain) {
#ifdef HIPLAR_WITH_MAGMA
    int i, m, n, *Adims, info, lwork;
    Rcomplex *work, tmp;
    double *rwork;
    SEXP val, nm, jpvt, tau, rank, A;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_zgeqrf");
#endif


    if (!(isMatrix(Ain) && isComplex(Ain))) {
		error("'a' must be a complex matrix");
	}

    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];

	jpvt = PROTECT(allocVector(INTSXP, n));
	for (i=0; i<n; i++) {
		INTEGER(jpvt)[i] = i+1;
	}

    tau = PROTECT(allocVector(CPLXSXP, m < n ? m : n));

    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgeqrf cpu interface");
#endif

        lwork = -1;
        magma_zgeqrf(m, n, (cuDoubleComplex *)COMPLEX(A), m,
            (cuDoubleComplex *)COMPLEX(tau), (cuDoubleComplex *)&tmp, lwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "magma_zgeqrf");
        }
        lwork = (int) tmp.r;
        work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
        magma_zgeqrf(m, n, (cuDoubleComplex *)COMPLEX(A), m,
            (cuDoubleComplex *)COMPLEX(tau), (cuDoubleComplex *)work, lwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "magma_zgeqrf");
        }

	} else {
		/* MAGMA gpu interface */
        cuDoubleComplex *d_A;
        int ldda  = ((m+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgeqrf gpu interface");
#endif

        magma_malloc((void**)&d_A, sizeof(cuDoubleComplex) * ldda * n);
        magma_zsetmatrix(m, n, (cuDoubleComplex *)COMPLEX(A), m, d_A, ldda);

        magma_zgeqrf2_gpu(m, n, d_A, ldda, (cuDoubleComplex *)COMPLEX(tau), &info);
        if (info < 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_zgeqrf2_gpu");
        }

        magma_zgetmatrix(m, n, d_A, ldda, (cuDoubleComplex *)COMPLEX(A), m);
        magma_free(d_A);

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
    SET_VECTOR_ELT(val, 2, tau);
    SET_VECTOR_ELT(val, 3, jpvt);

    setAttrib(val, ScalarString(mkChar("useMAGMA")), ScalarLogical(1));

    UNPROTECT(6);
    return val;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_coef_cmplx(SEXP Q, SEXP Bin) {
#ifdef HIPLAR_WITH_MAGMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims;
    SEXP B, qr=VECTOR_ELT(Q, 0), tau=VECTOR_ELT(Q, 2);
    Rcomplex *work, tmp;
    cuDoubleComplex *d_A, *d_C;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_coef_cmplx");
#endif


    k = LENGTH(tau);

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

    lwork = -1;

#if(0)
    magma_zunmqr('L', 'C', n, nrhs, k,
        (cuDoubleComplex *)COMPLEX(qr), n, (cuDoubleComplex *)COMPLEX(tau),
        (cuDoubleComplex *)COMPLEX(B), n, (cuDoubleComplex *)&tmp, lwork, &info );
#else
    F77_CALL(zunmqr)("L", "C", &n, &nrhs, &k,
             COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n,
             &tmp, &lwork, &info);
#endif

    if (info != 0) {
        error("error code %d from Lapack routine '%s'", info, "zunmqr");
    }
    lwork = (int) tmp.r;
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));

#if(0)
    magma_zunmqr('L', 'C', n, nrhs, k,
        (cuDoubleComplex *)COMPLEX(qr), n, (cuDoubleComplex *)COMPLEX(tau),
        (cuDoubleComplex *)COMPLEX(B), n, (cuDoubleComplex *)&tmp, lwork, &info );
#else
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_coef_cmplx fallback to Lapack");
#endif
    F77_CALL(zunmqr)("L", "C", &n, &nrhs, &k,
             COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n,
             work, &lwork, &info);
#endif
    if (info != 0) {
        error("error code %d from Lapack routine '%s'", info, "zunmqr");
    }

    /* ztrtrs checks for singularity */
    for (info=0; info<k; info++) {
        if ((COMPLEX(qr)[info + info*n].r == 0.0) && (COMPLEX(qr)[info + info*n].i == 0.0)) {
		    error("error code %d from singularity check", (info+1));
        }
    }

    magma_malloc((void**)&d_A, sizeof(cuDoubleComplex) * n * k);
    magma_malloc((void**)&d_C, sizeof(cuDoubleComplex) * n * nrhs);

    magma_zsetmatrix(n, k, (cuDoubleComplex *)COMPLEX(qr), n, d_A, n);
    magma_zsetmatrix(n, nrhs, (cuDoubleComplex *)COMPLEX(B), n, d_C, n);
    magma_ztrsm('L', 'U', 'N', 'N', k, nrhs, MAGMA_Z_ONE, d_A, n, d_C, n);
    magma_zgetmatrix(n, nrhs, d_C, n, (cuDoubleComplex *)COMPLEX(B), n);

    magma_free(d_A);
    magma_free(d_C);

    UNPROTECT(1);
    return B;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_qr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans) {
#ifdef HIPLAR_WITH_MAGMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims, tr;
    SEXP B, qr=VECTOR_ELT(Q, 0), tau=VECTOR_ELT(Q, 2);
    Rcomplex *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_qr_qy_cmplx");
#endif


    k = LENGTH(tau);

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

    lwork = -1;
#if(0)
    magma_zunmqr('L', tr ? 'C' : 'N', n, nrhs, k,
        (cuDoubleComplex *)COMPLEX(qr), n, (cuDoubleComplex *)COMPLEX(tau),
        (cuDoubleComplex *)COMPLEX(B), n, (cuDoubleComplex *)&tmp, lwork, &info );
#else
    F77_CALL(zunmqr)("L", tr ? "C" : "N", &n, &nrhs, &k,
             COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n,
             &tmp, &lwork, &info);
#endif
    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_zunmqr");
    }
    lwork = (int) tmp.r;
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
#if(0)
    magma_zunmqr('L', tr ? 'C' : 'N', n, nrhs, k,
        (cuDoubleComplex *)COMPLEX(qr), n, (cuDoubleComplex *)COMPLEX(tau),
        (cuDoubleComplex *)COMPLEX(B), n, (cuDoubleComplex *)work, lwork, &info );
#else
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_qr_qy_cmplx fallback to Lapack");
#endif
    F77_CALL(zunmqr)("L", tr ? "C" : "N", &n, &nrhs, &k,
             COMPLEX(qr), &n, COMPLEX(tau), COMPLEX(B), &n,
             work, &lwork, &info);
#endif
    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "zunmqr");
    }

    UNPROTECT(1);
    return B;

#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_dgeqrf(SEXP Ain) {
#ifdef HIPLAR_WITH_MAGMA
    int i, m, n, *Adims, info, lwork;
    double *work, tmp;
    SEXP val, nm, jpvt, tau, rank, A;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_dgeqrf");
#endif


    if (!(isMatrix(Ain) && isReal(Ain))) {
		error("'a' must be a numeric matrix");
	}

    PROTECT(A = duplicate(Ain));
    Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    m = Adims[0];
    n = Adims[1];

	jpvt = PROTECT(allocVector(INTSXP, n));
	for (i=0; i<n; i++) {
		INTEGER(jpvt)[i] = i+1;
	}
    tau = PROTECT(allocVector(REALSXP, m < n ? m : n)); 


    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgeqrf cpu interface");
#endif

    lwork = -1;
    magma_dgeqrf(m, n, REAL(A), m, REAL(tau), &tmp, lwork, &info);
    if (info < 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_dgeqrf");
    }
    lwork = (int) tmp; 
    work = (double *) R_alloc(lwork, sizeof(double));

    magma_dgeqrf(m, n, REAL(A), m, REAL(tau), work, lwork, &info);
    if (info < 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_dgeqrf");
    }

	} else {
		/* MAGMA gpu interface */
        double *d_A;
        int ldda  = ((m+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgeqrf gpu interface");
#endif

        magma_malloc((void**)&d_A, sizeof(double) * ldda * n);
        magma_dsetmatrix(m, n, REAL(A), m, d_A, ldda);

        magma_dgeqrf2_gpu(m, n, d_A, ldda, REAL(tau), &info);
        if (info < 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_dgeqrf2_gpu");
        }

        magma_dgetmatrix(m, n, d_A, ldda, REAL(A), m);
        magma_free(d_A);
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
    SET_VECTOR_ELT(val, 2, tau);
    SET_VECTOR_ELT(val, 3, jpvt);

    setAttrib(val, ScalarString(mkChar("useMAGMA")), ScalarLogical(1));

    UNPROTECT(6);
    return val;
#endif
}


SEXP magma_wrapper_coef_real(SEXP Q, SEXP Bin) {
#ifdef HIPLAR_WITH_MAGMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims;
    SEXP B, qr=VECTOR_ELT(Q, 0), tau=VECTOR_ELT(Q, 2);
    double *work, tmp;
    double *d_A, *d_C;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_coef_real");
#endif

    k = LENGTH(tau);

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

    lwork = -1;
#if(0)
    magma_dormqr('L', 'T', n, nrhs, k,
        REAL(qr), n, REAL(tau), REAL(B), n, &tmp, lwork, &info);
#else
    F77_CALL(dormqr)("L", "T", &n, &nrhs, &k,
        REAL(qr), &n, REAL(tau), REAL(B), &n, &tmp, &lwork, &info);
#endif

    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_dormqr");
    }
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

#if(0)
    magma_dormqr('L', 'T', n, nrhs, k,
        REAL(qr), n, REAL(tau), REAL(B), n, work, lwork, &info);
#else
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_coef_real fallback to Lapack");
#endif
    F77_CALL(dormqr)("L", "T", &n, &nrhs, &k,
        REAL(qr), &n, REAL(tau), REAL(B), &n, work, &lwork, &info);
#endif

    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_dormqr");
    }

    /* dtrtrs checks for singularity */
    for (info=0; info<k; info++) {
        if (REAL(qr)[info + info*n] == 0.0) {
		    error("error code %d from singularity check", (info+1));
        }
    }

    magma_malloc((void**)&d_A, sizeof(double) * n * k);
    magma_malloc((void**)&d_C, sizeof(double) * n * nrhs);

    magma_dsetmatrix(n, k, REAL(qr), n, d_A, n);
    magma_dsetmatrix(n, nrhs, REAL(B), n, d_C, n);
    magmablas_dtrsm('L', 'U', 'N', 'N', k, nrhs, 1.0, d_A, n, d_C, n);
    magma_dgetmatrix(n, nrhs, d_C, n, REAL(B), n);

    magma_free(d_A);
    magma_free(d_C);

    UNPROTECT(1);
    return B;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_qr_qy_real(SEXP Q, SEXP Bin, SEXP trans) {
#ifdef HIPLAR_WITH_MAGMA
    int n, nrhs, lwork, info, k, *Bdims, *Qdims, tr;
    SEXP B, qr=VECTOR_ELT(Q, 0), tau=VECTOR_ELT(Q, 2);
	SEXP rank=VECTOR_ELT(Q, 1);
    double *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_qr_qy_real");
#endif


    k = LENGTH(tau);

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

    lwork = -1;
#if(0)
    magma_dormqr('L',  tr ? 'T' : 'N', n, nrhs, k,
        REAL(qr), n, REAL(tau), REAL(B), n, &tmp, lwork, &info);
#else
    F77_CALL(dormqr)("L", tr ? "T" : "N", &n, &nrhs, &k,
        REAL(qr), &n, REAL(tau), REAL(B), &n, &tmp, &lwork, &info);
#endif
    if (info != 0) {
        error("error code %d from Lapack routine '%s'", info, "magma_dormqr");
    }
    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

#if(0)
    magma_dormqr('L',  tr ? 'T' : 'N', n, nrhs, k,
        REAL(qr), n, REAL(tau), REAL(B), n, work, lwork, &info);
#else
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_qr_qy_real fallback to Lapack");
#endif
    F77_CALL(dormqr)("L", tr ? "T" : "N", &n, &nrhs, &k,
        REAL(qr), &n, REAL(tau), REAL(B), &n, work, &lwork, &info);
#endif
    if (info != 0) {
        error("error code %d from Lapack routine '%s'", info, "magma_dormqr");
    }

    UNPROTECT(1);
    return B;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_chol(SEXP A, SEXP pivot) {
#ifdef HIPLAR_WITH_MAGMA
	int piv;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_chol");
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

            if (magma_interface == MAGMA_CPU_INTERFACE) {
				/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_chol cpu interface");
#endif
				magma_dpotrf('U', m, REAL(ans), m, &i);
			} else {
				/* MAGMA gpu interface */
                double *h_A, *h_R, *d_A;
                int ldda = ((n+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_chol gpu interface");
#endif

                //magma_malloc_cpu((void**)&h_A, sizeof(double) * m * m);
                //magma_malloc_pinned((void**)&h_R, sizeof(double) * m * m);
                h_A = REAL(ans);
                h_R = REAL(ans);
                magma_malloc((void**)&d_A, sizeof(double) * ldda * m);

                magma_dsetmatrix(m, m, h_A, m, d_A, ldda);
                magma_dpotrf_gpu('U', m, d_A, ldda, &i);
                magma_dgetmatrix(m, m, d_A, ldda, h_R, m);

                magma_free(d_A);
			}

			if (i != 0) {
			    if (i > 0) {
					error("the leading minor of order %d is not positive definite", i);
				}
				error("error code %d from MAGMA routine '%s'", i, "magma_dpotrf");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_chol2inv(SEXP A, SEXP size) {
#ifdef HIPLAR_WITH_MAGMA
    int sz = asInteger(size);

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_chol2inv");
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

        if (magma_interface == MAGMA_CPU_INTERFACE) {
	    	/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_chol2inv cpu interface");
#endif
	    	magma_dpotri('U', sz, REAL(ans), sz, &i);
	    } else {
	    	/* MAGMA gpu interface */
            double *h_A, *h_R, *d_A;
            int ldda = ((sz+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_chol2inv gpu interface");
#endif

            //magma_malloc_cpu((void**)&h_A, sizeof(double) * sz * sz);
            //magma_malloc_pinned((void**)&h_R, sizeof(double) * sz * sz);
            h_A = REAL(ans);
            h_R = REAL(ans);
            magma_malloc((void**)&d_A, sizeof(double) * ldda * sz);

            magma_dsetmatrix(sz, sz, h_A, sz, d_A, ldda);
            magma_dpotri_gpu('U', sz, d_A, ldda, &i);
            magma_dgetmatrix(sz, sz, d_A, ldda, h_R, sz);

            magma_free(d_A);

	    }

		if (i != 0) {
		    UNPROTECT(nprot);
		    if (i > 0) {
				error("element (%d, %d) is zero, so the inverse cannot be computed", i, i);
			}
			error("error code %d from MAGMA routine '%s'", i, "magma_dpotri");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_zgesv(SEXP A, SEXP Bin) {
#ifdef HIPLAR_WITH_MAGMA
    int n, p, info, *ipiv, *Adims, *Bdims;
    Rcomplex *avals;
    SEXP B, Ad, Bd;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_zgesv");
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

    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgesv cpu interface");
#endif
		magma_zgesv(n, p, (cuDoubleComplex *)avals, n, ipiv, (cuDoubleComplex *)COMPLEX(B), n, &info);
        if (info < 0) {
    		error("argument %d of MAGMA routine %s had invalid value", -info, "magma_zgesv");
    	}
	} else {
		/* MAGMA gpu interface */
        cuDoubleComplex *h_A, *h_B, *h_X;
        cuDoubleComplex *d_A, *d_B;
        int ldda = ((n+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgesv gpu interface");
#endif
    
        //magma_malloc_cpu((void**)h_A, sizeof(cuDoubleComplex) * n * n);  
        //magma_malloc_cpu((void**)h_B, sizeof(cuDoubleComplex) * n * p);
        //magma_malloc_cpu((void**)h_X, sizeof(cuDoubleComplex) * n * p);
        //magma_malloc_cpu((void**)ipiv, sizeof(int) * n);  
        h_A = (cuDoubleComplex *)avals;
        h_B = (cuDoubleComplex *)COMPLEX(B);
        h_X = (cuDoubleComplex *)COMPLEX(B);

        magma_malloc((void**)&d_A, sizeof(cuDoubleComplex) * ldda * n);  
        magma_malloc((void**)&d_B, sizeof(cuDoubleComplex) * ldda * p);

        magma_zsetmatrix(n, n, h_A, n, d_A, ldda);
        magma_zsetmatrix(n, p, h_B, n, d_B, ldda);
        magma_zgesv_gpu(n, p, d_A, ldda, ipiv, d_B, ldda, &info);
        if (info < 0) {
    		error("argument %d of MAGMA routine %s had invalid value", -info, "magma_zgesv_gpu");
    	}
        magma_zgetmatrix(n, p, d_B, ldda, h_X, n);

        magma_free(d_A);
        magma_free(d_B);
	}

    if (info > 0) {
		error("MAGMA routine magma_zgesv: system is exactly singular");
	}

    UNPROTECT(3);
    return B;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_dgesv(SEXP A, SEXP Bin, SEXP tolin) {
#ifdef HIPLAR_WITH_MAGMA
    int n, p, info, *ipiv, *Adims, *Bdims;
    double *avals, anorm, rcond, tol = asReal(tolin), *work, *work2;
    SEXP B, Ad, Bd;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_dgesv");
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

    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgesv cpu interface");
#endif
		magma_dgesv(n, p, avals, n, ipiv, REAL(B), n, &info);
        if (info < 0) {
		    error("argument %d of MAGMA routine %s had invalid value", -info, "magma_dgesv");
    	}
	} else {
		/* MAGMA gpu interface */
        double *h_A, *h_B, *h_X;
        double *d_A, *d_B;
        int ldda = ((n+31)/32)*32;
    
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgesv gpu interface");
#endif

        //magma_malloc_cpu((void**)&h_A, sizeof(double) * n * n);  
        //magma_malloc_cpu((void**)&h_B, sizeof(double) * n * p);
        //magma_malloc_cpu((void**)&h_X, sizeof(double) * n * p);
        //magma_malloc_cpu((void**)&ipiv, sizeof(int) * n);  
        h_A = avals;
        h_B = REAL(B);
        h_X = REAL(B);

        magma_malloc((void**)&d_A, sizeof(double) * ldda * n);  
        magma_malloc((void**)&d_B, sizeof(double) * ldda * p);

        magma_dsetmatrix(n, n, h_A, n, d_A, ldda);
        magma_dsetmatrix(n, p, h_B, n, d_B, ldda);
        magma_dgesv_gpu(n, p, d_A, ldda, ipiv, d_B, ldda, &info);
        if (info < 0) {
		    error("argument %d of MAGMA routine %s had invalid value", -info, "magma_dgesv_gpu");
    	}
        magma_dgetmatrix(n, p, d_B, ldda, h_X, n);

        magma_free(d_A);
        magma_free(d_B);
	}

    if (info < 0) {
		error("argument %d of MAGMA routine %s had invalid value", -info, "magma_dgesv");
	}
    if (info > 0) {
		error("MAGMA routine dgesv: system is exactly singular");
	}

	if (tol > 0) {
        anorm = F77_CALL(dlange)("1", &n, &n, REAL(A), &n, (double*) NULL);
    	work = (double *) R_alloc(4*n, sizeof(double));
    	F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
    	if (rcond < tol) {
            error("system is computationally singular: reciprocal condition number = %g", rcond);
		}
	}

    UNPROTECT(3);
    return B;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_dlange(SEXP A, SEXP type) {
#ifdef HIPLAR_WITH_MAGMA
	SEXP x, val;
	char tp;
	int m, n;
	int *xdims, nprot = 0;
	double *work;
    double *dA, *dW;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_dlange");
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

    val = PROTECT(allocVector(REALSXP, 1));
	nprot++;

#if(1)
	if (tp == 'I') {
		work = (double *) R_alloc(m, sizeof(double));
	} else {
		work = NULL;
	}

    REAL(val)[0] = F77_CALL(dlange)(&tp, &m, &n, REAL(x), &m, work);

#else
    /* MAGAM 1.3.0 only supports Inf norm and M & N must be mult. of 64 */

    if (MAGMA_SUCCESS != magma_dmalloc(&dA, m*n)) {
        error("magma_wrapper_dlange: not enough device memory\n");
        return R_NilValue;
    }
    magma_dsetmatrix(m, n, REAL(x), m, dA, m);

    dW = NULL;
    if (tp == 'I') {
        if (MAGMA_SUCCESS != magma_dmalloc(&dW, m)) {
            magma_free(dA);
            error("magma_wrapper_dlange: not enough device memory\n");
            return R_NilValue;
        }
    }

    REAL(val)[0] = magmablas_dlange(tp, m, n, dA, m, dW);
    magma_free(dA);
    magma_free(dW);
#endif

	UNPROTECT(nprot);
	return val;
#endif
	return R_NilValue; /* -Wall */
}


void magma_wrapper_bakslv(
double *t, int *ldt, int *n,
double *b, int *ldb, int *nb,
double *x, int *job, int *info) {
#ifdef HIPLAR_WITH_MAGMA
	char *side = "L", *uplo, *transa, *diag = "N";
	int i, ione = 1, j, nn = *n;
	double one = 1.0;
    double *dA, *dB;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_bakslv");
#endif

	*info = 0;
	for(i = 0; i < nn; i++) {   /* check for zeros on diagonal */
		if (t[i * (*ldt + 1)] == 0.0) {
			*info = i + 1;
			return;
		}
	}

#if(0)
// copy from b to dev, copy back from dev to x
	for(j = 0; j < *nb; j++) {  /* copy b to x */
		F77_CALL(dcopy)(n, &b[j * *ldb], &ione, &x[j * *ldb], &ione);
	}
#endif

	transa = ((*job) / 10) ? "T" : "N";
	uplo = ((*job) % 10) ? "U" : "L";
	if (*n > 0 && *nb > 0 && *ldt > 0 && *ldb > 0) {
        if (MAGMA_SUCCESS != magma_dmalloc(&dA, (*n)*(*ldt))) {
            error("magma_wrapper_bakslv: not enough device memory\n");
            return;
        }
        magma_dsetmatrix(*ldt, *n, t, *ldt, dA, *ldt);

        if (MAGMA_SUCCESS != magma_dmalloc(&dB, (*nb)*(*ldb) )) {
            magma_free(dA);
            error("magma_wrapper_bakslv: not enough device memory\n");
            return;
        }
        magma_dsetmatrix(*ldb, *nb, b, *ldb, dB, *ldb);

		magmablas_dtrsm(*side, *uplo, *transa, *diag, *n, *nb, one, dA, *ldt, dB, *ldb);
        magma_dgetmatrix(*ldb, *nb, dB, *ldb, x, *ldb);
        magma_free(dA);
        magma_free(dB);
	}

#endif
}


SEXP magma_wrapper_bakslv2(SEXP r, SEXP x, SEXP k, SEXP nb, SEXP utri, SEXP trans) {
#ifdef HIPLAR_WITH_MAGMA
	int info;
	int N, NRHS, LDA, LDB;
    int *rdims, *xdims;
	SEXP b;
	char u, t;
    double *dA, *dB;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_bakslv2");
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


    if (MAGMA_SUCCESS != magma_dmalloc( &dA, N*LDA )) {
        error("magma_wrapper_bakslv2: not enough device memory\n");
        return R_NilValue;
    }
    magma_dsetmatrix( LDA, N, REAL(r), LDA, dA, LDA);

    if (MAGMA_SUCCESS != magma_dmalloc( &dB, NRHS*LDB )) {
        magma_free(dA);
        error("magma_wrapper_bakslv2: not enough device memory\n");
        return R_NilValue;
    }
    magma_dsetmatrix( LDB, NRHS, REAL(b), LDB, dB, LDB);

    magmablas_dtrsm('L', u, t, 'N', N, NRHS, 1.0, dA, LDA, dB, LDB);
    magma_dgetmatrix( LDB, NRHS, dB, LDB, REAL(b), LDB);

    magma_free(dA);
    magma_free(dB);

	UNPROTECT(1);
	return b;
#endif
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v) {
#ifdef HIPLAR_WITH_MAGMA
    int *xdims, n, p, lwork, info = 0;
    double *xvals, *work, tmp;
    SEXP val, nm;
	int ldu, ldvt;
	char ju, jvt;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_svd");
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

    /* La.svd() uses dgesdd and jobz: adjust ju/jvt */
    if ((jvt != 'A') || (jvt != 'S') || (jvt != 'O') || (jvt != 'N')) {
        jvt = 'N';
        if ((ju == 'A') || (ju == 'S')) {
            jvt = ju;
        }
        if (ju == 'O') {
            if (n >= p) {
                ju  = 'O';
                jvt = 'A';
            } else {
                ju  = 'A';
                jvt = 'O';
            }
        }
    }


    lwork = -1;
    magma_dgesvd(ju, ju, n, p, xvals, n, REAL(s), REAL(u), ldu, REAL(v), ldvt, &tmp, lwork, &info);
    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_dgesvd");
    }
    lwork = (int)tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

    magma_dgesvd(ju, ju, n, p, xvals, n, REAL(s), REAL(u), ldu, REAL(v), ldvt, work, lwork, &info);

	if (info != 0) {
	    error("error code %d from MAGMA routine '%s'", info, "magma_dgesvd");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v) {
#ifdef HIPLAR_WITH_MAGMA
    int *xdims, n, p, info;
    Rcomplex *work, tmp;
    SEXP x, val, nm;
	char ju, jvt;
	int ldu, ldvt;
    int lwork;
    double *rwork;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_svd_cmplx");
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

    rwork = (double *) R_alloc(5*(n < p ? n:p), sizeof(double));
    /* ask for optimal size of work array */
    lwork = -1;
    magma_zgesvd(ju, jvt, n, p, (cuDoubleComplex *)COMPLEX(x), n, REAL(s),
        (cuDoubleComplex *)COMPLEX(u), ldu, (cuDoubleComplex *)COMPLEX(v),
        ldvt, (cuDoubleComplex *)&tmp, lwork, rwork, &info);
    if (info != 0) {
        error("error code %d from MAGMA routine '%s'", info, "magma_zgesvd");
    }
    lwork = (int) tmp.r;
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));

    magma_zgesvd(ju, jvt, n, p, (cuDoubleComplex *)COMPLEX(x), n, REAL(s),
        (cuDoubleComplex *)COMPLEX(u), ldu, (cuDoubleComplex *)COMPLEX(v),
        ldvt, (cuDoubleComplex *)work, lwork, rwork, &info);

    if (info != 0) {
		error("error code %d from MAGMA routine '%s'", info, "magma_zgesvd");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_rs(SEXP xin, SEXP only_values) {
#ifdef HIPLAR_WITH_MAGMA
    int *xdims, n, info = 0, ov;
    char jobv[1], uplo[1];
    SEXP values, ret, nm, x, z = R_NilValue;
    double *rz = NULL;
    int lwork, *iwork, liwork, itmp;
    double *work, tmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_rs");
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

//    if (!ov) {
/* with MAGMA 1.3.0 fallback to LAPACK */
//#ifdef HIPLAR_DBG
//R_ShowMessage("DBG: magma_wrapper_rs fallback to Lapack");
//#endif
        /* ask for optimal size of work arrays */
  /*      lwork = -1; liwork = -1;
        F77_CALL(dsyevd)(jobv, uplo, &n, REAL(x), &n, REAL(values), &tmp, &lwork, &itmp, &liwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "dsyevd");
        }
        lwork = (int) tmp;
        liwork = itmp; 
        work = (double *) R_alloc(lwork, sizeof(double));
        iwork = (int *) R_alloc(liwork, sizeof(int));
        F77_CALL(dsyevd)(jobv, uplo, &n, REAL(x), &n, REAL(values), work, &lwork, iwork, &liwork, &info);
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "dsyevd");
        }
*/
    //} else if (magma_interface == MAGMA_CPU_INTERFACE) {
	
	if (magma_interface == MAGMA_CPU_INTERFACE) {

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_rs cpu interface");
#endif
		/* MAGMA cpu interface */

        /* ask for optimal size of work arrays */
        lwork = -1; liwork = -1;
	
		magma_int_t num_gpus = magma_num_gpus();
		if ( num_gpus > 1 ) {
        	magma_dsyevd_m(num_gpus, jobv[0], uplo[0], n, REAL(x), n, REAL(values),
            &tmp, lwork, &itmp, liwork, &info);
		} else {
        	magma_dsyevd(jobv[0], uplo[0], n, REAL(x), n, REAL(values),
            &tmp, lwork, &itmp, liwork, &info);
		}
		if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_dsyevd");
        }

        lwork = (int) tmp;
#if(0)
/* workaround */
if (lwork < (magma_get_dsytrd_nb(n)+2)*n) {
    lwork = (magma_get_dsytrd_nb(n)+2)*n;
}
#endif
        liwork = itmp; 
        work = (double *) R_alloc(lwork, sizeof(double));
        iwork = (int *) R_alloc(liwork, sizeof(int));
		if(num_gpus > 1) {
        	magma_dsyevd_m(num_gpus, jobv[0], uplo[0], n, REAL(x), n, REAL(values),
           	work, lwork, iwork, liwork, &info);
		} else {
        	magma_dsyevd(jobv[0], uplo[0], n, REAL(x), n, REAL(values),
           	work, lwork, iwork, liwork, &info);
		}
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_dsyevd");
        }

	} else {
		/* MAGMA gpu interface */
        int ldda = ((n + 31)/32)*32;
        double *h_A, *h_R, *d_R, *h_work;
        double *w1;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_rs gpu interface");
#endif

// magma_malloc_cpu((void**)&h_A, sizeof(double) * n*n);
// magma_malloc_cpu((void**)&w1,  sizeof(double) * n  );

        h_A = REAL(x);
        w1 = REAL(values);

        magma_malloc_pinned((void**)&h_R, sizeof(double) * n*n);
        magma_malloc((void**)&d_R, sizeof(double) * n*ldda);

        /* ask for optimal size of work arrays */
        lwork = -1; liwork = -1;

        magma_dsyevd_gpu(jobv[0], uplo[0], n, d_R, ldda, w1,
            h_R, n, &tmp, lwork, &itmp, liwork, &info );
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_dsyevd_gpu");
        }

        lwork = (int) tmp;
        liwork = itmp; 

        magma_malloc_pinned((void**)&h_work, sizeof(double) * lwork);
        magma_malloc_cpu((void**)&iwork,  sizeof(int) * liwork);

        magma_dsetmatrix(n, n, h_A, n, d_R, ldda);

        magma_dsyevd_gpu(jobv[0], uplo[0], n, d_R, ldda, w1,
            h_R, n, h_work, lwork, iwork, liwork, &info );
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_dsyevd_gpu");
        }

//        magma_dgetmatrix(n, n, d_R, ldda, h_R, n);

        magma_free_pinned(h_R);
        magma_free(d_R);
        magma_free_pinned(h_work);
        magma_free_cpu(iwork);

    }


    if (!ov) {
		ret = PROTECT(allocVector(VECSXP, 2));
		nm = PROTECT(allocVector(STRSXP, 2));
		SET_STRING_ELT(nm, 1, mkChar("vectors"));
		SET_VECTOR_ELT(ret, 1, x);
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_rs_cmplx(SEXP xin, SEXP only_values) {
#ifdef HIPLAR_WITH_MAGMA
    int *xdims, n, lwork, lrwork, liwork, info, ov;
    char jobv[1], uplo[1];
    SEXP values, ret, nm, x, z = R_NilValue;
    Rcomplex *work, *rz = NULL, ctmp;

    int *iwork, itmp;
    double *rwork, rtmp;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_rs_cmplx");
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
/* with MAGMA 1.3.0 fallback to LAPACK */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_rs_cmplx fallback to Lapack");
#endif
        rwork = (double *) R_alloc((3*n-2) > 1 ? 3*n-2 : 1, sizeof(double));
        /* ask for optimal size of work array */
        lwork = -1;
        F77_CALL(zheev)(jobv, uplo, &n, COMPLEX(x), &n, REAL(values), &ctmp, &lwork, rwork,
               &info); 
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "zheev");
        }
        lwork = (int) ctmp.r;
        work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
        F77_CALL(zheev)(jobv, uplo, &n, COMPLEX(x), &n, REAL(values), work, &lwork, rwork,
            &info); 
        if (info != 0) {
            error("error code %d from Lapack routine '%s'", info, "zheev");
        }
    } else if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_rs_cmplx cpu interface");
#endif

        /* ask for optimal size of work array */
        lwork = -1;
        lrwork = -1;
        liwork = -1;

        magma_zheevd(jobv[0], uplo[0], n, (cuDoubleComplex *)(COMPLEX(x)), n, REAL(values),
            (cuDoubleComplex *)&ctmp, lwork, &rtmp, lrwork, &itmp, liwork, &info);
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_zheevd");
        }

        lwork = (int) ctmp.r;
        work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));
        lrwork = (int) rtmp;
        rwork = (double *) R_alloc(lrwork, sizeof(double));
        liwork = itmp;
        iwork = (int *) R_alloc(liwork, sizeof(int));

        magma_zheevd(jobv[0], uplo[0], n, (cuDoubleComplex *)(COMPLEX(x)), n, REAL(values),
            (cuDoubleComplex *)work, lwork, rwork, lrwork, iwork, liwork, &info);
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_zheevd");
        }

	} else {
		/* MAGMA gpu interface */
        cuDoubleComplex *h_A, *h_R, *d_R, *h_work;
        double *w1;
        int ldda = ((n + 31)/32)*32;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_rs_cmplx gpu interface");
#endif

        magma_zheevd_gpu(jobv[0],uplo[0], n, d_R, n, w1, h_R, n,
            (cuDoubleComplex *)&ctmp, -1, &rtmp, -1, &itmp, -1, &info);
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_zheevd_gpu");
        }

        lwork = (int) ctmp.r;
        lrwork = (int) rtmp;
        liwork = itmp;

//magma_malloc_cpu((void**)&h_A, sizeof(cuDoubleComplex) * n * n);
//magma_malloc_cpu((void**)&w1, sizeof(double) * n);  
        h_A = (cuDoubleComplex *)(COMPLEX(x));
        w1 = REAL(values);

        magma_malloc_pinned((void**)&h_R, sizeof(cuDoubleComplex) * n * n);
        magma_malloc((void**)&d_R, sizeof(cuDoubleComplex) * n * ldda );
        magma_malloc_pinned((void**)&h_work, sizeof(cuDoubleComplex) * lwork);  

//magma_malloc_cpu((void**)&rwork, sizeof(double) * lrwork);
        rwork = (double *) R_alloc(lrwork, sizeof(double));

        magma_malloc_cpu((void**)&iwork, sizeof(int) * liwork);

        magma_zsetmatrix(n, n, h_A, n, d_R, ldda);
        magma_zheevd_gpu(jobv[0], uplo[0], n, d_R, ldda, w1, h_R, n,
            h_work, lwork, rwork, lrwork, iwork, liwork, &info );
        if (info != 0) {
            error("error code %d from MAGMA routine '%s'", info, "magma_zheevd_gpu");
        }

        //magma_zgetmatrix(n, n, d_R, ldda, h_R, n);

        magma_free_pinned(h_R);
        magma_free(d_R);
        magma_free_pinned(h_work);
        magma_free_cpu(iwork);

    }


    if (!ov) {
		ret = PROTECT(allocVector(VECSXP, 2));
		nm = PROTECT(allocVector(STRSXP, 2));
		SET_STRING_ELT(nm, 1, mkChar("vectors"));
		SET_VECTOR_ELT(ret, 1, x);
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
	return R_NilValue; /* -Wall */
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
SEXP magma_wrapper_rg(SEXP x, SEXP only_values) {
#ifdef HIPLAR_WITH_MAGMA
    Rboolean vectors, complexValues;
    int i, n, lwork, info, *xdims, ov;
    double *work, *wR, *wI, *left, *right, *xvals, tmp;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, val;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_rg");
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


// only cpu interface in MAGMA 1.3.0
   
	printf("Here non summ \n");
	magma_int_t num_gpus = magma_num_gpus();
	if ( num_gpus > 1 ) {
		magma_dgeev_m(*jobVL, *jobVR, n, xvals, n, wR, wI, left, n, right, n, &tmp, lwork, &info);
	} else {
		magma_dgeev(*jobVL, *jobVR, n, xvals, n, wR, wI,
        	left, n, right, n, &tmp, lwork, &info);
	}
	if (info != 0) {
    	error("error code %d from MAGMA routine '%s'", info, "magma_dgeev");
    }

    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

// only cpu interface in MAGMA 1.3.0
    if(num_gpus > 1) {
		magma_dgeev_m(*jobVL, *jobVR, n, xvals, n, wR, wI,
        	left, n, right, n, work, lwork, &info);
	} else {
		magma_dgeev(*jobVL, *jobVR, n, xvals, n, wR, wI,
        	left, n, right, n, work, lwork, &info);
	}
	if (info != 0) {
    	error("error code %d from MAGMA routine '%s'", info, "magma_dgeev");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_rg_cmplx(SEXP x, SEXP only_values) {
#ifdef HIPLAR_WITH_MAGMA
    int  n, lwork, info, *xdims, ov;
    Rcomplex *work, *left, *right, *xvals, tmp;
    double *rwork;
    char jobVL[1], jobVR[1];
    SEXP ret, nm, values, val = R_NilValue;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_rg_cmplx");
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

// only cpu interface in MAGMA 1.3.0
    magma_zgeev(*jobVL, *jobVR, n, (cuDoubleComplex *)xvals, n,
        (cuDoubleComplex *)(COMPLEX(values)),
        (cuDoubleComplex *)left, n, (cuDoubleComplex *)right, n,
        (cuDoubleComplex *)&tmp, lwork, rwork, &info);
    if (info != 0) {
	    error("error code %d from MAGMA routine '%s'", info, "magma_zgeev");
    }

    lwork = (int) tmp.r;
    work = (Rcomplex *) R_alloc(lwork, sizeof(Rcomplex));

// only cpu interface in MAGMA 1.3.0
    magma_zgeev(*jobVL, *jobVR, n, (cuDoubleComplex *)xvals, n,
        (cuDoubleComplex *)(COMPLEX(values)),
        (cuDoubleComplex *)left, n, (cuDoubleComplex *)right, n,
        (cuDoubleComplex *)work, lwork, rwork, &info);
    if (info != 0) {
	    error("error code %d from MAGMA routine '%s'", info, "magma_zgeev");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_det_ge_real(SEXP Ain, SEXP logarithm) {
#ifdef HIPLAR_WITH_MAGMA
    int i, n, *Adims, info, *jpvt, sign, useLog;
    double modulus = 0.0; /* -Wall */
    SEXP val, nm, A;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_det_ge_real");
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

    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_det_ge_real cpu interface");
#endif
		magma_dgetrf(n, n, REAL(A), n, jpvt, &info);

	} else {
		/* MAGMA gpu interface */
        double *h_A, *h_R, *d_A;
        int ldda = ((n+31)/32)*32;
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_det_ge_real gpu interface");
#endif

        //magma_malloc_cpu((void**)&ipiv, sizeof(int) * n);
        //magma_malloc_cpu((void**)&h_A, sizeof(double) * n * n);  
        //magma_malloc_pinned((void**)&h_R, sizeof(double) * n * n);  
        h_A = REAL(A);
        h_R = REAL(A);

        magma_malloc((void**)&d_A, sizeof(double) * ldda* n);

        magma_dsetmatrix(n, n, h_R, n, d_A, ldda);
        magma_dgetrf_gpu(n, n, d_A, ldda, jpvt, &info);
        magma_dgetmatrix(n, n, d_A, ldda, h_A, n);

        magma_free(d_A);

	}


    sign = 1;
    if (info < 0) {
    	error("error code %d from MAGMA routine '%s'", info, "magma_dgetrf");
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_dgecon(SEXP A, SEXP norm) {
#ifdef HIPLAR_WITH_MAGMA
    SEXP x, val;
    int *xdims, m, n, info, *iwork;
    double anorm, *work, *work2;
    char typNorm[] = {'\0', '\0'};

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_dgecon");
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

    work = (double *) R_alloc((*typNorm == 'I' && m > 4*n) ? m : 4*n,
			      sizeof(double));
    anorm = F77_CALL(dlange)(typNorm, &m, &n, REAL(x), &m, work);

    iwork = (int *) R_alloc(m, sizeof(int));

    /* Compute the LU-decomposition and overwrite 'x' with result :*/
    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgecon cpu interface");
#endif
		magma_dgetrf(m, n, REAL(x), m, iwork, &info);

	} else {
		/* MAGMA gpu interface */
        double *h_A, *h_R, *d_A;
        int ldda = ((m+31)/32)*32;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_dgecon gpu interface");
#endif

        //magma_malloc_cpu((void**)&ipiv, sizeof(int) * m);
        //magma_malloc_cpu((void**)&h_A, sizeof(double) * m * n);  
        //magma_malloc_pinned((void**)&h_R, sizeof(double) * m * n);  
        h_A = REAL(x);
        h_R = REAL(x);

        magma_malloc((void**)&d_A, sizeof(double) * ldda * n);

        magma_dsetmatrix(m, n, h_R, m, d_A, ldda);
        magma_dgetrf_gpu(m, n, d_A, ldda, iwork, &info);
        magma_dgetmatrix(m, n, d_A, ldda, h_A, m);

        magma_free(d_A);

	}


    if (info) {
    	if (info < 0) {
    	    UNPROTECT(2);
    	    error("error [%d] from magma_dgetrf()", info);
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
	return R_NilValue; /* -Wall */
}


SEXP magma_wrapper_zgecon(SEXP A, SEXP norm) {
#ifdef HIPLAR_WITH_MAGMA
    SEXP val;
    Rcomplex *avals; /* copy of A, to be modified */
    int *dims, n, info;
    double anorm, *rwork;
    char typNorm[] = {'\0', '\0'};

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering magma_wrapper_zgecon");
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

    if (magma_interface == MAGMA_CPU_INTERFACE) {
		/* MAGMA cpu interface */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgecon cpu interface");
#endif
		magma_zgetrf(n, n, (cuDoubleComplex *)avals, n, (int *) R_alloc(n, sizeof(int)), &info);

	} else {
		/* MAGMA gpu interface */
        cuDoubleComplex *h_A, *h_R, *d_A;
        int *ipiv;
        int ldda = ((n+31)/32)*32;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: magma_wrapper_zgecon gpu interface");
#endif

        //magma_malloc_cpu((void**)&ipiv, sizeof(int) * n);
        //magma_malloc_cpu((void**)&h_A, sizeof(cuDoubleComplex) * n * n);  
        //magma_malloc_pinned((void**)&h_R, sizeof(cuDoubleComplex) * n * n);  
        ipiv = (int *) R_alloc(n, sizeof(int));
        h_A = (cuDoubleComplex *)avals;
        h_R = (cuDoubleComplex *)avals;

        magma_malloc((void**)&d_A, sizeof(cuDoubleComplex) * ldda * n);

        magma_zsetmatrix(n, n, h_R, n, d_A, ldda);
        magma_zgetrf_gpu(n, n, d_A, ldda, ipiv, &info);
        magma_zgetmatrix(n, n, d_A, ldda, h_A, n);

        magma_free(d_A);

	}


    if (info) {
    	UNPROTECT(1);
    	error("error [%d] from magma_zgetrf()'", info);
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
	return R_NilValue; /* -Wall */
}


