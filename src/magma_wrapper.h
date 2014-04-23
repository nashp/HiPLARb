#ifndef MAGMA_WRAPPER
#define MAGMA_WRAPPER

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


SEXP magma_wrapper_zgeqrf(SEXP Ain);
SEXP magma_wrapper_coef_cmplx(SEXP Q, SEXP Bin);
SEXP magma_wrapper_qr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans);
SEXP magma_wrapper_dgeqrf(SEXP Ain);
SEXP magma_wrapper_coef_real(SEXP Q, SEXP Bin);
SEXP magma_wrapper_qr_qy_real(SEXP Q, SEXP Bin, SEXP trans);

SEXP magma_wrapper_chol(SEXP A, SEXP pivot);
SEXP magma_wrapper_chol2inv(SEXP A, SEXP size);

SEXP magma_wrapper_zgesv(SEXP A, SEXP Bin);
SEXP magma_wrapper_dgesv(SEXP A, SEXP Bin, SEXP tolin);

SEXP magma_wrapper_dlange(SEXP A, SEXP type);

void magma_wrapper_bakslv(double *t, int *ldt, int *n,
    double *b, int *ldb, int *nb, double *x, int *job, int *info);
SEXP magma_wrapper_bakslv2(SEXP r, SEXP x, SEXP k, SEXP nb, SEXP utri, SEXP trans);

SEXP magma_wrapper_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v);
SEXP magma_wrapper_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v);

SEXP magma_wrapper_rs(SEXP xin, SEXP only_values);
SEXP magma_wrapper_rs_cmplx(SEXP xin, SEXP only_values);
SEXP magma_wrapper_rg(SEXP x, SEXP only_values);
SEXP magma_wrapper_rg_cmplx(SEXP x, SEXP only_values);

SEXP magma_wrapper_det_ge_real(SEXP Ain, SEXP logarithm);

SEXP magma_wrapper_dgecon(SEXP A, SEXP norm);
SEXP magma_wrapper_zgecon(SEXP A, SEXP norm);

#endif
