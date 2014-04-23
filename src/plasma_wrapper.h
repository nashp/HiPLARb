#ifndef PLASMA_WRAPPER
#define PLASMA_WRAPPER

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


SEXP plasma_wrapper_zgeqrf(SEXP Ain);
SEXP plasma_wrapper_coef_cmplx(SEXP Q, SEXP Bin);
SEXP plasma_wrapper_qr_qy_cmplx(SEXP Q, SEXP Bin, SEXP trans);
SEXP plasma_wrapper_dgeqrf(SEXP Ain);
SEXP plasma_wrapper_coef_real(SEXP Q, SEXP Bin);
SEXP plasma_wrapper_qr_qy_real(SEXP Q, SEXP Bin, SEXP trans);

SEXP plasma_wrapper_chol(SEXP A, SEXP pivot);
SEXP plasma_wrapper_chol2inv(SEXP A, SEXP size);

SEXP plasma_wrapper_zgesv(SEXP A, SEXP Bin);
SEXP plasma_wrapper_dgesv(SEXP A, SEXP Bin, SEXP tolin);

SEXP plasma_wrapper_dlange(SEXP A, SEXP type);

void plasma_wrapper_bakslv(double *t, int *ldt, int *n,
    double *b, int *ldb, int *nb, double *x, int *job, int *info);
SEXP plasma_wrapper_bakslv2(SEXP r, SEXP x, SEXP k, SEXP nb, SEXP utri, SEXP trans);

SEXP plasma_wrapper_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v);
SEXP plasma_wrapper_svd_cmplx(SEXP jobu, SEXP jobv, SEXP xin, SEXP s, SEXP u, SEXP v);

SEXP plasma_wrapper_rs(SEXP xin, SEXP only_values);
SEXP plasma_wrapper_rs_cmplx(SEXP xin, SEXP only_values);
SEXP plasma_wrapper_rg(SEXP x, SEXP only_values);
SEXP plasma_wrapper_rg_cmplx(SEXP x, SEXP only_values);

SEXP plasma_wrapper_det_ge_real(SEXP Ain, SEXP logarithm);

SEXP plasma_wrapper_dgecon(SEXP A, SEXP norm);
SEXP plasma_wrapper_zgecon(SEXP A, SEXP norm);

#endif
