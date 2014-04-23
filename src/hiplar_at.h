#ifndef HIPLAR_AT
#define HIPLAR_AT

#define HIPLAR_USE_PLASMA 1
#define HIPLAR_USE_MAGMA  2
#define HIPLAR_USE_AUTO   3

#define MAGMA_CPU_INTERFACE 0
#define MAGMA_GPU_INTERFACE 1


extern int hiplar_library;
extern int magma_interface;

extern int xover_zgeqrf;
extern int xover_dgeqrf;
extern int xover_chol;
extern int xover_chol2inv;
extern int xover_zgesv;
extern int xover_dgesv;
extern int xover_dlange;
extern int xover_bakslv;
extern int xover_svd;
extern int xover_svd_cmplx;
extern int xover_rs;
extern int xover_rs_cmplx;
extern int xover_rg;
extern int xover_rg_cmplx;
extern int xover_det_ge_real;
extern int xover_dgecon;
extern int xover_zgecon;
extern int xover_matprod;

extern int maxmagma_zgeqrf;
extern int maxmagma_dgeqrf;
extern int maxmagma_chol;
extern int maxmagma_chol2inv;
extern int maxmagma_zgesv;
extern int maxmagma_dgesv;
extern int maxmagma_dlange;
extern int maxmagma_bakslv;
extern int maxmagma_svd;
extern int maxmagma_svd_cmplx;
extern int maxmagma_rs;
extern int maxmagma_rs_cmplx;
extern int maxmagma_rg;
extern int maxmagma_rg_cmplx;
extern int maxmagma_det_ge_real;
extern int maxmagma_dgecon;
extern int maxmagma_zgecon;
extern int maxmagma_matprod;

#endif
