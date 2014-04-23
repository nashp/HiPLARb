#include <stdlib.h>
#include <string.h>

#ifdef HIPLAR_WITH_PLASMA
#include <plasma.h>

#include <R_ext/Complex.h>

#include "P.h"
 

int P_dpotrf(
const char *uplo,
int N,
double *A,
int LDA
) {
	int info;
	PLASMA_enum uo;

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}
	info = PLASMA_dpotrf(uo, N, A, LDA);

	return(info);
}


int P_dpotri(
const char *uplo,
int N,
double *A,
int LDA
) {
	int info;
	PLASMA_enum uo;

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}
	info = PLASMA_dpotri(uo, N, A, LDA);

	return(info);
}


int P_dpotrs(
const char *uplo,
int N,
int NRHS,
double *A,
int LDA,
double *B,
int LDB
) {
	PLASMA_enum uo;

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

    return PLASMA_dpotrs(uo, N, NRHS, A, LDA, B, LDB);
}


int P_dgeqrf(
int M,
int N,
double *A,
double *T
) {
	int info;

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_dgeqrf(M, N, A, M, T);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaRealDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);
/*
    printf("MB=%d NB=%d BSIZ=%d LM=%d LN=%d M=%d N=%d MT=%d NT=%d\n",
        descT->mb, descT->nb, descT->bsiz, descT->lm, descT->ln,
        descT->m, descT->n, descT->mt, descT->nt);
*/

	info = PLASMA_dgeqrf(M, N, A, M, descT);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);

}


int P_zgeqrf(
int M,
int N,
void *A,
void *T
) {
	int info;

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_zgeqrf(M, N, A, M, T);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_zgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaComplexDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

	info = PLASMA_zgeqrf(M, N, A, M, descT);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);

}


int factorQR_T_size(int M, int N) {
	int size;
	int NB, IB;
	int MT, NT;
    PLASMA_desc *descT;
    double *T;

    /* Get autotuned or set tile size; actual allocation of memory with R */
#if CHECK_VERSION_BEQ(2,4,5)
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &T);
#else
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &descT);
#endif
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
#if CHECK_VERSION_BEQ(2,4,5)
    free(T);
#else
    PLASMA_Dealloc_Handle_Tile(&descT);
#endif

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

	size = MT*NT*IB*NB;

	return(size);
}


int P_dormqr(
const char *side,
const char *trans,
int M,
int N,
int K,
double *A,
int LDA,
double *T,
double *B,
int LDB
) {
	PLASMA_enum s, t;
	int info;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*trans == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_dormqr(s, t, M, N, K, A, LDA, T, B, LDB);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaComplexDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

	info = PLASMA_dormqr(s, t, M, N, K, A, LDA, descT, B, LDB);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);
}


int P_zunmqr(
const char *side,
const char *trans,
int M,
int N,
int K,
void *A,
int LDA,
void *T,
void *B,
int LDB
) {
	PLASMA_enum s, t;
	int info;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*trans == 'C') {
		t = PlasmaConjTrans;
	} else {
		t = PlasmaNoTrans;
	}

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_zunmqr(s, t, M, N, K, A, LDA, T, B, LDB);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_zgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaComplexDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

	info = PLASMA_zunmqr(s, t, M, N, K, A, LDA, descT, B, LDB);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);
}


int P_dtrsm(
const char *side,
const char *uplo,
const char *transA,
const char *diag,
int N,
int NRHS,
double alpha,
double *A,
int LDA,
double *B,
int LDB
) {
	PLASMA_enum s, uo, t, d;
	int info;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	if (*transA == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}

	if (*diag == 'U') {
		d = PlasmaUnit;
	} else {
		d = PlasmaNonUnit;
	}

	info = PLASMA_dtrsm(s, uo, t, d, N, NRHS, alpha, A, LDA, B, LDB);

	return(info);
}


int P_ztrsm(
const char *side,
const char *uplo,
const char *transA,
const char *diag,
int N,
int NRHS,
Rcomplex alpha,
void *A,
int LDA,
void *B,
int LDB
) {
	PLASMA_enum s, uo, t, d;
	PLASMA_Complex64_t alpha1;
	int info;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	if (*transA == 'C') {
		t = PlasmaConjTrans;
	} else if (*transA == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}
	
	if (*diag == 'U') {
		d = PlasmaUnit;
	} else {
		d = PlasmaNonUnit;
	}

	alpha1 = alpha.r + alpha.i*I;

	info = PLASMA_ztrsm(s, uo, t, d, N, NRHS, alpha1, A, LDA, B, LDB);

	return(info);
}


int P_dtrtri(
const char *uplo,
const char *diag,
int N,
double *A,
int LDA
) {
	PLASMA_enum uo, d;

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	if (*diag == 'U') {
		d = PlasmaUnit;
	} else {
		d = PlasmaNonUnit;
	}

    return PLASMA_dtrtri(uo, d, N, A, LDA);
}


int P_dgesv(
int N,
int NRHS,
double *A,
int LDA,
int *IPIV,
double *B,
int LDB
) {
	int info;

	info = PLASMA_dgesv(N, NRHS, A, LDA, IPIV, B, LDB);

	return(info);
}


int P_zgesv(
int N,
int NRHS,
void *A,
int LDA,
int *IPIV,
void *B,
int LDB
) {
	int info;

	info = PLASMA_zgesv(N, NRHS, A, LDA, IPIV, B, LDB);

	return(info);
}


int P_dgemm(
const char *transA,
const char *transB,
int M,
int N,
int K,
double alpha,
double *A,
int LDA,
double *B,
int LDB,
double beta,
double *C,
int LDC
) {
	PLASMA_enum t1, t2;
	int info;

	if (*transA == 'T') {
		t1 = PlasmaTrans;
	} else {
		t1 = PlasmaNoTrans;
	}

	if (*transB == 'T') {
		t2 = PlasmaTrans;
	} else {
		t2 = PlasmaNoTrans;
	}

	info = PLASMA_dgemm(t1, t2, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC);

	return(info);
}


int P_zgemm(
const char *transA,
const char *transB,
int M,
int N,
int K,
Rcomplex alpha,
void *A,
int LDA,
void *B,
int LDB,
Rcomplex beta,
void *C,
int LDC
) {
	PLASMA_enum t1, t2;
	PLASMA_Complex64_t alpha1, beta1;
	int info;

	if (*transA == 'C') {
		t1 = PlasmaConjTrans;
	} else if (*transA == 'T') {
		t1 = PlasmaTrans;
	} else {
		t1 = PlasmaNoTrans;
	}

	if (*transB == 'C') {
		t2 = PlasmaConjTrans;
	} else if (*transB == 'T') {
		t2 = PlasmaTrans;
	} else {
		t2 = PlasmaNoTrans;
	}

	alpha1 = alpha.r + alpha.i*I;
	beta1 = beta.r + beta.i*I;

	info = PLASMA_zgemm(t1, t2, M, N, K, alpha1, A, LDA, B, LDB, beta1, C, LDC);

	return(info);
}


int P_dorgqr(
int M,
int N,
int K,
double *A,
int LDA,
double *T,
double *Q,
int LDQ
) {
	int info;

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_dorgqr(M, N, K, A, LDA, T, Q, LDQ);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaComplexDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

	info = PLASMA_dorgqr(M, N, K, A, LDA, descT, Q, LDQ);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);
}


int P_zungqr(
int M,
int N,
int K,
void *A,
int LDA,
void *T,
void *Q,
int LDQ
) {
	int info;

#if CHECK_VERSION_BEQ(2,4,5)
	info = PLASMA_zungqr(M, N, K, A, LDA, T, Q, LDQ);
#else
    PLASMA_desc *descT;
    int NB, IB;
    int MT, NT;

    /* Get autotuned or set tile size; T matrix allocated with R */
    PLASMA_Alloc_Workspace_dgeqrf(1, 1, &descT);
	PLASMA_Get(PLASMA_TILE_SIZE, &NB);
    PLASMA_Get(PLASMA_INNER_BLOCK_SIZE, &IB);
    PLASMA_Dealloc_Handle_Tile(&descT);

	MT = (M%NB==0) ? (M/NB) : (M/NB+1);
	NT = (N%NB==0) ? (N/NB) : (N/NB+1);

// possibly allocate space for descT in R and keep it in qr object instead
    info = PLASMA_Desc_Create(&descT, T, PlasmaComplexDouble,
         IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

	info = PLASMA_zungqr(M, N, K, A, LDA, descT, Q, LDQ);

    PLASMA_Desc_Destroy(&descT);
#endif

	return(info);
}


double P_dlange(
const char *type,
int M,
int N,
double *A,
int LDA,
double *WORK
) {
	PLASMA_enum norm;

	if (*type == 'M') {
		norm = PlasmaMaxNorm;
	} else if (*type == 'O') {
		norm = PlasmaOneNorm;
	} else if (*type == 'I') {
		norm = PlasmaInfNorm;
	} else if (*type == 'F') {
		norm = PlasmaFrobeniusNorm;
	} else {
		return 0.0;
	}

	return PLASMA_dlange(norm, M, N, A, LDA);
}


double P_dlansy(
const char *type,
const char *uplo,
int N,
double *A,
int LDA,
double *WORK
) {
	PLASMA_enum norm, uo;

	if (*type == 'M') {
		norm = PlasmaMaxNorm;
	} else if (*type == 'O') {
		norm = PlasmaOneNorm;
	} else if (*type == 'I') {
		norm = PlasmaInfNorm;
	} else if (*type == 'F') {
		norm = PlasmaFrobeniusNorm;
	} else {
		return 0.0;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

    return PLASMA_dlansy(norm, uo, N, A, LDA);
}


int P_dgesvd(
const char *jobu,
const char *jobvt,
int M,
int N,
double *A,
int LDA,
double *S,
double *U,
int LDU,
double *VT,
int LDVT
) {
	PLASMA_enum ju, jvt;
	PLASMA_desc *descT;
	int info;

/*
	if (*jobu != 'N') {
		return(-1);
	}
	if (*jobvt != 'N') {
		return(-2);
	}
*/

	ju = PlasmaNoVec;
	jvt = PlasmaNoVec;

	PLASMA_Alloc_Workspace_dgesvd(M, N, &descT);

	info = PLASMA_dgesvd(ju, jvt, M, N, A, LDA, S, descT, U, LDU, VT, LDVT);

	PLASMA_Dealloc_Handle_Tile(&descT);

	return(info);
}


int P_zgesvd(
const char *jobu,
const char *jobvt,
int M,
int N,
void *A,
int LDA,
double *S,
void *U,
int LDU,
void *VT,
int LDVT
) {
	PLASMA_enum ju, jvt;
	PLASMA_desc *descT;
	int info;

/*
	if (*jobu != 'N') {
		return(-1);
	}
	if (*jobvt != 'N') {
		return(-2);
	}
*/

	ju = PlasmaNoVec;
	jvt = PlasmaNoVec;

	PLASMA_Alloc_Workspace_zgesvd(M, N, &descT);

	info = PLASMA_zgesvd(ju, jvt, M, N, A, LDA, S, descT, U, LDU, VT, LDVT);

	PLASMA_Dealloc_Handle_Tile(&descT);

	return(info);
}


int P_dsyev(
const char *jobz,
const char *uplo,
int N,
double *A,
int LDA,
double *W,
double *Q,
int LDQ
) {
	PLASMA_enum jz, uo;
	PLASMA_desc *descT;
	int info;

	jz = PlasmaNoVec;
	if (*jobz != 'N') {
		jz = PlasmaVec;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	PLASMA_Alloc_Workspace_dsyev(N, N, &descT);

	info = PLASMA_dsyev(jz, uo, N, A, LDA, W, descT, Q, LDQ);

	PLASMA_Dealloc_Handle_Tile(&descT);

	return info;
}


int P_zheev(
const char *jobz,
const char *uplo,
int N,
void *A,
int LDA,
double *W,
void *Q,
int LDQ
) {
	PLASMA_enum jz, uo;
	PLASMA_desc *descT;
	int info;

	jz = PlasmaNoVec;
	if (*jobz != 'N') {
		jz = PlasmaVec;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	PLASMA_Alloc_Workspace_zheev(N, N, &descT);

	info = PLASMA_zheev(jz, uo, N, A, LDA, W, descT, Q, LDQ);

	PLASMA_Dealloc_Handle_Tile(&descT);

	return info;
}


int P_dgetrf(
int M,
int N,
double *A,
int LDA,
int *IPIV
) {
    return PLASMA_dgetrf(M, N, A, LDA, IPIV);
}


int P_zgetrf(
int M,
int N,
void *A,
int LDA,
int *IPIV
) {
    return PLASMA_zgetrf(M, N, A, LDA, IPIV);
}


int P_dgetri(
int N,
double *A,
int LDA,
int *IPIV
) {
    return PLASMA_dgetri(N, A, LDA, IPIV);
}


int P_dgetrs(
const char *trans,
int N,
int NRHS,
double *A,
int LDA,
int *IPIV,
double *B,
int LDB
) {
	PLASMA_enum t;

	if (*trans == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}

    return PLASMA_dgetrs(t, N, NRHS, A, LDA, IPIV, B, LDB);

}


int P_dsyrk(
const char *uplo,
const char *trans,
int N,
int K,
double alpha,
double *A,
int LDA,
double beta,
double *C,
int LDC
) {
	PLASMA_enum uo, t;

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	if (*trans == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}

    return PLASMA_dsyrk(uo, t, N, K, alpha, A, LDA, beta, C, LDC);

}

 
int P_dsymm(
const char *side,
const char *uplo,
int M,
int N,
double alpha,
double *A,
int LDA,
double *B,
int LDB,
double beta,
double *C,
int LDC
) {
	PLASMA_enum s, uo;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

    return PLASMA_dsymm(s, uo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);

}


int P_dtrmm(
const char *side,
const char *uplo,
const char *transA,
const char *diag,
int N,
int NRHS,
double alpha,
double *A,
int LDA,
double *B,
int LDB
) {
	PLASMA_enum s, uo, t, d;

	if (*side == 'L') {
		s = PlasmaLeft;
	} else {
		s = PlasmaRight;
	}

	if (*uplo == 'U') {
		uo = PlasmaUpper;
	} else {
		uo = PlasmaLower;
	}

	if (*transA == 'T') {
		t = PlasmaTrans;
	} else {
		t = PlasmaNoTrans;
	}

	if (*diag == 'U') {
		d = PlasmaUnit;
	} else {
		d = PlasmaNonUnit;
	}

    return PLASMA_dtrmm(s, uo, t, d, N, NRHS, alpha, A, LDA, B, LDB);

}


int P_dplrnt(
int M,
int N,
double *A,
int LDA,
unsigned long long int seed) {

    return PLASMA_dplrnt(M, N, A, LDA, seed);

}


int P_zplrnt(
int M,
int N,
void *A,
int LDA,
unsigned long long int seed) {

    return PLASMA_zplrnt(M, N, A, LDA, seed);

}


int P_dplgsy(
double bump,
int N,
double *A,
int LDA,
unsigned long long int seed) {

    return PLASMA_dplgsy(bump, N, A, LDA, seed);

}

#endif
