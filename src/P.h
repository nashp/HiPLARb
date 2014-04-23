
#define CHECK_VERSION_BEQ(MAJ,MIN,MIC) \
    ( (PLASMA_VERSION_MAJOR < MAJ) || \
      ((PLASMA_VERSION_MAJOR == MAJ) && (PLASMA_VERSION_MINOR < MIN)) || \
      ((PLASMA_VERSION_MAJOR == MAJ) && (PLASMA_VERSION_MINOR == MIN) && (PLASMA_VERSION_MICRO <= MIC)) \
    )

#define CHECK_VVERSION_BEQ(MAJ,MIN,MIC) \
    ( (R_PLASMA_MAJOR < MAJ) || \
      ((R_PLASMA_MAJOR == MAJ) && (R_PLASMA_MINOR < MIN)) || \
      ((R_PLASMA_MAJOR == MAJ) && (R_PLASMA_MINOR == MIN) && (R_PLASMA_MICRO <= MIC)) \
    )


int P_dpotrf(const char *uplo, int N, double *A, int LDA);
int P_dpotri(const char *uplo, int N, double *A, int LDA);
int P_dpotrs(const char *uplo, int N, int NRHS,
    double *A, int LDA, double *B, int LDB);

int P_zgeqrf(int M, int N, void *A, void *T);
int P_dgeqrf(int M, int N, double *A, double *T);

int factorQR_T_size(int M, int N);

int P_zunmqr(const char *side, const char *trans,
	int M, int N, int K,
	void *A, int LDA,
	void *T,
	void *B, int LDB);

int P_dormqr(const char *side, const char *trans,
	int M, int N, int K,
	double *A, int LDA,
	double *T,
	double *B, int LDB);

int P_dtrsm(const char *side, const char *uplo,
    const char *transA, const char *diag,
	int N, int NRHS, double alpha,
	double *A, int LDA,
	double *B, int LDB);

int P_ztrsm(const char *side, const char *uplo,
    const char *transA, const char *diag,
	int N, int NRHS, Rcomplex alpha,
	void *A, int LDA,
	void *B, int LDB);

int P_dtrtri(const char *uplo, const char *diag,
    int N, double *A, int LDA);

int P_dgesv(int N, int NRHS,
	double *A, int LDA, int *IPIV, double *B, int LDB);

int P_zgesv(int N, int NRHS,
	void *A, int LDA, int *IPIV, void *B, int LDB);

int P_dgemm(const char *transA, const char *transB,
	int M, int N, int K,
	double alpha, double *A, int LDA,
	double *B, int LDB, double beta,
	double *C, int LDC);

int P_zgemm(const char *transA, const char *transB,
	int M, int N, int K,
	Rcomplex alpha, void *A, int LDA,
	void *B, int LDB, Rcomplex beta,
	void *C, int LDC);

int P_dorgqr(int M, int N, int K,
	double *A, int LDA, double *T,
	double *Q, int LDQ);

int P_zungqr(int M, int N, int K,
	void *A, int LDA, void *T,
	void *Q, int LDQ);

double P_dlange(const char *type,
	int M, int N, double *A, int LDA, double *WORK);
double P_dlansy(const char *type, const char *uplo,
    int N, double *A, int LDA, double *WORK);

int P_dgesvd(const char *jobu, const char *jobvt,
	int M, int N,
	double *A, int LDA, double *S,
	double *U, int LDU,
	double *VT, int LDVT);

int P_zgesvd(const char *jobu, const char *jobvt,
	int M, int N,
	void *A, int LDA, double *S,
	void *U, int LDU,
	void *VT, int LDVT);

int P_dsyev(const char *jobz, const char *uplo,
	int N, double *A, int LDA,
	double *W, double *Q, int LDQ);

int P_zheev(const char *jobz, const char *uplo,
	int N, void *A, int LDA,
	double *W, void *Q, int LDQ);

int P_dgetrf(int M, int N, double *A, int LDA, int *IPIV);
int P_zgetrf(int M, int N, void *A, int LDA, int *IPIV);

int P_dgetri(int N, double *A, int LDA, int *IPIV);

int P_dgetrs(const char *trans, int N, int NRHS,
    double *A, int LDA, int *IPIV, double *B, int LDB);

int P_dsyrk(const char *uplo, const char *trans,
    int N, int K, double alpha, double *A, int LDA,
    double beta, double *C, int LDC);

int P_dsymm(const char *side, const char *uplo,
    int M, int N, double alpha, double *A, int LDA,
    double *B, int LDB, double beta, double *C, int LDC);

int P_dtrmm(const char *side, const char *uplo,
    const char *transA, const char *diag,
    int N, int NRHS, double alpha, double *A, int LDA,
    double *B, int LDB);

int P_dplrnt(int M, int N, double *A, int LDA,
    unsigned long long int seed);

int P_zplrnt(int M, int N, void *A, int LDA,
    unsigned long long int seed);

int P_dplgsy(double bump, int N,
    double *A, int LDA, unsigned long long int seed);

