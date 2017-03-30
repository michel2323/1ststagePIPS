#include "basics.h"
// declarations for LAPACK functions used to factor/solve:

// dsytrf_() factors a symmetric indefinite matrix A, see LAPACK 
// documentation for more details.
extern void FNAME(dgetrf)( 
	    int *m,
			int *n, 
			double A[], 
			int *lda, 
			int ipiv[], 
			int *info);

// dsytrs_() solves the system Ax = b using the factor obtained by dsytrf_().
extern void FNAME(dgetrs)(char *trans, 
			int *n, 
			int *nrhs, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double b[], 
			int *ldb,
			int *info);

int main(int argc, char *argv[]) {
  
  char trans = 'N';
  int info;
  int lwork=-1;
  int one=1;
  double lworkNew;
  int *ipiv;
  double *work;
  int i;
  
  init(argc,argv);
  mat M=openM();
  vec Sol=openSol();
  vec RHS=openRHS();
  diff(M, Sol, RHS);
  mat Mold=copyM(M);
  ipiv=(int*) calloc(M.n,sizeof(int));  
  double t0=MPI_Wtime();
  FNAME(dgetrf)( &M.n, &M.n, M.M[0], &M.n, ipiv, &info );
  double t1=MPI_Wtime();
  assign(Sol,RHS);
  double t2=MPI_Wtime();
  FNAME(dgetrs)( &trans, &M.n, &one,	M.M[0],	&M.n, ipiv, Sol.v,	&M.n,	&info);
  double t3=MPI_Wtime();
  diff(Mold, Sol, RHS);
  printtimes(t3-t2,t1-t2);
  
  closeSol(Sol);
  closeRHS(RHS);
  closeM(M);
}