#include "basics.h"

mat openM() {
  mat ret;
  int i=0;
  FILE *fp=fopen("./1ststageM.dmp","r");  
  if(fread(&ret.m, sizeof(int), 1, fp)==0) {
    printf("Error reading size m in 1ststageM.dmp %d\n", ret.m);
    exit(0);
  }
  if(fread(&ret.n, sizeof(int), 1, fp)==0) {
    printf("Error reading size n in 1ststageM.dmp %d\n", ret.n);
    exit(0);
  }
  printf("Matrix M is of size %dx%d.\n", ret.m, ret.n);
  ret.M=(double**) malloc(sizeof(double*)*ret.m);
  ret.M[0]=(double*) malloc(sizeof(double)*ret.m*ret.n);
  for(i=1;i<ret.m;i=i+1) {
    ret.M[i]=ret.M[0]+ret.n*i;
  }
  for(i=0;i<ret.m;i=i+1) {
    if(fread(ret.M[i], sizeof(double), ret.n, fp)==0) {
      printf("Error reading data in 1ststageM.dmp\n");
    }
  }
  fclose(fp);
  return ret;  
}
void closeM(mat M) {
  free(M.M[0]);
  free(M.M);
}

vec openSol() {
  vec ret;
  int tmp;
  FILE *fp=fopen("./1ststageSol.dmp","r");  
  if(fread(&ret.n, sizeof(int), 1, fp)==0) {
    printf("Error reading size n in 1ststageSol.dmp %d\n", ret.n);
    exit(0);
  }
  else {
    printf("Vector Sol is of size %d.\n", ret.n);
  }
  fread(&tmp, sizeof(int), 1, fp);
  ret.v=(double*) malloc(sizeof(double)*ret.n);
  if(fread(ret.v, sizeof(double), ret.n, fp)==0) {
      printf("Error reading data in 1ststageSol.dmp\n");
  }
  fclose(fp);
  return ret;  
  
}

void closeSol(vec v) {
  free(v.v);
}

vec openRHS() {
  vec ret;
  int tmp;
  FILE *fp=fopen("./1ststageRHS.dmp","r");  
  if(fread(&ret.n, sizeof(int), 1, fp)==0) {
    printf("Error reading size n in 1ststageRHS.dmp %d\n", ret.n);
    exit(0);
  }
  else {
    printf("Vector RHS is of size %d.\n", ret.n);
  }
  fread(&tmp, sizeof(int), 1, fp);
  ret.v=(double*) malloc(sizeof(double)*ret.n);
  if(fread(ret.v, sizeof(double), ret.n, fp)==0) {
      printf("Error reading data in 1ststageRHS.dmp\n");
  }
  fclose(fp);
  return ret;  
  
}

void closeRHS(vec v) {
  free(v.v);
}

double norm(vec v) {
  int i;
  double norm=0;
  for(i=0;i<v.n;i=i+1) {
    norm=norm+v.v[i]*v.v[i];  
  }
  norm=sqrt(norm);
  return norm;
}

void diff(mat M, vec Sol, vec RHS) {
  int i,j;
  vec diff;
  diff.v=(double*) malloc(sizeof(double)*M.n);
  diff.n=M.n;
  
  for(i=0;i<M.n;i=i+1) {
    diff.v[i]=0;
    for(j=0;j<M.n;j=j+1) {
      diff.v[i]=diff.v[i]+M.M[i][j]*Sol.v[j];
    }
  }
  printf("Norm of RHS: %e.\n", norm(RHS));
  printf("Norm of M*Sol: %e.\n", norm(diff));
  for(i=0;i<M.n;i=i+1) {
    diff.v[i]=RHS.v[i]-diff.v[i];
  }
  printf("Norm of RHS-M*Sol: %e.\n", norm(diff));
}

void assign(vec b, vec a) {
  int i;
  assert(b.n==a.n);
  for(i=0;i<a.n;i=i+1) b.v[i]=a.v[i];
}

mat copyM(mat A) {
  mat ret;
  int i,j;
  ret.n=A.n;
  ret.M=(double**) malloc(sizeof(double*)*ret.n);
  ret.M[0]=(double*) malloc(sizeof(double)*ret.n*ret.n);
  for(i=1;i<ret.n;i=i+1) {
    ret.M[i]=ret.M[0]+ret.n*i;
  }
  for(i=0;i<ret.n;i=i+1) {
    for(j=0;j<ret.n;j=j+1) {
      ret.M[i][j]=A.M[i][j];
    }
  }
  return ret;
}

void printtimes(double fact, double solv) {
  printf("Factorization: %f.\n", fact);
  printf("Solve: %f.\n", fact);
}

void init(int argc, char *argv[]) {
  MPI_Init(&argc,&argv);
}

void finalize() {
  MPI_Finalize();
}

void printM(mat M) {
  int i,j;
  for(i=0;i<M.n;i=i+1) {
    for(j=0;j<M.n;j=j+1) {
      printf("%e\n", M.M[j][i]);
    }
  }
  exit(0);
}