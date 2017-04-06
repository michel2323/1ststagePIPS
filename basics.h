#ifndef BASICS_H
#define BASICS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

typedef struct mat {
  double **M;
  int m;
  int n;
  
} mat;

typedef struct vec {
  double *v;
  int n;
  
} vec;

mat openM();
vec openSol();
vec openRHS();
void closeM();
void closeSol();
void closeRHS();
void diff(mat M, vec Sol, vec RHS);
void assign(vec b, vec a);
mat copyM(mat a);
void printtimes(double fact, double solv);
void init(int argc, char *argv[]);
void finalize();
void printM(mat M);
#endif