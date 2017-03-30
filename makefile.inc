

OMP_NUM_THREADS=4
NUMPROCS=4
CC=mpicc
CFLAGS=-g -O0
LFLAGS=-lm
LLAPACK=-llapack
LSCALAPACK=-L/home/michel/src/scalapack-2.0.2/build/lib -lscalapack -llapack -lblas -lgfortran