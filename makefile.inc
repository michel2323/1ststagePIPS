

OMP_NUM_THREADS=4
NUMPROCS=4
CC=mpicc
CXX=mpicxx
CFLAGS=-g -O0
CXXFLAGS=-g -O0 -std=c++11
LFLAGS=-lm
LLAPACK=-llapack
LSCALAPACK=-L/home/michel/src/scalapack-2.0.2/build/lib -lscalapack -llapack -lblas -lgfortran
ELEMENTAL_PATH=/home/michel/local
