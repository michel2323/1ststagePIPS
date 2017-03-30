
include makefile.inc


all: dsytr dgetr pdgetr
	
dsytr: dsytr.o basics.o	
	$(CC) -o $@ $^ $(LLAPACK) $(LFLAGS) 

dgetr: dgetr.o basics.o	
	$(CC) -o $@ $^ $(LLAPACK) $(LFLAGS) 
	
pdgetr: pdgetr.o basics.o	
	$(CC) -o $@ $^ $(LSCALAPACK) $(LFLAGS) 
	
dsytr.o: dsytr.c
	$(CC) $(CFLAGS) -c $<
	
dgetr.o: dgetr.c
	$(CC) $(CFLAGS) -c $<
	
pdgetr.o: pdgetr.c
	$(CC) $(CFLAGS) -c $<
	
basics.o: basics.c
	$(CC) $(CFLAGS) -c $<
	
.PHONY: clean run
	
run:
	./dsytr
	./dgetr
	mpiexec -n 4 ./pdgetr
	
clean: 
	rm *.o dsytr dgetr pdgetr
	