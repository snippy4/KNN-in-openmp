
CFLAGS = -std=c99 -march=native -mavx2 -O3 

all: gccserial gccomp

gccnearly: k-folds.c
	gcc $(CFLAGS) -fopenmp -o k-folds-gcc k-folds.c

gcccomplete: mpi-k-folds.c
	mpicc $(CFLAGS) -fopenmp -o k-folds-complete-gcc mpi-k-folds.c
	
iccnearly: k-folds.c
	icc $(CFLAGS) -qopenmp -o k-folds-gcc k-folds.c

icccomplete: mpi-k-folds.c
	mpiicc $(CFLAGS) -qopenmp -o k-folds-complete-gcc mpi-k-folds.c

