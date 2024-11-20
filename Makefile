CC = gcc

CFLAGS = -std=c99 -march=native -mavx2 -O3
OMP_FLAGS = -fopenmp

all: gccserial gccomp

gccserial: knn.c
	$(CC) $(CFLAGS) -o knn-gcc knn.c

gccomp: knnomp.c
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o knnomp-gcc knnomp.c
	
iccserial: knn.c
	icc $(CFLAGS) -o knn-icc knn.c

iccomp: knnomp.c
	icc $(CFLAGS) -qopenmp -o knnomp-icc knnomp.c

gccaknn: aknn.c
	$(CC) $(CFLAGS) $(OMP_FLAGS) -o aknn-gcc aknn.c

clean:
	rm -f knn knnomp

