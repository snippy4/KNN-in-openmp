#!/bin/bash -l
#Written by Dr Maryam Abo Tabik
# Specify the current working directory as the location for executables/files
# This is the default setting.
#SBATCH -D ./

# Export the current environment to the compute node
# This is the default setting.
#SBATCH --export=ALL

# Specific course queue, exclusive use (for timings), max 1 min wallclock time
#SBATCH -p lowpriority
#SBATCH -t 30:00
#SBATCH -N 2
#SBATCH -n 80
export OMP_NUM_THREADS=8
# load modules
module load compilers/intel/2019u5
module load mpi/intel-mpi/2019u5/bin
module load libs/nvidia-cuda/12.4.0/bin
# GNU no-opt
mpicc -fopenmp -march=native -lm -std=c99 -O3 k-folds-mpi.c -o knnmpi.out
time mpirun -np 1 ./knnmpi.out "data/asteroids.csv" "out.csv" 3 10

echo '-------'