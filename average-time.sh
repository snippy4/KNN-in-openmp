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
#SBATCH -t 5:00
#SBATCH -N 1
#SBATCH -n 40

# load modules
module load compilers/intel/2019u5

# number of threads to run on
export OMP_NUM_THREADS=40
# The program to run (replace with your actual executable)

# Number of times to run the program
NUM_RUNS=5

# Variable to accumulate total time
total_time=0
gcc -fopenmp -march=native -mavx2 -O3 knnomp.c -std=c99 -o knn.out
# Loop to run the program NUM_RUNS times
for ((i = 1; i <= NUM_RUNS; i++)); do
    # Measure the time taken to run the program
    start_time=$(date +%s.%N) # Get start time in seconds
    ./knn.out "data/asteroids_train.csv" "data/asteroids_test.csv" "out.csv" 3
    end_time=$(date +%s.%N)   # Get end time in seconds

    # Calculate elapsed time
    elapsed_time=$(echo "$end_time - $start_time" | bc)

    # Add elapsed time to total
    total_time=$(echo "$total_time + $elapsed_time" | bc)

    echo "Run $i: $elapsed_time seconds"
done

# Calculate average time
average_time=$(echo "scale=2; $total_time / $NUM_RUNS" | bc)
echo "Average time over $NUM_RUNS runs: $average_time seconds"
