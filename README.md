# Optimised K-Nearest Neighbours 

This repo is a C implementation of the K-Nearest Neighbours algorithm, and the K-Folds algorithm implemented in C using OpenMP for multi-threading.
Different optimisation techniques were used in order to achieve the lowest possible runtime on large data sets with multiple cores.

| Cores   | Runtime (s) |
| -------- | ------- |
| 1 | 24.84   |
| 2 | 12.41 |
| 4    | 6.59    |
| 8    | 3.63    |
| 16    | 2.27    |
| 32    |  1.59   |

benchmarks done using asteroids_train, asteroids_test with k value 3 running the knnomp version for exact K-Nearest Neighbours

### Prerequisites
- gcc compiler

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/snippy4/KNN-in-openmp
   ```
2. Compile using the command
   ```bash
   gcc -o knnomp -march=native -mavx2 -fopenmp -O3 -std=c99 knnompmain.c
   ```
3. Run the program using the following syntax
   ```bash
   ./knnomp <train data> <test data> <output file> <k value>
   ```
