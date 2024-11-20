# K-Nearest Neighbours implemented in C

This repo is a C implementation of the K-Nearest Neighbours algorithm implemented in C using OpenMP to parallelise the solution.
Different optimisation techniques were used in order to achieve the lowest possible runtime on large data sets.

| Cores   | Runtime (s) |
| -------- | ------- |
| 1 | 24.84   |
| 2 | 12.41 |
| 4    | 6.59    |
| 8    | 3.63    |
| 16    | 2.27    |
| 32    |  1.59   |

### Prerequisites
- gcc Compiler

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/snippy4/KNN-in-openmp
   ```
2. Compile using the makefile
   ```bash
   make gccomp
   ```
3. Run the program using the following syntax
   ```bash
   ./knnomp-gcc <train data> <test data> <output file> <k value>
   ```
