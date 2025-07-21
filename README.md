# Matrix Multiplication: Algorithm Variants & Performance Analysis

## Project description
Matrix multiplication is a key computing operation affected by memory hierarchy and CPU performance. This project examines single-core (standard, element-line, and block multiplication) and multi-core (two parallel element-line versions) implementations. Performance is evaluated using execution time, cache misses, cache hits, speedup, and effi ciency via PAPI.

#### Grade: 19/20

### 1. Standard Matrix Multiplication
The standard matrix multiplication algorithm follows a triple-nested loop structure with O(n³) complexity, multiplying each row of the first matrix by each column of the second. It was tested in C++ and Python for matrix sizes from 600×600 to 3000×3000 (increments of 400), collecting execution times to assess performance differences.

### 2. Element-by-Line Multiplication
This variant modifi es the standard approach by multiplying each element of matrix A by a full row of B, maintaining O(n³) complexity while optimizing memory access for cache efficiency.

To improve computational efficiency, two OpenMP parallel implementations were developed. Parallel Version 1 parallelizes the outer loop (matrix A’s rows), assigning each thread a full row of A to compute with all columns of B, using `#pragma omp parallel for private(i,j,k)`. Parallel Version 2 parallelizes the inner loop (matrix B’s columns) using `#pragma omp for` inside a shared parallel block.

### 3. Block Matrix Multiplication
The block matrix multiplication algorithm improves cache effi ciency by dividing matrices into bkSize × bkSize submatrices, reducing cache misses while maintaining O(n³) complexity.

## Conclusions
This study highlights the often-overlooked impact of memory management on program efficiency. While parallel computing offers significant speedups, with Parallel 1 (outer for loop version) offering the highest performance gains, particularly on smaller matrix sizes, even sequential programs can achieve substantial improvements through strategic memory utilization.

As expected, the Python implementation performs the worst, particularly for larger matrices. It serves as a baseline comparison, highlighting the inefficiencies of interpreted languages for computationally intensive tasks like matrix multiplication.
