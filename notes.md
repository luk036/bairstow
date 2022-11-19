# Polynomial root-finding for Filter design
用于 FIR 滤波器设计的多项式寻根法


## Motivation (test)

Auto-correlation polynomial
Spectral factorization

1. FFT
2. Polynomial root-finding
   1. Convert to eigenvalue problem
      Matrix: O(n^2) storage
      Stable
   2. Based on Newton's method 
      1. Aberth's method
      2. Bairstow's method

Stability
Sensitivity: Round-off error
Concurrency

LSP:

## Newton's method based algorithms

Single vs. Pacipicate

- Mathematic
  - Newton's method based
  - Local convergence

- Numerical algorithm and analysis
  - Stability
  - Round-off error
  - Sensibility

- Leja ordering

- Programming
  - Multi-threading
  - Round-robin


