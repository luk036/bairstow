[![Documentation Status](https://readthedocs.org/projects/bairstow/badge/?version=latest)](https://bairstow.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/luk036/bairstow/branch/main/graph/badge.svg?token=hbjnTDpNeg)](https://codecov.io/gh/luk036/bairstow)
[![Coverage Status](https://coveralls.io/repos/github/luk036/bairstow/badge.svg?branch=main)](https://coveralls.io/github/luk036/bairstow?branch=main)

# bairstow

> Bairstow's Rootfinding Algorithm Python Code

The Aberth-Ehrlich (AE) method introduced by Ehrlich (1967); Aberth (1973) combines the elements of Newton’s method with an implicit deflation strategy, which allows for the computation of all roots of a polynomial simultaneously and converges cubically. This method is considered an improvement over the Durand-Kerner method, another simultaneous root solver method which converges quadratically and is 10-100 times slower (see, e.g., Ghidouche et al. 2017). The facts that AE is extremely fast for various degrees of polynomials, its ability to find all the roots at once (unlike other iterative root-solving methods such as Laguerre’s) and its root polishing procedure, which is inside the main iteration loop and can be controlled for a chosen accuracy.

Bairstow's method is an iterative method that is used to find complex roots of a polynomial. This method is based on synthetic division and can be used to find all roots of a polynomial.

Parallel Bairstow's method refers to a parallel algorithm based on Bairstow's method. This algorithm uses parallel computation to accelerate the process of finding complex roots of polynomials.

## ✨ Features

- Pure Python. No `numpy` required
- O(N) storage requirement, where N = degree of polynomial.
- Preserve structure, e.g. auto correlation function

## Used by

- [multiplierless](https://github.com/luk036/multiplierless)
