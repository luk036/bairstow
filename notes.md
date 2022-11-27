## FIR Filter Design
FIR 滤波器设计

input: specification
output: filter coefficients (or a Verilog file)
method: Convex optimization via Spectral factorization

Convex optimization:
    Ellipsoid method + parallel cuts
Spectral factorization: 
    FFT
    Polynomial root-finding 
        Auto-correlation function
        Parallel Bairstow's method
        Aberth's method
        [ ] Leja ordering

---

## Multiplier-less FIR Filter Design

low-cost low-power
Discrete version of Ellipsoid method
Canonical Signed Digit (CSD) 
input: number of non-zeros (nnz)
output: filter coefficients in CSD form
        [ ] Common Sub-expression Extraction/Elimination
        [ ] a Verilog file

---

## Possible contributions

Write more test cases
C++ Porting
Code Clean-up
