# Native contiguous band factorization controls

Configure this standalone CMake fixture in Debug and Release and run CTest.
The fixture compiles the owning solver's actual factorization, triangular
solves and embedded BLAS routines. No external BLAS or new dependency is used.
A retained pre-change factorization supplies a supplementary bitwise check.
The independent oracle is a dense Gaussian solve with row pivoting plus the
original dense matrix times the native solution. It exercises positive and
indefinite matrices, bandwidth exceeding dimension and padded storage.
A separate control covers empty, singleton and exact-zero singular pivots.

`band_timing N M` compares CPU times for the same manufactured diagonally
dominant band matrix. It requires identical factors; it is a kernel benchmark,
not an equilibrium run or an accuracy-matched cross-solver speed comparison.
Use the same compiler flags and hardware for both routines. The native global
bandwidth, pivot order and singularity threshold for nonzero pivots are retained.

The loop saves the original pivot row, normalizes it, and visits the triangular
update by contiguous destination rows instead of stride MP-1. Every entry gets
the same product at every pivot. The scratch row is O(M). Singular NaN inputs
are outside this test's contract; finite factors are compared exactly.

Trailing32-pivot updates now use128-entry cache tiles in original pivot order.
The near-threshold oracle covers accepted/rejected pivots across the panel
boundary. Default precision/build conventions are unchanged; no dependency.
