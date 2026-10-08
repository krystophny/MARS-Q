# Optional folded angular ordering

CHEASE_BAND_ORDER=folded packs active bulk variables and retains a three-value
center Schur border. Native assembly/reconstruction remain unchanged; native
ordering remains default and is retained when it gives a narrower band. The
module is embedded in the existing monolithic source, preserving Make profiles
and compiler precision conventions. No external dependency is added.

The independent cell-clique known-solution oracle includes periodic seam,
center and boundary derivatives. Invalid/inactive solve guards fail closed.
The native successful-status legacy is initialized at call sites; separate
factorization PR32 repairs the routine contract and supplies tiled updates.
