# Native profile bracket search

PPSPLN previously scanned every knot across every query, recomputing its
radius before each unreached knot. Scan queries first and calculate the
radius once for each bracket search; stop at the same first qualifying knot.
Preserve clamping, knot ties, endpoint/exterior interval selection and all
coefficient arithmetic. The evaluation stage retains its native radius
calculation. No sorted-query assumption or new dependency is introduced.

CMake/CTest Debug and Release compare values/derivatives against an exact cubic
on nonuniform knots, both signed axis fluxes, axis clipping, knot neighbourhoods
and exterior queries. Outputs are byte-equal to the frozen native parent in
48 controls. Public derivatives here are with respect to rho (KOPT2); the
special KOPT1 axis convention is not an independent smooth-axis-limit test.

The optional `profile_oracle timing` performs200000 native calls with32 queries
and64 intervals, alternating parent/current three times. Ryzen9 single-core,
Release native-vector flags with FMA disabled: parent medians0.877/0.877s,
current0.159/0.157s (public/MARS). Full-equilibrium controls must assess the
combined stack separately. Tests use native precision/includes and keep the
project Make build authoritative.
