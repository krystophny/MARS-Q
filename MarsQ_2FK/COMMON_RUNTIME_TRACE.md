# Common runtime response trace

The `WRITEKJPCOMMONTRACE` hook is default-off.  It is reached only in the
trapped `KJPCOEFF` lane, after `KCHI`, `KEQUILK`, `KBTIME`, `KPHI`, `KG`, and
`KH` have run for the current `(JS,KGRID,lambda)`.  Enable it by placing the
same selected-surface/ell rows used by the existing action traces in
`ELL_TRACE.REQUEST`; the requested integer must equal the executed
`KNTVELL` (the assembled `999` selector is intentionally rejected by the
existing selector contract).

For each selected `(JS,KGRID,ell)` the writer appends
`ELL_TRACE_JS####_G{1,2}_COMMON.OUT`.  Rows are emitted for every stability
harmonic `m` and every native trapped quadrature sample.  Interior rows expose
the exact complex phase and endpoint-subtracted G/H integrands used in
`KG`/`KH`; `endpoint_flag=-1/+1` rows carry only the analytic lower/upper
singular add-back, respectively.  Summing the two endpoint rows therefore
reproduces the single `(FL+FU)` add-back in each native factor.
`g_normalization` and `h_normalization` are the native factors that must be
applied after the quadrature sum.

The chart fields are MARS-native: `chi=RCHIK`, `phi=RPHIK`, and `tau=RTK` on
the ordered half-bounce.  `rho_pol` is `CS` on the full mesh or `CSM` on the
half mesh.  No `s_tor` or cylindrical `(R,Z)` value is fabricated; a common
comparison must join the accepted radial/physical map and retain its hash.

This packet is a producer audit, not a torque correction.  It does not include
the energy-integrated `I_ell` drive, resonance coarea/Jacobian, pressure
recovery, radial folding, toroidal covector, or final work assembly.  Those
remain separate gates in the controller `PLAN.md`.
