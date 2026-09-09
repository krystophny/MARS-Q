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
half mesh.  Schema v2 additionally emits `dpsids`, `hchi`, the normalized
signed state `vpar_state`, `orientation_state`, and physical
`orientation_vpar`.  The source identity is
`hchi = B.grad(chi)/B = DPSIDS/(J*B) = dpsids/jb`; on the native increasing-
`chi` leg, `vpar_state = sign(hchi)*(v_parallel/v)` and
`orientation_vpar = sign(vpar_state*hchi)`.  Exact turning-point rows carry
zero state/orientation; their open-leg limits are represented by interior
rows.  No `s_tor` or cylindrical `(R,Z)` value is fabricated; a common
comparison must join the accepted radial/physical map and retain its hash.

The same request also writes `ELL_TRACE_JS####_G{1,2}_ORBIT.OUT` in the
accepted zero-FOW trapped lane.  This is an ordered source-native full bounce:
the first leg is the exact `RCHIK/RPHIK/RTK` lower-to-upper sequence, and the
second leg reverses the open samples while time and `ell*bounce_angle`
continue to `2*RTK(end)` and `2*pi*ell`.  The upper turning point is emitted
and owned once.  The last lower-turn row is explicitly an unowned periodic
closure of the owned initial row.  `delta_tau_measure` is the positive
native-`RTK` measure of the preceding cell and sums to `2*RTK(end)`.
`vpar_state` and both orientation fields change sign on the reflected leg.

This exact reflection does not provide a finite-orbit toroidal position or a
leg-resolved complex amplitude.  Native `KG` and `KH` have already combined
the two velocity signs into `cos(ell*OMEGAB*RTK)` on the half bounce, so their
individual complex leg amplitudes cannot be recovered algebraically from the
folded factor.  Those remain required external/source extensions before a
common physical phase or work map can consume this geometry record.

This packet is a producer audit, not a torque correction.  It does not include
the energy-integrated `I_ell` drive, resonance coarea/Jacobian, pressure
recovery, radial folding, toroidal covector, or final work assembly.  Those
remain separate gates in the controller `PLAN.md`.
