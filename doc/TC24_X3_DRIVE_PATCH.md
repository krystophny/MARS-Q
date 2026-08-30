# Diagnostic `gg3 X3 d|B0|/dchi` drive term for `KNTV=21`

## The candidate discrepancy

`CALCDWKCOMP` builds the kinetic drive from `X1U, X2U, B1U, B2U, B3U`. There
is no `X3` column anywhere in `kinetic.f`. The advective drive is

    dB = X1 d|B0|/ds + xi^chi d|B0|/dchi
    xi^chi = gg1 X1 + gg2 X2 + gg3 X3

so `KNTV=21` does not carry an explicit `gg3 X3 d|B0|/dchi` column. MARS's own `KNTV=11` path
includes it — `ntv_pre.f:346-350` assembles all three `gg` terms — so the two
MARS routes disagree on the drive. ARES measures the omitted term on the
converter-side control (phiI010) at core `7.07e-4`, mid `0.248` and edge
`1.157` of the full drive: at the edge the missing term is larger than the
drive it is missing from, in the region carrying about three quarters of the
volume-integrated torque.

## The coefficient, derived

`ntv_pre.f:589-592` gives the geometry factors:

    gg1 = -(dpsids/(J B))^2 g_schi
    gg2 = F / B^2
    gg3 = dpsids / J

`kinetic.f` carries the `xi^2` channel as (`KEQUILK`, both meshes)

    RW2(J) = d(B0/B)/dchi                                  ! kinetic.f:6321
    RX2(J) = -RW2(J) * RJA(JS,J) * T(JS) / DPSIDS(JS) / B0K ! kinetic.f:6376

Since `RW2 = d(B0/B)/dchi = -(B0/B^2) dB/dchi`,

    RX2 = (F/B^2) (dB/dchi) * [ B0 J / (DPSIDS B0K) ]
        = gg2 (dB/dchi) C ,   C = B0 J / (DPSIDS B0K)

The same common factor `C` gives the `xi^3` channel, and `J` and `dpsids`
cancel exactly:

    RX3 = gg3 (dB/dchi) C = (dpsids/J)(dB/dchi) B0 J/(dpsids B0K)
        = (dB/dchi) B0 / B0K
        = -RW2(J) B^2 / B0K
        = -RW2(J) / HK(JS,J,KGRID)**2 / B0K        (B0 = 1, HK = B0/B)

so the new coefficient is the existing line with `RJA*T/DPSIDS` replaced by
`1/HK**2`:

```fortran
C     COEFFICIENTS FOR XI^3 IN H-FACTOR
      IF (KGRID.EQ.1) THEN
         DO J=1,NCHI
            RX3(J)=-RW2(J)/HK(JS,J,1)**2/B0K
         ENDDO
      ELSEIF (KGRID.EQ.2) THEN
         DO J=1,NCHI
            RX3(J)=-RW2(J)/HK(JS,J,2)**2/B0K
         ENDDO
      ENDIF
      RX3(NCHI+1)=RX3(1)
```

`RX3` then splines to `RX3K` beside `RX2K` (`kinetic.f:6528`), and the pitch
weight is the one `RVAK2` already uses, because `X3` enters `xi^chi` with the
same structure as `X2`:

```fortran
      RVAK8(J) = PHASE*(2-LAM/RHK(J))*RX3K(J)
      VX3(M,L) = VX3(M,L) + EPHASE*RVAK8(J)
```

## What is still required

The coefficient above is self-contained and checkable. Wiring `X3` through to
the torque is not, and touches four files:

1. `kineticm.f` — declare `VX3, VX31, VX3PARA/PERP(+M,+C)` beside the `VX2`
   family, including the `THREADPRIVATE` list.
2. `kinetic.f` — allocate and zero them at both mesh sizes
   (`:83`, `:96`, `:191`, `:204`); add `X3PARA/X3PERP` accumulation beside
   `X2PARA/X2PERP` in every `ICASE` branch of the convolution
   (`:4521, 4542, 4585, 4604, 4653`) and the `VX3PARA` reduction at `:4753`;
   add the half-mesh copies at `:952` and `:997`.
3. `mpienv.f` — the `VX2PARA = 0.` family at `:52`/`:63` and the
   `BUFFER_DATA` slot layout at `kinetic.f:552`/`:711`. **The buffer index map
   is positional**, so adding a slot renumbers it; this is where a careless
   patch silently corrupts a different channel.
4. `specmat.inc` / `newrun.inc` — a `KYX3` column key, the `DSUBM(KYPPARA+…,
   KYX3+…)` and `DSUBM(KYPPERP+…, KYX3+…)` fills in `FILLMATDWKCOMP`
   (`kinetic.f:9311` is the `KYX2` template), the `X3U` terms in
   `CALCPRECOMP` (`:9383`, `:9391`), and the input guard at `:9608-9610`
   which currently omits `X3U`.

## Validation contract

This is a physics change to a frozen production path, so it is promoted only
against an independent oracle, never against its own plausibility:

1. With `X3` forced to zero the branch must reproduce the accepted
   `tc24-marsk-passive-keeptfun-20260827-1dee54df` `TORQUENTV.OUT`
   byte-for-byte. That isolates the wiring from the physics.
2. The resulting `Delta|B|` must match ARES's `advective_drive` — which
   already computes all three `gg` terms and is tested — on the same frozen
   inputs, before any torque is compared.
3. Only then is the torque difference attributable to `S1`.

Until 1 and 2 pass, no run from this branch is a TC24 result.

## Coefficient re-derivation (2026-08-29)

The first campaign (`tc24-marsk-x3-drive-20260829`) passed its inertness gate
— `prod_x3off` reproduced the accepted `phiI010-l5` torque byte for byte — and
`KX3DRIVE=1` then changed the radial structure exactly as the diagnosis
predicted: outer-radius sign changes fell from 4 to 1, the `q=8/3` neighbour
ratio from 59.4 to 3.08 and `q=3` from 140 to 2.21, reproduced independently
on the fast and production decks.

But the magnitude was absurd. Converted to SI against the canonical NEO-RT
MARS-field NTVTOK control, the peak `|dT/drho|` went from 8.93e4 N m (within
12% of NEO-RT's 1.02e5) to 2.45e7, a factor of 274 above NEO-RT. A correction
term cannot be 274 times the thing it corrects, so the coefficient was
re-derived from scratch rather than defended.

The original derivation asserted a common factor `B0*J/(DPSIDS*B0K)` by
inspection. That was wrong. Deriving it instead of guessing:

`HK = B0K/B` (`kinetic.f:2955`), and `RW2 = dHK/dchi`, so

```
dB/dchi = -RW2*B0K/HK**2
```

Write the known-good `RX2 = -RW2*RJA*T/DPSIDS/B0K` as `gg2*(dB/dchi)*C` with
`gg2 = F/B**2` and `T = F`, using `B**2 = B0K**2/HK**2`:

```
gg2*(dB/dchi)*C = (F*HK**2/B0K**2)*(-RW2*B0K/HK**2)*C = -F*RW2*C/B0K
```

Matching against `-RW2*F*RJA/DPSIDS/B0K` gives

```
C = RJA/DPSIDS
```

with no leftover field factor at all. That same `C` reproduces the `gg1` term
of `RX1B` exactly — `gg1 = -(DPSIDS/(RJA*B))**2*g_schi` yields
`RW2*DPSIDS*G12L/RJA/B0K`, which is the second term of `RX1B` verbatim — so
`C` is pinned by two independent coefficients, not fitted to one.

Then `gg3 = DPSIDS/RJA` cancels `C` outright:

```
RX3 = gg3*(dB/dchi)*C = (DPSIDS/RJA)*(-RW2*B0K/HK**2)*(RJA/DPSIDS)
    = -RW2*B0K/HK**2
```

The committed coefficient was `-RW2/HK**2/B0K`, whereas this derivation gives
`-RW2*B0K/HK**2`. For the TC24 deck, however, MARS prints `B0K=1.0`.
Consequently the two expressions are numerically identical and the proposed
`B0K**2` explanation does not apply to this case.

The corrected-coefficient campaign `tc24-marsk-x3-fix-20260829` completed all
four jobs. Every corrected `TORQUENTV.OUT` is byte-identical to the
corresponding first-campaign output because TC24 has `B0K=1`. In particular,
`prod_x3on` still peaks at `2.4471e7 N m`, about 240 times the NEO-RT peak.

That 240-times response was subsequently traced to a second defect in the
diagnostic patch. `KH` multiplies `VX1`, `VX2`, `VQ1`, `VQ2`, `VQ3`, and
`VDP` by the common bounce-average factor
`RCHIHK*OMEGAB/(4*pi)` after integrating the orbit. The patch omitted this
factor from `VX3`. The old result is therefore an unnormalized sensitivity
test and cannot reject the X3 channel. Branch
`fix/tc24-x3-orbit-average-normalization-20260830` applies the same factor to
`VX3`. Its exact-input A/B must pass the byte-identical X3-off gate before the
X3-on result is interpreted. The complete Park/PENTRC action identity remains
open even if that A/B improves the torque comparison.
