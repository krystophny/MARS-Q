# DWK contraction guard: probe every surface, not just NR/2

`CALCDWKCOMP` refuses to build the NTV torque from zeroed inputs:

```fortran
IF (CACHEMAX.LE.0.0.OR.FIELDMAX.LE.0.0.OR.
&    OPPARAMAX.LE.0.0.OR.OPPERPMAX.LE.0.0.OR.
&    MAX(PPARAMAX,PPERPMAX).LE.0.0)
&   STOP 'INVALID ZERO DWK CONTRACTION INPUT'
```

`FIELDMAX`, `OPPARAMAX`, `OPPERPMAX`, `PPARAMAX` and `PPERPMAX` are all
maxima over the whole radius. `CACHEMAX` was not: it sampled a single
hard-coded surface,

```fortran
IF (IS.EQ.MAX(1,NR/2)) CACHEMAX = MAX(...)
```

`DELRATS` is a rational-surface exclusion width, so for `DELRATS > 0` the
surface `NR/2` can legitimately fall inside an excluded band and cache exactly
zero while every other surface is healthy. That aborted a valid run in
`tc24-marsk-discriminator-scan-20260829`, with the diagnostic line showing a
zero cache beside parallel and perpendicular pressure maxima byte-identical to
the healthy baseline.

This patch makes `CACHEMAX` a running maximum over the same `IS` loop, so it
matches the other five probes. A genuine failure zeroes the cache on every
surface — as the large-`NEPK` failures in the same campaign do, where the
pressures are zero as well — so the guard still fires on those.

The change is diagnostic-only. `CACHEMAX` is used solely to decide whether to
`STOP`; it enters no computed quantity, so no accepted torque number can move.
A rerun under this build stays comparable to the accepted campaign.
