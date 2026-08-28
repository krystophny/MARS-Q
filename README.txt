MarsQ: source code for MARS-F/K/Q, run "make" under Linux to compile, update makefile for proper compiler
CheaseMerge: source code for CHEASE, solve fixed boundary Grad-Shafranov equation to provide input to MARS-*
EXAMPLE: various test examples
RZplot: post-processing Matlab scripts
Manual: MARS-F/K/Q user manual




Portable local builds
=====================

The `MarsQ_2FK/makefile` provides explicit compiler targets which leave the
historical tracked `marsq.x` untouched and write executables under `build/`:

    make -C MarsQ_2FK gnu
    source /opt/intel/oneapi/setvars.sh
    make -C MarsQ_2FK ifx
    export PATH=/path/to/nvhpc/comm_libs/hpcx/bin:/path/to/nvhpc/compilers/bin:$PATH
    make -C MarsQ_2FK nvhpc

The GNU target includes the compatibility flags required by this legacy source.
The NVIDIA target intentionally uses `-O1`; nvfortran 26.5 crashes internally
while compiling `kinetic.f` at `-O2`.
The equivalent-current (`KKF=-3`) Biot-Savart kernel uses OpenMP; set
`OMP_NUM_THREADS` to the desired physical-core count at runtime.

CHEASE has out-of-tree `gnu` and `ifx` targets in `CheaseMerge/makefile`.
Both retain the zero-initialization and large-memory-model assumptions of the
legacy source; omitting zero initialization can crash before `EQDATA` is read.

The `gnu`, `ifx`, and `nvhpc` MARS build targets write
`<executable>.provenance.json` beside every executable.  The manifest records
the exact compiler wrapper and backend, versions, flags and commands, Git
commit and dirty-tree fingerprint, linked libraries, and binary SHA-256.  The
wrapper builds to a temporary path, verifies that the source tree did not
change during compilation, installs the executable atomically, and validates
the manifest/artifact pair.  For a production GNU build from a clean checkout:

```
python3 tools/build_with_provenance.py --profile gnu \
    --build-dir /mnt/storage/codex-mars/build/production --require-clean
```

Verify the artifact later without executing MARS:

```
python3 tools/build_with_provenance.py \
    --verify-manifest /path/to/marsq-gnu.x.provenance.json
```


Perturbative MARS-K NTV on a frozen MARS-F response
===================================================

Use `INCKIN=1`, `IPERTURB=1`, `KNTV=20` or `21`, and `ODWKCOM=.true.`.
The perturbative kinetic terms are evaluated from the fluid response without
feeding back into that response.  `IPERTURB=1` retains rotation in the fluid
part; `IPERTURB=2` excludes it.

NTV is quadratic in the perturbation.  For multiple coil rows, combine their
complex currents/responses before evaluating MARS-K NTV.  Do not add torque
profiles from independent row runs.  `KEYTORQ=2` cross-coupling is currently
implemented only for the Shaing-formula path (`KNTV=10/11`); MARS now rejects
that unsupported combination with `KNTV=20/21` instead of returning empty
cross terms.

To compute the MARS-K kinetic-work NTV torque on an already calculated
MARS-F perturbation, set `KPERTREAD=1` and provide `BPLASMA_INPUT` and
`XPLASMA_INPUT` in native MARS output format.  The importer runs after the
fluid solve and before the kinetic-work diagnostic, validates dimensions,
toroidal and poloidal modes, rejects trailing data, and assigns neither field
unless both inputs pass.  XPLASMA's duplicated equilibrium-profile columns
are validated; their relative mismatch against the active run equilibrium is
reported because the active equilibrium remains authoritative.  This is
supported with `INCKIN>0`, `IPERTURB>0`, and `KNTV=20/21`, or for the Shaing
`KNTV=10/11` diagnostic.  `INCFEED=4` may be used as a vacuum-only carrier:
the external import explicitly enables the NTV output path and skips unrelated
`FEEDOUT` feedback diagnostics before installing the supplied perturbation, so
no discarded plasma response has to be solved and no carrier field can leak
into NTV.  This also avoids requiring feedback-only workspaces in a frozen-field
postprocessing run.
JxB, Reynolds, and ergodic torque diagnostics are deliberately skipped in
this mode because an external B/X pair does not supply their consistent J/V.
The full post-output `ENERGYMAT` diagnostic is skipped for the same reason; it
assumes and mutates a self-consistent MHD carrier.  Before NTV contraction,
however, MARS must assemble the reciprocal pressure-to-displacement operator
blocks omitted by `IPERTURB=1` (`KPBKEY=0`).  For `KNTV=21`, the converged-
eigenvalue passive assembly must occur before `OUTPUT` forms torque for every
native or reloaded response; the later `ENERGYMAT` diagnostic is too late to
supply the torque operator even in an ordinary native `IPERTURB=0` run.  MARS
therefore performs one pre-output assembly with `KJPKEY=0,KPBKEY=1`.  A native
field remains in place, while the frozen-field path restores the validated
external B/X arrays immediately afterward.  The native path performs its
historical one-time KJP workspace release during this assembly; the validated
cache path has no live orbit workspace to release.  Neither path recomputes
the kinetic response coefficients or feeds pressure back into the MARS-F
field.
With `ODWKCOM=.true.`, MARS serializes the thread-local kinetic-work
components for every radial surface before assembling the radial NTV torque
profile.  This production path is independent of the optional surface
fraction distribution diagnostic (`INCSFD=0`).  Serialization only preserves
already computed kinetic-work coefficients for the later radial assembly; it
does not modify the imported `B1U/B2U/B3U` or `X1U/X2U/X3U` arrays, rerun the
force solve, or feed the kinetic response back into the frozen MARS-F field.
Native vector output uses 16 significant digits so an imported perturbation
can be audited component by component without an eight-digit formatting floor.

For `KNTV=21`, the final kinetic-work diagnostic reuses the DWK component
workspace retained by the master KJP assembly.  The allocator checks the
complete component map and every array shape before reuse; an inconsistent or
partially allocated workspace is a fatal error.  This avoids the historical
second-allocation failure at the end of otherwise successful kinetic runs.
