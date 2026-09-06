# MARS regression tests

The fast standard-library test suite covers the complete local delta from the
upstream baseline `8824bb1`:

| Change | Fast coverage | System coverage |
|---|---|---|
| portable namelist errors | source contract and black-box malformed input | compiler matrix |
| strict `INPUT_BNM` header | parser contract | MAST-U `KKF=-3` run |
| OpenMP equivalent-current kernel | scoping and arithmetic contract | 1-thread/parallel `CURHARMO` comparison |
| GNU/ifx/NVHPC targets | target/flag contract | clean build matrix |
| NVHPC `-O1` workaround | exact flag contract | NVHPC build |
| CHEASE legacy initialization | GNU/ifx flag contract | ITER EQDSK-to-OUTRMAR run |
| perturbative MARS-K NTV | input contracts and unsupported cross-term rejection | paired MARS-F/MARS-K field comparison |
| external frozen B/X import | strict/atomic parser contract and invalid-mode rejection | exact BPLASMA/XPLASMA round trip plus MARS-K torque |
| perturbative passive-operator ordering | every final-sweep `KNTV=21` field, including ordinary native `IPERTURB=0`, assembles passive `KPBKEY=1` blocks once before torque; only imported B/X is restored | paired native/reload B/X and DWK equality plus matching nonzero `TORQUENTV.OUT` |
| MARS-K matrix call contract | all `CALCDWKCOMP` callers pass eight live matrices | debug cached continuation reaches all DWK and torque outputs |
| validated DWK cache recovery | narrow namelist contract and KJP bypass | failed-after-cache ITER runs resume without coefficient regeneration |
| deterministic frequency diagnostic | `RUU2` is thread-private with its per-surface consumers | repeated OpenMP `FREQUENCIES.OUT` comparison |
| two-species frequency diagnostic | hot-ion `SLAM0(:,3)` is read only when species three exists | GNU `-fcheck=all` MAST-U run |
| retained DWK component workspace | idempotent allocation, component-map, and shape contracts | KNTV=21 reaches `TORQUENTV.OUT` after kinetic assembly |
| default-off trapped selected-`ell` trace | explicit `JS KGRID ELL` request gate, executed `KIA_TRAP` path, requested/served selector check, selected-surface replay, and cache-write suppression | accepted full-harmonic cache replay reproduces torque while emitting only requested pitch-energy traces |
| independent DWK density check | numerical angular-quadrature oracle and source normalization contract | opt-in direct check closes against the component density without an extra `4*pi^2` factor |
| five-term DWK drive ledger | independent complex linearity/sign oracle and request-file contract | opt-in `X1`, `X2`, `B1`, `B2`, `B3` sum closes against the unchanged pre-edge production density |
| 2x5 DWK bilinear ledger | independent integer/half-mesh combination oracle and request-file contract | each pressure drive is separated into X1- and X2-work rows whose ten terms close against the unchanged pre-edge production density |
| executable build provenance | profile/flag and manifest hash/tamper contracts | clean-tree build plus independent manifest verification |

Run the tests with:

```sh
make -C MarsQ_2FK test
```

The source-contract tests need no compiler and always run.  The black-box
executable tests need a built MARS, which is resolved in this order:

1. `MARS_EXE`, when set.  It always wins, and a path that does not exist is a
   hard error rather than a silent fallback, so the explicit workflow below
   still fails loudly on a typo.
2. otherwise the in-tree default build `build/marsq-gnu.x`, relative to the
   repository root, when it exists, is executable, and its provenance says it
   was built from a commit this checkout contains, with a clean worktree.
3. otherwise the executable tests skip, naming which of those it was.

So a working tree that has been built from itself runs the runtime tier by
default; a tree with no build, or with one left behind by another branch, falls
back to the source-contract tier and says so.

Point the runtime tier at a different executable with:

```sh
make -C MarsQ_2FK test-runtime MARS_EXE="$PWD/build/marsq-ifx.x"
```

The default build is whatever was last written to `build/`, and a build
directory outlives the branch it was made on.  That comparison is no longer
left to the reader: the runtime tier reads
`build/marsq-gnu.x.provenance.json` and refuses a default build whose
`source.commit` is not an ancestor of `HEAD`, or that was made from a modified
worktree, skipping with the recorded commit and branch in the message.  Rebuild
with `python3 tools/build_with_provenance.py --profile gnu` and the tier
re-enables itself.

An explicit `MARS_EXE` is never subject to that check: choosing a binary by
hand is the caller's decision.

The long MAST-U and ITER system runs use external/private fixtures and are
recorded in the consuming project rather than copied into this public source
repository.
