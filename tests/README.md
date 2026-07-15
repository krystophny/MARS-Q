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
| retained DWK component workspace | idempotent allocation, component-map, and shape contracts | KNTV=21 reaches `TORQUENTV.OUT` after kinetic assembly |

Run fast tests with:

```sh
make -C MarsQ_2FK test
```

Add executable-level tests with:

```sh
make -C MarsQ_2FK test-runtime MARS_EXE="$PWD/build/marsq-gnu.x"
```

The long MAST-U and ITER system runs use external/private fixtures and are
recorded in the consuming project rather than copied into this public source
repository.
