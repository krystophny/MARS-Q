# Prescribed-q continuation on missing inner surfaces

- `PROFILE` fills untraced inner surfaces when its first traced index `IP` is greater than1.
- It initializes axis values only for `NSTTP<=2` and `NSTTP=3`. Prescribed-q `NSTTP=4` nevertheless passes those undefined values to cubic interpolation.
- `CINT` defines the prescribed-q quantity as `CID2=1/CIDQ`, and `CID0` is the raw coarea integral. Continue raw `CID0` with the existing four-point interpolation; preserve the reciprocal identity after continuing `CIDQ`. Leave legacy selectors unchanged.

## Independent native oracle

Run from the ordinary repository root. Requires Python3 and gfortran. No external data, solver installation or BLAS library is used.

```sh
PARENT_DIR=$(mktemp -d)
git worktree add --detach "$PARENT_DIR/source" 8824bb18e1514b4a27f6357d5fe859d4fe690542
python3 tests/q_inner_coarea/run_oracle.py "$PARENT_DIR/source/CheaseMerge" CheaseMerge "$PARENT_DIR/evidence"
```

The test compiles the whole native `PROFILE`, with an independent analytic coarea supplier for circular and elongated quadratic physical flux. It deliberately selects `IP=2`. The parent fails under signaling-NaN initialization; the candidate recovers the known coarea and reciprocal identity. Native output bytes for `NSTTP=1/3` are unchanged. All12 controls must pass.

These manufactured fields qualify the profile continuation; they are not claimed finite-R force-balanced equilibria. No PDE is solved by this oracle.

## Complete generated native inputs

```sh
python3 tests/q_inner_coarea/generate_inputs.py "$PARENT_DIR/native-inputs" --variant mars
```

This creates circular and elongated tokamaks with R0=6.2m, a=0.62m, elongation1/1.7, zero pressure, prescribed q=1.5, NS32/NT32. Build with the repository's native Make configuration and run the absolute binary path in either generated directory. For example:

```sh
cd "$PARENT_DIR/native-inputs/circular"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 180 /absolute/path/to/chease < chease_namelist > native.log 2>&1
```

This branch is independent of [current-source sign PR39](https://github.com/krystophny/MARS-Q/pull/39) and [profile-primitive PR30](https://github.com/krystophny/MARS-Q/pull/30). Both are separate requirements for full inverse-q physical validation. A current-sign-only native control can converge its first solve and then stop during the profile update. An exit code0 is insufficient: require native convergence at the requested final mesh, nonempty finite equilibrium outputs, and independent q/flux/current/force checks.
