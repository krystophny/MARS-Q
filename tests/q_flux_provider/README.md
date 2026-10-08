# Prescribed-q finite source and primitive

The `NSTTP=4`, `NPROFZ=0` current uses a cubic in rho and divides its derivative
by rho. That cubic can have a nonzero derivative at rho=0, giving an infinite
axis/under-axis current. The primitive is obtained by a separate interpolation
and integration, so it is not the same function whose derivative drives current.

The candidate constructs one cubic F² in x=rho² from the native q/coarea knots.
Its derivative gives FFprime, its value gives the primitive at profile knots, and the
current uses jphi=-R*pprime-FFprime/R. It preserves knot values and changes the
numerical interpolation between knots. This is an explicit source discretization
change; it does not establish full inverse-q convergence or admit nonsmooth
physical q laws. Native CHECK rejects NPROFZ=1 before this selector runs;
that unsupported combination is not a claimed preserved equilibrium lane.

The direct QFLUXSOURCE value/derivative API is consistent. Downstream PREMAP,
MAPPIN and GUESS still interpolate TMF values in rho, and ISOFUN retains a
rho spline of TMF for other consumers. Those continuous field reconstructions
are not the same function as sqrt(F²(x)). No exact off-knot native field or
export derivative/primitive consistency is claimed by this patch.

Requirements: Python3, NumPy, SciPy, gfortran. No external equilibrium data.
The solver gains no dependency. From an ordinary clone with this branch:

```sh
git worktree add /tmp/mars-q-parent 8824bb18e1514b4a27f6357d5fe859d4fe690542
python tests/q_flux_provider/check.py /tmp/mars-q-parent/CheaseMerge CheaseMerge /tmp/mars-q-source-oracle
python tests/q_flux_provider/check_full_profile.py /tmp/mars-q-parent/CheaseMerge CheaseMerge /tmp/mars-q-profile-oracle --mesh-axis
```

- The first command compiles whole native ISOFUN/CURENT and native spline
  routines. Twenty-four parent cases fail axis/under-axis finiteness;
  twenty-four repaired cases recover exact linear/quadratic/cubic F² and
  current for both q signs and finite pressure derivatives. Four unselected
  selector1/3 packets remain byte-identical.
- The second command adds whole PROFILE and exact analytic shifted
  circular/elliptical coarea. Thirty-two cases cover constant/linear-in-x q,
  signed q and finite pressure derivatives. The mesh origin is at the axis
  to isolate source interpolation from missing-inner contour tracing.
- These are mathematical source controls, not force-balanced equilibria.
  F² tests qualify magnitudes; they do not independently qualify the global
  signed-field convention for negative-q equilibria.

The generator supplies complete circular and elongated elliptical native inputs:

```sh
python tests/q_flux_provider/generate_inputs.py /tmp/mars-q-native-cases --variant mars
make -B -C CheaseMerge -j2 F95=gfortran F95FLAGS='-O3 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch -mcmodel=medium' LDFLAGS=-mcmodel=medium
```

The remaining full-iteration gap has a separate complete reproducer. It creates
a detached parent worktree, downloads public review PRs30/39/41/42, applies
only their source patches, generates both inputs and builds with explicit
gfortran flags. It excludes this unrun source-provider candidate:

```sh
bash tests/q_flux_provider/reproduce_q4_iteration.sh /tmp/mars-q4-prepared --prepare-only
bash tests/q_flux_provider/reproduce_q4_iteration.sh /tmp/mars-q4-executed
```

The first command tests patch assembly, inputs and build without executing a
solver. The second uses one thread and a180-second timeout per case, retains
exit status/logs and checks requested final-grid, nonempty finite NOUT. Previously
observed numerical values belong to a separately sealed integration build;
the minimal assembled build is not claimed to reproduce identical numbers.

To run this candidate itself later, first apply the separate primitive,
current-sign, SMOOTH-ordering and missing-inner-contour fixes (review
PRs30/39/41/42), then use its absolute executable path:

```sh
cd /tmp/mars-q-native-cases/circular
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 180 /absolute/clone/CheaseMerge/chease.x < chease_namelist > chease.log 2>&1
cd /tmp/mars-q-native-cases/shaped
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 180 /absolute/clone/CheaseMerge/chease.x < chease_namelist > chease.log 2>&1
```

A native STOP may return0; require the requested final mesh,
nonempty finite equilibrium files and independent q/flux/current/force checks.
The same generated circular problem still fails the full native iteration after
the contour/ordering fixes. This provider has not been tested as a full native
equilibrium repair; its finite source contract is the scope of this PR.

An optional read-only conditioning control requires mpmath:

```sh
python tests/q_flux_provider/diagnose_chart.py /tmp/mars-q-chart.json
```

It interpolates an exact shifted physical quadratic with analytically consistent
native chart jets. The exact infinitesimal coarea source is FFprime=-45; the
cubic chart interpolant gives about+3.60e6 at its local axis. Ordinary angular
refinement improves the Hessian but does not converge that fourth-derivative
limit. This is a representation/ordering-of-limits test, not an equilibrium or
a claim that all finite-level contour derivatives diverge.
