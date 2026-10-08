# radial endpoint safety

Radial trapezoid request1 evaluates an existing polar-axis GS term at sigma0 and divides by zero. Stop early for that radial request while preserving general GAUSS1 and angular trapezoid semantics.

## Independent native reproducer

Run from the ordinary repository root. Requires Python3 and gfortran; No external equilibrium files or BLAS library are needed. The script compiles actual repository routines and generates all control inputs. No equilibrium solve is launched by the oracle.

```sh
PARENT_DIR=$(mktemp -d)
git worktree add --detach "$PARENT_DIR/source" 8824bb18e1514b4a27f6357d5fe859d4fe690542
python3 tests/radial-endpoint-safety/run_oracle.py "$PARENT_DIR/source/CheaseMerge" "$PARENT_DIR/parent-evidence" --mars
python3 tests/radial-endpoint-safety/run_oracle.py CheaseMerge "$PARENT_DIR/fixed-evidence" --mars --guard
```

The receipt records native source, fixture and script hashes, parent failure and repaired behavior. These tests establish the stated kernel/API behavior. They do not certify all equilibrium selectors, physical continuum errors or nonlinear convergence.
