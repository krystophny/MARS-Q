# MAPPIN K1 axis-field interpolation

The harness extracts the actual two-line T0 assignment from native MAPPIN and
executes it with the unchanged native CUCCCC statement function. It compares
against exact constant, even quadratic, cubic and reversed-polarity polynomial
axis values on two nonuniform radial grids. Parent passes only2/8 controls;
CSM(4) passes8/8. This is an isolated native-call oracle, not a complete MAPPIN
or downstream equilibrium test.

Affected assignment requires NSURF!=1 and K1 mapping: NIDEAL1/2 regardless NRFP,
or NIDEAL6 with NRFP0. Constant-F cases hide the defect. The mapped F-axis feeds
GEQDSK and JSOLVER; native q-axis is computed separately, and the pre-solve
normalization precedes this assignment. No PDE/psi impact is established.

Run CMake/CTest with Debug or Release; CHEASE_SOURCE can select a frozen parent
without editing the test. Downstream MARS Pair A/B delivery remains open.
