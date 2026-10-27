# Regularity conditions at the magnetic axis

This note derives the behaviour of every MARS-F/K/Q unknown at the magnetic axis
and states the conditions that `BOUNDC` and `AXISY` (`MarsQ_2FK/marsq.f`)
impose on node 1 and half-cell 1 when the radial grid reaches the axis
(`SLEFT = 0`, where `FIXORI` sets the pseudo-axis `CS(1) = 1e-3 CS(2)`).

Labels: **[V]** derived here and checked against the source and a numerical run,
**[H]** derived but not checked by an independent run.

## 1. Representation

MARS uses `(s, chi, phi)` with `s = sqrt(psi/psi_a)`, an equal-arc poloidal
angle `chi` and the geometric toroidal angle `phi`, and the Jacobian
`J = [grad s . (grad chi x grad phi)]^(-1)` (Liu et al. 2000, Sec. II). A
perturbation is `sum_m f_m(s) exp(i m chi + i n phi)`.

The variables are defined by the code itself **[V]**:

* `B1, B2, B3` and `J1, J2U, J3` are Jacobian-scaled contravariant components,
  `B^i = J b . grad x^i`. The equilibrium parts added in the field diagnostic
  are `B2 = psi'(s)` and `B3 = J F/R^2` (`marsq.f`, field reconstruction near
  `HB2 = ... + DPSIDS(I)`), and `|b|^2 = g_ij B^i B^j / J^2`.
* Displacement and velocity are
  `xi = X1 (e_s - B_s B/B^2) + X2 J (B x grad s)/B^2 + X3 B`
  (`docs/ntv_hamada_gauge.md`, Sec. 3), and `v` is expanded in the same way with
  `V1, V2, V3` (`GETXYMORE`: `X1 = V1/(lambda + i n Omega)`, likewise for `X2`, `X3`).
  `W1, W2, W3` use the same basis **[H]**.
* `J2L, B2L, B3L` are covariant components `A . e_chi`, `A . e_phi`.
* `PR, PE, PP, PPARA, PPERP, RHOP, PD, PED, DPHI, JRE` are scalars.

## 2. Local analysis

Near the axis `r ~ s` is the distance from the axis, `J = R (R_s Z_chi - R_chi Z_s) ~ R0 r r_s ~ s`,
`|grad s| ~ 1/r_s` and `|grad chi| ~ 1/r`. A smooth perturbation has smooth
Cartesian components. Hence a scalar harmonic behaves as `f_m ~ s^|m|`
(times a series in `s^2`), and for a smooth vector field `a` the polar components obey
`a_r, a_theta ~ s^(|m|-1)` for `|m| >= 1`, `a_r, a_theta ~ s` for `m = 0`, and
`a_phi ~ s^|m|`. For `|m| = 1` the axis values are related by
`a_theta = i m a_r`, which is the statement that the field at the axis is a
constant Cartesian vector. This is the standard analyticity argument for
Fourier expansions about a polar origin. In a shaped equilibrium the equal-arc
angle differs from the geometric polar angle at `O(1)`, so harmonics of the same
parity mix, but the statement that `f_m(0) = 0` for every `m != 0` and the
leading powers below remain those of the lowest harmonic present. The table
uses the circular-limit powers, which MARS itself assumes in its
interpolation exponents (`IEXV2 = -1`, `IEXJ1 = 1`, `IEXB3 = 1`).

Multiplying by the metric factors gives the leading power `p` in `f ~ s^p`:

| unknown | mesh | definition | `m = 0` | `|m| = 1` | `|m| >= 2` |
|---|---|---|---|---|---|
| `V1, X1, W1` | node | `a . grad s ~ a_r` | 1 | **0** (finite) | `|m|-1` |
| `B1` | node | `J b . grad s ~ s b_r` | 2 | 1 | `|m|` |
| `J2U` | node | `J j . grad chi ~ j_theta` | 1 | **0** (finite) | `|m|-1` |
| `J3` | node | `J j . grad phi ~ s j_phi` | 1 | 2 | `|m|+1` |
| `J2L, B2L` | node | `a . e_chi ~ s a_theta` | 2 | 1 | `|m|` |
| `B3L, PD, PED, DPHI, JRE` | node | scalar-like | **0** | 1 | `|m|` |
| `B2` | half | `J b . grad chi ~ b_theta` | 1 | **0** (finite) | `|m|-1` |
| `B3` | half | `J b . grad phi ~ s b_phi` | 1 | 2 | `|m|+1` |
| `J1` | half | `J j . grad s ~ s j_r` | 2 | 1 | `|m|` |
| `V2, X2, W2` | half | `xi_b B/(J |grad s|) ~ a_theta/s` | **0** | **-1** (1/s) | `|m|-2` (finite for 2) |
| `V3, X3, W3` | half | `a_par/B` | **0** | 1 | `|m|` |
| `PR, PE, PP, PPARA, PPERP, RHOP` | half | scalar | **0** | 1 | `|m|` |

Bold entries do not vanish at the axis. `JRE2`, `JRE3`, `JRE2L` follow `J2U`,
`J3`, `J2L`, and `JRE1` follows `J1`.

Two exact relations hold on the axis cell. The divergence of `b` in these
variables is `d_s B1 + i m B2 + i n B3 = 0` for every harmonic, so for `m != 0`

    B2 = -(d_s B1 + i n B3)/(i m),   B3 = O(s^(|m|+1)),

and for `|m| = 1` this is the regularity relation `b_theta = i m b_r` above.

## 3. What the code did before

`BOUNDC` (A. Bondeson 1991) annihilated the node-1 X rows and set `X(1) = 0`,
except `V1(1) = V1(2)` for `|m| = 1`. For the half-cell-1 Y rows the `ANNIHY` call
was commented out (already in the upstream import 8824bb1), so only the diagonal
was overwritten by 1 while every coupling of the Galerkin row was kept. Such a
row is neither the equation nor a boundary condition **[V]**. For `B2` the
original diagonal is `lambda` (the time derivative of the induction equation),
so the row returns about `lambda` times the correct `B2(1)`, and the `B3` row,
which is the divergence relation, then puts the whole `d_s B1` into `B3(1)`.
In the RMP example (`EXAMPLE/RmpLin`, plasma run) `|B3(1)|` for `m = -1` was
`3.0e6` times the power-law trend of half-cells 2 and 3 and numerically equal to
the regular `|B2(1)|`; the TC24 NCHI240 field shows the same factor `3.4e6`
**[V]**.

Keeping the plain Galerkin row on half-cell 1 is not a remedy. With the
`FIXORI` pseudo-metric it gives `|B2(1)| ~ 1e2` and `|X2(1)| ~ 1e6` against
`1e-2` and `1e2` in the Kink example **[V]**.

## 4. Conditions now imposed (grid reaching the axis)

Node 1 (X rows): `X(1) = 0` where `p > 0`, and `X(1) = X(2)` where `p = 0`.
The latter is the zero-slope condition of a function even in `s`, with an error
of `O(CS(2)^2)`. Compared with before, this adds `J2U` (and `X1`, `W1`,
`JRE2`) for `|m| = 1` and the scalar-like node unknowns for `m = 0`.

Half-cell 1 (Y rows):

* `B2`, `m != 0`: the divergence-free row
  `(B1(2) - B1(1))/CSH(1) + i m B2(1) + i n B3(1) = 0`. It converges at second
  order to `B2(CSM(1))` for a manufactured divergence-free field
  (`tests/test_axis_regularity.py`) **[V]**.
* every other component: `Y(1) = 0` by `ANNIHY` and a unit diagonal, which is
  exact at the axis for `p > 0`.
* components with `p <= 0` (`V2, X2, W2` for `|m| <= 2`, and `V3, X3, W3` and
  the scalars for `m = 0`) are rebuilt after the solve by `AXISY` as
  `Y(1) = Y(2) (CSM(1)/CSM(2))^p`. The block-tridiagonal storage gives a Y row
  access only to `X(I)`, `X(I+1)` and `Y(I)`, so this extrapolation cannot be a
  matrix row. It is consistent because half-cell 1 enters the rest of the
  system only through the node-2 X rows, which are insensitive to it: changing
  `Y(1)` by six orders of magnitude moved the outer solution by less than
  `3e-5` (Sec. 5). The reconstruction runs before `GETXY`, so it never touches
  an imported field (`KPERTREAD = 1`).

`OUTPUT_RECTRZ` used to overwrite the axis cell after the solve with
`B3U(1) = B3U(2)` and `B2U(1)` from the divergence (`B2U(1) = B2U(3)` for `m = 0`).
It changes the global arrays, so with `CTJXB = 0` (no later `GETXY`) the copy
reached `BPLASMA.OUT` and the kinetic NTV contraction. `B3U(1) = B3U(2)` contradicts
`B3 ~ s^(|m|+1)`, so this override is now skipped when the grid reaches the axis **[V]**.

An inner boundary (`SLEFT > 0`) keeps the earlier conditions unchanged. The
`n = 0, m = 0` free-boundary branch for `J2U`/`J3` (`K_BC_N0 = 20..23`) is kept
as it was. It compares `ARM` (which is 2 for `m = 0`) with zero and is therefore
never taken, which is a pre-existing issue not changed here.

## 5. Checks

* `tests/test_axis_regularity.py` compiles the production decks and checks the
  node-1 and half-cell-1 rows for `m = -3..3` against the table above, the
  unchanged inner-boundary rule, the convergence of the divergence row on a
  manufactured field and the `AXISY` power laws **[V]**.
* `EXAMPLE/Kink` at `CFBAL = 2.8` (no wall): the growth rate is
  `5.64616e-3` before and after (stored 2014 reference `5.64874e-3`). `B1..B3`
  for `s >= 0.5` and in the vacuum change by at most `6e-6` (relative L2).
  Half-cell 1 becomes regular: `B2(1)/B2(2) = 1.01` for `m = -1`, `X2 ~ 1/s`
  for `|m| = 1` and `B3(1) = 0` **[V]**.
* `EXAMPLE/RmpLin` plasma response (rotation, resistivity): `B1..B3` for
  `s >= 0.5` and in the vacuum change by at most `2.7e-5` **[V]**.
* Whether the TC24 NCHI240 own-solved field and kinetic torque behave as
  expected is the subject of a separate run (see the pull request) **[H]**.
  The kinetic side of the axis cell is treated in Sec. 6.

## 6. Kinetic drive at the axis node

The kinetic pressures and the same-drive work of a cell (`FILLMATDWKCOMP`,
`CALCDWKCOMP` in `MarsQ_2FK/kinetic.f`) evaluate the KH drive
`H = VX1 X1 + VX2 X2 + VQ1 B1 + VQ2 B2 + VQ3 B3` at the lower node, the centre and
the upper node of the cell (weights `PTRAPH/2`, `1 - PTRAPH`, `PTRAPH/2`). The
half-mesh unknowns are carried from the centre to a node with their leading
power, `X2(node) = X2(CSM) (CSM/CS)` (`IEXV2 = -1`) and `B3(node) = B3(CSM) (CS/CSM)`
(`IEXB3 = 1`), and multiplied by the coefficients of that node.

The coefficients carry the inverse powers **[V, source]**. `RX2 = -dB/dchi J T/(psi' B0)`
with `dB/dchi ~ s`, `J ~ s` and `psi' ~ s`, so `VX2 ~ s` (at least), and
`RQ3 = T/(psi' B0) ~ 1/s`. The products `VX2 X2` and `VQ3 B3` are therefore
finite at the axis, as the drive of a regular field must be.

The axis node has no coefficients of its own. `KJPCOEFF` is not run on node 1 of
the integer grid, and the cell-centre coefficients stand in for it (the copy
`VX1PARA(:,:,1) = VX1PARAM(:,:,1)` etc. in `KINETIC`, `READ_SURFACE_QUANTITIES(1,2)`
for the lower node in `CALCDWKCOMP`, `HLO = HMID` for the same-drive blocks)
**[V, source]**. Scaling `X2` to the node while keeping the centre coefficient
evaluates `VX2(CSM) X2(CSM) CSM/CS(1)`: the `1/s` growth of `X2` is applied, but the
compensating `s` of `VX2` is not. With the `FIXORI` pseudo-axis `CS(1) = 1e-3 CS(2)`
the factor is `CSM(1)/CS(1) = 500`, and the `X2-X2` part of the quadratic same-drive
work grows like `(CSM(1)/CS(1))^2`, which diverges as `CS(1) -> 0`. `B3` is
suppressed by `1/500` for the same reason. Nothing cancels this in the code. The
term is the only singular one: `X1(1)` and `B1(1)` are node values fixed by
`BOUNDC`, `B2` is not rescaled, and the centre and upper-node terms use their own
coefficients.

The consistent node drive with the stand-in coefficients is
`VX2(CS) X2(CS) = VX2(CSM) (CS/CSM) X2(CSM) (CSM/CS) = VX2(CSM) X2(CSM)` and
likewise `VQ3 B3`, so at the axis node (`SLEFT = 0`) `X2` and `B3` now enter
unscaled (`ZV2M = ZB3M = 1` for `IS = 1` in `FILLMATDWKCOMP` and in the same-drive
loop of `CALCDWKCOMP`). The remaining error of the stand-in coefficients is
`O(CSM(1))`, the same order as the rest of the three-point rule on the axis cell.

The exact axis limit of the same-drive integrand is zero. The cell work is the
integral of `psi'(s)/2 Q(s)` with `Q` the quadratic form of the drive, which is
finite for a regular field, and `psi'(0) = 0` (`DPSIDS(1) = 0`). The code uses
`DPSIDSM(1)` at the axis node, i.e. it takes the centre value of the integrand
there. This is a first-order quadrature error, not a singularity; on the TC24
NCHI240 field it changes the axis-cell ion torque by 3 % (below), and it is kept
so that the same-drive, weak and dissipative-kernel forms keep agreeing under
ideal induction.

Checks:

* `tests/test_hermitian_ntv_torque.py` compiles the production routines. On a
  grid that reaches the axis the torque of every cell, including the axis cell,
  equals the analytic same-drive work, the weak and dissipative-kernel forms agree
  with it, and it is independent of the pseudo-axis (`CS(1) = 0.05` and `1e-5`
  agree to 12 digits). With `SLEFT > 0` the node scaling is kept **[V]**.
* Offline evaluation of `DWHERM(1)` with the HCOEF blocks of the TC24 NCHI240
  supplier run and the own-solved axis-fixed fluid field (same grid): the ion
  work of the axis cell is `12.7` times that of cell 2 before (the lower-node
  `X2` term is 99.6 % of it), and `0.051` times after; with the exact limit
  `psi'(0) = 0` it would be `0.050`. The integrand of this field rises like
  `s^4..s^5` over cells 1 and 2, for which the three-point rule gives `0.04..0.05`,
  so the axis cell is on trend **[V for the evaluation, H for the kinetic run,
  whose blocks use NCHI0 = 120]**.
* The pressure rows of a self-consistent kinetic solve (`IPERTURB = 0`, `FILLMAT`
  and its variants in `kinetic.f`) use the same stand-in coefficients with the
  scaled `X2`. They are not changed here and need their own validation **[H]**.

## References

* Y. Q. Liu, A. Bondeson, C. M. Fransson, B. Lennartson and C. Breitholtz,
  Phys. Plasmas 7, 3681 (2000), Sec. II: coordinates and Jacobian of MARS-F.
* A. Bondeson, G. Vlad and H. Lütjens, Phys. Fluids B 4, 1889 (1992): the MARS
  code (bibliographic reference; the axis treatment is not described there in
  any text we could check).
* The axis analyticity of Fourier harmonics is the standard polar-coordinate
  argument (for example G. M. Lewis and P. M. Bellan, J. Math. Phys. 31, 2592
  (1990)); the table above is derived here and checked against MARS's own
  interpolation exponents and the solved power laws.
