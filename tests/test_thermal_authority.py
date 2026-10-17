"""Dimensional thermal inputs through native KEQPROF and DFFFDPSI.

Manufactured pressures are affine in s**2, so the native interpolation and
differentiation have exact independent gradients, including the full-mesh
edge.  Compile the committed legacy routine too: opting out must retain its
numerical behavior, including its failure to honor the three species temperatures.
"""
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(os.environ.get("MARSQ_ROOT") or Path(__file__).resolve().parents[1])
NATIVE = ROOT / "MarsQ_2FK"
BASELINE = "a602c861711e14a310ea22d352e29a1e6b71bdd2"
NR = 12
KAPPA = 1.7
EV_J = 1.602176634e-19
MU0 = 4.e-7 * math.pi


def routine(source, name):
    match = re.search(
        rf"^ +SUBROUTINE +{name}\b.*?^ +END[ \t]*$",
        source, re.MULTILINE | re.DOTALL | re.IGNORECASE,
    )
    if match is None:
        raise ValueError("Missing production routine " + name)
    return match.group() + "\n"


DRIVER = r"""
program thermal_authority
    use dimensim
    use globalm
    use rcomdm
    implicit none
    integer :: mode, shape, baseline, j, grid, k, readbacks
    real*8 :: amplitude, s, x, density, carbon
    character(len=32) :: arg
    common /thermal_test_readback/ readbacks

    call get_command_argument(1, arg)
    read(arg, *) nspecies
    call get_command_argument(2, arg)
    read(arg, *) mode
    call get_command_argument(3, arg)
    read(arg, *) shape
    call get_command_argument(4, arg)
    read(arg, *) amplitude
    call get_command_argument(5, arg)
    read(arg, *) baseline
    call get_command_argument(6, arg)
    read(arg, *) fracpth
    nr = 12
    nrp1 = nr + 1
    kproftauth = mode
    nexpv = 1
    nprofie = 4
    ispecies_f0(1:nspecies) = 0
    ifow = 0
    nksmoothr = 0
    zne0 = 3.2d19 * amplitude
    zti0 = 2100.d0 * amplitude
    zte0 = 1300.d0 * (0.75d0 + 0.25d0 * amplitude)
    b0exp = 2.5d0
    alphap = 0.5d0
    especies_z(1:3) = [1.d0, -1.d0, 6.d0]
    readbacks = 0
    allocate(cs(nrp1), csm(nrp1), rho(nrp1), rhom(nrp1), &
             peq(nrp1), peqm(nrp1), tempi(nrp1), tempim(nrp1), &
             tempe(nrp1), tempem(nrp1), dpsids(nrp1), dpsidsm(nrp1), &
             fff(nrp1,2), dfff(nrp1,2))
    allocate(especies_denf(nrp1,2,nspecies), &
             especies_pref(nrp1,2,nspecies), especies_den(nrp1,2,nspecies), &
             especies_pre(nrp1,2,nspecies), especies_tem(nrp1,2,nspecies), &
             especies_prep(nrp1,2,nspecies))
    do j = 1, nrp1
        cs(j) = (real(j-1,8)/nr)**1.3d0
    end do
    csm(1:nr) = (cs(1:nr)+cs(2:nrp1))/2.d0
    csm(nrp1) = cs(nrp1)
    especies_denf = 0.d0
    especies_pref = 0.d0
    dfff = -999.d0
    if (nspecies == 3) especies_pref(:,:,3) = 0.15d0
    do grid = 1, 2
        do j = 1, nrp1
            s = cs(j)
            if (grid == 2) s = csm(j)
            x = s*s
            density = 1.d0 - 0.25d0*x
            carbon = 1.d0/30.d0
            if (shape == 1) then
                density = 1.d0
                carbon = carbon + x/60.d0
            end if
            if (nspecies == 3) especies_denf(j,grid,3) = carbon
            if (grid == 1) then
                rho(j) = density
                peq(j) = 0.018d0*(1.d0-0.31d0*x)+0.002d0*x*x
                tempi(j) = (1.d0-0.4d0*x)/density
                tempe(j) = (1.d0-0.22d0*x)/density
                if (shape == 1) then
                    tempi(j) = 1.d0
                    tempe(j) = 1.d0+0.3d0*x
                end if
                dpsids(j) = 2.d0*1.7d0*s
            else
                rhom(j) = density
                peqm(j) = 0.018d0*(1.d0-0.31d0*x)+0.002d0*x*x
                tempim(j) = (1.d0-0.4d0*x)/density
                tempem(j) = (1.d0-0.22d0*x)/density
                if (shape == 1) then
                    tempim(j) = 1.d0
                    tempem(j) = 1.d0+0.3d0*x
                end if
                dpsidsm(j) = 2.d0*1.7d0*s
            end if
        end do
    end do
    if (baseline == 1) then
        call keqprof_baseline
    else
        call keqprof
    end if
    write(*,'(A,2ES26.17,I5)') 'META ', alphap, fracpth, readbacks
    do grid = 1, 2
        do j = 1, nrp1
            if (grid == 1) then
                write(*,'(A,2I5,ES26.17)') 'EQ ', grid, j, peq(j)
            else
                write(*,'(A,2I5,ES26.17)') 'EQ ', grid, j, peqm(j)
            end if
            do k = 1, nspecies
                write(*,'(A,3I5,5ES26.17)') 'ROW ', grid, j, k, &
                    especies_den(j,grid,k), especies_tem(j,grid,k), &
                    especies_pre(j,grid,k), especies_prep(j,grid,k), &
                    especies_pref(j,grid,k)
            end do
        end do
    end do
end program

subroutine kthermreadback
    implicit none
    integer :: readbacks
    common /thermal_test_readback/ readbacks
    readbacks = readbacks + 1
end subroutine
"""


class ThermalAuthorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.work = Path(cls.tmp.name)
        source = (NATIVE / "kinetic.f").read_text()
        baseline = subprocess.run(
            ["git", "show", BASELINE + ":MarsQ_2FK/kinetic.f"],
            cwd=ROOT, capture_output=True, text=True, timeout=30,
        )
        if baseline.returncode:
            raise RuntimeError("Missing committed legacy oracle: " + baseline.stderr)
        legacy = routine(baseline.stdout, "KEQPROF").replace(
            "SUBROUTINE KEQPROF", "SUBROUTINE KEQPROF_BASELINE", 1,
        )
        (cls.work / "kernels.f").write_text(
            routine(source, "KEQPROF") + legacy + routine(source, "DFFFDPSI"),
        )
        (cls.work / "driver.f90").write_text(DRIVER)
        cls.exe = cls.work / "thermal_authority"
        command = [
            os.environ.get("FC", "gfortran"), "-O0", "-fcheck=all",
            "-ffpe-trap=invalid,zero,overflow", "-std=legacy",
            "-fdefault-real-8", "-fdefault-double-8",
            "-ffixed-line-length-none", "-ffree-line-length-none",
            "-I" + str(NATIVE),
            *[str(NATIVE / name) for name in ("dimensim.f", "globalm.f", "rcomdm.f")],
            "kernels.f", str(NATIVE / "splines.f"), "driver.f90", "-o", str(cls.exe),
        ]
        compiled = subprocess.run(
            command, cwd=cls.work, capture_output=True, text=True, timeout=120,
        )
        if compiled.returncode:
            raise RuntimeError(compiled.stdout + compiled.stderr)

    def run_native(self, species, mode, shape=0, amplitude=1., *,
                   baseline=False, thermal_fraction=.73):
        result = subprocess.run(
            [str(self.exe), str(species), str(mode), str(shape), str(amplitude),
             str(int(baseline)), str(thermal_fraction)],
            cwd=self.work, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows, equilibrium, meta = {}, {}, None
        for line in result.stdout.splitlines():
            fields = line.split()
            if fields[0] == "ROW":
                rows[tuple(map(int, fields[1:4]))] = tuple(map(float, fields[4:]))
            elif fields[0] == "EQ":
                equilibrium[tuple(map(int, fields[1:3]))] = float(fields[3])
            elif fields[0] == "META":
                meta = tuple(map(float, fields[1:]))
        self.assertEqual(len(rows), 2*(NR+1)*species, result.stdout)
        self.assertEqual(len(equilibrium), 2*(NR+1), result.stdout)
        self.assertIsNotNone(meta, result.stdout)
        self.assertTrue(all(math.isfinite(v) for row in rows.values() for v in row))
        return rows, equilibrium, meta

    def close(self, actual, expected, label):
        self.assertAlmostEqual(
            actual, expected, delta=3.e-11*abs(expected)+3.e-18, msg=label,
        )

    @staticmethod
    def point(grid, j):
        nodes = [(i/NR)**1.3 for i in range(NR+1)]
        if grid == 1 or j == NR+1:
            return nodes[j-1]
        return (nodes[j-1]+nodes[j])/2.

    def check_dimensional(self, species, shape, amplitude):
        rows, equilibrium, meta = self.run_native(species, 2, shape, amplitude)
        unit = 3.2e19*amplitude*EV_J*MU0/2.5**2
        ti0 = 2100.*amplitude
        te0 = 1300.*(.75+.25*amplitude)
        for (grid, j), value in equilibrium.items():
            x = self.point(grid, j)**2
            self.close(value, .018*(1.-.31*x)+.002*x*x, "PEQ changed")
            density = 1.-.25*x if shape == 0 else 1.
            carbon = 1./30. + (x/60. if shape == 1 else 0.)
            fractions = [1., 1.] if species == 2 else [1.-6.*carbon, 1., carbon]
            temperatures = [ti0*(1.-.4*x)/density, te0*(1.-.22*x)/density]
            if shape == 1:
                temperatures = [ti0, te0*(1.+.3*x)]
            if species == 3:
                temperatures.append(temperatures[0])
            pressures = [density*f*t*unit for f, t in zip(fractions, temperatures)]
            charge = 0.
            for k in range(1, species+1):
                n, temperature, pressure, gradient, fraction = rows[grid, j, k]
                label = str((species, shape, amplitude, grid, j, k))
                self.close(n, density*fractions[k-1], "density " + label)
                self.close(temperature/unit, temperatures[k-1], "eV " + label)
                self.close(pressure, pressures[k-1], "pressure " + label)
                self.close(fraction, pressure/sum(pressures), "fraction " + label)
                charge += n*[1., -1., 6.][k-1]
                if grid == 1 or j <= NR:
                    if shape == 0:
                        slope = [-.4*ti0*fractions[0], -.22*te0]
                        if species == 3:
                            slope.append(-.4*ti0/30.)
                    else:
                        slope = [0. if species == 2 else -.1*ti0, .3*te0]
                        if species == 3:
                            slope.append(ti0/60.)
                    exact = slope[k-1]*unit/KAPPA
                    # Differencing a constant pressure has cancellation error
                    # on the uniform psi grid, even though its exact slope is zero.
                    roundoff = 64.*math.ulp(max(pressures))*NR/KAPPA
                    self.assertAlmostEqual(
                        gradient, exact, delta=3.e-11*abs(exact)+roundoff,
                        msg="gradient " + label,
                    )
            self.assertLess(abs(charge), 3.e-15, "charge neutrality")
        axis = [rows[1, 1, k][2] for k in range(1, species+1)]
        self.close(meta[0], 1.-axis[1]/sum(axis), "aggregate ion fraction")
        self.close(meta[1], sum(axis)/equilibrium[1, 1], "thermal PEQ fraction")
        self.assertEqual(meta[2], 1., "active branch must reach readback")
        self.assertNotEqual(rows[1, NR+1, 2][1], rows[1, NR, 2][1], "edge flattened")

    def test_dimensional_two_species_amplitudes_edges_and_gradients(self):
        for shape in (0, 1):
            for amplitude in (1., 1.7):
                with self.subTest(shape=shape, amplitude=amplitude):
                    self.check_dimensional(2, shape, amplitude)

    def test_dimensional_three_species_amplitudes_edges_and_gradients(self):
        for shape in (0, 1):
            for amplitude in (1., 1.7):
                with self.subTest(shape=shape, amplitude=amplitude):
                    self.check_dimensional(3, shape, amplitude)

    def test_legacy_two_species_matches_committed_native_baseline(self):
        for shape in (0, 1):
            with self.subTest(shape=shape):
                live = self.run_native(2, 0, shape)
                self.assertEqual(live, self.run_native(2, 0, shape, baseline=True))
                self.assertEqual(live[2][2], 0.)

    def test_accepted_mode_one_matches_committed_two_species_baseline(self):
        for shape in (0, 1):
            for fraction in (1., .73):
                with self.subTest(shape=shape, thermal_fraction=fraction):
                    live = self.run_native(2, 1, shape, thermal_fraction=fraction)
                    baseline = self.run_native(
                        2, 1, shape, baseline=True, thermal_fraction=fraction,
                    )
                    self.assertEqual(live, baseline)
                    self.assertEqual(live[2][2], 0.)

    def test_legacy_three_species_is_an_unmatched_temperature_control(self):
        legacy = self.run_native(3, 0)
        self.assertEqual(legacy, self.run_native(3, 0, baseline=True))
        active = self.run_native(3, 2)
        unit = 3.2e19*EV_J*MU0/2.5**2
        for k, supplied in ((1, 2100.), (2, 1300.), (3, 2100.)):
            self.close(active[0][1, 1, k][1]/unit, supplied, "active axis eV")
            self.assertGreater(abs(legacy[0][1, 1, k][1]/unit/supplied-1.), .05)
        self.assertEqual(legacy[2][2], 0.)


if __name__ == "__main__":
    unittest.main()
