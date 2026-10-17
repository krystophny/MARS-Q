"""Native diamagnetic frequencies for a constant-density, linear-pressure plasma.

For p(psi)=p0+b*psi and constant density, both diamagnetic frequencies
are constant, have opposite signs, and have zero radial derivatives.
This independently fixes their boundary values on any radial mesh.
"""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


ROOT = Path(os.environ.get("MARSQ_ROOT") or Path(__file__).resolve().parents[1])

DECLARATIONS = """      SUBROUTINE NATIVE_DIAMAG(NRP1,NR,CS,CSM,RHO,RHOM,
     & PEQ,PEQM,PPEQ,PPEQM,DPSIDS,DPSIDSM,B0K,OMEGACI0,
     & ALPHAP,ETASTI,ETASTE,OMEGASI,OMEGASE,OMEGASIM,OMEGASEM,
     & DOMEGASI,DOMEGASE,DOMEGASIM,DOMEGASEM,DLNRHO,DLNRHOM)
      IMPLICIT NONE
      INTEGER NRP1,NR,J
      REAL CS(NRP1),CSM(NRP1),RHO(NRP1),RHOM(NRP1),
     & PEQ(NRP1),PEQM(NRP1),PPEQ(NRP1),PPEQM(NRP1),
     & DPSIDS(NRP1),DPSIDSM(NRP1),B0K,OMEGACI0,ALPHAP,
     & ETASTI,ETASTE,OMEGASI(NRP1),OMEGASE(NRP1),
     & OMEGASIM(NRP1),OMEGASEM(NRP1),DOMEGASI(NRP1),
     & DOMEGASE(NRP1),DOMEGASIM(NRP1),DOMEGASEM(NRP1),
     & DLNRHO(NRP1),DLNRHOM(NRP1),H1,H2,ZTEMP
"""

DRIVER = """program physical_diamagnetic_profile
  implicit none
  integer, parameter :: np=6, nr=np-1
  integer :: mesh, initial
  real :: cs(np), csm(np), rho(np), rhom(np), p(np), pm(np)
  real :: pp(np), ppm(np), dpsi(np), dpsim(np), oi(np), oe(np)
  real :: oim(np), oem(np), doi(np), doe(np), doim(np), doem(np)
  real :: dn(np), dnm(np), bfield, cyclotron, alpha, slope
  real :: expected_i, expected_e, frequency_error, derivative_error
  bfield=2.
  cyclotron=5.
  alpha=.37
  slope=-3.
  expected_i=-alpha*bfield*slope/cyclotron
  expected_e=(1.-alpha)*bfield*slope/cyclotron
  do mesh=1,2
    if (mesh == 1) then
      cs=[0., .03, .11, .29, .63, 1.]
    else
      cs=[0., .08, .25, .51, .9, 1.]
    end if
    csm(1:nr)=(cs(1:nr)+cs(2:np))/2.
    csm(np)=1.
    rho=1.
    rhom=1.
    p=4.+slope*cs**2
    pm=4.+slope*csm**2
    pp=slope
    ppm=slope
    dpsi=2.*cs
    dpsim=2.*csm
    do initial=1,2
      oi=0.
      oe=0.
      if (initial == 2) then
        oi=71.
        oe=-93.
      end if
      oim=0.
      oem=0.
      doi=0.
      doe=0.
      doim=0.
      doem=0.
      dn=0.
      dnm=0.
      call native_diamag(np,nr,cs,csm,rho,rhom,p,pm,pp,ppm, &
        dpsi,dpsim,bfield,cyclotron,alpha,1.,1.,oi,oe,oim,oem, &
        doi,doe,doim,doem,dn,dnm)
      frequency_error=max(maxval(abs(oi-expected_i)), &
        maxval(abs(oe-expected_e)), &
        maxval(abs(oim(1:nr)-expected_i)), &
        maxval(abs(oem(1:nr)-expected_e)))
      derivative_error=max(maxval(abs(doi)),maxval(abs(doe)), &
        maxval(abs(doim)),maxval(abs(doem)))
      print *, 'PHYSICAL_DIAMAG',mesh,initial,frequency_error, &
        derivative_error,oe(np),expected_e
      if (.not.(frequency_error < 5.e-14)) stop 1
      if (.not.(derivative_error < 5.e-12)) stop 2
    end do
  end do
end program
"""


def run_native(source_file):
    text = Path(source_file).read_text()
    start = text.index("C     COMPUTE ADDITIONAL EQUILIBRIUM PROFILES")
    end = text.index("      IF (ISMPIRUN.EQ.0", start)
    block = DECLARATIONS + text[start:end] + "      END\n"
    with tempfile.TemporaryDirectory(prefix="mars-diamagnetic-") as directory:
        root = Path(directory)
        (root / "native.f").write_text(block)
        (root / "oracle.f90").write_text(DRIVER)
        compile_result = subprocess.run(
            [os.environ.get("FC", "gfortran"), "-O0", "-g", "-fcheck=all",
             "-fbacktrace", "-fdefault-real-8", "-fdefault-double-8",
             "-ffixed-line-length-none", "-std=legacy",
             *shlex.split(os.environ.get("MARS_TEST_FCFLAGS", "")), str(root / "native.f"),
             str(root / "oracle.f90"), "-o", str(root / "oracle")],
            capture_output=True, text=True,
        )
        if compile_result.returncode:
            raise RuntimeError(compile_result.stdout + compile_result.stderr)
        return subprocess.run([str(root / "oracle")], capture_output=True, text=True)


class DiamagneticEdge(unittest.TestCase):
    def test_constant_physical_frequencies_and_zero_gradients(self):
        result = run_native(ROOT / "MarsQ_2FK/marsq.f")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.count("PHYSICAL_DIAMAG"), 4)
        print(result.stdout, end="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-file", type=Path)
    args, remaining = parser.parse_known_args()
    if args.source_file:
        result = run_native(args.source_file)
        print(result.stdout, end="")
        print(result.stderr, end="")
        raise SystemExit(result.returncode)
    unittest.main(argv=[__file__, *remaining])
