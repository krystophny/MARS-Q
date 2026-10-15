"""Execute the native opt-in DIRECT routine against a flat-field work oracle.

No solver replay: all radial/angular coefficients are constant, so both the
continuum integral and native linear finite element integral are exact.
The dependencies provide only this geometry and file IO; the tested routine
is extracted unchanged from kinetic.f. Scratch uses TMPDIR or the system's
disk-backed /var/tmp, never the shared RAM /tmp.
"""
from pathlib import Path
import math
import os
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = """
module dimensim
  implicit none
  integer, parameter :: nr=2,nrp1=3,nchi=8,msmax=1
end module
module globalm
  use dimensim
  implicit none
  real(8), parameter :: pi=acos(-1.d0)
  complex(8), parameter :: ci=(0.d0,1.d0)
  real(8) :: rm(msmax,2)=0, csh(nr)=1, csm(nr)=[.5d0,1.5d0]
  complex(8) :: b1u(nrp1,msmax)=0,b2u(nrp1,msmax)=0
  complex(8) :: b3u(nrp1,msmax)=0,x1u(nrp1,msmax)=0,x2u(nrp1,msmax)=0
end module
module rcomdm
  use dimensim
  implicit none
  real(8) :: t(nrp1)=1,tm(nr)=1,tsave(nrp1)=2,tmsave(nr)=2
  real(8) :: dpsids(nrp1)=0,dpsidsm(nr)=0,ppeqm(nr)=0
  real(8) :: g22l(nrp1,nchi)=1,g22lm(nr,nchi)=1,g12lm(nr,nchi)=0
  real(8) :: rja(nrp1,nchi)=1,rjam(nr,nchi)=1
  real(8) :: req(nrp1,nchi)=1,reqm(nr,nchi)=1
end module
module kineticm
end module
module feedbackm
  implicit none
  integer :: ktrest=1
end module
module toolbox
contains
  integer function assignfreefileunit()
    assignfreefileunit=77
  end function
end module
! Constant field geometry has exactly zero angular derivative.
subroutine derchi(a,b,n,m)
  implicit none
  integer :: n,m
  real(8) :: a(n,8),b(m,8)
  b=0
end subroutine
program fixture
  use dimensim
  use globalm
  use rcomdm
  use feedbackm
  implicit none
  integer :: k
  real(8), parameter :: alpha=.3d0
  complex(8), parameter :: pressure=(.7d0,-.2d0)
  complex(8) :: pa(nrp1,msmax,1),pe(nrp1,msmax,1)
  complex(8) :: ya(nrp1,1),ye(nrp1,1)
  character(8) :: arg
  call get_command_argument(1,arg)
  read(arg,*) k
  ! x_s=alpha*s, B_phi=-F*alpha, P_perp=pressure, P_parallel=0.
  ! J=R=1 and B^2=F^2 give integral -2*pi^2*alpha*pressure.
  ! The passive matrix and linear FE integral are exact for these constants.
  x1u(:,1)=alpha*[0.d0,1.d0,2.d0]
  b3u=-2*alpha
  pa=0
  pe=pressure
  ya=0
  ye=2*pi*pi*alpha*pressure
  if (k == 0) then
    ktrest=0
    t=2
    tm=2
    tsave=97
    tmsave=97
  endif
  call calcdwkdirectcheck(pa,pe,ya,ye,1)
  if (k == 0) then
    if (any(t /= 2.d0) .or. any(tm /= 2.d0)) stop 2
    if (ktrest /= 0) stop 3
  else
    if (any(t /= 1.d0) .or. any(tm /= 1.d0)) stop 4
    if (ktrest /= 1) stop 5
  endif
end program
"""


@unittest.skipUnless(shutil.which("gfortran"), "native fixture needs gfortran")
class DirectSavedField(unittest.TestCase):
    def test_native_direct_saved_and_current_field(self):
        source = (ROOT / "MarsQ_2FK/kinetic.f").read_text()
        start = source.index("      SUBROUTINE CALCDWKDIRECTCHECK(")
        end = source.index("      END SUBROUTINE CALCDWKDIRECTCHECK", start)
        end += len("      END SUBROUTINE CALCDWKDIRECTCHECK")
        scratch = os.environ.get("TMPDIR", "/var/tmp")
        with tempfile.TemporaryDirectory(prefix="mars-direct-field-", dir=scratch) as d:
            p = Path(d)
            (p / "fixture.f90").write_text(FIXTURE)
            (p / "native.f").write_text(source[start:end] + "\n")
            (p / "comioc.inc").write_text("")
            subprocess.run(["gfortran", "-fcheck=all", "-ffixed-line-length-none",
                            "fixture.f90", "native.f", "-o", "check.x"],
                           cwd=p, check=True, capture_output=True, text=True)
            expected = -2 * math.pi**2 * .3 * complex(.7, -.2)
            for k in (1, 0):
                subprocess.run([str(p / "check.x"), str(k)], cwd=p,
                               check=True, capture_output=True, text=True)
                rows = [line.split() for line in (p / "DWK_DIRECT_CHECK.OUT").read_text().splitlines()
                        if line.strip() and not line.lstrip().startswith("%")]
                self.assertEqual(len(rows), 2)
                for row in rows:
                    direct = complex(float(row[3]), float(row[4]))
                    actual = complex(float(row[5]), float(row[6]))
                    self.assertLess(abs(direct / expected - 1), 5e-10)
                    self.assertLess(abs(actual / expected - 1), 5e-10)
                    self.assertLess(abs(direct / actual - 1), 5e-10)


if __name__ == "__main__":
    unittest.main()
