"""Execute the production passive-carrier lifecycle in a minimal harness.

The regression suite's source-contract test catches accidental deletion of the
save/restore statements.  This test extracts the actual
``PREPAREKINETICENERGYMAT`` routine, supplies only its module/common state and
a ``LINEAR`` spy, and checks the values visible inside and after the call.
The same harness is compiled once with the post-call restore block removed;
that mutation must fail the behavioral assertions.
"""

from __future__ import annotations

import re
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
KINETIC = ROOT / "MarsQ_2FK" / "kinetic.f"


PREAMBLE = r"""
      MODULE RCOMDM
      END MODULE RCOMDM
      MODULE DIMENSIM
      INTEGER, PARAMETER :: MXMAX=2, MYMAX=2
      END MODULE DIMENSIM
      MODULE GLOBALM
      INTEGER :: NRP1, NR
      REAL*8 :: T(4), TSAVE(4), TM(4), TMSAVE(4)
      END MODULE GLOBALM
      MODULE KINETICM
      INTEGER :: KJPKEY, KPBKEY
      END MODULE KINETICM
      MODULE FEEDBACKM
      INTEGER :: KTREST
      END MODULE FEEDBACKM
      MODULE SPY
      INTEGER :: CALLS, SEEN_KTREST
      REAL*8 :: SEEN_T(4), SEEN_TM(4)
      END MODULE SPY

      SUBROUTINE LINEAR(ASUBM,BSUBM,CSUBM,DSUBM,ESUBM,FSUBM,GSUBM,HSUBM)
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE FEEDBACKM
      USE SPY
      IMPLICIT NONE
      COMPLEX*16 ASUBM(MXMAX,MXMAX,*), BSUBM(MXMAX,MXMAX,*),
     &           CSUBM(MXMAX,MXMAX,*), DSUBM(MYMAX,MYMAX,*),
     &           ESUBM(MXMAX,MYMAX,*), FSUBM(MYMAX,MXMAX,*),
     &           GSUBM(MYMAX,MXMAX,*), HSUBM(MXMAX,MYMAX,*)
      CALLS = CALLS + 1
      SEEN_T = T
      SEEN_TM = TM
      SEEN_KTREST = KTREST
      IF (KJPKEY.NE.0 .OR. KPBKEY.NE.1) ERROR STOP 11
      END SUBROUTINE LINEAR
"""


MAIN = r"""
      PROGRAM CHECK_PASSIVE_CARRIER
      USE GLOBALM
      USE KINETICM
      USE FEEDBACKM
      USE SPY
      IMPLICIT NONE
      COMPLEX*16 AL0, ALNORM
      COMMON /PAMARG/ AL0, ALNORM
      COMPLEX*16 A(2,2,1),B(2,2,1),C(2,2,1),D(2,2,1),
     &           E(2,2,1),F(2,2,1),G(2,2,1),H(2,2,1)
      NRP1=3
      NR=2
      AL0=(3.0,4.0)
      ALNORM=(9.0,2.0)
      T=(/1.0,1.0,1.0,1.0/)
      TSAVE=(/2.0,3.0,5.0,7.0/)
      TM=(/1.0,1.0,1.0,1.0/)
      TMSAVE=(/11.0,13.0,17.0,19.0/)
      KTREST=1
      CALLS=0
      CALL PREPAREKINETICENERGYMAT(A,B,C,D,E,F,G,H)
      IF (CALLS.NE.1 .OR. SEEN_KTREST.NE.0) ERROR STOP 21
      IF (MAXVAL(ABS(SEEN_T(1:3)-TSAVE(1:3))).GT.1D-12) ERROR STOP 22
      IF (MAXVAL(ABS(SEEN_TM(1:2)-TMSAVE(1:2))).GT.1D-12) ERROR STOP 23
      IF (MAXVAL(ABS(T(1:3)-1.0D0)).GT.1D-12) ERROR STOP 24
      IF (MAXVAL(ABS(TM(1:2)-1.0D0)).GT.1D-12) ERROR STOP 25
      IF (KTREST.NE.1) ERROR STOP 26
      IF (ABS(AL0-(3.0,4.0)).GT.1D-12) ERROR STOP 27

      T=TSAVE
      TM=TMSAVE
      KTREST=0
      CALLS=0
      AL0=(3.0,4.0)
      CALL PREPAREKINETICENERGYMAT(A,B,C,D,E,F,G,H)
      IF (CALLS.NE.1 .OR. SEEN_KTREST.NE.0) ERROR STOP 31
      IF (MAXVAL(ABS(SEEN_T(1:3)-TSAVE(1:3))).GT.1D-12) ERROR STOP 32
      IF (MAXVAL(ABS(SEEN_TM(1:2)-TMSAVE(1:2))).GT.1D-12) ERROR STOP 33
      IF (MAXVAL(ABS(T(1:3)-TSAVE(1:3))).GT.1D-12) ERROR STOP 34
      IF (MAXVAL(ABS(TM(1:2)-TMSAVE(1:2))).GT.1D-12) ERROR STOP 35
      IF (KTREST.NE.0) ERROR STOP 36
      IF (ABS(AL0-(3.0,4.0)).GT.1D-12) ERROR STOP 37
      WRITE(*,'(A)') 'PASS passive carrier lifecycle'
      END PROGRAM CHECK_PASSIVE_CARRIER
"""


def routine(source: str) -> str:
    start = source.index("      SUBROUTINE PREPAREKINETICENERGYMAT(")
    end_match = re.search(r"^      END$", source[start:], re.MULTILINE)
    if end_match is None:
        raise AssertionError("production routine terminator not found")
    end = start + end_match.end()
    body = source[start:end]
    body = body.replace("      USE RCOMDM\n", "")
    body = body.replace("      INCLUDE 'specmat.inc'", """      COMPLEX*16 ASUBM(2,2,*),BSUBM(2,2,*),CSUBM(2,2,*),
     &           DSUBM(2,2,*),ESUBM(2,2,*),FSUBM(2,2,*),GSUBM(2,2,*),
     &           HSUBM(2,2,*)""")
    body = body.replace("      INCLUDE 'compam.inc'", """      COMPLEX*16 AL0,ALNORM
      COMMON /PAMARG/ AL0,ALNORM""")
    body = body.replace("      INCLUDE 'comioc.inc'\n", "")
    return body


def remove_post_restore(body: str) -> str:
    block = re.compile(
        r"      IF \(KTRESTSAVE\.NE\.0\) THEN\n"
        r"(?:     &?.*\n|      .*\n)*?"
        r"      ENDIF\n"
    )
    matches = list(block.finditer(body))
    if len(matches) != 2:
        raise AssertionError(f"expected two lifecycle blocks, got {len(matches)}")
    match = matches[1]
    return body[: match.start()] + "      CONTINUE\n" + body[match.end() :]


class PassiveCarrierLifecycleTests(unittest.TestCase):
    def compile_and_run(self, body: str, directory: Path) -> subprocess.CompletedProcess[str]:
        source = directory / "lifecycle.f"
        executable = directory / "lifecycle.x"
        source.write_text(PREAMBLE + body + MAIN)
        subprocess.run(
            ["gfortran", "-ffixed-form", "-ffixed-line-length-none", "-O0", str(source), "-o", str(executable)],
            cwd=directory,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return subprocess.run([str(executable)], cwd=directory, check=False, text=True, capture_output=True)

    def test_actual_routine_lifecycle_and_restore_mutation(self) -> None:
        body = routine(KINETIC.read_text())
        with tempfile.TemporaryDirectory(prefix="mars-passive-carrier-") as tmp:
            directory = Path(tmp)
            normal = self.compile_and_run(body, directory)
            self.assertEqual(normal.returncode, 0, normal.stdout + normal.stderr)
            self.assertIn("PASS passive carrier lifecycle", normal.stdout)

            mutated = self.compile_and_run(remove_post_restore(body), directory)
            self.assertNotEqual(mutated.returncode, 0)


if __name__ == "__main__":
    unittest.main()
