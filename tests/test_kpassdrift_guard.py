"""Input guard of KPASSDRIFT in marsq.f, executed on the production statements.

IFOW is referenced only by the KORBIT orbit tracer (IORBIT > 0).  The
orbit-width response is switched by PSPECIES_FOWP/FOWT, whose KIA_PASS terms
(KOPT = 1..3) share the denominator that carries the KPASSDRIFT drift.  The
guard must therefore reject only KFASTRUN /= 0, whose analytic passing energy
integral has no drift term.
"""
import re

from support import ROOT, NativeTest

SRC = ROOT / "MarsQ_2FK" / "marsq.f"


def guard():
    source = SRC.read_text()
    start = source.index("      IF (KPASSDRIFT.NE.0.AND.KPASSDRIFT.NE.1)")
    end = source.index("      READ(CHNAME,QLIN", start)
    return "\n".join(line for line in source[start:end].splitlines()
                     if not re.match(r"^[Cc*!]", line)) + "\n"


PROGRAM = """      PROGRAM KPDGUARD
      INTEGER KPASSDRIFT,KFASTRUN,IFOW,INCKIN
      REAL*8 PSPECIES_NP(2)
      CHARACTER(16) ARG
      INCKIN = 1
      PSPECIES_NP = 1.D0
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) KPASSDRIFT
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) KFASTRUN
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      READ(ARG,*) IFOW
@GUARD@
      WRITE(*,'(A)') 'ACCEPTED'
      END
"""


class KpassdriftGuard(NativeTest):
    def run_guard(self, exe, *flags):
        result = self.execute(exe, *flags)
        return "ACCEPTED" in result.stdout, result.stdout + result.stderr

    def test_orbit_tracer_flag_accepted_fast_run_rejected(self):
        exe = self.compile({"guard.f": PROGRAM.replace("@GUARD@\n", guard())})
        cases = {(1, 0, 0): True, (1, 0, 1): True, (0, 1, 1): True, (1, 1, 0): False}
        for flags, accepted in cases.items():
            ok, text = self.run_guard(exe, *flags)
            self.assertEqual(ok, accepted, (flags, text))


if __name__ == "__main__":
    import unittest
    unittest.main()
