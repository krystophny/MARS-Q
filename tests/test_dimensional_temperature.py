"""Native compatibility profiles retain dimensional temperatures for Z=1,2.

The independent oracle converts the supplied eV and electron density to the
code pressure unit. Changing ion charge changes the neutralizing ion density,
and must not change either supplied temperature amplitude.
"""
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DimensionalTemperature(unittest.TestCase):
    def test_dimensional_input_requires_native_amplitudes(self):
        from test_mars_input import minimal_input, run_with_input
        for dimensional in (0, 1):
            text = minimal_input(kinetic='KPROFTAUTH=1')
            text = text.replace('&BASIC',
                                f'&BASIC\n NPROFIE=4, NEXPV={dimensional}')
            result = run_with_input(text)
            if dimensional:
                self.assertIn('&BASIC', result.stdout)
                self.assertNotIn('KPROFTAUTH=1 REQUIRES', result.stdout)
            else:
                self.assertIn('KPROFTAUTH=1 REQUIRES NPROFIE=4,NEXPV=1,2 MAXWELLIANS',
                              result.stdout)

    def test_native_temperature_authority_and_legacy_control(self):
        source = Path(os.environ.get('MARS_ORACLE_SOURCE', ROOT))/'MarsQ_2FK'
        text = (source/'kinetic.f').read_text()
        begin = text.index('      SUBROUTINE KEQPROF\n')
        end = text.index('C     COMPUTE DP/DPSI FOR EACH PARTICLE SPECIES', begin)
        # Execute the complete native density/pressure/temperature stage.
        # The downstream derivative stage cannot affect those amplitudes.
        kernel = text[begin:end]+' 100  CONTINUE\n      RETURN\n      END\n'
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', ROOT)) as td:
            work = Path(td)
            (work/'native.f').write_text(kernel)
            (work/'driver.f90').write_text(DRIVER)
            command = [os.environ.get('FC', 'gfortran'), '-O0', '-fcheck=all',
                       '-fdefault-real-8', '-fdefault-double-8',
                       '-ffixed-line-length-none', '-ffree-line-length-none',
                       str(source/'globalm.f'), str(source/'dimensim.f'),
                       'native.f', 'driver.f90', '-o', 'oracle']
            result = subprocess.run(command, cwd=work, text=True, capture_output=True,
                                    timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(work/'oracle')], cwd=work, text=True,
                                    capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            rows = [list(map(float, line.split()[1:])) for line in result.stdout.splitlines()
                    if line.startswith('TEMPERATURE ')]
            self.assertEqual(len(rows), 16)
            for charge, authority, grid, js, ion, electron, ni, ne in rows:
                if authority:
                    scale = 5e19*1.6021917e-19*(4e-7*math.pi)/2.3**2
                    expected_i = 1300*scale*(1+.1*(js-1)+.03*(grid-1))
                    expected_e = 700*scale*(1+.2*(js-1)+.05*(grid-1))
                else:
                    # With authority disabled the native path repartitions PEQ.
                    p = 1.7 if grid == 1 else 1.9
                    ratio_i = (1300/2000)*(1+.1*(js-1)+.03*(grid-1))
                    ratio_e = (700/2000)*(1+.2*(js-1)+.05*(grid-1))
                    ratio = ratio_i/(ratio_i+ratio_e)
                    expected_i = p*ratio/ni
                    expected_e = p*(1-ratio)/ne
                self.assertTrue(math.isclose(ion, expected_i, rel_tol=3e-13),
                                (charge, authority, grid, js, ion, expected_i))
                self.assertTrue(math.isclose(electron, expected_e, rel_tol=3e-13),
                                (charge, authority, grid, js, electron, expected_e))


DRIVER = '''program temperature_authority
use globalm
use dimensim
implicit none
integer :: iz, mode, j, g
NR=2; NRP1=3; NSPECIES=2; KPROFSPEC=0; NPROFIE=4; NEXPV=1
B0EXP=2.3d0; ZNE0=5d19; ZTI0=1300d0; ZTE0=700d0
ALPHAP=ZTI0/(ZTI0+ZTE0)
allocate(RHO(3),RHOM(3),PEQ(3),PEQM(3),TEMPI(3),TEMPIM(3),TEMPE(3),TEMPEM(3))
allocate(ESPECIES_DENF(3,2,2),ESPECIES_DEN(3,2,2),ESPECIES_PREF(3,2,2))
allocate(ESPECIES_PRE(3,2,2),ESPECIES_TEM(3,2,2))
RHO=1d0; RHOM=.9d0; PEQ=1.7d0; PEQM=1.9d0
do j=1,3
 TEMPI(j)=1d0+.1d0*(j-1); TEMPIM(j)=TEMPI(j)+.03d0
 TEMPE(j)=1d0+.2d0*(j-1); TEMPEM(j)=TEMPE(j)+.05d0
enddo
FRACPTH=(ZTI0+ZTE0)*ZNE0*1.6021917d-19*(4d-7*PI)/B0EXP**2/PEQ(1)
do iz=1,2
 ESPECIES_Z(1)=iz; ESPECIES_Z(2)=-1d0
 do mode=0,1
  KPROFTAUTH=mode
  call KEQPROF
  do g=1,2
   do j=1,2
    write(*,'(A,4I3,4ES25.16)') 'TEMPERATURE ',iz,mode,g,j, &
       ESPECIES_TEM(j,g,1),ESPECIES_TEM(j,g,2), &
       ESPECIES_DEN(j,g,1),ESPECIES_DEN(j,g,2)
   enddo
  enddo
 enddo
enddo
end program
subroutine KREADEXPLICITPROFILES
error stop 'unexpected explicit-profile branch'
end subroutine
'''


if __name__ == '__main__':
    unittest.main()
