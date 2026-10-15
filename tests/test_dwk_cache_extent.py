"""Execute the production cache guard with poisoned unpopulated storage."""
from pathlib import Path
import re,shutil,subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]

@unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
class CacheExtent(unittest.TestCase):
 def test_only_serialized_coefficients_are_validated(self):
  source=(ROOT/'MarsQ_2FK/kinetic.f').read_text()
  routine=source[source.index('      SUBROUTINE CALCDWKCOMP'):]
  match=re.search(r'         OCACHEFINITE = OCACHEFINITE.AND\.(.*?)         DO INDX=',routine,re.S)
  block=match.group(0).removesuffix('         DO INDX=')
  # These executable statements are the production verifier. The behavioral
  # oracle is that serialization populates ten coefficient slots; spare
  # capacity must not affect finiteness, magnitude or a ten-channel contraction.
  fixture=DRIVER.replace('C     PRODUCTION_GUARD',block)
  with tempfile.TemporaryDirectory(dir=ROOT,prefix='.cache-extent-') as temp:
   work=Path(temp);(work/'fixture.f').write_text(fixture)
   subprocess.run(['gfortran','-O2','-ffixed-line-length-none','-fcheck=all','fixture.f','-o','fixture'],cwd=work,check=True,capture_output=True,text=True)
   output=subprocess.check_output([str(work/'fixture')],text=True)
  rows=[line.split() for line in output.splitlines()]
  # Baseline and inactive NaN/Inf/huge values must yield identical admitted data.
  for row in rows[:4]:
   self.assertEqual(row[1],'T',row)
   self.assertAlmostEqual(float(row[2]),10.0)
   self.assertAlmostEqual(float(row[3]),55.0)
  # A NaN in the real or imaginary part of ANY populated channel must reject.
  self.assertTrue(all(row[1]=='F' for row in rows[4:]),output)
  print('64 cases passed: four finite controls; 60 rejected active NaNs')

DRIVER='''      program fixture
      use, intrinsic :: ieee_arithmetic
      implicit none
      complex*16 buffer_i(1,1,1,30),bufferm_i(1,1,1,30),
     &           buffer_i1(1,1,1,30)
      real*8 cachemax,nan,inf,contraction
      logical ocachefinite
      integer mode,j
      nan=ieee_value(0d0,ieee_quiet_nan)
      inf=ieee_value(0d0,ieee_positive_inf)
      do mode=0,63
         buffer_i=0d0
         bufferm_i=0d0
         buffer_i1=0d0
         do j=1,10
            buffer_i(1,1,1,j)=dcmplx(dble(j),0d0)
            bufferm_i(1,1,1,j)=dcmplx(dble(j),0d0)
            buffer_i1(1,1,1,j)=dcmplx(dble(j),0d0)
         enddo
         if (mode.eq.1) buffer_i(:,:,:,11:30)=dcmplx(nan,nan)
         if (mode.eq.2) bufferm_i(:,:,:,11:30)=dcmplx(inf,0d0)
         if (mode.eq.3) buffer_i1(:,:,:,11:30)=dcmplx(1d100,0d0)
         if (mode.ge.4) then
            j=mod(mode-4,10)+1
            select case ((mode-4)/10)
            case (0)
               buffer_i(1,1,1,j)=dcmplx(nan,0d0)
            case (1)
               bufferm_i(1,1,1,j)=dcmplx(nan,0d0)
            case (2)
               buffer_i1(1,1,1,j)=dcmplx(nan,0d0)
            case (3)
               buffer_i(1,1,1,j)=dcmplx(0d0,nan)
            case (4)
               bufferm_i(1,1,1,j)=dcmplx(0d0,nan)
            case (5)
               buffer_i1(1,1,1,j)=dcmplx(0d0,nan)
            end select
         endif
         ocachefinite=.true.
         cachemax=0d0
C     PRODUCTION_GUARD
         contraction=sum(real(buffer_i(:,:,:,1:10)))
         print *, mode,ocachefinite,cachemax,contraction
      enddo
      end program
'''
if __name__=='__main__':unittest.main()
