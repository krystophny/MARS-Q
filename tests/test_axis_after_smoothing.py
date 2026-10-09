"""The next nonlinear source must use the minimum of the smoothed field."""
from pathlib import Path
import subprocess
from support import NATIVE, NativeTest

def build(source,native,out):
 out.mkdir(parents=True,exist_ok=False);public=(source/'nonlin.f90').exists()
 if public:
  (out/'nonlin.f90').write_text(native)
  globals='''module globals
implicit none
integer,parameter::rkind=kind(1.d0)
integer::NINSCA=2,N4NSNT=4,NSMOOTH,NSURF=6,NVERBOSE=0,NCON=0,index_out=1
real(rkind)::CPSICL(4),CPSIO(4),RELAX=0,RESIDU=0,CEPS=1e-10,SPSIM
 type code_type
  integer::output_flag=0
  character(180)::output_diag(3)=''
 end type
 type out_type
  type(code_type)::codeparam
 end type
 type(out_type)::eqchease_out(1)
end module
''';(out/'globals.f90').write_text(globals);nativefile='nonlin.f90';globalsfile='globals.f90';flags=['-ffree-line-length-none']
 else:
  start=native.index('         SUBROUTINE NONLIN\n');end=native.index('\nC*DECK',start);(out/'nonlin.f').write_text(native[start:end]+'\n');globals='''      MODULE GLOBALS
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBLA.inc'
      INCLUDE 'COMBAL.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMSOL.inc'
      INCLUDE 'COMPHY.inc'
      END MODULE
''';(out/'globals.f').write_text(globals);nativefile='nonlin.f';globalsfile='globals.f';flags=['-fdefault-real-8','-fdefault-double-8','-ffixed-line-length-none','-mcmodel=medium','-I'+str(source.resolve())]
 fixture='''module oracle_state
implicit none
integer::source_calls=0,conver_calls=0
real(8)::bubble
end module
program physical_stage_oracle
use globals
use oracle_state
implicit none
character(30)::arg
call get_command_argument(1,arg)
read(arg,*)NSMOOTH
call get_command_argument(2,arg)
read(arg,*)bubble
NINSCA=2
N4NSNT=4
NSURF=6
RELAX=0.d0
CPSICL(1:4)=[-1.d0,0.d0,0.d0,0.d0]
SPSIM=-1.d0
call nonlin
if(source_calls/=2)error stop 8
open(17,file='physical_state.dat',form='unformatted')
write(17)CPSICL(1:4),SPSIM
close(17)
print *,'Native NONLIN field/axis PASS'
end program
subroutine setupb
use globals
use oracle_state
implicit none
real(8)::stage_beta,u,psi,minimum
source_calls=source_calls+1
stage_beta=CPSICL(2)
! Independent known polynomial: psi(u,z)=u²-1+stage_beta*u*(1-u)+z²/kappa².
! Smoothing changes its derivative jet stage_beta and preserves endpoint values.
minimum=-1.d0-stage_beta*stage_beta/(4.d0*(1.d0-stage_beta))
u=-stage_beta/(2.d0*(1.d0-stage_beta))
psi=.02d0**2-1.d0+stage_beta*.02d0*(1.d0-.02d0)
! The flux source requires the saved axis minimum of this actual field.
if(abs(SPSIM-minimum)>2.d-14)error stop 7
if(1.d0-psi/SPSIM<=0.d0)error stop 9
end subroutine
subroutine solvit
use globals
CPSICL(2)=0.d0
end subroutine
subroutine smooth
use globals
use oracle_state
CPSICL(2)=bubble
end subroutine
subroutine magaxe
use globals
implicit none
real(8)::stage_beta
stage_beta=CPSICL(2)
SPSIM=-1.d0-stage_beta*stage_beta/(4.d0*(1.d0-stage_beta))
end subroutine
subroutine output(n)
integer n
end subroutine
subroutine errorch(x,y)
real(8)::x(*),y(*)
end subroutine
subroutine error1(x,y)
real(8)::x(*),y(*)
end subroutine
subroutine conver(n,flag)
use oracle_state
integer n,flag
conver_calls=conver_calls+1
flag=0
if(conver_calls==2)flag=1
end subroutine
subroutine dcopy(n,x,ix,y,iy)
integer n,ix,iy,j
real(8)::x(*),y(*)
do j=1,n
 y(1+(j-1)*iy)=x(1+(j-1)*ix)
enddo
end subroutine
subroutine scopy(n,x,ix,y,iy)
integer n,ix,iy
real(8)::x(*),y(*)
call dcopy(n,x,ix,y,iy)
end subroutine
subroutine scopyr(r,n,x,ix,y,iy)
integer n,ix,iy,j
real(8)::x(*),y(*),r
do j=1,n
 y(1+(j-1)*iy)=r*x(1+(j-1)*ix)+(1-r)*y(1+(j-1)*iy)
enddo
end subroutine
''';(out/'fixture.f90').write_text(fixture)
 # Compile each format with its own flags, link the actual whole NONLIN.
 for file in [globalsfile,nativefile]:
  r=subprocess.run(['gfortran','-O0','-g','-std=legacy','-fcheck=all',*flags,'-c',file],cwd=out,capture_output=True,text=True);assert r.returncode==0,r.stderr[-2000:]
 r=subprocess.run(['gfortran','-O0','-g','-std=legacy','-fcheck=all','-ffree-line-length-none','-mcmodel=medium','-c','fixture.f90'],cwd=out,capture_output=True,text=True);assert r.returncode==0,r.stderr[-2000:]
 r=subprocess.run(['gfortran','-mcmodel=medium','globals.o','nonlin.o','fixture.o','-o','oracle'],cwd=out,capture_output=True,text=True);assert r.returncode==0,r.stderr[-2000:]


class AxisAfterSmoothing(NativeTest):
    def test_axis_matches_actual_field_in_all_source_calls(self):
        folder=self.work/'build'
        build(NATIVE,(NATIVE/'chease.f').read_text(),folder)
        for smoothing in (0,1):
            for bubble in (-.1,-.2):
                with self.subTest(smoothing=smoothing,bubble=bubble):
                    result=subprocess.run([str(folder/'oracle'),str(smoothing),str(bubble)],
                        cwd=self.work,capture_output=True,text=True,timeout=30)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
