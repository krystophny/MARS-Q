"""Inner coareas against closed circle/ellipse integrals and independent roots."""
from pathlib import Path
import re
import subprocess
import math
import struct
from support import NATIVE, NativeTest

def routine(text, name):
    start = text.index('         SUBROUTINE '+name+'(')
    end = start+re.search(r'^         END\s*$', text[start:], re.MULTILINE).end()
    return text[start:end]+'\n'


def build(source, folder):
    folder.mkdir(parents=True, exist_ok=False)
    text = (source/'chease.f').read_text()
    code = routine(text, 'PROFILE')+routine(text, 'GAUSS')
    if '         SUBROUTINE QAXISCOAREA(' in text:
        code += routine(text, 'QAXISCOAREA')+routine(text, 'QAXISVALUE')
    (folder/'native.f').write_text(code)
    (folder/'oracle.f').write_text('''      PROGRAM ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBLA.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMISO.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMMAP.inc'
      INCLUDE 'COMESH.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMSOL.inc'
      COMMON /MODEL/ KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ,IMODE
      REAL*8 KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ
      INTEGER IMODE
      CHARACTER*64 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NSTTP
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) KAPPA
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      READ(ARG,*) IMODE
      CALL GET_COMMAND_ARGUMENT(4,ARG)
      READ(ARG,*) NT
      CALL GET_COMMAND_ARGUMENT(5,ARG)
      READ(ARG,*) QSIGN
      CPI=ACOS(-1.D0)
      R0=1.D0
      RZ0=0.D0
      RMAG=1.00125D0
      RZMAG=.00025D0
      IF (IMODE.EQ.8) THEN
         RMAG=.99875D0
         RZMAG=0.D0
      ENDIF
      AXR=RMAG
      AXZ=RZMAG
      SPSIM=-.1D0
      NISO=100
      NMGAUS=4
      IF (IMODE.EQ.8) NMGAUS=1
      NT1=NT+1
      A0=.1D0
      WIDTH=1.D-8
      BUMP=0.D0
      IF (IMODE.EQ.1) BUMP=.0075D0
      CC=(.1D0-BUMP*(1.D0-EXP(-A0*A0/WIDTH)))/(A0*A0)
      DO J=1,NT1
         CT(J)=2.D0*CPI*(J-1)/NT
      ENDDO
      DO J=1,NISO
         PSIISO(J)=SPSIM*(1.D0-((J-.5D0)/NISO)**2)
      ENDDO
      CALL EVLATE(1,R0,RZ0,D1,D2,CPSICL(1))
      IF (IMODE.EQ.2) SPSIM=.1D0
      IF (IMODE.EQ.3) RMAG=0.D0
      IF (IMODE.EQ.4) PSIISO(1)=SPSIM
      CALL PROFILE(NISO)
      OPEN(17,FILE='coarea.dat',FORM='UNFORMATTED')
      WRITE(17) CID0(1:NISO),CIDR(1:NISO),CIDQ(1:NISO),CID2(1:NISO)
      CLOSE(17)
      END
      SUBROUTINE EVLATE(KCASE,R,Z,DR,DZ,P)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER KCASE,IMODE
      COMMON /MODEL/ KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ,IMODE
      DOUBLE PRECISION KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ
      X=R-AXR
      Y=Z-AXZ
      S=X*X+Y*Y/(KAPPA*KAPPA)
      P=-.1D0+CC*S+BUMP*(1.D0-EXP(-S/WIDTH))
      DER=CC+BUMP/WIDTH*EXP(-S/WIDTH)
      DR=2.D0*DER*X
      DZ=2.D0*DER*Y/(KAPPA*KAPPA)
      IF (IMODE.EQ.5) THEN
         P=-.1D0+CC*(X*X-Y*Y)
         DR=2.D0*CC*X
         DZ=-2.D0*CC*Y
      ENDIF
      IF (IMODE.EQ.6) P=SQRT(-ABS(R))
      END
      SUBROUTINE BOUND(N,T,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N,IMODE
      DOUBLE PRECISION T(N),B(N)
      COMMON /MODEL/ KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ,IMODE
      DOUBLE PRECISION KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ
      DO J=1,N
         CO=COS(T(J))
         SI=SIN(T(J))
         AA=CO*CO+SI*SI/(KAPPA*KAPPA)
         BB=2.D0*((1.D0-AXR)*CO-AXZ*SI/(KAPPA*KAPPA))
         DD=(1.D0-AXR)**2+AXZ**2/(KAPPA*KAPPA)-A0*A0
         B(J)=(-BB+SQRT(BB*BB-4.D0*AA*DD))/(2.D0*AA)
         IF (IMODE.EQ.7) B(J)=1.D-10
      ENDDO
      END
      SUBROUTINE ISOFIND(I1,I2,S,T,W,AX,EDGE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER I1,I2
      DOUBLE PRECISION S(*),T(*),W(*),AX,EDGE
      END
      SUBROUTINE CINT(K,S,T,W)
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMISO.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMSOL.inc'
      COMMON /MODEL/ KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ,IMODE
      REAL*8 KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ
      INTEGER IMODE
      REAL*8 S(*),T(*),W(*),LO,HI,MID,VAL,DER,SROOT
      LO=0.D0
      HI=A0*A0
      DO J=1,100
         MID=.5D0*(LO+HI)
         VAL=-.1D0+CC*MID+BUMP*(1.D0-EXP(-MID/WIDTH))
         IF (VAL.LT.PSIISO(K)) THEN
            LO=MID
         ELSE
            HI=MID
         ENDIF
      ENDDO
      SROOT=.5D0*(LO+HI)
      DER=CC+BUMP/WIDTH*EXP(-SROOT/WIDTH)
      CID0(K)=CPI*KAPPA/DER
      CIDR(K)=AXR*CID0(K)
      CIDQ(K)=CID0(K)/SQRT(AXR**2-SROOT)
      CID2(K)=1.D0/CIDQ(K)
      IF (NSTTP.LE.2) THEN
         CID0(K)=1.D0
         CID2(K)=1.D0
      ELSE IF (NSTTP.EQ.3) THEN
         CID0(K)=1.D0
         CID2(K)=0.D0
      ENDIF
      END
      SUBROUTINE ISOFUN(KN)
      INTEGER KN
      END
      INTEGER FUNCTION ISRCHFGE(N,A,INC,TARGET)
      INTEGER N,INC,J
      DOUBLE PRECISION A(*),TARGET
      ISRCHFGE=N+1
      DO J=1,N
         IF (A(1+(J-1)*INC).GE.TARGET) THEN
            ISRCHFGE=J
            RETURN
         ENDIF
      ENDDO
      END
''')
    cmd = ['gfortran', '-O0', '-g', '-fdefault-real-8', '-fdefault-double-8',
           '-ffixed-line-length-none', '-std=legacy', '-mcmodel=medium', '-fcheck=all',
           '-I'+str(source.resolve()), 'native.f', 'oracle.f', '-o', 'oracle']
    result = subprocess.run(cmd, cwd=folder, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2500:]


def expected(kappa, model):
    bump = .0075 if model == 1 else 0.
    width = 1e-8
    cc = (.1-bump)/.01
    rows = [[], [], [], []]
    for i in range(100):
        target = .1*((i+.5)/100)**2
        lo, hi = 0., .01
        for _ in range(90):
            mid = (lo+hi)/2
            if cc*mid+bump*(-math.expm1(-mid/width)) < target:
                lo = mid
            else:
                hi = mid
        root = (lo+hi)/2
        derivative = cc+bump/width*math.exp(-root/width)
        c0 = math.pi*kappa/derivative
        axis_r = .99875 if model == 8 else 1.00125
        cq = c0/math.sqrt(axis_r**2-root)
        for row, value in zip(rows, (c0,axis_r*c0,cq,1/cq)):
            row.append(value)
    return rows


class InnerCoarea(NativeTest):
    def test_shifted_contours_and_invalid_state_rejection(self):
        folder=self.work/'build'
        build(NATIVE,folder)
        cases=[(k,m,nt,q) for k in (1.,1.7) for m in (0,1)
               for nt in (16,32,64) for q in (-1.5,1.5)]
        cases += [(1.,m,32,1.5) for m in range(2,9)]
        for k,m,nt,q in cases:
            with self.subTest(kappa=k,model=m,panels=nt,q=q):
                result=subprocess.run([str(folder/'oracle'),*map(str,(4,k,m,nt,q))],
                    cwd=self.work,capture_output=True,text=True,timeout=30)
                if 2 <= m <= 7:
                    self.assertNotEqual(result.returncode,0)
                    continue
                self.assertEqual(result.returncode,0,result.stderr)
                measured=struct.unpack('<400d',(self.work/'coarea.dat').read_bytes()[4:-4])
                truth=expected(k,m)
                for row in range(4):
                    for column in range(27):
                        self.assertLessEqual(abs(measured[100*row+column]/truth[row][column]-1),
                                             1e-6 if nt==16 else 3e-8)
