"""Whole native PROFILE/axis-coarea against analytic shifted circle/ellipse fields."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.optimize import brentq


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
    x = ((np.arange(100)+.5)/100)**2
    root = np.array([brentq(lambda s: cc*s+bump*(-np.expm1(-s/width))-.1*v,
                           0., .01, xtol=1e-25, rtol=1e-14) for v in x])
    derivative = cc+bump/width*np.exp(-root/width)
    c0 = np.pi*kappa/derivative
    axis_R = .99875 if model == 8 else 1.00125
    cq = c0/np.sqrt(axis_R**2-root)
    return np.stack([c0, axis_R*c0, cq, 1/cq]), x


def main():
    p = argparse.ArgumentParser()
    p.add_argument('parent', type=Path); p.add_argument('fixed', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    args.parent = args.parent.resolve(); args.fixed = args.fixed.resolve()
    args.output = args.output.resolve(); args.output.mkdir(parents=True, exist_ok=False)
    rows, parity = [], {}
    for variant, source in [('parent', args.parent), ('fixed', args.fixed)]:
        folder = args.output/variant; build(source, folder)
        cases = [(4, k, m, nt, q) for k in [1., 1.7] for m in [0, 1]
                 for nt in [16, 32, 64] for q in [1.5, -1.5]]
        cases += [(4, 1., 8, 32, 1.5)]
        cases += [(s, k, 0, 32, 1.5) for s in [1, 3] for k in [1., 1.7]]
        if variant == 'fixed': cases += [(4, 1., m, 32, 1.5) for m in range(2, 8)]
        for i, (sel, k, m, nt, q) in enumerate(cases):
            run = folder/f'case{i}'; run.mkdir()
            result = subprocess.run([str(folder/'oracle'), *map(str, [sel, k, m, nt, q])],
                                    cwd=run, capture_output=True, text=True)
            (run/'stdout.txt').write_text(result.stdout); (run/'stderr.txt').write_text(result.stderr)
            if 2 <= m <= 7:
                assert result.returncode != 0 and not (run/'coarea.dat').exists()
                rows.append(dict(variant=variant, invalid_model=m, returncode=result.returncode)); continue
            assert result.returncode == 0, (variant, k, m, nt, result.stderr[-1500:])
            payload = (run/'coarea.dat').read_bytes()
            if sel != 4: parity[variant, sel, k] = payload; continue
            measured = np.frombuffer(payload[4:-4], '<f8').reshape(4, 100)
            truth, x = expected(k, m)
            err = float(np.max(abs(measured[:, :27]/truth[:, :27]-1)))
            if variant == 'fixed':
                assert err < (1e-6 if nt == 16 else 3e-8), (k, m, nt, err)
            rows.append(dict(variant=variant, kappa=k, model=m, panels=nt, signed_q=q,
                             inner27_max_relative_coarea_error=err,
                             reconstructed_q_max_error=float(np.max(abs((2*np.pi*q/measured[2])*truth[2]/(2*np.pi)-q)))))
    for sel in [1, 3]:
        for k in [1., 1.7]: assert parity['parent', sel, k] == parity['fixed', sel, k]
    assert max(r.get('inner27_max_relative_coarea_error', 0) for r in rows if r['variant']=='parent') > .1
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
    report = dict(status='PASS', scope='Whole native PROFILE and new helper/GAUSS, exact manufactured shifted circular/elliptical regular fields supplied by EVLATE and exact outer coareas supplied by CINT. No native PDE. Independent scalar roots and closed coarea integrals; localized smooth axis-curvature control demonstrates parent extrapolation error. Two signed q reconstructions, three angular refinements, four legacy byte guards, six invalid controls.',
                  controls=rows, source_sha256={str(s): sha(s/'chease.f') for s in [args.parent, args.fixed]},
                  oracle_sha256=sha(__file__))
    (args.output/'receipt.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', controls=len(rows), legacy_byte_guards=4)))


if __name__ == '__main__':
    main()
