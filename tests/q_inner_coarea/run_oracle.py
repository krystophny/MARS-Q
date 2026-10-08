"""Whole native PROFILE missing-inner q coarea oracle, no equilibrium solve."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def fixture(source,folder):
    folder.mkdir(parents=True,exist_ok=False)
    native=(source/'chease.f').read_text()
    start=native.index('         SUBROUTINE PROFILE(KN)')
    end=native.index('\nC*DECK',start)
    (folder/'profile.f').write_text(native[start:end]+'\n')
    main='''      PROGRAM COAREA_ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMBLA.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMISO.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMMAP.inc'
      DOUBLE PRECISION KAPPA,COAREA,RHO2,EXPECTED
      CHARACTER*32 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NSTTP
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) KAPPA
      CPI=ACOS(-1.D0)
      SPSIM=-.1D0
      RMAG=1.D0
      R0=1.D0
      RZ0=0.D0
      RZMAG=0.D0
      CPSICL(1)=-.099D0
      PSIISO(1:6)=[-.0995D0,-.098D0,-.095D0,-.09D0,-.08D0,-.06D0]
      COAREA=CPI*KAPPA*.01D0/.1D0
      BSFRAC=KAPPA
      CALL PROFILE(6)
      IF (NSTTP.EQ.4) THEN
         RHO2=1.D0-PSIISO(1)/SPSIM
         EXPECTED=SQRT(1.D0-.01D0*RHO2)/COAREA
         IF (ABS(CID0(1)/COAREA-1.D0).GT.2.D-13) STOP 7
         IF (ABS(CID2(1)*CIDQ(1)-1.D0).GT.2.D-13) STOP 8
         IF (ABS(CID2(1)/EXPECTED-1.D0).GT.2.D-11) STOP 9
      ENDIF
      OPEN(17,FILE='coarea.dat',FORM='UNFORMATTED')
      WRITE(17) CID0(1:6),CIDR(1:6),CIDQ(1:6),CID2(1:6)
      CLOSE(17)
      PRINT *,'Whole native PROFILE coarea PASS',NSTTP,KAPPA
      END
      SUBROUTINE ISOFIND(K1,K2,S,T,W,AX,EDGE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER K1,K2
      DOUBLE PRECISION S(*),T(*),W(*),AX,EDGE
      RETURN
      END
      SUBROUTINE CINT(K,S,T,W)
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMISO.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMCON.inc'
      DOUBLE PRECISION S(*),T(*),W(*),COAREA,RHO2
      COAREA=ACOS(-1.D0)*BSFRAC*.01D0/.1D0
      RHO2=1.D0-PSIISO(K)/SPSIM
      CIDR(K)=COAREA
      CIDQ(K)=COAREA/SQRT(1.D0-.01D0*RHO2)
      IF (NSTTP.EQ.4) THEN
         CID0(K)=COAREA
         CID2(K)=1.D0/CIDQ(K)
      ELSE IF (NSTTP.LE.2) THEN
         CID0(K)=1.D0
         CID2(K)=1.D0
      ELSE
         CID0(K)=1.D0
         CID2(K)=0.D0
      ENDIF
      RETURN
      END
      SUBROUTINE ISOFUN(KN)
      INTEGER KN
      RETURN
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
'''
    (folder/'oracle.f').write_text(main)
    command=['gfortran','-O0','-g','-fdefault-real-8','-fdefault-double-8',
             '-ffixed-line-length-none','-std=legacy','-mcmodel=medium',
             '-fcheck=all','-finit-real=snan','-ffpe-trap=invalid',
             '-I'+str(source),'profile.f','oracle.f','-o','oracle']
    result=subprocess.run(command,cwd=folder,capture_output=True,text=True)
    if result.returncode:raise RuntimeError(result.stderr[-2000:])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parent',type=Path);parser.add_argument('fixed',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for name,source in [('parent',args.parent),('fixed',args.fixed)]:
        build=args.output/name;fixture(source,build)
        for selector in [1,3,4]:
            for kappa in [1.,1.7]:
                folder=build/f'{selector}_{kappa}';folder.mkdir()
                result=subprocess.run([str(build/'oracle'),str(selector),str(kappa)],
                    cwd=folder,capture_output=True,text=True)
                if name=='parent'and selector==4:
                    assert result.returncode!=0 and 'SIGFPE'in result.stderr,result.stderr
                else:assert result.returncode==0,(name,selector,result.stderr)
                rows.append(dict(variant=name,selector=selector,kappa=kappa,
                    returncode=result.returncode,stdout=result.stdout.strip(),stderr=result.stderr.strip()))
    for selector in [1,3]:
        for kappa in [1.,1.7]:
            name=f'{selector}_{kappa}/coarea.dat'
            assert(args.output/'parent'/name).read_bytes()==(args.output/'fixed'/name).read_bytes()
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    report=dict(status='PASS',scope='Whole actual PROFILE IP2 with independently analytic circular/elliptical quadratic coarea data supplied by CINT fixture. No native PDE. Undefined q4 axis locals cause parent SIGFPE; candidate retains raw coarea continuation and CID2=1/CIDQ. Legacy1/3 outputs byte-identical.',
        controls=rows,source_sha256={k:sha(v/'chease.f')for k,v in [('parent',args.parent),('fixed',args.fixed)]},
        oracle_sha256=sha(__file__),fixture_sha256={k:sha(args.output/k/'oracle.f')for k in ['parent','fixed']})
    (args.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',controls=len(rows))))


if __name__=='__main__':main()
