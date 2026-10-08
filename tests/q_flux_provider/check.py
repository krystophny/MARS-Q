"""Whole native ISOFUN/CURENT against exact smooth signed-flux F² controls."""
import argparse,hashlib,json,subprocess,re
from pathlib import Path
import numpy as np

def extract(text,name):
 a=text.index('         SUBROUTINE '+name+'(');m=re.search(r'^      +END *$',text[a:],re.M)
 return text[a:a+m.end()]+'\n'

def build(source,out):
 out.mkdir(parents=True,exist_ok=False);text=(source/'chease.f').read_text();names=['ISOFUN','CURENT','SPLINE','TRIDAG','SCOPY','SSCAL','RESETI','VZERO'];code=''.join(extract(text,n)for n in names)
 if '         SUBROUTINE QFLUXSOURCE('in text:code+=extract(text,'QFLUXSOURCE')
 (out/'native.f').write_text(code)
 fixture='''      PROGRAM SOURCE_ORACLE
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMCON.inc'
      INCLUDE 'COMESH.inc'
      INCLUDE 'COMPHY.inc'
      INCLUDE 'COMNUM.inc'
      INCLUDE 'COMSUR.inc'
      INCLUDE 'COMSOL.inc'
      REAL*8 P(5),R(5),CUR(5),D(8),W(8),V(8),X,C,Q,EPS
      CHARACTER*64 ARG
      COMMON /ORACLEMODEL/ QORACLE,PP0,KAPPA,MODEL
      INTEGER MODEL
      REAL*8 QORACLE,PP0,KAPPA
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NSTTP
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) MODEL
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      READ(ARG,*) QORACLE
      CALL GET_COMMAND_ARGUMENT(4,ARG)
      READ(ARG,*) PP0
      CALL GET_COMMAND_ARGUMENT(5,ARG)
      READ(ARG,*) KAPPA
      SPSIM=-.1D0
      CPI=ACOS(-1.D0)
      NISO=8
      NSURF=6
      NPROFZ=0
      NTMF0=0
      NRFP=0
      PREDGE=0.D0
      C=-SPSIM/(CPI*KAPPA*.01D0)
      DO J=1,NISO
         CSIPR(J)=.005D0+(J-1)*.995D0/7.D0
         X=CSIPR(J)**2
         PSIISO(J)=SPSIM*(1.D0-X)
         IF (MODEL.EQ.3) THEN
            CID2(J)=C*(1.D0+.05D0*X)
         ELSE
            CID2(J)=C*SQRT(1.D0-.01D0*X)
         ENDIF
         CIDQ(J)=1.D0/CID2(J)
         CID0(J)=1.D0/C
         TMF(J)=1.D0
      ENDDO
      CALL ISOFUN(NISO)
      P=(/SPSIM*(1.D0+.001D0),SPSIM,
     +    SPSIM*(1.D0-1.D-8),SPSIM*.91D0,SPSIM*.36D0/)
      R=(/1.001D0,1.D0,1.002D0,1.03D0,1.08D0/)
      CALL CURENT(5,P,R,CUR)
      OPEN(17,FILE='profiles.dat',FORM='UNFORMATTED')
      WRITE(17) TTP(1:NISO),TMF(1:NISO),CUR
      CLOSE(17)
      DO J=1,NISO
         WRITE(*,'(A,I3,3ES26.17)') 'knot',J,
     +     CSIPR(J)**2,TTP(J),TMF(J)**2
      ENDDO
      DO J=1,5
         WRITE(*,'(A,I3,3ES26.17)') 'query',J,P(J),R(J),CUR(J)
      ENDDO
      END
      SUBROUTINE PPRIME(N,P,V)
      INTEGER N,J,MODEL
      REAL*8 P(N),V(N),QORACLE,PP0,KAPPA
      COMMON /ORACLEMODEL/ QORACLE,PP0,KAPPA,MODEL
      V=PP0
      END
      SUBROUTINE PRFUNC(N,P,V,KD)
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMSOL.inc'
      INCLUDE 'COMPHY.inc'
      INTEGER N,KD,MODEL,J
      REAL*8 P(N),V(N),QORACLE,PP0,KAPPA,X
      COMMON /ORACLEMODEL/ QORACLE,PP0,KAPPA,MODEL
      DO J=1,N
         X=MAX(0.D0,1.D0-P(J)/SPSIM)
         V(J)=QORACLE
         IF (MODEL.EQ.2) V(J)=QORACLE*(1.D0+.1D0*X)
         IF (KD.EQ.1) THEN
            V(J)=0.D0
            IF (MODEL.EQ.2) V(J)=2.D0*QORACLE*.1D0*SQRT(X)
         ENDIF
      ENDDO
      END
      SUBROUTINE POLYNM(N,A,V,K)
      INTEGER N,K
      REAL*8 A(*),V(N)
      STOP 'unselected thermal branch'
      END
'''
 (out/'fixture.f').write_text(fixture)
 command=['gfortran','-O0','-g','-std=legacy','-fdefault-real-8','-fdefault-double-8','-ffixed-line-length-none','-mcmodel=medium','-fcheck=all','-I'+str(source.resolve()),'native.f','fixture.f','-o','oracle'];r=subprocess.run(command,cwd=out,capture_output=True,text=True)
 assert r.returncode==0,r.stderr[-2500:]

def expected(model,q,kappa,x):
 c=.1/(np.pi*kappa*.01)
 if model==3:
  f2=4*np.pi**2*q*q*c*c*(1+.05*x)**2;d=4*np.pi**2*q*q*c*c*.1*(1+.05*x)
 else:
  factor=1+.1*x if model==2 else np.ones_like(x);fp=.1 if model==2 else 0.
  f2=4*np.pi**2*q*q*c*c*factor**2*(1-.01*x);d=4*np.pi**2*q*q*c*c*(2*factor*fp*(1-.01*x)-.01*factor**2)
 return f2,-d/(2*(-.1))

def main():
 p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('fixed',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);rows=[];parity={}
 for variant,source in [('parent',a.parent),('fixed',a.fixed)]:
  out=a.output/variant;build(source,out)
  cases=[(4,m,q,pp,k)for m in [1,2,3]for q in [1.5,-1.5]for pp in [0.,.13]for k in [1.,1.7]]+[(sel,1,1.5,pp,1.)for sel in [1,3]for pp in [0.,.13]]
  for i,(sel,m,q,pp,k)in enumerate(cases):
   folder=out/f'case{i}';folder.mkdir();r=subprocess.run([str(out/'oracle'),*map(str,[sel,m,q,pp,k])],cwd=folder,capture_output=True,text=True);assert r.returncode==0,(variant,i,r.stderr)
   lines=r.stdout.splitlines();knots=np.array([[float(x)for x in row.split()[2:]]for row in lines if row.startswith('knot')]);queries=np.array([[float(x)for x in row.split()[2:]]for row in lines if row.startswith('query')]);binary=(folder/'profiles.dat').read_bytes()
   if sel!=4:parity[variant,i]=binary;continue
   f2,ff=expected(m,q,k,knots[:,0]);x=1-queries[:,0]/(-.1);_,sourceff=expected(m,q,k,x);j=-queries[:,1]*pp-sourceff/queries[:,1]
   if variant=='fixed':
    np.testing.assert_allclose(knots[:,1],ff,rtol=2e-10,atol=2e-10);np.testing.assert_allclose(knots[:,2],f2,rtol=2e-12,atol=2e-12);np.testing.assert_allclose(queries[:,2],j,rtol=2e-10,atol=2e-10)
   else:assert not np.isfinite(queries[:2,2]).all()
   rows.append(dict(variant=variant,model=m,q=q,pprime=pp,kappa=k,query_current=queries[:,2].tolist(),expected=j.tolist(),max_finite_query_relative_error=float(np.max(abs(queries[2:,2]-j[2:])/np.maximum(1,abs(j[2:]))))))
 for key in parity:
  if key[0]=='parent':assert parity[key]==parity['fixed',key[1]]
 sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest();report=dict(status='PASS',scope='Whole native ISOFUN and CURENT, existing native SPLINE/TRIDAG:24 circular/elliptical coarea or even-polynomial controls per parent/fixed, q±1.5,pprime0/.13. Actual profile function supplied as independent exact even-polynomial fixture. Axis and negative-x extrapolation current finite/exact after repair; parent nonfinite. Four legacy1/3 native output packets byte-identical. No PDE.',controls=rows,legacy_byteparity_controls=4,source_sha256={str(x):sha(x/'chease.f')for x in [a.parent,a.fixed]},oracle_sha256=sha(__file__))
 (a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status='PASS',q_controls=len(rows),legacy_byteparity_controls=4)))
if __name__=='__main__':main()
