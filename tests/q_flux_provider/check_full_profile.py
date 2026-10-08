"""Whole PROFILE/ISOFUN/CURENT with shifted analytic circle/ellipse coarea."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.io import FortranFile
from chease_q_axis_coarea_oracle import build as coarea_build, routine


PROVIDERS = '''      SUBROUTINE PPRIME(N,P,V)
      INTEGER N,J,QM
      REAL*8 P(N),V(N),PP
      COMMON /QMODEL/ PP,QM
      V=PP
      END
      SUBROUTINE PRFUNC(N,P,V,KD)
      INCLUDE 'DECLAR.inc'
      INCLUDE 'COMDIM.inc'
      INCLUDE 'COMSOL.inc'
      INCLUDE 'COMPHY.inc'
      INTEGER N,KD,QM,J,IMODE
      REAL*8 P(N),V(N),PP,X,KAPPA,A0,BUMP,WIDTH,CC,
     & QSIGN,AXR,AXZ
      COMMON /MODEL/ KAPPA,A0,BUMP,WIDTH,CC,QSIGN,AXR,AXZ,IMODE
      COMMON /QMODEL/ PP,QM
      DO J=1,N
         X=1.D0-P(J)/SPSIM
         V(J)=QSIGN
         IF(QM.EQ.1) V(J)=QSIGN*(1.D0+.1D0*X)
         IF(KD.EQ.1) THEN
            V(J)=0.D0
            IF(QM.EQ.1) V(J)=2.D0*.1D0*QSIGN*SQRT(MAX(X,0.D0))
         ENDIF
      ENDDO
      END
      SUBROUTINE POLYNM(N,P)
      STOP 88
      END
'''


def build(source, output, mesh_axis=False):
    coarea_build(source, output)
    text = (source/'chease.f').read_text()
    names = ['ISOFUN', 'CURENT', 'SPLINE', 'TRIDAG', 'SCOPY', 'SSCAL', 'RESETI']
    if '         SUBROUTINE QFLUXSOURCE(' in text:
        names += ['QFLUXSOURCE']
    with (output/'native.f').open('a') as f:
        f.write(''.join(routine(text, n) for n in names))
    fixture = (output/'oracle.f').read_text()
    if mesh_axis:
        fixture=fixture.replace('      AXZ=RZMAG',
            '      AXZ=RZMAG\n      R0=RMAG\n      RZ0=RZMAG')
    fixture = fixture.replace('      SUBROUTINE ISOFUN(KN)\n      INTEGER KN\n      END\n', '')
    fixture = fixture.replace('      CHARACTER*64 ARG', '''      REAL*8 PP,QX(5),QP(5),QR(5),QJ(5)
      INTEGER QM
      COMMON /QMODEL/ PP,QM
      CHARACTER*64 ARG''')
    fixture = fixture.replace('      CPI=ACOS(-1.D0)', '''      CALL GET_COMMAND_ARGUMENT(6,ARG)
      READ(ARG,*) QM
      CALL GET_COMMAND_ARGUMENT(7,ARG)
      READ(ARG,*) PP
      NSURF=6
      NPROFZ=0
      NTMF0=0
      NRFP=0
      PREDGE=0.D0
      CPI=ACOS(-1.D0)''')
    fixture = fixture.replace('      CALL PROFILE(NISO)', '''      DO J=1,NISO
         CSIPR(J)=SQRT(1.D0-PSIISO(J)/SPSIM)
      ENDDO
      CALL PROFILE(NISO)''')
    fixture = fixture.replace('      CLOSE(17)\n      END', '''      WRITE(17) TTP(1:NISO),TMF(1:NISO)
      QX=(/-1.D-6,0.D0,1.D-8,.09D0,.64D0/)
      QP=SPSIM*(1.D0-QX)
      QR=(/1.001D0,1.D0,1.002D0,1.03D0,1.08D0/)
      CALL CURENT(5,QP,QR,QJ)
      WRITE(17) QX,QR,QJ
      CLOSE(17)
      END''', 1)
    (output/'oracle.f').write_text(fixture+PROVIDERS)
    cmd = ['gfortran','-O0','-g','-fdefault-real-8','-fdefault-double-8',
           '-ffixed-line-length-none','-std=legacy','-mcmodel=medium','-fcheck=all',
           '-I'+str(source.resolve()),'native.f','oracle.f','-o','oracle']
    result = subprocess.run(cmd, cwd=output, capture_output=True, text=True)
    (output/'build.stderr').write_text(result.stderr)
    assert result.returncode == 0, result.stderr[-1500:]


def expected(kappa, q, qm, x):
    # psi=10[(R-Ra)^2+(Z-Za)^2/kappa²-a²]. Analytic coarea.
    c = 10/(np.pi*kappa)
    factor = 1+.1*x if qm else np.ones_like(x)
    d = .1 if qm else 0.
    f2 = (2*np.pi*q*c)**2*factor**2*(1.00125**2-.01*x)
    dx = (2*np.pi*q*c)**2*(2*factor*d*(1.00125**2-.01*x)-.01*factor**2)
    return f2, dx/.2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('parent', type=Path);p.add_argument('fixed', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--mesh-axis', action='store_true',
        help='Place the manufactured mesh origin at the axis; isolate provider from missing-inner tracing.')
    a = p.parse_args();a.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for lane, source in [('parent',a.parent),('fixed',a.fixed)]:
        folder = a.output/lane;build(source,folder,a.mesh_axis)
        for kappa in [1.,1.7]:
            for q in [1.5,-1.5]:
                for qm in [0,1]:
                    for pp in [0.,.13]:
                        run=folder/f'k{kappa}_q{q}_m{qm}_p{pp}';run.mkdir()
                        r=subprocess.run([str((folder/'oracle').resolve()),
                            *map(str,[4,kappa,0,64,q,qm,pp])],cwd=run,
                            capture_output=True,text=True)
                        (run/'stdout.txt').write_text(r.stdout)
                        (run/'stderr.txt').write_text(r.stderr)
                        assert r.returncode==0 and(run/'coarea.dat').exists()
                        with FortranFile(run/'coarea.dat')as f:
                            coarea=f.read_record('<f8').reshape(4,100)
                            source_data=f.read_record('<f8').reshape(2,100)
                            query=f.read_record('<f8').reshape(3,5)
                        x=((np.arange(100)+.5)/100)**2
                        f2,ff=expected(kappa,q,qm,x)
                        _,qff=expected(kappa,q,qm,query[0])
                        current=-query[1]*pp-qff/query[1]
                        if lane=='fixed':
                            np.testing.assert_allclose(source_data[0],ff,rtol=2e-7,atol=2e-5)
                            np.testing.assert_allclose(source_data[1]**2,f2,rtol=2e-11)
                            np.testing.assert_allclose(query[2],current,rtol=2e-7,atol=2e-5)
                        else:
                            assert not np.isfinite(query[2,:2]).all()
                        rows.append(dict(lane=lane,kappa=kappa,signed_q=q,
                            q_model='linear_x'if qm else'constant',pprime=pp,
                            current_finite=bool(np.isfinite(query[2]).all()),
                            max_FFprime_error=float(np.max(abs(source_data[0]-ff))),
                            primitive_F2_relative_error=float(np.max(abs(source_data[1]**2/f2-1)))))
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    report=dict(status='PASS',controls=rows,script_sha256=sha(__file__),
        source_sha256={k:sha(s/'chease.f')for k,s in [('parent',a.parent),('fixed',a.fixed)]},
        mesh_origin_at_axis=a.mesh_axis,
        scope='Whole native PROFILE/ISOFUN/CURENT and GAUSS/SPLINE. Analytic shifted quadratic circle/ellipse and exact outer CINT fixture; signed q and finite pressure derivatives test the source contract, not force-balanced equilibria. With mesh-axis, all contours are traced by the exact CINT adapter and the missing-inner helper is unused. Native F² value and derivative follow one polynomial in x. Parent has axis/underaxis nonfinite. No PDE or full inverse-map convergence claim.')
    (a.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','controls':len(rows),'receipt_sha256':sha(a.output/'receipt.json')}))


if __name__=='__main__':main()
