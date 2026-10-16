#!/usr/bin/env python3
"""Native PAMS regression: stable solves of well-conditioned block systems.

The test links unmodified pams.f and its bundled smlpck.f. It requires gfortran
and BLAS. No equilibrium data or solver-specific physical assumptions enter the
manufactured system. An uncorrected restricted-pivot implementation fails this
test. Source revisions may be supplied to reproduce a historical failure.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

EPSILONS = (1e-8, 1e-10, 1e-12, 1e-14, 1e-16)
FLAGS = ['-O2', '-fdefault-real-8', '-fdefault-double-8',
         '-ffixed-line-length-none', '-ffree-line-length-none',
         '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections']


def driver(source):
    """Use the native dummy declarations and verbatim numerical procedures."""
    start = source.index('         SUBROUTINE REDUCY(')
    end = source.index('      REAL*8            RNORM', start)
    header = source[start:end]
    signature = ''.join(line for line in header.splitlines()
                        if not line.startswith('C'))
    args = signature[signature.index('(') + 1:signature.index(')')]
    args = [arg.strip() for arg in args.replace('$', '').split(',')]
    declarations = header[header.index('      IMPLICIT LOGICAL'):]
    declarations = '\n'.join(line for line in declarations.splitlines()
                             if not line.startswith('C'))
    declarations = declarations.replace(
        '         INTEGER   MD,ND,MDY,NXC,NYC,NCASE,NITMAX',
        'INTEGER,PARAMETER::MD=1,ND=1,MDY=1,NXC=1,NYC=1,NCASE=2,NITMAX=1')

    def call(name):
        return ('call ' + name + '( &\n' + ', &\n'.join(
            ', '.join(args[i:i + 8]) for i in range(0, len(args), 8)) + ')\n')

    return '''program native_pams_manufactured
''' + declarations + '''
INTEGER NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG,NEV,M,N,MY
COMMON /AUXINT/ NEG,NEGEIG,NIT,NXDEV,NYDEV,NXPNEG,NYPNEG,NEV,M,N,MY
INTEGER INFO,IPIV(3)
CHARACTER*100 ARG
REAL*8 EPS,BERR
COMPLEX*16 FULL(3,3),COPY(3,3),RHS(3),EXACT(3),GOT(3),SOL(3,1)
INTEGER,PARAMETER::QP=SELECTED_REAL_KIND(30)
COMPLEX(KIND=QP) AQ(3,3),BQ(3),RQ(3)
call get_command_argument(1,arg)
read(arg,*)eps
m=1
my=1
n=1
full(1,:)=[(1d0,0d0),(1d0,0d0),(0d0,0d0)]
full(2,:)=[(1d0,0d0),cmplx(eps,0d0,kind=8),(1d0,0d0)]
full(3,:)=[(0d0,0d0),(1d0,0d0),(1d0,0d0)]
exact=[(1d0,0d0),(2d0,0d0),(3d0,0d0)]
rhs=matmul(full,exact)
aq=cmplx(real(full,kind=8),aimag(full),kind=qp)
bq=cmplx(real(rhs,kind=8),aimag(rhs),kind=qp)
copy=full
sol(:,1)=rhs
call zgetrf(3,3,copy,3,ipiv,info)
if(info/=0)stop 11
call zgetrs('N',3,1,copy,3,ipiv,sol,3,info)
if(info/=0)stop 12
write(*,'(a,es24.16)')'GLOBAL_LU_MAXERR ',maxval(abs(sol(:,1)-exact))
a=0
b=0
c=0
d=0
e=0
f=0
g=0
h=0
b(1,1,1:2)=1
d(1,1,1)=eps
e(1,1,1)=1
f(1,1,1)=1
g(1,1,1)=1
h(1,1,1)=1
! Zero all inverse workspaces, including upstream's unused extra RHS.
sc1=0
scy1=0
epsdet=1d-199
''' + call('reducy') + call('reducx') + '''
dx(1,1)=rhs(1)
dy(1,1)=rhs(2)
dx(1,2)=rhs(3)
''' + call('redy') + call('redx') + call('solvex') + call('solvey') + '''
got=[x(1,1),y(1,1),x(1,2)]
rq=bq-matmul(aq,cmplx(real(got,kind=8),aimag(got),kind=qp))
berr=maxval(abs(rq)/(matmul(abs(aq),abs(cmplx(real(got,kind=8),aimag(got),kind=qp)))+abs(bq)))
write(*,'(a,3es24.16)')'PAMS_RESULT ',maxval(abs(got-exact)),maxval(abs(rq)),berr
write(*,*)'COMPLETED_NATIVE_FIXTURE'
end program
'''


def checked(command, cwd):
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True,
                            timeout=60)
    if result.returncode:
        raise RuntimeError(' '.join(map(str, command)) + '\n' + result.stderr)
    return result.stdout


def evaluate(source, revision, work):
    provenance = {}
    for name in ('pams.f', 'smlpck.f'):
        if revision:
            data = subprocess.check_output(
                ['git', '-C', str(source), 'show', revision + ':MarsQ_2FK/' + name])
        else:
            data = (source / 'MarsQ_2FK' / name).read_bytes()
        (work / name).write_bytes(data)
        provenance[name] = hashlib.sha256(data).hexdigest()
    (work / 'fixture.f90').write_text(driver((work / 'pams.f').read_text()))
    checked(['gfortran', *FLAGS, '-c', 'pams.f', 'smlpck.f'], work)
    checked(['gfortran', *FLAGS, 'fixture.f90', 'pams.o', 'smlpck.o',
             '-Wl,--gc-sections', '-lblas', '-o', 'fixture'], work)
    cases = []
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    for eps in EPSILONS:
        result = subprocess.run([str(work / 'fixture'), str(eps)], cwd=work,
                                env=environment, text=True, capture_output=True,
                                timeout=10)
        eigenvalues = [1.0, (1 + eps + math.sqrt((1 - eps)**2 + 8)) / 2,
                       (1 + eps - math.sqrt((1 - eps)**2 + 8)) / 2]
        condition = max(map(abs, eigenvalues)) / min(map(abs, eigenvalues))
        lines = result.stdout.splitlines()
        oracle = [float(line.split()[1]) for line in lines
                  if line.startswith('GLOBAL_LU_MAXERR ')]
        native = [list(map(float, line.split()[1:])) for line in lines
                  if line.startswith('PAMS_RESULT ')]
        complete = 'COMPLETED_NATIVE_FIXTURE' in result.stdout
        errors = []
        if condition >= 2.000001 or len(oracle) != 1 or oracle[0] > 1e-14:
            errors.append('independent well-conditioned global-LU oracle failed')
        # Fortran STOP can return status0; require the completion marker too.
        if result.returncode != 0 or not complete or len(native) != 1:
            errors.append('native solve did not complete')
        elif not all(math.isfinite(x) for x in native[0]):
            errors.append('native solve returned nonfinite values')
        elif native[0][0] > 1e-11 or native[0][2] > 1e-12:
            errors.append('native forward/backward accuracy gate failed')
        cases.append(dict(eps=eps,condition2=condition,global_lu_maxerr=oracle,
                          pams_result=native,completed=complete,errors=errors,
                          returncode=result.returncode,stdout=result.stdout,
                          stderr=result.stderr))
    pin = checked(['git', '-C', str(source), 'rev-parse', revision or 'HEAD'],work).strip()
    return dict(revision=pin,source_sha256=provenance,cases=cases,
                passed=all(not c['errors'] for c in cases))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--revision', action='append', help='git revision; repeatable')
    parser.add_argument('--build-dir', type=Path,
                        default=Path.home()/'.cache'/'mars-pams-stability')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    args.build_dir.mkdir(parents=True,exist_ok=True)
    reports = []
    for revision in args.revision or [None]:
        with tempfile.TemporaryDirectory(prefix='native-', dir=args.build_dir) as root:
            reports.append(evaluate(args.source.resolve(),revision,Path(root)))
    receipt = dict(compiler=checked(['gfortran','--version'],args.source).splitlines()[0],
                   flags=FLAGS,threads=1,reports=reports)
    if args.report:
        args.report.write_text(json.dumps(receipt,indent=2)+'\n')
    for report in reports:
        print(report['revision'])
        for case in report['cases']:
            print(case['eps'], 'PASS' if not case['errors'] else 'FAIL',
                  '; '.join(case['errors']), case['pams_result'])
    return 0 if all(r['passed'] for r in reports) else 1


if __name__ == '__main__':
    raise SystemExit(main())
