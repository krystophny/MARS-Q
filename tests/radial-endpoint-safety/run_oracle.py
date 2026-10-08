"""Native GAUSS/BASIS2/SETUPA circular-axis endpoint arithmetic, no PDE solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import re


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--mars',action='store_true')
    p.add_argument('--guard',action='store_true')
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    source=args.source
    if args.mars:
        text=(source/'chease.f').read_text()
        for routine in ['GAUSS','BASIS2']:
            begin=text.index('         SUBROUTINE '+routine+'(')
            end=text.index('\nC*DECK',begin)
            (args.output/(routine.lower()+'.f')).write_text(text[begin:end])
        for f in ['HERMIT.inc','DECLAR.inc']:
            shutil.copyfile(source/f,args.output/f)
        setup=text
        init=text[text.index('         SUBROUTINE INITIA'):]
        kernels=['gauss.f','basis2.f']
    else:
        for f in ['gauss.f90','basis2.f90','HERMIT.inc']:
            shutil.copyfile(source/f,args.output/f)
        setup=(source/'setupa.f90').read_text()
        init=(source/'initia.f90').read_text()
        kernels=['gauss.f90','basis2.f90']
    (args.output/'globals.f90').write_text('module globals\ninteger,parameter::rkind=kind(1.d0)\nend module\n')
    arithmetic=next(line.strip()for line in setup.splitlines()if line.strip().startswith('ZV(J8, 1) ='))
    code='''program circular_axis_rule
use globals
use,intrinsic::ieee_arithmetic
implicit none
integer::k,nsgaus,j8,j22,j12
real(rkind)::xr(20),wr(20),s1(1),s2(1),t1(1),t2(1),s(1),t(1)
real(rkind)::zdbds(1,16),zdbdt(1,16),rsint(1,1),zv(1,16),zfrac
character(20)::arg
call get_command_argument(1,arg)
read(arg,*)nsgaus
GUARDPLACE
call gauss(nsgaus,xr,wr)
s1=0
s2=.1_rkind
t1=0
t2=.2_rkind
t=.1_rkind
do k=1,nsgaus
 s=(xr(k)+1)/2*s2
 rsint(1,1)=s(1)
 call basis2(1,1,s1,s2,t1,t2,s,t,zdbds,zdbdt)
 zfrac=0
 j8=1
 j22=1
 j12=1
 ARITHMETICPLACE
 if(.not.ieee_is_finite(zv(1,1)))error stop 3
end do
print *, 'Native SETUPA axis radial-rule finite arithmetic PASS',nsgaus
end program
'''.replace('ARITHMETICPLACE',arithmetic).replace('GUARDPLACE','call axis_guard(nsgaus)'if args.guard else '')
    guard=[]
    if args.guard:
        begin=re.search(r'^.*IF\s*\(NSGAUS\s*\.EQ\.\s*1\)',init,re.M).start()
        end=re.search(r'^.*END\s*IF',init[begin:],re.M).end()+begin
        block=init[begin:end]
        suffix='.f'if args.mars else '.f90'
        guard_file=args.output/('axis_guard'+suffix)
        guard_file.write_text('      SUBROUTINE AXIS_GUARD(NSGAUS)\n      IMPLICIT NONE\n      INTEGER NSGAUS\n'+block+'\n      END\n')
        guard=[guard_file.name]
    (args.output/'oracle.f90').write_text(code)
    build=subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-none',
        '-ffixed-line-length-none','-fdefault-real-8','-fdefault-double-8','-std=legacy',
        'globals.f90',*kernels,*guard,'oracle.f90','-o','oracle'],
        cwd=args.output,capture_output=True,text=True)
    if build.returncode:raise RuntimeError(build.stderr[-2000:])
    rows=[]
    for order,expected in [(1,1 if args.guard else 3),(2,0),(4,0),(8,0)]:
        result=subprocess.run(['./oracle',str(order)],cwd=args.output,capture_output=True,text=True)
        assert result.returncode==expected,(order,result)
        rows.append(dict(radial_request=order,returncode=result.returncode,expected=expected,
            stdout=result.stdout.strip(),stderr=result.stderr.strip()))
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    receipt=dict(scope='Native exact BASIS2 derivatives and SETUPA angular-gradient division evaluated on circular first radial element; radial trapezoid endpoint sigma0 is nonfinite. No equilibrium solve or fitted geometry.',
        controls=rows,guarded=args.guard,mars=args.mars,
        source_sha256={f:sha(source/f)for f in (['chease.f','HERMIT.inc','DECLAR.inc']if args.mars else ['gauss.f90','basis2.f90','HERMIT.inc','setupa.f90','initia.f90'])},
        adapter_sha256=sha(__file__),fixture_sha256=sha(args.output/'oracle.f90'))
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(native_radial_rule1_axis_nonfinite=not args.guard,guarded=args.guard,mars=args.mars,interior_rules_pass=[2,4,8])))


if __name__=='__main__':main()
