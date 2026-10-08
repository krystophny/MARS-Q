"""Generate self-contained circular/shaped prescribed-q native CHEASE inputs."""
import argparse
import math
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--variant',choices=['public','mars'],required=True)
    args=parser.parse_args();args.output.mkdir(exist_ok=False,parents=True)
    selector=5 if args.variant=='public'else 4
    for name,kappa in [('circular',1.),('shaped',1.7)]:
        folder=args.output/name;folder.mkdir()
        lines=['0.1','0','0','513']
        for i in range(513):
            t=2*math.pi*i/512
            lines.append(f'{1+.1*math.cos(t):.16e} {kappa*.1*math.sin(t):.16e}')
        lines.append('101 4')
        if args.variant=='public':lines.append('5')
        lines.extend(f'{i/100:.16e}'for i in range(101))
        lines.extend('0.0000000000000000e+00'for _ in range(101))
        lines.extend('1.5000000000000000e+00'for _ in range(101))
        (folder/'EXPEQ').write_text('\n'.join(lines)+'\n')
        deck=f'''&EQDATA
 NEQDSK=0, NSURF=6, NTCASE=0, NIDEAL=6,
 NBLOPT=0, NBSOPT=0, CPRESS=1., CFBAL=1.,
 NCSCAL=4, NTMF0=0, NSTTP={selector}, NFUNC=4, NPPFUN=4,
 NIPR=1, NISO=100, NPP=1, NPPR=30, NSOUR=2, NPROPT=2,
 NS=32, NT=32, NPSI=128, NCHI=128,
 NV=160, REXT=6., NVEXP=3, R0W=.9, RZ0W=0.,
 NMESHA=0, NEGP=-1, NER=1, EPSLON=1.E-10,
 NOPT=0, NPLOT=1, NBAL=0, B0EXP=5.3, R0EXP=6.2,
 RELAX=.3, NINMAP=100, NINSCA=100, NRBOX=129, NZBOX=129,
'''
        if args.variant=='public':
            deck+=' COCOS_IN=2, COCOS_OUT=2, SIGNB0XP=-1, SIGNIPXP=1,\n TENSBND=0., TENSPROF=0., NFUNRHO=0, NVERBOSE=2,\n &END\n***\n*** Generated prescribed-q circular/shaped control\n***\n***\n'
        else:deck+=' &END\n &NEWRUN\n /\n'
        (folder/'chease_namelist').write_text(deck)
    print('Prepared circular and elongated prescribed-q inputs; no solver executed.')


if __name__=='__main__':
    main()
