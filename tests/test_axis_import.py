#!/usr/bin/env python3
"""Axis slots of a perturbation imported with KPERTREAD=1 (docs/axis_regularity.md).

READPERTURB, AXISIMPORT, IAXPWX and IAXPWY are compiled unchanged against small
fixture modules. The imported files carry a field that is a regular power of s on
node 2 and half-cell 2 outward, but node 1 and half-cell 1 follow another axis
convention (B3 there holds the B2 value, X2 and B2 vanish), as in the TC24
supplier files. After the import the axis slots must be those BOUNDC and AXISY
give a solved field: the power law at CSM(1), zero for positive node powers, and
B2 from div b = 0 on the axis cell. The powers come from the independent tables
of test_axis_regularity.
"""
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

from test_axis_regularity import PX, PY, deck, power

DECKS = ('READPERTURB', 'AXISIMPORT', 'IAXPWX', 'IAXPWY')
NR, NV, N = 8, 3, -3.0
MODES = list(range(-3, 4))

MODULES = '''module dimensim
implicit none
integer NR,NRP1,NTP1,MSMAX
end module
module globalm
implicit none
integer KXV1,KXB1,KXJ2U,KXJ3,KXJ2L,KXX1,KXW1,KXJRE2,KXJRE3,KXJRE2L,KXB2L
integer KYV2,KYB2,KYB3,KYJ1,KYPR,KYX2,KYX3,KYW2,KYJRE1
real*8 RNTOR,SLEFT
real*8,allocatable :: RM(:,:),CS(:),CSM(:),CSH(:)
complex*16,allocatable :: B1U(:,:),B2U(:,:),B3U(:,:),X1U(:,:),X2U(:,:),X3U(:,:)
end module
module rcomdm
implicit none
real*8,allocatable :: DPSIDS(:),T(:)
end module
'''

DRIVER = '''program axis_import
use dimensim
use globalm
use rcomdm
implicit none
integer :: i,ms
real*8 :: s(@NRP1@)
NR=@NR@;NRP1=NR+1;NTP1=NRP1+@NV@;MSMAX=@MSMAX@;RNTOR=@N@;SLEFT=0
KXV1=1;KXB1=2;KXJ2U=3;KXJ3=4;KXJ2L=5;KXX1=6;KXW1=-1;KXJRE2=-1;KXJRE3=-1;KXJRE2L=-1;KXB2L=-1
KYV2=1;KYB2=2;KYB3=3;KYJ1=4;KYPR=5;KYX2=6;KYX3=7;KYW2=-1;KYJRE1=-1
allocate(RM(MSMAX,2),CS(NRP1),CSM(NR),CSH(NR),DPSIDS(NRP1),T(NRP1))
allocate(B1U(NTP1,MSMAX),B2U(NTP1,MSMAX),B3U(NTP1,MSMAX))
allocate(X1U(NRP1,MSMAX),X2U(NRP1,MSMAX),X3U(NRP1,MSMAX))
s=[@CS@]
CS=s;CSM=0.5d0*(CS(1:NR)+CS(2:NRP1));CSH=CS(2:NRP1)-CS(1:NR)
do ms=1,MSMAX
RM(ms,:)=ms-1-(MSMAX-1)/2
end do
DPSIDS=1;T=1
call READPERTURB
do ms=1,MSMAX
do i=1,2
write(*,'(12es25.16)') B1U(i,ms),B2U(i,ms),B3U(i,ms),X1U(i,ms),X2U(i,ms),X3U(i,ms)
end do
end do
end program
'''


def grids():
    cs = (np.arange(NR + 1) / NR) ** 1.3
    return cs, 0.5 * (cs[:-1] + cs[1:]), cs[1:] - cs[:-1]


def regular(m):
    """Coefficients and powers of a field regular at the axis, for harmonic m."""
    rng = np.random.default_rng(10 + m)
    c = rng.normal(size=6) + 1j * rng.normal(size=6)
    p = [power(PX, 'B1', m), power(PY, 'B2', m), power(PY, 'B3', m),
         power(PX, 'X1', m), power(PY, 'X2', m), power(PY, 'X3', m)]
    return c, p


def write_inputs(work):
    cs, csm, _ = grids()
    xb = np.concatenate([cs, cs[-1] + np.arange(1, NV + 1) * (cs[-1] - cs[-2])])
    hb = np.concatenate([csm, xb[NR + 1:] - 0.5 * (cs[-1] - cs[-2]), [xb[-1]]])
    header = [f'{len(MODES)} {{}} {N} 0 0 0'] + [' '.join([str(float(m))] * 6) for m in MODES]
    brows, xrows = [], []
    for m in MODES:
        c, p = regular(m)
        b = np.array([c[0] * xb ** p[0], c[1] * hb ** p[1], c[2] * hb ** p[2]])
        x = np.array([c[3] * cs ** p[3], c[4] * np.append(csm, 1.0) ** p[4],
                      c[5] * np.append(csm, 1.0) ** p[5]])
        # supplier axis convention: B3 slot carries the B2 value, B2 and X2 vanish
        b[:, 0] = [0, 0, b[1, 1]]
        x[:, 0] = [0, 0, 0]
        brows += [' '.join(f'{v.real:.17e} {v.imag:.17e}' for v in b[:, i]) for i in range(b.shape[1])]
        xrows += [' '.join(f'{v.real:.17e} {v.imag:.17e}' for v in x[:, i]) for i in range(x.shape[1])]
    (work / 'BPLASMA_INPUT').write_text('\n'.join([header[0].format(NR + 1 + NV)] + header[1:] + brows) + '\n')
    profile = ['1 1 1 1 1 1'] * (NR + 1)
    (work / 'XPLASMA_INPUT').write_text('\n'.join([header[0].format(NR + 1)] + header[1:] + profile + xrows) + '\n')


def expected(m):
    """Axis slots of a solved field: BOUNDC (nodes) and AXISY (half-cell 1)."""
    cs, csm, csh = grids()
    c, p = regular(m)
    node2 = [c[0] * cs[1] ** p[0], c[3] * cs[1] ** p[3]]
    b1, x1 = [0 if q > 0 else v for q, v in zip((p[0], p[3]), node2)]
    b3 = c[2] * csm[0] ** p[2]
    b2 = (-((node2[0] - b1) / csh[0] + 1j * N * b3) / (1j * m) if m else c[1] * csm[0] ** p[1])
    return np.array([b1, b2, b3, x1, c[4] * csm[0] ** p[4], c[5] * csm[0] ** p[5]])


def run(source, build_parent):
    src = source.resolve() / 'MarsQ_2FK'
    text = (src / 'marsq.f').read_text()
    native = ''.join(deck(text, name) for name in DECKS)
    cs, _, _ = grids()
    driver = (DRIVER.replace('@NRP1@', str(NR + 1)).replace('@NR@', str(NR))
              .replace('@NV@', str(NV)).replace('@MSMAX@', str(len(MODES)))
              .replace('@N@', f'{N}d0').replace('@CS@', ','.join(f'{v:.17e}' for v in cs)))
    with tempfile.TemporaryDirectory(prefix='axis-import-', dir=build_parent) as directory:
        work = Path(directory)
        (work / 'modules.f90').write_text(MODULES)
        (work / 'native.f').write_text(native)
        (work / 'driver.f90').write_text(driver)
        command = ['gfortran', '-O0', '-ffixed-line-length-none', '-ffree-line-length-none',
                   '-fdefault-real-8', '-fdefault-double-8', '-fcheck=bounds',
                   '-ffpe-trap=zero,invalid', 'modules.f90', 'native.f', 'driver.f90', '-o', 'test']
        compiled = subprocess.run(command, cwd=work, capture_output=True, text=True)
        if compiled.returncode:
            raise RuntimeError((compiled.stdout + compiled.stderr)[-8000:])
        write_inputs(work)
        result = subprocess.run([str(work / 'test')], cwd=work, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stdout[-4000:] + result.stderr)
    rows = [line for line in result.stdout.splitlines() if 'MARS-K' not in line and 'L2' not in line]
    v = np.array(' '.join(rows).split(), float).reshape(len(MODES), 2, 6, 2)
    return v[..., 0] + 1j * v[..., 1]


class TestAxisImport(unittest.TestCase):
    def test_imported_axis_slots(self):
        got = run(Path(__file__).resolve().parents[1], Path.cwd())
        names = ['B1', 'B2', 'B3', 'X1', 'X2', 'X3']
        for k, m in enumerate(MODES):
            want = expected(m)
            scale = np.maximum(np.abs(got[k, 1]), np.abs(want))
            for j, name in enumerate(names):
                with self.subTest(m=m, field=name):
                    self.assertLessEqual(abs(got[k, 0, j] - want[j]), 1e-12 * scale[j] + 1e-300)


if __name__ == '__main__':
    unittest.main()
