"""Reference spectrum for tests/test_ntv_hamada_gauge.py case (b).

Usage (rmp_torque checkout at c3614334 or later, which provides
rmp_torque.boozer_gauge from PR #55)::

    PYTHONPATH=$RMP_TORQUE/src python3 tools/ntv_hamada_gauge/rmp_torque_reference.py

On a Hamada base (J_H constant, B.grad = D/J_H with D = iota d_theta + d_phi)
the Hamada gauge div xi = 0 gives D w = -d for w = lam/J_H and
dB_H = dB_L(xi_perp) + w D(B0).  This is rmp_torque's delta_b_boozer kernel
(Eqs. 4.6-4.8) with source T = d, so the reference is

    delta_b_boozer(dB_L, d, B0, n, iota, sigma=|iota| sigma_q)

with the fixture of CASE_B: b.deltaB = B0 C2 (dm6 = -dm7**2, gg1 = gg2 = 0),
dB_L = B0 C2 + xi_perp.grad B0 and d = -dB_L/B0 - kappa.xi_perp
= -C2 - 2 xi_perp.grad B0/B0 - cp X1/B0**2.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from rmp_torque import boozer_gauge as bg

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
from test_ntv_hamada_gauge import CASE_B, FIX  # noqa: E402


def main() -> int:
    c = CASE_B
    s = np.array(c["s"])
    i = c["surface"] - 1
    modes = np.array(c["modes"])
    a = np.array([complex(*v) for v in c["a"]])
    b = np.array([complex(*v) for v in c["b"]])
    x1 = bg.Spectrum(modes, a * s[i] + b)
    c2 = bg.Spectrum(modes, c["fr"] * a.imag * s[i] - b.real
                     + 1j * (c["fi"] * a.real + b.imag * s[i]))
    b0 = bg.Spectrum(np.array([-1, 0, 1]),
                     np.array([-0.5, 1.0, -0.5]) * np.array([c["eps"], 1, c["eps"]]) * c["b00"])
    dbds = bg.Spectrum(np.array([-1, 0, 1]), np.array([0.5 * c["c1"], c["c0"], 0.5 * c["c1"]]))
    top, pts = 40, 2048
    bc2 = bg.grid_product([b0, c2], [1, 1], top, pts)
    gb = bg.grid_product([x1, dbds], [1, 1], top, pts)
    gbb = bg.grid_product([x1, dbds, b0], [1, 1, -1], top, pts)
    x1b2 = bg.grid_product([x1, b0], [1, -2], top, pts)
    c2t = bg.grid_product([c2], [1], top, pts)
    dbl = bg.Spectrum(bc2.modes, bc2.values + gb.values)
    div = bg.Spectrum(bc2.modes, -c2t.values - 2 * gbb.values - c["cp"] * x1b2.values)
    iota = c["iota"]
    sigma_q = abs(c["n"] * c["dqds"]) * (s[i + 1] - s[i - 1]) * 0.5
    sigma = bg.park_sigma_in_iota_units(sigma_q, iota)
    out = bg.delta_b_boozer(dbl, div, b0, c["n"], iota, sigma=sigma)
    k = np.arange(-c["m_max"], c["m_max"] + 1)
    vals = out.on(k)
    src = Path(bg.__file__).resolve().parents[2]
    rev = subprocess.run(["git", "-C", str(src), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    FIX.write_text(json.dumps({
        "generator": "tools/ntv_hamada_gauge/rmp_torque_reference.py",
        "rmp_torque_revision": rev, "case": c, "sigma_q": sigma_q,
        "delta_b_hamada": {str(int(m)): [v.real, v.imag] for m, v in zip(k, vals)}},
        indent=1) + "\n")
    print("wrote", FIX, "max|dB_H|", np.abs(vals).max())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
