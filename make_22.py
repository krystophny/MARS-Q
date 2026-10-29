"""Circular benchmark: MARS-K NTV torque with master and corrected coefficient vs PENTRC and MARS Shaing.

The master curve is the corrected run times 2 (the coefficient is the only difference).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STORE = Path("/mnt/storage/codex-ntv-20261014")
OUT = Path(__file__).resolve().parent / "22.png"
SURF = ("0.04", "0.09", "0.16", "0.2", "0.36")

mk = np.array([(d["eps"], d["ratio"], d["ratio_upstream_coeff"]) for d in (
    json.load(open(STORE / f"ntv_s1_frozen_marsk_20261014_02/out/res_s{s}/analyze.log")) for s in SURF)])
sh = np.array([(d["eps"], d["ratio"]) for d in (
    json.load(open(STORE / f"ntv_s1_frozen_mars_shaing_20261014_02/out/sh_s{s}/analyze.log")) for s in SURF)])
pe = np.array([(r["eps"], r["ratio"]) for r in json.load(open(
    STORE / "ntv_s1_frozen_pentrc_radial_20261014_01/out/sbp0_256/verdict.json"))["rows"] if r["applicable"]])
print("master", mk[:, 2].round(3), "fixed", mk[:, 1].round(3), "shaing", sh[:, 1].round(3), "pentrc", pe[:, 1].round(3))

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
fig, a = plt.subplots(figsize=(5, 3))
a.axhline(1.0, color="k", lw=0.8)
a.plot(mk[:, 0], mk[:, 2], "v-", color="#D55E00", label="MARS-K, master (-2n Im DWK)")
a.plot(mk[:, 0], mk[:, 1], "^-", color="#0072B2", label="MARS-K, this PR (-n Im DWK)")
a.plot(sh[:, 0], sh[:, 1], "s--", color="#009E73", ms=4, label="MARS Shaing formula (KNTV=10)")
a.plot(pe[:, 0], pe[:, 1], "D:", color="#CC79A7", ms=3.5, label="PENTRC")
a.set_xlabel(r"inverse aspect ratio $\epsilon$ (dimensionless)")
a.set_ylabel("NTV torque / analytic reference")
a.set_ylim(0.6, 2.3)
a.set_title("Circular tokamak, n=3, superbanana-plateau regime", fontsize=9)
a.legend(fontsize=7, loc="center right")
fig.tight_layout()
fig.savefig(OUT, dpi=150)
print("wrote", OUT)
