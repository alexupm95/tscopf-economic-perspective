"""
Section IV-A - angles relative to the COI of the three generators of the 9-bus TSCOPF run 112948
(short circuit at bus 7, branch 5-7 opened, 5 ms step, 100 ms fault).

Writes figures/angles_tscopf_9bus.{pdf,svg}.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent
RUN = ROOT.parent / "results" / "9bus" / "Results - 2026-10-08 112948"
ang = pd.read_csv(RUN / "Transient_Stability" / "CSV" / "angle_rel_COI.csv", sep=";")
COLORS = {"G1": "#9671bd", "G2": "#378d94", "G3": "#C51212"}

plt.rcParams.update({
    "font.family": "Courier New",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Courier New",
    "mathtext.it": "Courier New:italic",
    "mathtext.bf": "Courier New:bold",
    "mathtext.default": "regular",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.size": 24,
})

fig, ax = plt.subplots(figsize=(10.0, 6.0), layout="constrained")
ax.axvspan(0.01, 0.11, color="0.85", zorder=0)               # fault: 0.010 to 0.105 s, line open from 0.110 s
for lim in (-100, 100):
    ax.axhline(lim, color="0.35", linewidth=1.2, linestyle=(0, (3, 2)), zorder=1)
for g, c in COLORS.items():
    ax.plot(ang.t, ang[g], color=c, linewidth=2.2, label=g, zorder=3)
ax.set_xlim(0, 5)
ax.set_ylim(-125, 125)
ax.set_yticks([-100, -50, 0, 50, 100])
ax.set_xlabel("time [s]", fontsize=24)
ax.set_ylabel("angles\nw.r.t. COI [deg]", fontsize=24)
ax.tick_params(labelsize=22)
ax.grid(axis="y", linewidth=0.5, alpha=0.5)
ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.legend(frameon=False, fontsize=22, ncol=3, loc="outside upper center", handlelength=1.3, columnspacing=1.0,
           handletextpad=0.5)

# the limit is reached by G3 only, at the two instants the text quotes
hit = ang.t[ang["G3"] > 100 - 1e-4].round(3).tolist()   # 1e-4: at 0.285 s G3 is at 99.9997 deg with a zero dual
assert ang[["G1", "G2"]].abs().max().max() < 99 and hit == [0.29, 0.765], hit
print("G3 at the upper limit at t =", hit)

FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)
for ext in ("svg", "pdf"):
    fig.savefig(FIGDIR / f"angles_tscopf_9bus.{ext}", bbox_inches="tight")
print("written:", FIGDIR / "angles_tscopf_9bus.pdf")
