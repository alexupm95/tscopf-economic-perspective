"""Bar plot of the nodal energy component lambda_k at every bus, co-located case (run 113323, 5 ms step, 100 ms fault).

Plots the printed values (0.1 $/MWh, as in tables/tab_colocated_bus1_9bus.tex; run tables_9bus.py first).
Writes figures/lambda_colocated_bus1_9bus.{pdf,svg}.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "results" / "9bus" / "Results - 2026-10-08 113323"

plt.rcParams.update({
    "font.family": "Courier New",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Courier New",
    "mathtext.it": "Courier New:italic",
    "mathtext.bf": "Courier New:bold",
    "mathtext.default": "regular",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.size": 24, "axes.labelsize": 26, "xtick.labelsize": 21, "ytick.labelsize": 21,
})

# this run has no kkt_Pd_prices.csv: lambda_k straight from the dual export
dual = pd.read_csv(RUN / "Dispatch" / "CSV_duals" / "dual_P_balance.csv", sep=";")
df = pd.DataFrame({"bus": dual["Bus_ID"], "energy": dual.set_index("Bus_ID").iloc[:, 0].to_numpy()})
assert (df["energy"] > 0).all()
lam = (df["energy"] * 10 + 0.5).astype(int) / 10  # half-up to 0.1, as printed

# check against the printed table
table = (HERE / "tables" / "tab_colocated_bus1_9bus.tex").read_text(encoding="utf-8")
printed = {float(x) for x in re.findall(r"\\rowhl\{(\d+\.\d0)\}|& (\d+\.\d0) & \$", table) for x in x if x}
for bus in (1, 2, 3, 5, 6, 8):
    assert lam[df["bus"] == bus].iloc[0] in printed, f"bus {bus} not as printed in the table"

fig, ax = plt.subplots(figsize=(10, 3.6))  # short panel: every bar sits in 40-50, zero baseline kept
bars = ax.bar(df["bus"].astype(str), lam, color="#378d94", width=0.65)
ax.bar_label(bars, fmt="%.1f", fontsize=22, padding=3)
ax.set_xlabel("bus")
ax.set_ylabel(r"$\lambda_k$ [\$/MWh]")
ax.set_ylim(0, 60)
ax.set_yticks([0, 20, 40])
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", alpha=0.3)
ax.set_axisbelow(True)
fig.tight_layout()

out = HERE / "figures" / "lambda_colocated_bus1_9bus"
out.parent.mkdir(exist_ok=True)
fig.savefig(out.with_suffix(".pdf"))
fig.savefig(out.with_suffix(".svg"))
print(lam.tolist())
