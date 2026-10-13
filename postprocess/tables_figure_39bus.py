"""
Section IV-B on the single run with the converter P/Q overload loops on (100 ms fault).

Reads `results/39bus/Results - 2026-10-07 124340` (nothing is solved here) and the PSS/E replay of its
dispatch (`psse/aligned.csv`), and writes
  tables/tab_prices_39bus_gfm.tex            generators, ACOPF against TSCOPF
  tables/tab_prices_profit_39bus_gfm.tex     same, with the profit of each generator
  tables/tab_settlement_39bus_gfm.tex        settlement at the shadow prices
  figures/angles_psse_39bus_gfm.{pdf,svg}    angle to the COI of the ten units, TSCOPF and PSS/E

The prices are rebuilt from the dual exports: psi_g = minus the sum of the dyn rows of kkt_Pg_terms.csv,
and psi_d with `prices()` of prices.py.
"""

from decimal import Decimal as Dec
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from prices import prices
from printed_values import q, r1, tex_money, tex_mw

ROOT = Path(__file__).resolve().parent
RUN = ROOT.parent / "results" / "39bus" / "Results - 2026-10-07 124340"
PSSE = ROOT / "psse" / "aligned.csv"
SG, GFM = [1, 2, 7, 8, 9, 10], [3, 4, 5, 6]
HL = {6, 7}
TABS = ROOT / "tables"
TABS.mkdir(exist_ok=True)
TAB_A, TAB_B = TABS / "tab_prices_39bus_gfm.tex", TABS / "tab_settlement_39bus_gfm.tex"
FIGDIR = ROOT / "figures"

rd = lambda run, *p: pd.read_csv(run.joinpath(*p), sep=";")
r2 = lambda v: Dec(repr(float(v))).quantize(Dec("0.01"), rounding="ROUND_HALF_UP") + Dec("0")
k2 = lambda v: q(Dec(repr(float(v))) / 1000)
money = lambda x, **kw: tex_money(x, **kw).removesuffix(".00")   # two decimals unless they are zero

# ---------------------------------------------------------------- prices from the dual exports
gen, load, _, binding = prices(RUN)
terms = rd(RUN, "Dispatch", "CSV_duals", "kkt_Pg_terms.csv")
psi_g = -terms[terms.store == "dyn"].groupby("gen_id").contribution.sum()
assert (gen["psi"] - psi_g.reindex(gen.index)).abs().max() < 1e-6      # stationarity closes
gen["psi"] = psi_g.reindex(gen.index)
assert list(binding) == [("Transient_Stability", "Gen_7", "upper")] and len(binding[list(binding)[0]]) == 1
load = load[load["P"] > 0]

# ACOPF: the warm start of this run
lam_ac = rd(RUN, "Dispatch_WarmStart", "CSV_duals", "dual_P_balance.csv").set_index("Bus_ID").iloc[:, 0]
gen_ac = rd(RUN, "Dispatch_WarmStart", "CSV", "generators_report.csv").set_index("ID")

A = {}
for gid, r in gen.iterrows():
    lam, psi = r2(r.lam), r2(r.psi)
    A[gid] = dict(c=r2(r.c), P_ac=r1(gen_ac.P_MW[gid]) + Dec("0"), lam_ac=r2(lam_ac[int(r.bus)]),
                  P=r1(r.P) + Dec("0"), lam=lam, psi=psi, pi=q(lam + psi))
    assert abs(float(A[gid]["pi"]) - r.pi) < 0.011, gid
    A[gid]["profit"] = q((A[gid]["pi"] - A[gid]["c"]) * A[gid]["P"] / 1000)     # k$/h, from the printed factors


def row(gid):
    d = A[gid]
    cells = [f"G{gid}", tex_mw(d["c"]), tex_mw(d["P_ac"]), money(d["lam_ac"]), tex_mw(d["P"]),
             money(d["lam"]), money(d["psi"], signed=True), money(d["pi"])]
    return " & ".join(rf"\rowhl{{{c}}}" if gid in HL else c for c in cells) + r" \\"


def group(label, ids):
    return [rf"\multirow{{{len(ids)}}}{{*}}{{{label}}} & " + row(ids[0])] + [" & " + row(g) for g in ids[1:]]


def write(path, lines):
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("written:", path)


CAP = (r"\captionsetup{justification=centering, textfont={sc,footnotesize}, labelfont=footnotesize, "
       r"labelsep=newline}")
U = lambda x: r"{\scriptsize " + x + "}"
write(TAB_A, [
    r"\begin{table}[!tb]", CAP, r"\centering",
    r"\label{tab:prices_39bus_gfm}", r"\footnotesize", r"\setlength{\tabcolsep}{2pt}",
    r"\newcommand{\rowhl}[1]{\cellcolor{gray!25}#1}",
    r"\begin{tabular}{llr>{\hspace{5pt}}rr>{\hspace{5pt}}rrrr}", r"\toprule",
    r" & & & \multicolumn{2}{c}{\acs{ACOPF}} & \multicolumn{4}{c}{\acs{TSCOPF}} \\",
    r"\cmidrule(lr){4-5}\cmidrule(lr){6-9}",
    r" & Unit & $\mathrm{c}_g$ & $P$ & $\lambda_k$ & $P$ & $\lambda_k$ & $\psi_{g(k)}$ & $\pi_{g(k)}$ \\",
    " & & " + " & ".join(U(x) for x in [r"[\$/MWh]", "[MW]", r"[\$/MWh]", "[MW]"] + [r"[\$/MWh]"] * 3) + r" \\",
    r"\midrule", *group(r"\acs{SG}", SG), r"\midrule", *group(r"\acs{IBG}", GFM),
    r"\bottomrule", r"\end{tabular}", r"\end{table}",
])

# Same table with the profit at pi_g and without the SG/IBG column (the rule separates the two groups)
TAB_C = TABS / "tab_prices_profit_39bus_gfm.tex"
rowp = lambda gid: row(gid).removesuffix(r" \\") + " & " + (
    rf"\rowhl{{{money(A[gid]['profit'])}}}" if gid in HL else money(A[gid]["profit"])) + r" \\"
TAB_C.write_text("\n".join([
    r"\begin{table}[!tb]", CAP, r"\centering",
    r"\label{tab:prices_39bus_gfm}", r"\footnotesize", r"\setlength{\tabcolsep}{2pt}",
    r"\newcommand{\rowhl}[1]{\cellcolor{gray!25}#1}",
    r"\begin{tabular}{lr>{\hspace{5pt}}rr>{\hspace{5pt}}rrrrr}", r"\toprule",
    r" & & \multicolumn{2}{c}{\acs{ACOPF}} & \multicolumn{5}{c}{\acs{TSCOPF}} \\",
    r"\cmidrule(lr){3-4}\cmidrule(lr){5-9}",
    r"Unit & $\mathrm{c}_g$ & $P$ & $\lambda_k$ & $P$ & $\lambda_k$ & $\psi_{g(k)}$ & $\pi_{g(k)}$ & Profit \\",
    " & " + " & ".join(U(x) for x in [r"[\$/MWh]", "[MW]", r"[\$/MWh]", "[MW]"] + [r"[\$/MWh]"] * 3
                       + [r"[k\$/h]"]) + r" \\",
    r"\midrule", *map(rowp, SG), r"\midrule", *map(rowp, GFM),
    r"\bottomrule", r"\end{tabular}", r"\end{table}",
]) + "\n", encoding="utf-8")
print("written:", TAB_C)

# ---------------------------------------------------------------- settlement in k$/h
pay_tsc = (load["pi"] * load["P"]).sum()
B = {
    "ACOPF": dict(cost=sum(d["c"] * d["P_ac"] for d in A.values()),
                  gen=sum(q(d["lam_ac"] * d["P_ac"]) for d in A.values()),
                  load=(lam_ac.reindex(load.index) * load["P"]).sum()),
    "TSCOPF": dict(cost=sum(d["c"] * d["P"] for d in A.values()),
                   gen=sum(q(d["pi"] * d["P"]) for d in A.values()),
                   load=pay_tsc),
}
for d in B.values():
    d.update(cost=q(d["cost"] / 1000), gen=q(d["gen"] / 1000), load=k2(d["load"]))
    d["R"] = q(d["load"] - d["gen"])
rev_tsc = (gen["pi"] * gen["P"]).sum()
assert abs(B["TSCOPF"]["gen"] - k2(rev_tsc)) <= Dec("0.01")
assert abs(B["TSCOPF"]["cost"] - k2((gen["c"] * gen["P"]).sum())) <= Dec("0.01")
assert abs(B["TSCOPF"]["R"] - k2(pay_tsc - rev_tsc)) <= Dec("0.01")
# the printed profits add up to printed revenues minus printed cost
assert sum(d["profit"] for d in A.values()) == B["TSCOPF"]["gen"] - B["TSCOPF"]["cost"]

write(TAB_B, [
    r"\begin{table}[!tb]", CAP, r"\centering",
    r"\label{tab:settlement_39bus_gfm}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
    r"\begin{tabular}{lrrrr}", r"\toprule",
    r" & Cost & Load payments & Generator revenues & $R$ \\",
    r" & [k\$/h] & [k\$/h] & [k\$/h] & [k\$/h] \\", r"\midrule",
    *[rf"\acs{{{n}}} & {tex_money(d['cost'])} & {tex_money(d['load'])} & {tex_money(d['gen'])} & "
      rf"{tex_money(d['R'])} \\" for n, d in B.items()],
    r"\bottomrule", r"\end{tabular}", r"\end{table}",
])

# ---------------------------------------------------------------- figure: angle to the COI, TSCOPF and PSS/E
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
al = pd.read_csv(PSSE)
ang = rd(RUN, "Transient_Stability", "CSV", "angle_rel_COI.csv")
assert len(al) == len(ang) and np.allclose(al.t, ang.t)
for g in SG + GFM:                                           # the replay was aligned to this run
    assert np.allclose(al[f"G{g}_angle_rel_coi_deg_tsc"], ang[f"G{g}"], atol=1e-6), g

# Whole 5 s window, both limits, TSCOPF as lines and PSS/E as dots.
RED, TEAL, PURPLE = "#C51212", "#378d94", "#9671bd"
colour = {**{g: PURPLE for g in SG}, **{g: TEAL for g in GFM}, 7: RED}
fig, ax = plt.subplots(figsize=(10.0, 6.0), layout="constrained")
ax.axvspan(0.01, 0.01 + 0.100, color="0.85", zorder=0)        # fault on for 100 ms
for lim in (-100, 100):
    ax.axhline(lim, color="0.35", linewidth=1.2, linestyle=(0, (3, 2)), zorder=1)
for g in SG + GFM:
    z = 4 if g == 7 else 3
    ax.plot(al.t, al[f"G{g}_angle_rel_coi_deg_tsc"], color=colour[g], linewidth=2.0, zorder=z)
    ax.plot(al.t, al[f"G{g}_angle_rel_coi_deg_psse"], color=colour[g], linestyle="none",
            marker="o", markersize=5, markevery=(2 * g, 25), zorder=z)
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
line = lambda c, **kw: plt.Line2D([], [], color=c, **kw)
fig.legend(
    [line(RED, linewidth=2), line(PURPLE, linewidth=2), line(TEAL, linewidth=2),
     line("0.15", linewidth=2), line("0.15", linestyle="none", marker="o", markersize=7)],
    ["G7", "other SGs", "IBGs", "TSCOPF", "PSS/E"],
    frameon=False, fontsize=22, ncol=5, loc="outside upper center",
    handlelength=1.3, columnspacing=1.0, handletextpad=0.5)

FIGDIR.mkdir(exist_ok=True)
for ext in ("svg", "pdf"):
    fig.savefig(FIGDIR / f"angles_psse_39bus_gfm.{ext}", bbox_inches="tight")
print("written:", FIGDIR / "angles_psse_39bus_gfm.pdf")

# ---------------------------------------------------------------- numbers and checks for the prose
print("\nTable A (printed):")
print(pd.DataFrame(A).T.loc[SG + GFM].to_string())
print("\nTable B (printed, k$/h):")
print(pd.DataFrame(B).T.to_string())
cost_ac, cost = B["ACOPF"]["cost"], B["TSCOPF"]["cost"]
print(f"\ncost increase {float((cost - cost_ac) / cost_ac * 100):.1f} %")
print("G7 at the limit at t =", binding[list(binding)[0]], "s; peak", ang["G7"].max())
print("IBG peak |angle| [deg]:", {g: round(ang[f"G{g}"].abs().max(), 2) for g in GFM})
print("largest |psi_g| besides G7:", float(max(abs(d["psi"]) for g, d in A.items() if g != 7)))
mu = gen[["mu_up", "mu_lo"]].abs().max(axis=1)
print("units within their output limits:", [f"G{g}" for g in gen.index if mu[g] < 1e-6])
print(f"psi_d: min {r2(load.psi.min())} (bus {load.psi.idxmin()}), max {r2(load.psi.max())} (bus {load.psi.idxmax()})")
stab_load = (load["psi"] * load["P"]).sum()
stab_gen = (gen["psi"] * gen["P"]).sum()
print(f"sum psi_d P_d {stab_load:.2f} $/h = {stab_load / pay_tsc * 100:.3f} % of the load payments; "
      f"sum psi_g P_g {stab_gen:.2f} $/h; R_lambda {pay_tsc - rev_tsc - stab_load + stab_gen:.2f}")
for tag in ("Dispatch_WarmStart", "Dispatch"):
    c = rd(RUN, tag, "CSV", "circuits_report.csv")
    lam = rd(RUN, tag, "CSV_duals", "dual_P_balance.csv").set_index("Bus_ID").iloc[:, 0]
    c["rent"] = -(lam.reindex(c.FROM_BUS).values * c.Pik_MW + lam.reindex(c.TO_BUS).values * c.Pki_MW)
    r = c[(c.FROM_BUS == 21) & (c.TO_BUS == 22)].iloc[0]
    print(f"{tag}: line 21-22 loading {r.Loading:.3f}, rent {r.rent:.0f} $/h; at limit:",
          [f"{int(a)}-{int(b)}" for a, b in c[c.Loading > 0.985][["FROM_BUS", "TO_BUS"]].values],
          f"; total branch rent {c.rent.sum():.0f}")
d = (al.filter(like="_angle_rel_coi_deg_psse").to_numpy() - al.filter(like="_angle_rel_coi_deg_tsc").to_numpy())
print("PSS/E against TSCOPF, largest angle difference [deg]:", round(float(np.abs(d).max()), 2), )
