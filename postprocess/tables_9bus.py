"""Printed values of the 9-bus runs at a 5 ms step (100 ms fault), Section IV-A.

Inputs are the raw exports of the two runs:
  lambda_k, P, c, mu : prices() of prices.py
  psi_g              : minus the store = dyn rows of Dispatch/CSV_duals/kkt_Pg_terms.csv (never sigma + chi alone)
  psi_d per load     : -sum_t lambda_hat^t_k a^t, a = aZ v_t^2 + aI v_t V   (ZIP data of INPUT_FILES/<case>/load_data.csv)

Money is in k$/h: every cell is the rounded product of the printed P and pi, and the cost, the totals, R and
the profits are sums and differences of printed cells.

Writes tables/ and tables/prose_numbers_9bus.txt here. The table header and footer come from templates/.

    python tables_9bus.py
"""
from decimal import ROUND_HALF_UP, Decimal as Dec
from pathlib import Path

import pandas as pd

from prices import prices, rd
from printed_values import q, tex_mw

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
RUNS = {"linear": REPO / "results" / "9bus" / "Results - 2026-10-08 112948",
        "coloc": REPO / "results" / "9bus" / "Results - 2026-10-08 113323"}
CFG = {
    "linear": dict(case="9bus_linear_H24_6_3", step=Dec("0.01"), label={1: "G1", 2: "G2", 3: "G3"},
                   table="tab_prices_9bus.tex", at_limit=set()),
    "coloc": dict(case="9bus_colocated_bus1_H24_6_3", step=Dec("0.1"), label={1: "G1", 2: "G2", 3: "G3", 4: "G1b"},
                  table="tab_colocated_bus1_9bus.tex", at_limit={4}),
}
ACOPF = {"linear": Dec("8775"), "coloc": Dec("7775")}       # $/h, unchanged benchmarks
K = lambda x: q(x / 1000)                                   # k$/h, two decimals


def raw(run, case):
    """Full-precision generators and loads of a run, from the raw exports."""
    gen, bus_load, _, binding = prices(run)
    terms = rd(run, "Dispatch", "CSV_duals", "kkt_Pg_terms.csv")
    gen["psi"] = -terms[terms["store"] == "dyn"].groupby("gen_id")["contribution"].sum()
    assert (gen["psi"] - (gen["c"] + gen["mu_up"] - gen["mu_lo"] - gen["lam"])).abs().max() < 1e-6

    buses = rd(run, "Dispatch", "CSV", "buses_report.csv").set_index("BUS")
    volt = rd(run, "Transient_Stability", "CSV", "bus_voltage_magnitude.csv")
    dual = rd(run, "Transient_Stability", "CSV_duals", "dual_Pbalance.csv")
    f = REPO / "INPUT_FILES" / case / "load_data.csv"
    if f.exists():
        loads = pd.read_csv(f, sep=";").set_index("id")
    else:                                                   # no per-load file: one constant-impedance load per bus
        pd_bus = buses["PD_MW"][buses["PD_MW"] > 0]
        loads = pd.DataFrame({"bus": pd_bus.index, "p_d": pd_bus.to_numpy(), "aZ_p": 1, "aI_p": 0, "aP_p": 0},
                             index=range(1, len(pd_bus) + 1))
    assert (loads["aP_p"] == 0).all()                       # constant-power loads are prohibited
    assert (loads.groupby("bus")["p_d"].sum() - buses["PD_MW"][buses["PD_MW"] > 0]).abs().max() < 1e-9
    psi = {}
    for l, r in loads.iterrows():
        j = buses.index.get_loc(r["bus"])
        v, V = volt.drop(columns="t").iloc[:, j], buses.at[r["bus"], "V_pu"]
        psi[l] = -(dual.iloc[:, j] * (r["aZ_p"] * v ** 2 + r["aI_p"] * v * V)).sum()
    loads = loads.assign(psi=pd.Series(psi), lam=bus_load["lam"].reindex(loads["bus"]).to_numpy(),
                         kind=["Z" if z == 1 else "I" for z in loads["aZ_p"]])
    for b, grp in loads.groupby("bus"):                     # one Z load at a bus: the bus value of prices()
        if len(grp) == 1:
            assert abs(grp["psi"].iloc[0] - bus_load.at[b, "psi"]) < 1e-9, b

    # statements of the text that do not depend on the printed digits
    assert set(binding) == {("Transient_Stability", "Gen_3", "upper")}, binding
    t_hit = binding["Transient_Stability", "Gen_3", "upper"]
    win = (volt["t"] > 0.01 + 1e-9) & (volt["t"] <= t_hit[-1] + 1e-9)
    j5 = buses.index.get_loc(5)
    assert (volt.drop(columns="t").iloc[:, j5][win] < buses.at[5, "V_pu"]).all(), "V5 above its pre-fault value"
    return gen, loads, t_hit


def printed(key, run):
    cfg = CFG[key]
    p2 = lambda v: Dec(repr(float(v))).quantize(cfg["step"], ROUND_HALF_UP)
    gen, loads, t_hit = raw(run, cfg["case"])
    at_limit = set(gen.index[(gen["mu_up"] > 1e-6) | (gen["mu_lo"] > 1e-6)])
    assert at_limit == cfg["at_limit"], at_limit

    Pd = {l: Dec(repr(float(loads.at[l, "p_d"]))).normalize() for l in loads.index}
    Pg = {g: Dec(repr(float(gen.at[g, "P"]))).quantize(Dec("0.1"), ROUND_HALF_UP) for g in gen.index}
    assert sum(Pg.values()) == sum(Pd.values()), (Pg, "lossless: sum P_g = sum P_d on the page")

    G, L = {}, {}
    for g in gen.index:
        lam, c = p2(gen.at[g, "lam"]), Dec(str(gen.at[g, "c"]))
        psi = p2(gen.at[g, "psi"]) if g in at_limit else c - lam      # within its limits: pi = c_g on the page
        assert abs(float(psi) - gen.at[g, "psi"]) <= float(cfg["step"]) + 1e-7, (g, psi)
        G[g] = dict(name=cfg["label"][g], bus=int(gen.at[g, "bus"]), P=Pg[g], lam=lam, psi=psi, pi=lam + psi, c=c,
                    pay=K(Pg[g] * (lam + psi)), cost=K(Pg[g] * c), loss_lam=max(Dec(0), K(Pg[g] * (c - lam))))
    for l in loads.index:
        lam, psi = p2(loads.at[l, "lam"]), p2(loads.at[l, "psi"])
        L[l] = dict(name=f"{loads.at[l, 'kind']}$_{loads.at[l, 'bus']}$", bus=int(loads.at[l, "bus"]),
                    kind=loads.at[l, "kind"], P=Pd[l], lam=lam, psi=psi, pi=lam + psi, pay=K(Pd[l] * (lam + psi)))
    for d, full in ((G, gen), (L, loads)):                  # no sign flips at the printed resolution
        assert all((v["psi"] > 0) == (full.at[i, "psi"] > 0) and v["psi"] != 0 for i, v in d.items())
    pay, rev = sum(v["pay"] for v in L.values()), sum(v["pay"] for v in G.values())
    cost = sum(v["cost"] for v in G.values())
    return G, L, dict(pay=pay, rev=rev, R=pay - rev, cost=cost, t_hit=t_hit,
                      pct=((cost / K(ACOPF[key]) - 1) * 100).quantize(Dec("0.1"), ROUND_HALF_UP))


def rows(key, G, L):
    """Data rows in the layout of the paper table: one list of cell strings per asset."""
    coloc = key == "coloc"
    sgn = lambda x: ("$+$" if x > 0 else "$-$") + f"{abs(x):.2f}"
    pi = lambda x: f"{x:.0f}" if x == x.to_integral_value() else f"{x:.2f}"
    out = []
    for kind, d in (("Gen", G), ("Load", L)):
        order = sorted(d, key=lambda i: (d[i]["bus"], d[i]["name"]) if kind == "Gen" else (d[i]["bus"], d[i].get("kind") != "Z"))
        for n, i in enumerate(order):
            v = d[i]
            cells = [v["name"] if coloc else str(v["bus"]), tex_mw(v["P"]), f"{v['lam']:.2f}", sgn(v["psi"]), pi(v["pi"]),
                     f"{v['pay']:.2f}"]
            if coloc:
                cells.append(f"{v['pay'] - v['cost']:.2f}" if kind == "Gen" else "---")
            # shaded: the co-located pairs (IV-A.3); the assets that improve the stability margin (IV-A.1)
            shade = sum(w["bus"] == v["bus"] for w in d.values()) > 1 if coloc else (v["psi"] > 0) == (kind == "Gen")
            cells = [rf"\rowhl{{{c}}}" for c in cells] if shade else cells
            out.append([rf"\multirow{{{len(d)}}}{{*}}{{{kind}}}" if n == 0 else ""] + cells)
    return out


def write_table(key, body):
    old = (ROOT / "templates" / CFG[key]["table"]).read_text(encoding="utf-8")
    head, foot = old.split(r"\midrule", 1)[0], old[old.index(r"\bottomrule"):]
    n = len(body) // 2
    lines = [" & ".join(r) + r" \\" for r in body]
    text = head + "\\midrule\n" + "\n".join(lines[:n]) + "\n\\midrule\n" + "\n".join(lines[n:]) + "\n" + foot
    (ROOT / "tables" / CFG[key]["table"]).write_text(text, encoding="utf-8", newline="\n")


def notes(key, run, G, L, S):
    out = [f"{'IV-A.1' if key == 'linear' else 'IV-A.3'} ({run.name})",
           f"  cost {S['cost']} k$/h (ACOPF {K(ACOPF[key])}; +{S['pct']} %)",
           f"  G3 at the upper angle limit at t = {', '.join(f'{t:g}' for t in S['t_hit'])} s",
           *[f"  {v['name']}: P {tex_mw(v['P'])}  lam {v['lam']}  psi {v['psi']:+}  pi {v['pi']}  pi-c {v['pi'] - v['c']}  "
             f"receives {v['pay']}  profit {v['pay'] - v['cost']}  loss if paid lam alone {v['loss_lam']}" for v in G.values()],
           *[f"  load {v['kind']}{v['bus']}: P {tex_mw(v['P'])}  lam {v['lam']}  psi {v['psi']:+}  pi {v['pi']}  pays {v['pay']}"
             for v in L.values()],
           f"  lambda_k at generator and load buses: {min(v['lam'] for v in (G | L).values())} to "
           f"{max(v['lam'] for v in G.values())};  lambda_3 / c_3 = {max(v['lam'] for v in G.values()) / 10:.2f}",
           f"  loads pay {S['pay']}  generators receive {S['rev']}  R {S['R']}  [k$/h]"]
    b5 = {v["kind"]: v["psi"] for v in L.values() if v["bus"] == 5}
    if len(b5) == 2:
        out.append(f"  bus 5: constant-current load pays {b5['I'] - b5['Z']:.2f} $/MWh more than constant-impedance")
    return out


def main():
    (ROOT / "tables").mkdir(exist_ok=True)
    text = []
    for key, run in RUNS.items():
        G, L, S = printed(key, run)
        write_table(key, rows(key, G, L))
        text += notes(key, run, G, L, S) + [""]
    (ROOT / "tables" / "prose_numbers_9bus.txt").write_text("\n".join(text), encoding="utf-8")
    print("\n".join(text))


if __name__ == "__main__":
    main()
