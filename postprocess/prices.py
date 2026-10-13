"""Shadow prices of a TSCOPF run (constant-impedance loads), from its raw exports.

    from prices import prices

lambda_k : Dispatch/CSV_duals/dual_P_balance.csv
pi_g     : c_g + mu_upper - mu_lower (stationarity);  psi_g = pi_g - lambda_k
psi_d    : sum over contingency folders and time steps of  lambda_hat^t_k * (V^t_k / V_k)^2
Contingency folders: Transient_Stability/ and Transient_Stability_c<id>/.
"""
import tomllib
from pathlib import Path

import pandas as pd


def rd(*parts):
    return pd.read_csv(Path(*parts), sep=";")


def prices(run):
    run = Path(run)
    buses = rd(run, "Dispatch", "CSV", "buses_report.csv").set_index("BUS")
    lam = rd(run, "Dispatch", "CSV_duals", "dual_P_balance.csv").set_index("Bus_ID").iloc[:, 0]

    blocks = sorted(p for p in run.glob("Transient_Stability*") if (p / "CSV_duals").is_dir())
    psi_d = pd.DataFrame(index=lam.index)
    binding = {}
    for b in blocks:
        dual = rd(b, "CSV_duals", "dual_Pbalance.csv")
        volt = rd(b, "CSV", "bus_voltage_magnitude.csv").drop(columns="t")
        assert len(dual) == len(volt), b
        # Coefficient and sign as in kkt/pd_coeffs.py (_fullbus_poly), ZIP (1,0,0):
        # v^2 on load buses, v^2 / V^2 on zero-load buses (their balance row is written unscaled).
        scale = (buses["V_pu"] ** 2).where(buses["PD_MW"].abs() < 1e-12, 1.0).to_numpy()
        a = volt.to_numpy() ** 2 / scale
        psi_d[b.name] = -(dual.to_numpy() * a).sum(axis=0)
        t = rd(b, "CSV", "bus_voltage_magnitude.csv")["t"]
        for side in ("upper", "lower"):
            d = rd(b, "CSV_duals", f"dual_delta_COI_{side}.csv").filter(like="Gen_")
            for g in d.columns:
                steps = d.index[d[g].abs() > 1e-4]
                if len(steps):
                    binding[(b.name, g, side)] = [float(t[s]) for s in steps]

    man = tomllib.loads((run / "run_manifest.toml").read_text(encoding="utf-8"))["resolved"]
    gens = rd(run, "Dispatch", "CSV", "generators_report.csv")
    mu_up = rd(run, "Dispatch", "CSV_duals", "dual_Pg_upper.csv").iloc[:, 1].abs()
    mu_lo = rd(run, "Dispatch", "CSV_duals", "dual_Pg_lower.csv").iloc[:, 1].abs()
    gen = pd.DataFrame({
        "bus": man["gen_bus"], "P": gens["P_MW"].to_numpy(), "c": man["c1"],
        "mu_up": mu_up.to_numpy(), "mu_lo": mu_lo.to_numpy(),
    }, index=man["gen_ids"])
    gen["lam"] = lam.loc[gen["bus"]].to_numpy()
    gen["pi"] = gen["c"] + gen["mu_up"] - gen["mu_lo"]
    gen["psi"] = gen["pi"] - gen["lam"]

    load = pd.DataFrame({"P": buses["PD_MW"], "lam": lam, "psi": psi_d.sum(axis=1)})
    load["pi"] = load["lam"] + load["psi"]
    return gen, load, psi_d, binding

