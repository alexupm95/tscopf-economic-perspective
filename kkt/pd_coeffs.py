"""Coefficients of the load parameter: ``∂expr_r / ∂P_{d,k}`` for every model row that
contains demand.

Why this file exists
--------------------
``P_g`` is a JuMP variable, so ``F_P`` is a free slice of the exported Jacobian and the
Tier-2 notebook computes no physics at all. ``P_d`` is a ``Float64`` inlined into the
builders — MOI has no column for it (see ``param_jacobian.jl``) — so its column has to be
*computed*. Two consumers need the same coefficients:

* **Tier 1-D** (``KKT_Pd_Duals.ipynb``) contracts them with duals into the nodal price,
* **Tier 2-D** (``S_Sensitivity_Pd.ipynb``) keeps them as the matrix ``F_Pd`` and feeds it
  to ``lu.solve``.

Keeping one implementation is the point: a second copy is how the family lists in this
project went stale before.

Sign convention
---------------
Everything here returns ``∂expr/∂P_d`` with ``expr`` exactly as the model stores it
(``expr {==,<=} 0``). Both consumers then apply the same MOI rule as the generator side,

    contribution = −dual × ∂expr/∂P_d

The two paths carry opposite signs, and that is physical, not a stray negation: the
transient nodal balance is oriented the other way round from the Kron ``Pe`` definition.
Verified against a model-rebuild finite difference at 1.4e-14 (FULL_BUS) and 5.1e-7
(Kron) when Tier 1-D was built.

Traps preserved from the original implementation
------------------------------------------------
* ``angle_abs.csv`` is in **degrees**; reading it as radians turns a 26° swing into 1487°.
* ``angle_abs.csv`` carries both a ``t`` column and a ``COI`` column — ``COI`` is an
  inertia-weighted aggregate, not a machine, and must be dropped along with ``t``.
* ``generator_internal_voltage.csv`` is **long** format (``generator;E_pu``). An earlier
  prototype refitted ``E`` by least squares because it did not know this is exported.
* At a bus whose base load is zero the FULL_BUS builder writes the row *without* its
  ``V_k²`` factor, so its dual carries the opposite scaling and the coefficient is
  ``Z/V_k²`` rather than ``Z``.
* ``zip_load_p`` is ``(Z, I, P)`` order. The default ``(1, 0, 0)`` is constant impedance.
* A run with a ``load_data.csv`` splits a bus's demand into loads with different ZIP
  shares. ``zip_p`` then has to be the per-bus mean (``zip_by_bus``), not the config
  tuple, or the bus coefficient silently describes a load model the run did not use.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["coeff_arrays", "f_pd_matrix", "PD_FAMILIES",
           "read_load_table", "zip_by_bus", "load_coeff_arrays"]

# Which row families carry P_d on each path. Used to slice `jacobian_rows.csv`; the
# builders' own family names, not a private vocabulary.
PD_FAMILIES = {
    "FULL_BUS": ("eq_const_Pbalance_tf", "eq_const_Pbalance_tpf"),
    "KRON_REDUCED": ("eq_const_Pe_tf", "eq_const_Pe_tpf"),
}


def _rd(run_dir, *parts):
    return pd.read_csv(run_dir.joinpath(*parts), sep=";")


def _schur_factors(Y: np.ndarray, n_bus: int):
    """Rank-one factors of ``dY_red/dy_k``.

    Folding a bus shunt into the reduced matrix goes through the Schur complement, so one
    bus's load reaches *every* generator's ``Pe`` row:

        dY_red/dy_k = u_k v_kᵀ ,   u_k = (Y_IB Y_BB⁻¹)[:, k] ,   v_kᵀ = (Y_BB⁻¹ Y_BI)[k, :]
    """
    Ybb, Ybi, Yib = Y[:n_bus, :n_bus], Y[:n_bus, n_bus:], Y[n_bus:, :n_bus]
    W = np.linalg.solve(Ybb, Ybi)           # Y_BB⁻¹ Y_BI
    U = np.linalg.solve(Ybb.T, Yib.T).T     # Y_IB Y_BB⁻¹
    return U, W


def _load_Y(run_dir, name: str) -> np.ndarray:
    g = pd.read_excel(run_dir / "Bus_Matrices" / name, sheet_name="Real_G", header=0)
    b = pd.read_excel(run_dir / "Bus_Matrices" / name, sheet_name="Imaginary_B", header=0)
    return g.to_numpy(float) + 1j * b.to_numpy(float)


def _fullbus_poly(v, V, zip3, *, zero_load: bool):
    """``∂expr/∂P_d`` for one FULL_BUS nodal-balance row, given a (Z, I, P) triple.

    The one place the ZIP polynomial is written down. Per-bus (`coeff_arrays`) and
    per-load (`load_coeff_arrays`) both come through here, so the two cannot drift.
    """
    aZ, aI, aP = zip3
    Z = aZ * v**2 + aI * v * V + aP * V**2
    # Zero-load buses: the builder wrote the row unscaled, so its dual carries
    # 1/V² and the coefficient must be divided by V² to match.
    return Z / (V**2 if zero_load else 1.0)


def _zip_at(zip_p, bus: int):
    """Accept either one system-wide triple or a ``{bus: triple}`` mapping."""
    if isinstance(zip_p, dict):
        return zip_p[bus]
    return zip_p


def read_load_table(run_dir, *, load_factor: float = 1.0):
    """The run's archived ``Inputs/load_data.csv``, or ``None`` when the case had none.

    Returns the table with `bus` translated to the model's internal 1..n_bus numbering
    and `p_d` / `q_d` scaled by `load_factor` — the archived copy is raw input, exactly
    like `Inputs/bus_data.csv`, and the internal id of a bus is its row position in that
    same file (`Mapping_Buses_Labels`).
    """
    path = run_dir / "Inputs" / "load_data.csv"
    if not path.is_file():
        return None
    loads = pd.read_csv(path, sep=";")
    bus_csv = run_dir / "Inputs" / "bus_data.csv"
    if not bus_csv.is_file():
        # MATPOWER-mode runs archive the .m instead, and the label -> index map lives
        # nowhere else in the folder. Say so rather than guessing an identity map.
        raise FileNotFoundError(
            f"{bus_csv} is missing, so load_data.csv bus labels cannot be mapped to the "
            "model's internal numbering. This run was read from a MATPOWER .m file.")
    buses = pd.read_csv(bus_csv, sep=";")
    label_to_index = {int(b): i + 1 for i, b in enumerate(buses.iloc[:, 0])}
    loads["bus"] = [label_to_index[int(b)] for b in loads["bus"]]
    loads["p_d"] = loads["p_d"].astype(float) * load_factor
    loads["q_d"] = loads["q_d"].astype(float) * load_factor
    return loads


def zip_by_bus(loads, zip_p, *, channel: str = "p"):
    """Per-bus (Z, I, P) triples: the demand-weighted mean of the loads at each bus.

    Mirrors `zip_by_bus` in `functions_4_TS_fullbus_helpers.jl` — the model builds its
    balance rows from exactly these numbers, so reading a run means rebuilding them.
    Returns `zip_p` unchanged when the run had no load table.
    """
    if loads is None or len(loads) == 0:
        return zip_p
    w = "p_d" if channel == "p" else "q_d"
    cols = [f"aZ_{channel}", f"aI_{channel}", f"aP_{channel}"]
    out = {}
    for bus, grp in loads.groupby("bus"):
        total = float(grp[w].sum())
        if abs(total) < 1e-12:
            continue                      # falls back to the system-wide split below
        shares = (grp[cols].to_numpy(float) * grp[w].to_numpy(float)[:, None]).sum(0)
        out[int(bus)] = tuple(shares / total)
    return _DefaultZip(out, zip_p)


class _DefaultZip(dict):
    """`{bus: triple}` that answers the system-wide split for any bus it does not hold."""

    def __init__(self, mapping, default):
        super().__init__(mapping)
        self._default = tuple(default)

    def __missing__(self, key):
        return self._default


def load_coeff_arrays(run_dir, *, loads, n_tf: int, n_tpf: int, channel: str = "p"):
    """``∂expr/∂p_{d,l}`` per individual load — the per-load twin of `coeff_arrays`.

    FULL_BUS only. A load enters exactly one row, its own bus's nodal balance, with its
    *own* ZIP polynomial; the bus dual is shared, the coefficient is not. That is what
    lets two loads on the same bus carry different stability prices while the LMP they
    pay for energy is identical.

    Returns ``{load_id: array of shape (n_steps,)}``, aligned with `coeff_arrays` rows.
    The identity that ties the two together, per bus:

        Σ_l p_d,l · coeff_l  ==  p_d,bus · coeff_bus

    exact by construction, because the per-bus row is the demand-weighted mean.
    """
    if loads is None or len(loads) == 0:
        return {}
    buses = _rd(run_dir, "Dispatch", "CSV", "buses_report.csv")
    bus_ids = buses["BUS"].astype(int).tolist()
    V_ss = dict(zip(bus_ids, buses["V_pu"].astype(float)))
    PD_MW = dict(zip(bus_ids, buses["PD_MW"].astype(float)))
    Vt = _rd(run_dir, "Transient_Stability", "CSV", "bus_voltage_magnitude.csv")
    cols = [f"aZ_{channel}", f"aI_{channel}", f"aP_{channel}"]

    out: dict[int, np.ndarray] = {}
    for row in loads.itertuples(index=False):
        b = int(row.bus)
        col = f"Bus_{b}"
        if col not in Vt.columns or b not in V_ss:
            continue
        zip3 = tuple(float(getattr(row, c)) for c in cols)
        out[int(row.id)] = _fullbus_poly(Vt[col].to_numpy(float), V_ss[b], zip3,
                                         zero_load=abs(PD_MW[b]) < 1e-12)
    return out


def coeff_arrays(run_dir, *, network: str, zip_p, n_tf: int, n_tpf: int):
    """``∂expr/∂P_{d,k}`` for every dynamic row, per perturbed bus.

    The manifest-derived values are passed in rather than re-read here, so a caller's
    existing parse stays authoritative and this refactor cannot move a number.

    Returns ``(bus_ids, coeff)`` where ``coeff[bus]`` is

    * **FULL_BUS** — shape ``(n_steps,)``. Load at bus *k* appears only in bus *k*'s own
      nodal balance, so the entry lands on ``eq_const_Pbalance_*`` of that same bus.
    * **KRON_REDUCED** — shape ``(n_steps, n_gen)``. The Schur complement spreads one
      bus's shunt over the whole reduced matrix, so the entry lands on **every**
      generator's ``eq_const_Pe_*`` row.

    Row ``i`` of either array is trajectory step ``i``: ``i < n_tf`` is the fault window,
    the rest post-fault. That is the ordering the trajectory CSVs already use.
    """
    buses = _rd(run_dir, "Dispatch", "CSV", "buses_report.csv")
    bus_ids = buses["BUS"].astype(int).tolist()
    V_ss = dict(zip(bus_ids, buses["V_pu"].astype(float)))
    PD_MW = dict(zip(bus_ids, buses["PD_MW"].astype(float)))

    coeff: dict[int, np.ndarray] = {}

    if network == "FULL_BUS":
        Vt = _rd(run_dir, "Transient_Stability", "CSV", "bus_voltage_magnitude.csv")
        for b in bus_ids:
            col = f"Bus_{b}"
            if col not in Vt.columns:
                continue
            coeff[b] = _fullbus_poly(Vt[col].to_numpy(float), V_ss[b],
                                     _zip_at(zip_p, b), zero_load=abs(PD_MW[b]) < 1e-12)

    elif network == "KRON_REDUCED":
        ang = _rd(run_dir, "Transient_Stability", "CSV", "angle_abs.csv")
        gen_cols = [c for c in ang.columns if c not in ("t", "COI")]
        dlt = np.deg2rad(ang[gen_cols].to_numpy(float))          # degrees on disk

        Eg = _rd(run_dir, "Transient_Stability", "CSV",
                 "generator_internal_voltage.csv")["E_pu"].to_numpy(float)
        n_gen = len(Eg)
        assert dlt.shape[1] == n_gen, f"{dlt.shape[1]} angle columns vs {n_gen} machines"
        EE = np.outer(Eg, Eg)

        Y_F, Y_P = _load_Y(run_dir, "Ybus_fault.xlsx"), _load_Y(run_dir, "Ybus_postfault.xlsx")
        n_bus_red = Y_F.shape[0] - n_gen
        U_F, W_F = _schur_factors(Y_F, n_bus_red)
        U_P, W_P = _schur_factors(Y_P, n_bus_red)

        n_steps = len(dlt)
        for b in bus_ids:
            k = b - 1                       # bus ids are 1..n_bus and index the B partition
            if k >= n_bus_red:
                continue
            out = np.zeros((n_steps, n_gen))
            for lo, hi, U, W in ((0, n_tf, U_F, W_F),
                                 (n_tf, n_tf + n_tpf, U_P, W_P)):
                dY = np.outer(U[:, k], W[k, :])
                for t in range(lo, min(hi, n_steps)):
                    dd = dlt[t][:, None] - dlt[t][None, :]
                    M = EE * (dY.real * np.cos(dd) + dY.imag * np.sin(dd))
                    # expr_Pe(g,t) = Pe − E_g Σ_i E_i(G cos Δδ + B sin Δδ)
                    out[t] = -M.sum(axis=1)
            coeff[b] = out

        # The faulted bus keeps no load shunt during the fault window: the builder
        # overwrites that diagonal with FAULT_BUS_SHUNT, so a ~0 fault-window
        # coefficient there is expected rather than a missing term.

    else:
        raise ValueError(f"network_form {network!r} has no P_d coefficient rule here.")

    return bus_ids, coeff


def f_pd_matrix(run_dir, rows_df: pd.DataFrame, *, network: str, zip_p,
                n_tf: int, n_tpf: int):
    """Assemble ``F_Pd``: ``(n_F, n_bus)``, aligned row-for-row with the ``F`` block of
    ``jacobian_rows.csv`` and column-for-column with the returned bus ids.

    This is the piece Tier 2-D adds. Tier 1-D never needed it because it contracts in
    physics space; here the coefficients have to be placed at their Jacobian row
    positions, which is the only genuinely new step and the one that can go wrong
    silently — a missed entry shows up as a smaller charge, never as an error. Every
    intended placement is therefore asserted.
    """
    fams = PD_FAMILIES.get(network)
    if fams is None:
        raise ValueError(f"network_form {network!r} has no P_d row families here.")

    F = rows_df[rows_df.block == "F"].reset_index(drop=True)
    # (family, key1, key2) -> position within the F block, i.e. the row index of F_Pd.
    pos = {(f, int(k1), int(k2)): i
           for i, (f, k1, k2) in enumerate(zip(F.family, F.key1, F.key2))
           if f in fams and pd.notna(k1) and pd.notna(k2)}

    bus_ids, coeff = coeff_arrays(run_dir, network=network, zip_p=zip_p,
                                  n_tf=n_tf, n_tpf=n_tpf)
    fam_tf, fam_tpf = fams

    def row_of(key1: int, step: int):
        """Trajectory step (0-based, fault window first) -> F-block row position.

        `key2` restarts at 1 in each window, so the post-fault index is rebased. Getting
        this wrong is the failure this function exists to prevent.
        """
        if step < n_tf:
            return pos.get((fam_tf, key1, step + 1))
        return pos.get((fam_tpf, key1, step - n_tf + 1))

    F_Pd = np.zeros((len(F), len(bus_ids)))
    placed = missed = 0

    for c, b in enumerate(bus_ids):
        arr = coeff.get(b)
        if arr is None:
            continue
        if network == "FULL_BUS":
            for t, val in enumerate(arr):                 # own bus row only
                r = row_of(b, t)
                if r is None:
                    missed += 1
                    continue
                F_Pd[r, c] = val
                placed += 1
        else:
            for t in range(arr.shape[0]):                 # every machine's Pe row
                for g in range(arr.shape[1]):
                    r = row_of(g + 1, t)
                    if r is None:
                        missed += 1
                        continue
                    F_Pd[r, c] = arr[t, g]
                    placed += 1

    assert missed == 0, (
        f"{missed} coefficient(s) found no matching F row — the (family, key1, key2) "
        "join or the fault/post-fault step rebase is wrong, and the charge would come "
        "out too small rather than failing.")
    assert placed > 0, "no P_d coefficient was placed; check network_form and families."

    return F_Pd, bus_ids, placed
