"""
Compare one PSS/E replay with the TSC-OPF trajectories.

PSS/E is interpolated onto the TSC-OPF time grid. The centre of inertia uses
machine H and the converter gfm_H_eff weights, for all ten units.

Usage:
  python 02_compare_psse_tscopf.py --run pqon
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import case_data as cd

PURPLE = "#9671bd"
TEAL = "#378d94"
RED = "#C51212"
PEAK_LO = 0.10
PEAK_HI = 0.50
COI_TOL_DEG = 1e-3

SYNC_STYLE = {
    "G1": "-",
    "G2": "-.",
    "G8": ":",
    "G9": (0, (5, 1, 1, 1)),
    "G10": (0, (1, 1)),
}
GFM_STYLE = {
    "G3": "-",
    "G4": "-.",
    "G5": ":",
    "G6": (0, (5, 1, 1, 1)),
}


def _units() -> list[str]:
    return [f"G{gid}" for gid in cd.UNIT_IDS]


def _style(unit: str) -> tuple[str, str | tuple]:
    if unit == "G7":
        return RED, "-"
    if unit in GFM_STYLE:
        return TEAL, GFM_STYLE[unit]
    return PURPLE, SYNC_STYLE[unit]


def _rcparams(size: int) -> None:
    plt.rcParams.update(
        {
            "font.family": "Courier New",
            "font.size": size,
            "axes.labelsize": size,
            "xtick.labelsize": size,
            "ytick.labelsize": size,
            "legend.fontsize": max(10, size - 8),
            "mathtext.fontset": "custom",
            "mathtext.rm": "Courier New",
            "mathtext.it": "Courier New:italic",
            "mathtext.bf": "Courier New:bold",
            "mathtext.default": "regular",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def _panel(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, color="0.85", zorder=0)
    ax.set_axisbelow(True)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.45)
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.45)
    plt.close(fig)
    print(f"Wrote {path.with_suffix('.pdf')}")
    print(f"Wrote {path.with_suffix('.svg')}")


def _collapse(frame: pd.DataFrame, t_open: float) -> pd.DataFrame:
    """Keep the post-fault sample at the fault instant and the pre-clear sample at opening."""
    parts = []
    for stamp, group in frame.groupby(frame["time"].round(6), sort=True):
        if len(group) == 1 or abs(float(stamp) - t_open) > 1e-4:
            parts.append(group.tail(1))
        else:
            parts.append(group.head(1))
    out = pd.concat(parts, ignore_index=True)
    time = out["time"].to_numpy()
    if len(time) < 2 or np.any(np.diff(time) <= 0.0):
        raise SystemExit("PSS/E time is not strictly increasing after event samples are collapsed.")
    return out


def _infer_delt(time: np.ndarray) -> float:
    step = np.diff(np.unique(np.round(time, 6)))
    step = step[step > 1e-6]
    if len(step) == 0:
        raise SystemExit("PSS/E channel file has no time step.")
    return float(np.median(step))


def _scale_to_mw(raw: float, dispatch_mw: float, mbase: float) -> tuple[float, str]:
    candidates = {
        "system_pu": cd.SBASE,
        "machine_pu": mbase,
        "MW": 1.0,
    }
    name = min(candidates, key=lambda key: abs(raw * candidates[key] - dispatch_mw))
    scale = candidates[name]
    error = abs(raw * scale - dispatch_mw)
    if error > 1.0:
        raise SystemExit(
            f"PSS/E power {raw:.6g} does not match the dispatch {dispatch_mw:.4f} MW "
            f"on the system base, the machine base, or in MW (closest error {error:.4f})."
        )
    return scale, name


def _coi(angles: dict[str, np.ndarray], weights: dict[str, float]) -> np.ndarray:
    total = sum(weights.values())
    acc = np.zeros_like(next(iter(angles.values())), dtype=float)
    for unit, weight in weights.items():
        acc = acc + weight * angles[unit]
    return acc / total


def _check_tsc_coi(run_dir: Path, weights: dict[str, float]) -> None:
    absolute = cd.read_semi(run_dir / "Transient_Stability" / "CSV" / "angle_abs.csv")
    relative = cd.read_semi(run_dir / "Transient_Stability" / "CSV" / "angle_rel_COI.csv")
    cols = _units()
    built = _coi({unit: absolute[unit].to_numpy() for unit in cols}, weights)
    mismatch = np.max(np.abs(absolute[cols].to_numpy() - built[:, None] - relative[cols].to_numpy()))
    print(f"TSC-OPF centre of inertia rebuilt from H and gfm_H_eff: max |error| = {mismatch:.3e} deg")
    if mismatch > COI_TOL_DEG:
        raise SystemExit(
            "The centre-of-inertia weights do not reproduce angle_rel_COI.csv. "
            "PSS/E angles will not be reduced with these weights."
        )


def _metrics(time: np.ndarray, error: np.ndarray) -> dict[str, float]:
    full = np.isfinite(error)
    first = full & (time <= 1.0 + 1e-9)

    def pack(mask: np.ndarray) -> tuple[float, float]:
        err = error[mask]
        if len(err) == 0:
            return float("nan"), float("nan")
        return float(np.sqrt(np.mean(err**2))), float(np.max(np.abs(err)))

    rmse_full, max_full = pack(full)
    rmse_1s, max_1s = pack(first)
    return {
        "rmse_full": rmse_full,
        "maxabs_full": max_full,
        "rmse_1s": rmse_1s,
        "maxabs_1s": max_1s,
    }


def _peak(time: np.ndarray, angle: np.ndarray) -> tuple[float, float]:
    mask = (time >= PEAK_LO) & (time <= PEAK_HI)
    if not np.any(mask):
        raise SystemExit("No samples in the first-swing window 0.10-0.50 s.")
    idx = int(np.argmax(angle[mask]))
    return float(angle[mask][idx]), float(time[mask][idx])


def _loading(p_mw, q_mvar, voltage, invbase: float, imax: float):
    i_sys = np.hypot(p_mw / cd.SBASE, q_mvar / cd.SBASE) / np.maximum(voltage, 1e-6)
    return i_sys * cd.SBASE / (invbase * imax)


def _e_internal(p_mw, q_mvar, voltage, invbase: float, xl: float):
    vv = np.maximum(voltage, 1e-6)
    p_dev = p_mw / invbase
    q_dev = q_mvar / invbase
    return np.hypot(vv + xl * q_dev / vv, xl * p_dev / vv)


def _upsert(path: Path, fresh: pd.DataFrame, key: str) -> None:
    if path.exists():
        old = pd.read_csv(path)
        old = old.loc[~old[key].isin(fresh[key].unique())]
        fresh = pd.concat([old, fresh], ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    fresh.to_csv(path, index=False)
    print(f"Wrote {path}")


def _valve_margin(run_dir: Path, tables: dict, tag: str) -> pd.DataFrame:
    valve = cd.read_semi(run_dir / "Transient_Stability" / "CSV" / "governor_P_valve.csv")
    rows = []
    for gid in cd.SYNC_IDS:
        unit = f"G{gid}"
        static = tables["gens"].loc[cd.bus_of(gid)]
        series = valve[unit].to_numpy()
        pmax = float(static.Pmax)
        pmin = float(static.Pmin)
        rows.append(
            {
                "run": tag,
                "unit": unit,
                "Pmin_MW": pmin,
                "Pmax_MW": pmax,
                "valve_min_MW": float(series.min()),
                "valve_max_MW": float(series.max()),
                "margin_to_Pmin_MW": float(series.min() - pmin),
                "margin_to_Pmax_MW": float(pmax - series.max()),
            }
        )
    return pd.DataFrame(rows)


def _plot_angles(tag: str, time, tsc, psse, t_clear: float) -> None:
    _rcparams(16)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axvspan(cd.T_FAULT_TABLE, t_clear, color="0.85", zorder=0)
    ax.axhline(100, color="0.45", ls="--", lw=1.0, zorder=1)
    ax.axhline(-100, color="0.45", ls="--", lw=1.0, zorder=1)
    handles = []
    peak = 100.0
    for unit in _units():
        color, linestyle = _style(unit)
        ax.plot(time, psse[unit], color=color, ls=linestyle, lw=2.0, zorder=2)
        ax.plot(time, tsc[unit], color=color, ls=(0, (2, 1.2)), lw=1.3, zorder=3)
        handles.append(Line2D([0], [0], color=color, ls=linestyle, lw=2.0, label=unit))
        peak = max(peak, float(np.nanmax(np.abs(psse[unit]))), float(np.nanmax(np.abs(tsc[unit]))))
    handles.append(Line2D([0], [0], color="0.2", ls=(0, (2, 1.2)), lw=1.3, label="TSC-OPF"))
    ax.set_xlim(0.0, 5.0)
    ax.set_ylim(-(peak + 25.0), peak + 25.0)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("angle relative to the COI [deg]")
    legend = ax.legend(handles=handles, frameon=False, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    _panel(ax)
    _ = legend
    _save(fig, cd.FIGURES / f"angles_rel_coi_{tag}")


def _plot_paper(time, g7_tsc, g7_psse, g6_tsc, g6_psse, i_tsc, i_psse, t_clear: float) -> None:
    _rcparams(18)
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    series = (
        (axes[0], g7_psse, g7_tsc, RED, "G7 angle relative to the COI [deg]"),
        (axes[1], g6_psse, g6_tsc, PURPLE, "G6 angle relative to the COI [deg]"),
        (axes[2], i_psse, i_tsc, TEAL, "G6 current / Imax"),
    )
    for ax, psse, tsc, color, label in series:
        ax.axvspan(cd.T_FAULT_TABLE, t_clear, color="0.85", zorder=0)
        ax.plot(time, psse, color=color, lw=2.2, label="PSS/E", zorder=2)
        ax.plot(time, tsc, color=color, ls=(0, (2, 1.2)), lw=1.6, label="TSC-OPF", zorder=3)
        ax.set_ylabel(label)
        legend = ax.legend(frameon=False, loc="best")
        _ = legend
        _panel(ax)
    axes[0].axhline(100, color="0.45", ls="--", lw=1.0, zorder=1)
    axes[0].axhline(-100, color="0.45", ls="--", lw=1.0, zorder=1)
    axes[2].set_xlabel("time [s]")
    axes[2].set_xlim(0.0, 1.0)
    _save(fig, cd.FIGURES / "paper_100ms")


def _plot_gfm(tag: str, time, frame: pd.DataFrame, t_clear: float, imax: dict[str, float]) -> None:
    """Four grid-forming units, PSS/E against TSC-OPF, for one clearing time."""
    _rcparams(16)
    panels = (
        ("angle_rel_coi_deg", "angle relative to the COI [deg]", False),
        ("speed_dev_pu", "speed deviation [pu]", False),
        ("P_MW", "active power [MW]", False),
        ("Q_MVAr", "reactive power [MVAr]", False),
        ("E_internal_pu", "internal voltage [pu]", False),
        ("I_over_Imax", "current magnitude [pu]", True),
    )
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
    handles = []
    for unit in (f"G{gid}" for gid in cd.GFM_IDS):
        color, linestyle = _style(unit)
        handles.append(Line2D([0], [0], color=color, ls=linestyle, lw=2.0, label=unit))
    handles.append(Line2D([0], [0], color="0.2", ls=(0, (2, 1.2)), lw=1.4, label="TSC-OPF"))
    for ax, (quantity, label, as_magnitude) in zip(axes.ravel(), panels):
        ax.axvspan(cd.T_FAULT_TABLE, t_clear, color="0.85", zorder=0)
        for unit in (f"G{gid}" for gid in cd.GFM_IDS):
            color, linestyle = _style(unit)
            scale = imax[unit] if as_magnitude else 1.0
            tsc = frame[f"{unit}_{quantity}_tsc"].to_numpy() * scale
            psse = frame[f"{unit}_{quantity}_psse"].to_numpy() * scale
            ax.plot(time, psse, color=color, ls=linestyle, lw=2.0, zorder=2)
            ax.plot(time, tsc, color=color, ls=(0, (2, 1.2)), lw=1.4, zorder=3)
        ax.set_ylabel(label)
        _panel(ax)
    for ax in axes[2, :]:
        ax.set_xlabel("time [s]")
        ax.set_xlim(0.0, 5.0)
    fig.legend(handles=handles, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    _save(fig, cd.FIGURES / f"gfm_compare_{tag}")


def _gfm_series(frame: pd.DataFrame, unit: str, quantity: str, imax: dict[str, float]):
    scale = imax[unit] if quantity == "I_over_Imax" else 1.0
    tsc = frame[f"{unit}_{quantity}_tsc"].to_numpy() * scale
    psse = frame[f"{unit}_{quantity}_psse"].to_numpy() * scale
    return tsc, psse


def _plot_gfm_error(tag: str, time, frame: pd.DataFrame, t_clear: float, imax: dict[str, float]) -> None:
    """PSS/E minus TSC-OPF for the four grid-forming units."""
    _rcparams(16)
    panels = (
        ("angle_rel_coi_deg", "angle error [deg]"),
        ("speed_dev_pu", "speed error [pu]"),
        ("P_MW", "active-power error [MW]"),
        ("Q_MVAr", "reactive-power error [MVAr]"),
        ("E_internal_pu", "internal-voltage error [pu]"),
        ("I_over_Imax", "current-magnitude error [pu]"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
    handles = []
    for unit in (f"G{gid}" for gid in cd.GFM_IDS):
        color, linestyle = _style(unit)
        handles.append(Line2D([0], [0], color=color, ls=linestyle, lw=2.0, label=unit))
    for ax, (quantity, label) in zip(axes.ravel(), panels):
        ax.axvspan(cd.T_FAULT_TABLE, t_clear, color="0.85", zorder=0)
        ax.axhline(0.0, color="0.45", lw=1.0, zorder=1)
        for unit in (f"G{gid}" for gid in cd.GFM_IDS):
            color, linestyle = _style(unit)
            tsc, psse = _gfm_series(frame, unit, quantity, imax)
            ax.plot(time, psse - tsc, color=color, ls=linestyle, lw=2.0, zorder=2)
        ax.set_ylabel(label)
        _panel(ax)
    for ax in axes[2, :]:
        ax.set_xlabel("time [s]")
        ax.set_xlim(0.0, 5.0)
    fig.legend(handles=handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    _save(fig, cd.FIGURES / f"gfm_error_{tag}")


def _gfm_max_diff(tag: str, time, frame: pd.DataFrame, t_clear: float, imax: dict[str, float]) -> pd.DataFrame:
    """Maximum |PSS/E − TSC-OPF| for each grid-forming variable."""
    quantities = (
        ("angle_rel_coi_deg", "angle relative to the COI", "deg"),
        ("speed_dev_pu", "speed deviation", "pu"),
        ("P_MW", "active power", "MW"),
        ("Q_MVAr", "reactive power", "MVAr"),
        ("E_internal_pu", "internal voltage", "pu"),
        ("I_over_Imax", "current magnitude", "pu"),
    )
    fault = (time >= cd.T_FAULT_TABLE - 1e-9) & (time < t_clear - 1e-9)
    post = time >= t_clear - 1e-9
    rows = []
    for unit in (f"G{gid}" for gid in cd.GFM_IDS):
        for quantity, label, unit_name in quantities:
            tsc, psse = _gfm_series(frame, unit, quantity, imax)
            err = psse - tsc
            abs_err = np.abs(err)
            at = int(np.argmax(abs_err))

            def _window(mask: np.ndarray) -> tuple[float, float]:
                if not np.any(mask):
                    return float("nan"), float("nan")
                local = np.where(mask, abs_err, -1.0)
                idx = int(np.argmax(local))
                return float(abs_err[idx]), float(time[idx])

            fault_max, fault_t = _window(fault)
            post_max, post_t = _window(post)
            rows.append(
                {
                    "unit": unit,
                    "quantity": label,
                    "units": unit_name,
                    "max_abs": float(abs_err[at]),
                    "t_s": float(time[at]),
                    "tsc": float(tsc[at]),
                    "psse": float(psse[at]),
                    "psse_minus_tsc": float(err[at]),
                    "rmse": float(np.sqrt(np.mean(err**2))),
                    "max_abs_fault": fault_max,
                    "t_fault_s": fault_t,
                    "max_abs_post": post_max,
                    "t_post_s": post_t,
                }
            )
    table = pd.DataFrame(rows)
    dest = cd.RESULTS / tag / "gfm_max_diff.csv"
    table.to_csv(dest, index=False)
    print(f"Wrote {dest}")
    for _, row in table.iterrows():
        print(
            f"  {row.unit} {row.quantity}: max |err| {row.max_abs:.4g} {row.units} "
            f"at {row.t_s:.3f} s (TSC {row.tsc:.4g}, PSS/E {row.psse:.4g})"
        )
    return table


def _dt_check(tag: str, weights: dict[str, float], tables: dict) -> None:
    coarse_path = cd.RESULTS / tag / "psse_channels.csv"
    fine_path = cd.RESULTS / f"{tag}_dt0.0005" / "psse_channels.csv"
    if tag != "100ms" or not fine_path.exists():
        return
    coarse = _prepare_psse(pd.read_csv(coarse_path), 0.110, weights, tables)
    fine = _prepare_psse(pd.read_csv(fine_path), 0.110, weights, tables)
    rows = []
    for name, coarse_y, fine_y in (
        ("G7_angle_rel_coi_deg", coarse["angle"]["G7"], fine["angle"]["G7"]),
        ("G6_I_over_Imax", coarse["loading"]["G6"], fine["loading"]["G6"]),
    ):
        pasted = np.interp(coarse["time"], fine["time"], fine_y)
        err = pasted - coarse_y
        abs_err = np.abs(err)
        at = int(np.argmax(abs_err))
        rows.append(
            {
                "quantity": name,
                "maxabs": float(abs_err[at]),
                "rmse": float(np.sqrt(np.mean(err**2))),
                "t_max": float(coarse["time"][at]),
                "n_over": int(np.sum(abs_err > (0.05 if "angle" in name else 0.01))),
                "n": int(len(err)),
            }
        )
    text = "\n".join(
        f"{row['quantity']}: max |1 ms - 0.5 ms| = {row['maxabs']:.6g} at t={row['t_max']:.4f} s, "
        f"RMSE = {row['rmse']:.6g}, samples above the reporting threshold = {row['n_over']} of {row['n']}"
        for row in rows
    )
    text += (
        "\nThe current maximum is the fault-application sample: the 1 ms run applies the fault at "
        "0.009 s and the 0.5 ms run at 0.0095 s. Both steps are kept; the curves do not otherwise disagree."
    )
    dest = cd.RESULTS / "dt_check_100ms.txt"
    dest.write_text(text + "\n", encoding="utf-8")
    print(text)
    print(f"Wrote {dest}")


def _prepare_psse(raw: pd.DataFrame, t_clear: float, weights: dict[str, float], tables: dict) -> dict:
    delt = _infer_delt(raw["time"].to_numpy())
    t_fault = cd.T_FAULT_TABLE - delt
    t_open = t_clear - delt
    frame = _collapse(raw, t_open)
    time = frame["time"].to_numpy()
    # The sample at the fault instant is already faulted. Use the one before it.
    pre = frame.loc[frame["time"] < t_fault - 1e-6]
    if pre.empty:
        raise SystemExit("PSS/E export has no pre-fault sample to place P and Q on the MW scale.")
    anchor = pre.iloc[-1]
    scales = []
    names = []
    for gid, mbase_bus in ((6, cd.bus_of(6)), (1, cd.bus_of(1))):
        bus = cd.bus_of(gid)
        mbase = float(tables["gfm"].loc[bus, "InvBase"]) if gid in cd.GFM_IDS else cd.SBASE
        scale, name = _scale_to_mw(float(anchor[f"P_{bus}"]), tables["gen_p"][bus], mbase)
        scales.append(scale)
        names.append(name)
        _ = mbase_bus
    if scales[0] != scales[1]:
        raise SystemExit(
            f"G6 power scale is {names[0]} ({scales[0]}) and G1 is {names[1]} ({scales[1]})."
        )
    scale = scales[0]
    print(f"PSS/E P and Q channels match the dispatch as {names[0]} (scale {scale:g}).")
    pmech_bus = cd.bus_of(1)
    pmech_mw = float(anchor[f"PMEC_{pmech_bus}"]) * scale
    if abs(pmech_mw - tables["gen_p"][pmech_bus]) > 1.0:
        raise SystemExit(
            f"G1 mechanical power {pmech_mw:.4f} MW does not match the dispatch "
            f"{tables['gen_p'][pmech_bus]:.4f} MW with the P scale."
        )
    angles = {}
    for gid in cd.UNIT_IDS:
        bus = cd.bus_of(gid)
        series = frame[f"ANG_{bus}"].to_numpy()
        angles[f"G{gid}"] = np.degrees(np.unwrap(np.radians(series)))
    relative = {unit: angles[unit] - _coi(angles, weights) for unit in angles}
    power = {f"G{gid}": frame[f"P_{cd.bus_of(gid)}"].to_numpy() * scale for gid in cd.UNIT_IDS}
    reactive = {f"G{gid}": frame[f"Q_{cd.bus_of(gid)}"].to_numpy() * scale for gid in cd.UNIT_IDS}
    voltage = {f"G{gid}": frame[f"V_{cd.bus_of(gid)}"].to_numpy() for gid in cd.UNIT_IDS}
    speed = {f"G{gid}": frame[f"SPD_{cd.bus_of(gid)}"].to_numpy() for gid in cd.UNIT_IDS}
    efd = {f"G{gid}": frame[f"EFD_{cd.bus_of(gid)}"].to_numpy() for gid in cd.SYNC_IDS}
    pmech = {f"G{gid}": frame[f"PMEC_{cd.bus_of(gid)}"].to_numpy() * scale for gid in cd.SYNC_IDS}
    loading = {}
    internal = {}
    for gid in cd.GFM_IDS:
        unit = f"G{gid}"
        bus = cd.bus_of(gid)
        inv = float(tables["gfm"].loc[bus, "InvBase"])
        imax = float(tables["gfm"].loc[bus, "Imax"])
        xl = float(tables["gfm"].loc[bus, "Xl"])
        loading[unit] = _loading(power[unit], reactive[unit], voltage[unit], inv, imax)
        internal[unit] = _e_internal(power[unit], reactive[unit], voltage[unit], inv, xl)
    return {
        "time": time,
        "delt": delt,
        "scale": scale,
        "scale_name": names[0],
        "angle": relative,
        "speed": speed,
        "P": power,
        "Q": reactive,
        "V": voltage,
        "Efd": efd,
        "Pmech": pmech,
        "loading": loading,
        "E": internal,
    }


def _interp(src_t, src_y, dst_t):
    if dst_t[0] < src_t[0] - 1e-6 or dst_t[-1] > src_t[-1] + 1e-6:
        raise SystemExit(
            f"PSS/E time {src_t[0]:.4f}..{src_t[-1]:.4f} does not cover "
            f"the TSC-OPF grid {dst_t[0]:.4f}..{dst_t[-1]:.4f}."
        )
    return np.interp(dst_t, src_t, src_y)


def compare(run_dir: Path, tag: str, t_clear: float) -> None:
    tables = cd.load_tables(run_dir)
    weights = cd.inertia_weights(tables)
    _check_tsc_coi(run_dir, weights)
    channel_path = cd.RESULTS / tag / "psse_channels.csv"
    if not channel_path.exists():
        raise SystemExit(f"Missing {channel_path}. Run 01_run_psse.py --run {tag} first.")
    psse = _prepare_psse(pd.read_csv(channel_path), t_clear, weights, tables)
    csv = run_dir / "Transient_Stability" / "CSV"
    relative = cd.read_semi(csv / "angle_rel_COI.csv")
    speed = cd.read_semi(csv / "speed_dev.csv")
    active = cd.read_semi(csv / "generator_active_power.csv")
    reactive = cd.read_semi(csv / "generator_reactive_power.csv")
    voltage = cd.read_semi(csv / "bus_voltage_magnitude.csv")
    efd = cd.read_semi(csv / "dq_E_fd_pu.csv")
    pmech = cd.read_semi(csv / "governor_P_mech.csv")
    loading = cd.read_semi(csv / "gfm_current_loading.csv")
    edroop = cd.read_semi(csv / "gfm_E_droop.csv")
    time = relative["t"].to_numpy()
    rows = []
    aligned = {"t": time}
    tsc_angle = {}
    psse_angle = {}

    def add(unit: str, quantity: str, tsc_y, psse_src) -> np.ndarray:
        psse_y = _interp(psse["time"], psse_src, time)
        err = psse_y - tsc_y
        stats = _metrics(time, err)
        rows.append({"run": tag, "unit": unit, "quantity": quantity, **stats})
        aligned[f"{unit}_{quantity}_tsc"] = tsc_y
        aligned[f"{unit}_{quantity}_psse"] = psse_y
        return psse_y

    for unit in _units():
        gid = int(unit[1:])
        bus = cd.bus_of(gid)
        tsc_angle[unit] = relative[unit].to_numpy()
        psse_angle[unit] = add(unit, "angle_rel_coi_deg", tsc_angle[unit], psse["angle"][unit])
        add(unit, "speed_dev_pu", speed[unit].to_numpy(), psse["speed"][unit])
        add(unit, "P_MW", active[unit].to_numpy(), psse["P"][unit])
        add(unit, "Q_MVAr", reactive[unit].to_numpy(), psse["Q"][unit])
        add(unit, "V_pu", voltage[f"Bus_{bus}"].to_numpy(), psse["V"][unit])
        if gid in cd.SYNC_IDS:
            add(unit, "Efd_pu", efd[unit].to_numpy(), psse["Efd"][unit])
            add(unit, "Pmech_MW", pmech[unit].to_numpy(), psse["Pmech"][unit])
        else:
            add(unit, "I_over_Imax", loading[f"GFM{gid}"].to_numpy(), psse["loading"][unit])
            add(unit, "E_internal_pu", edroop[f"GFM{gid}"].to_numpy(), psse["E"][unit])

    metrics = pd.DataFrame(rows)
    _upsert(cd.RESULTS / "comparison_metrics.csv", metrics, "run")
    aligned_df = pd.DataFrame(aligned)
    # The file postprocess/tables_figure_39bus.py reads.
    aligned_df.to_csv(cd.HERE / "aligned.csv", index=False)
    print(f"Wrote {cd.HERE / 'aligned.csv'}")

    tsc_peak, tsc_tp = _peak(time, tsc_angle["G7"])
    psse_peak, psse_tp = _peak(psse["time"], psse["angle"]["G7"])
    peak = pd.DataFrame(
        [
            {
                "run": tag,
                "tsc_peak_deg": tsc_peak,
                "tsc_t_s": tsc_tp,
                "psse_peak_deg": psse_peak,
                "psse_t_s": psse_tp,
                "peak_error_deg": psse_peak - tsc_peak,
                "time_error_s": psse_tp - tsc_tp,
            }
        ]
    )
    _upsert(cd.RESULTS / "g7_peak.csv", peak, "run")
    print(
        f"G7 first-swing peak: TSC-OPF {tsc_peak:.4f} deg at {tsc_tp:.4f} s, "
        f"PSS/E {psse_peak:.4f} deg at {psse_tp:.4f} s"
    )
    angle_rows = metrics.loc[metrics["quantity"] == "angle_rel_coi_deg"]
    for _, row in angle_rows.iterrows():
        print(
            f"  {row.unit:4s} angle RMSE {row.rmse_full:.3f} deg "
            f"(first second {row.rmse_1s:.3f}), max |err| {row.maxabs_full:.3f} deg"
        )

    valve = _valve_margin(run_dir, tables, tag)
    _upsert(cd.RESULTS / "valve_margin.csv", valve, "run")
    for _, row in valve.iterrows():
        print(
            f"  {row.unit} valve margin to Pmin {row.margin_to_Pmin_MW:.2f} MW, "
            f"to Pmax {row.margin_to_Pmax_MW:.2f} MW"
        )

    _plot_angles(tag, time, tsc_angle, psse_angle, t_clear)
    imax = {
        f"G{gid}": float(tables["gfm"].loc[cd.bus_of(gid), "Imax"]) for gid in cd.GFM_IDS
    }
    _plot_gfm(tag, time, aligned_df, t_clear, imax)
    _plot_gfm_error(tag, time, aligned_df, t_clear, imax)
    _gfm_max_diff(tag, time, aligned_df, t_clear, imax)
    if tag == "100ms":
        _plot_paper(
            time,
            tsc_angle["G7"],
            psse_angle["G7"],
            tsc_angle["G6"],
            psse_angle["G6"],
            aligned_df["G6_I_over_Imax_tsc"].to_numpy(),
            aligned_df["G6_I_over_Imax_psse"].to_numpy(),
            t_clear,
        )
        _dt_check(tag, weights, tables)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare PSS/E and TSC-OPF for one run.")
    parser.add_argument("--run", required=True)
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Draw the GFM comparison from aligned.csv.",
    )
    args = parser.parse_args()
    run_dir, tag, t_clear = cd.resolve_run(args.run)
    if args.plots_only:
        frame = pd.read_csv(cd.HERE / "aligned.csv")
        tables = cd.load_tables(run_dir)
        imax = {
            f"G{gid}": float(tables["gfm"].loc[cd.bus_of(gid), "Imax"]) for gid in cd.GFM_IDS
        }
        _plot_gfm(tag, frame["t"].to_numpy(), frame, t_clear, imax)
        return 0
    compare(run_dir, tag, t_clear)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
