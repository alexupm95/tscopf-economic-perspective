"""
PSS/E dynamic replay of one 39-bus GFM TSC-OPF run.

--flat runs 5 s with no disturbance. Otherwise a three-phase fault is applied
at bus 23 one time step before 0.010 s and, one time step before the first
post-fault TSC-OPF sample, the fault is cleared and line 23-24 is opened.

The fault admittance is [0, -1e6] pu, so the PSS/E call uses [G, B] = [0, -1e8]
MVA. If that does not collapse bus 23, the opposite B sign is tried once.

Usage:
  python 01_run_psse.py --run pqon --flat
  python 01_run_psse.py --run pqon
  python 01_run_psse.py --run pqon --delt 0.0005
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

import case_data as cd

# PSS/E 36.6 install; edit both paths for another install.
PSSE_BIN = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSBIN")
PSSE_PY = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSPY313")
MACH_ID = "1 "
SYNC_BUSES = tuple(cd.bus_of(i) for i in cd.SYNC_IDS)
GEN_BUSES = tuple(cd.bus_of(i) for i in cd.UNIT_IDS)
COLLAPSE_PU = 0.02
_PSSE_READY = False


class CollapseMiss(RuntimeError):
    """Bus 23 stayed up, so this fault susceptance is not the inductive short."""


def _require_psspy():
    if (sys.version_info.major, sys.version_info.minor) != (3, 13):
        raise SystemExit(
            f"PSS/E 36.6 bindings need Python 3.13. This interpreter is {sys.version.split()[0]}."
        )
    os.environ["PATH"] = str(PSSE_BIN) + os.pathsep + os.environ.get("PATH", "")
    sys.path.insert(0, str(PSSE_PY))
    import psse36  # noqa: F401
    import psspy
    import dyntools

    return psspy, dyntools


def _check(ierr, what: str) -> None:
    code = ierr[0] if isinstance(ierr, (tuple, list)) else ierr
    if code != 0:
        raise RuntimeError(f"{what} failed with ierr={ierr}")


def result_dir(tag: str, delt: float, flat: bool, probe: bool, fault_b: float | None) -> Path:
    if flat:
        name = f"flat_{tag}"
    elif abs(delt - 0.001) < 1e-12:
        name = tag
    else:
        name = f"{tag}_dt{delt:.4g}"
    if probe:
        name += f"_probe_{fault_b:.0e}"
    return cd.RESULTS / name


def confirm_tsc_window(run_dir: Path, t_clear: float) -> None:
    """The table instant is the first open-line sample, not manifest clearing_time."""
    voltage = cd.read_semi(run_dir / "Transient_Stability" / "CSV" / "bus_voltage_magnitude.csv")
    time = voltage["t"].to_numpy()
    bus = voltage["Bus_23"].to_numpy()

    def at(instant: float) -> float:
        idx = int(abs(time - instant).argmin())
        if abs(float(time[idx]) - instant) > 1e-6:
            raise SystemExit(f"Bus 23 voltage has no sample at t = {instant}.")
        return float(bus[idx])

    v_on = at(cd.T_FAULT_TABLE)
    v_last = at(t_clear - 0.005)
    v_off = at(t_clear)
    print(
        f"TSC-OPF bus 23 voltage: t=0.010 -> {v_on:.3e}, "
        f"t={t_clear - 0.005:.3f} -> {v_last:.3e}, t={t_clear:.3f} -> {v_off:.4f}"
    )
    if v_on > 1e-3 or v_last > 1e-3 or v_off < 0.5:
        raise SystemExit(
            "Bus 23 in the TSC-OPF file does not match the fault window. "
            "PSS/E events were not applied."
        )


def _channels(psspy) -> list[tuple[int, str]]:
    mapping: list[tuple[int, str]] = []
    channel = 1
    for bus in GEN_BUSES:
        for code, name in ((1, "ANG"), (2, "P"), (3, "Q"), (7, "SPD")):
            _check(
                psspy.machine_array_channel([channel, code, bus], MACH_ID, ""),
                f"{name} bus {bus}",
            )
            mapping.append((channel, f"{name}_{bus}"))
            channel += 1
        _check(psspy.voltage_channel([channel, -1, -1, bus], ""), f"V bus {bus}")
        mapping.append((channel, f"V_{bus}"))
        channel += 1
    for bus in SYNC_BUSES:
        for code, name in ((5, "EFD"), (6, "PMEC")):
            _check(
                psspy.machine_array_channel([channel, code, bus], MACH_ID, ""),
                f"{name} bus {bus}",
            )
            mapping.append((channel, f"{name}_{bus}"))
            channel += 1
    _check(psspy.voltage_channel([channel, -1, -1, cd.FAULT_BUS], ""), "V bus 23")
    mapping.append((channel, "V_23"))
    return mapping


def _export(dyntools, out_path: Path, mapping: list[tuple[int, str]], csv_path: Path) -> pd.DataFrame:
    _title, _chanid, chandata = dyntools.CHNF(str(out_path)).get_data()
    if "time" in chandata:
        time_key = "time"
    elif 0 in chandata:
        time_key = 0
    else:
        raise RuntimeError(f"Channel file has no time key. Keys start {list(chandata)[:8]}")
    frame = {"time": [float(t) for t in chandata[time_key]]}
    for number, name in mapping:
        if number not in chandata:
            raise RuntimeError(f"Channel {number} ({name}) is missing from the output file.")
        frame[name] = [float(v) for v in chandata[number]]
    out = pd.DataFrame(frame)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path}  ({len(out)} rows, t={out['time'].iloc[0]:.4f}..{out['time'].iloc[-1]:.4f})")
    return out


def _log_hits(log: Path) -> list[str]:
    if not log.exists():
        return []
    text = log.read_text(encoding="utf-8", errors="replace")
    needles = (
        "initial condition",
        "out of limit",
        "suspect",
        "illegal",
        "invalid",
        "fatal",
        "not less than",
        "blown up",
        "iteration limit",
        "network not converged",
    )
    hits = []
    for line in text.splitlines():
        low = line.lower()
        if any(needle in low for needle in needles):
            hits.append(line.strip())
    return hits


def _prepare(psspy, sav: Path, dyr: Path, log: Path, out_path: Path, delt: float):
    global _PSSE_READY
    # The previous initialisation keeps the progress file open. Release it
    # before a second psseinit in the same process. A fresh process has no session.
    if _PSSE_READY:
        for _closer in (
            lambda: psspy.progress_output(6, "", [0, 0]),
            lambda: psspy.alert_output(6, "", [0, 0]),
            lambda: psspy.report_output(6, "", [0, 0]),
        ):
            try:
                _closer()
            except Exception:
                pass
    log.parent.mkdir(parents=True, exist_ok=True)
    if log.exists():
        log.unlink()
    _check(psspy.psseinit(200), "psseinit")
    _PSSE_READY = True
    psspy.progress_output(2, str(log), [0, 0])
    psspy.alert_output(2, str(log), [0, 0])
    psspy.prompt_output(6, "", [0, 0])
    psspy.report_output(6, "", [0, 0])
    _check(psspy.case(str(sav)), "case")
    _check(psspy.cong(0), "cong")
    zload = [0.0, 100.0, 0.0, 100.0]
    for apiopt in (1, 2, 3):
        returned = psspy.conl(0, 1, apiopt, [0, 0], zload)
        ierr = returned[0] if isinstance(returned, (tuple, list)) else returned
        _check(ierr, f"conl {apiopt}")
    _check(psspy.ordr(0), "ordr")
    _check(psspy.fact(), "fact")
    _check(psspy.tysl(0), "tysl")
    _check(psspy.dyre_new([1, 1, 1, 1], str(dyr), "", "", ""), "dyre_new")
    mapping = _channels(psspy)
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    _check(
        psspy.dynamics_solution_param_2(
            [50, _i, _i, _i, _i, _i, _i, _i],
            [_f, _f, delt, _f, _f, _f, _f, _f],
        ),
        "dynamics_solution_param_2",
    )
    if out_path.exists():
        out_path.unlink()
    _check(psspy.strt_2([0, 1], str(out_path)), "strt_2")
    return mapping


def _write_init_notes(log: Path, dest: Path, extra: str) -> list[str]:
    hits = _log_hits(log)
    body = extra.rstrip() + "\n"
    if hits:
        body += "\n".join(hits[:200]) + "\n"
    else:
        body += "No initial-condition, limit, or convergence messages matched.\n"
    dest.write_text(body, encoding="utf-8")
    print(f"Wrote {dest}")
    return hits


def _pmech_lifts(log_text: str) -> dict[int, float]:
    """VMAX must sit above the Pm GENROU initialises, which includes I^2 Ra."""
    bus = None
    lifts = {}
    bus_re = re.compile(r"Model TGOV1 Bus\s+(\d+)")
    pmech_re = re.compile(r"PMECH \(pu on Trate\)=\s+([0-9.]+)\s+VMAX=\s+([0-9.]+)")
    for line in log_text.splitlines():
        found = bus_re.search(line)
        if found:
            bus = int(found.group(1))
        found = pmech_re.search(line)
        if found and bus is not None:
            pmech = float(found.group(1))
            vmax = float(found.group(2))
            if pmech > vmax + 1e-6:
                lifts[bus] = pmech + 0.001
    return lifts


def _apply_vmax_lift(tables: dict, log: Path) -> bool:
    if cd.VMAX_PATH.exists():
        return False
    text = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    lifts = _pmech_lifts(text)
    if not lifts:
        return False
    lines = [
        "TGOV1 VMAX was raised once, to 0.001 pu above the mechanical power",
        "PSS/E initialises. That power is the scheduled electrical power plus",
        "stator I^2 Ra, so a unit dispatched at Pmax starts above Pmax/MBASE.",
        "The TSC-OPF valve limit stays at Pmax. This is the only change.",
    ]
    for bus, vmax in sorted(lifts.items()):
        lines.append(f"  bus {bus} ({cd.unit_of(bus)}): VMAX = {vmax:.4f} pu on 100 MVA")
    note = cd.CASE_DIR / "model_adjustments.txt"
    previous = note.read_text(encoding="utf-8") if note.exists() else ""
    note.write_text(previous + "\n" + "\n".join(lines) + "\n", encoding="utf-8")
    cd.VMAX_PATH.write_text(
        "".join(f"{bus},{vmax:.6f}\n" for bus, vmax in sorted(lifts.items())),
        encoding="ascii",
    )
    xl, gap = cd.genrou_adjustment()
    cd.write_dyr(cd.DYR_PATH, tables, xl=xl, xpp_gap=gap)
    print("\n".join(lines))
    return True


def _fatal_init(hits: list[str]) -> str | None:
    for line in hits:
        low = line.lower()
        if "suspect" in low or "initialization failed" in low or "fatal" in low:
            return line
        if "illegal" in low or "invalid" in low or "not less than" in low:
            return line
    return None


def _bus23_during(frame: pd.DataFrame, t0: float, t1: float) -> float:
    part = frame.loc[(frame["time"] >= t0) & (frame["time"] <= t1), "V_23"]
    if part.empty:
        raise RuntimeError(f"No bus 23 voltage samples between {t0} and {t1}.")
    return float(part.min())


def _start(psspy, tables, sav, log, out_path, delt, folder, *, allow_fallback: bool):
    try:
        return _prepare(psspy, sav, cd.DYR_PATH, log, out_path, delt)
    except RuntimeError as exc:
        if not allow_fallback:
            raise
        print(exc)
        _write_fallback(log, folder)
        cd.write_dyr(cd.DYR_PATH, tables, xl=cd.XL_FALLBACK, xpp_gap=cd.XPP_GAP)
        cd.GENROU_PATH.write_text(f"{cd.XL_FALLBACK} {cd.XPP_GAP}\n", encoding="ascii")
        print("Retrying initialisation with Xl = 1e-5 and X\" = X' - 1e-6.")
        try:
            mapping = _start(
                psspy, tables, sav, log, out_path, delt, folder, allow_fallback=False
            )
        except RuntimeError as retry_exc:
            note = cd.CASE_DIR / "model_adjustments.txt"
            note.write_text(
                note.read_text(encoding="utf-8")
                + f"\nFallback initialisation also failed:\n{retry_exc}\n",
                encoding="utf-8",
            )
            raise SystemExit(
                "PSS/E dynamic initialisation failed after the one allowed GENROU change. "
                "See PSSE/case/model_adjustments.txt."
            ) from retry_exc
        note = cd.CASE_DIR / "model_adjustments.txt"
        note.write_text(
            note.read_text(encoding="utf-8") + "\nFallback initialisation was accepted.\n",
            encoding="utf-8",
        )
        return mapping


def _write_fallback(log: Path, folder: Path) -> None:
    fresh = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    note = cd.CASE_DIR / "model_adjustments.txt"
    note.write_text(
        "PSS/E rejected GENROU with Xl = 0 and X\" = X'.\n"
        "One data-check change, not iterated: Xl = 1e-5, X\" = X' - 1e-6.\n"
        f"Run folder: {folder}\n"
        "Progress excerpt:\n"
        f"{fresh[-3000:]}\n",
        encoding="utf-8",
    )
    print(note.read_text(encoding="utf-8"))


def _accept(psspy, tables, sav, log, out_path, delt, folder):
    mapping = _start(psspy, tables, sav, log, out_path, delt, folder, allow_fallback=True)
    if _apply_vmax_lift(tables, log):
        mapping = _start(psspy, tables, sav, log, out_path, delt, folder, allow_fallback=True)
    return mapping


def run_flat(psspy, dyntools, tables: dict, sav: Path, delt: float, folder: Path) -> int:
    log = folder / "psse.log"
    out_path = folder / "psse.out"
    mapping = _accept(psspy, tables, sav, log, out_path, delt, folder)
    hits = _write_init_notes(log, folder / "init_notes.txt", "Flat run, no disturbance.")
    fatal = _fatal_init(hits)
    if fatal:
        raise SystemExit(f"Dynamic initialisation is not usable: {fatal}")
    nprt = max(1, int(round(0.05 / delt)))
    _check(psspy.run(0, cd.T_END, nprt, 1, 0), "run to 5 s")
    frame = _export(dyntools, out_path, mapping, folder / "psse_channels.csv")
    ranges = _drift(frame)
    lines = ["Range over the flat 5 s run (max - min):"]
    failed = False
    # Machine power channels are pu on the system base in the earlier REGFMA1
    # runs (hundreds of MW appear as single-digit pu). MW channels are larger.
    power_is_mw = float(frame["P_30"].abs().max()) > 30.0
    p_limit = 0.05 if power_is_mw else 5e-4
    lines.append(f"Power channels treated as {'MW' if power_is_mw else 'pu on the system base'}.")
    limits = {"ANG": 0.01, "V": 5e-4, "SPD": 1e-4, "P": p_limit, "Q": p_limit, "EFD": 5e-4, "PMEC": p_limit}
    for name, span in ranges.items():
        kind = name.split("_", 1)[0]
        limit = limits[kind]
        flag = "FAIL" if span > limit else "ok"
        if span > limit:
            failed = True
        lines.append(f"  {name:10s}  {span:.3e}  limit {limit:.0e}  {flag}")
    report = "\n".join(lines) + "\n"
    (folder / "flat_report.txt").write_text(report, encoding="utf-8")
    print(report)
    if failed:
        raise SystemExit("The flat run drifted. The fault study will not be run.")
    print("Flat run: nothing moved beyond the tolerance.")
    return 0


def _drift(frame: pd.DataFrame) -> dict[str, float]:
    ranges = {}
    for column in frame.columns:
        if column == "time":
            continue
        kind = column.split("_", 1)[0]
        if kind not in ("ANG", "V", "SPD", "P", "Q", "EFD", "PMEC"):
            continue
        series = frame[column]
        ranges[column] = float(series.max() - series.min())
    return ranges


def run_fault(
    psspy,
    dyntools,
    tables: dict,
    sav: Path,
    delt: float,
    t_clear: float,
    fault_b: float,
    folder: Path,
    probe: bool,
) -> pd.DataFrame:
    t_fault = cd.T_FAULT_TABLE - delt
    t_open = t_clear - delt
    t_end = t_fault + 0.03 if probe else cd.T_END
    log = folder / "psse.log"
    out_path = folder / "psse.out"
    print(
        f"Events: fault at {t_fault:.6f} s, clear and trip line "
        f"{cd.TRIP[0]}-{cd.TRIP[1]} at {t_open:.6f} s, end {t_end:.4f} s, "
        f"delt={delt:.4g} s, B={fault_b:.0e} MVA"
    )
    mapping = _accept(psspy, tables, sav, log, out_path, delt, folder)
    hits = _write_init_notes(
        log,
        folder / "init_notes.txt",
        f"Fault B = {fault_b:.6g} MVA. probe={probe}.",
    )
    fatal = _fatal_init(hits)
    if fatal:
        raise SystemExit(f"Dynamic initialisation is not usable: {fatal}")
    nprt = max(1, int(round(0.05 / delt)))
    _check(psspy.run(0, t_fault, nprt, 1, 0), "run to fault")
    _check(psspy.dist_bus_fault(cd.FAULT_BUS, 1, 0.0, [0.0, fault_b]), "dist_bus_fault")
    if not probe:
        _check(psspy.run(0, t_open, nprt, 1, 0), "run to clear")
        _check(psspy.dist_clear_fault(1), "dist_clear_fault")
        _check(psspy.dist_branch_trip(cd.TRIP[0], cd.TRIP[1], "1 "), "dist_branch_trip")
    _check(psspy.run(0, t_end, nprt, 1, 0), "run to end")
    frame = _export(dyntools, out_path, mapping, folder / "psse_channels.csv")
    v_fault = _bus23_during(frame, t_fault + delt, t_fault + 0.025)
    (folder / "bus23_fault_voltage.txt").write_text(
        f"min V_23 between {t_fault + delt:.6f} and {t_fault + 0.025:.6f} s = {v_fault:.6e} pu\n"
        f"fault_b_mva = {fault_b:.6g}\n",
        encoding="utf-8",
    )
    print(f"Bus 23 voltage during the fault: {v_fault:.6e} pu")
    if v_fault > COLLAPSE_PU:
        raise CollapseMiss(
            f"Bus 23 did not collapse (V = {v_fault:.4f} pu) with B = {fault_b:.0e} MVA."
        )
    return frame


def _probe_sign(psspy, dyntools, tables, sav, delt, t_clear, tag) -> float:
    if cd.FAULT_SIGN_PATH.exists():
        sign = float(cd.FAULT_SIGN_PATH.read_text(encoding="ascii").strip())
        print(f"Using the fault susceptance already accepted: B = {sign:.0e} MVA")
        return sign
    # Inductive short first: Y = -j 1e6 pu -> B = -1e8 MVA. The opposite sign
    # is tried only if bus 23 stays up.
    for sign in (-1.0e8, 1.0e8):
        folder = result_dir(tag, delt, flat=False, probe=True, fault_b=sign)
        print(f"Probe B = {sign:.0e} MVA")
        try:
            run_fault(psspy, dyntools, tables, sav, delt, t_clear, sign, folder, probe=True)
        except CollapseMiss as exc:
            print(exc)
            continue
        cd.FAULT_SIGN_PATH.write_text(f"{sign:.6g}\n", encoding="ascii")
        print(f"Bus 23 collapses with B = {sign:.0e} MVA. That sign is kept.")
        return sign
    raise SystemExit(
        "Neither sign of B collapsed bus 23. The fault was not applied, "
        "and the study stops here."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one PSS/E dynamic replay.")
    parser.add_argument("--run", required=True)
    parser.add_argument("--flat", action="store_true")
    parser.add_argument("--delt", type=float, default=0.001)
    parser.add_argument("--probe", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--fault-b", type=float, default=None)
    args = parser.parse_args()
    if args.delt <= 0.0 or args.delt > 0.005:
        raise SystemExit("--delt must be positive and no coarser than the 5 ms TSC-OPF step.")
    run_dir, tag, t_clear = cd.resolve_run(args.run)
    tables = cd.load_tables(run_dir)
    sav = cd.CASE_DIR / tag / "ieee39.sav"
    dyr = cd.CASE_DIR / tag / "ieee39.dyr"
    if not sav.exists() or not dyr.exists():
        raise SystemExit(f"Run 00_build_case.py --run {tag} first.")
    # Keep a data-check rewrite inside this case. The shared paper record stays.
    cd.DYR_PATH = dyr
    if not args.flat:
        confirm_tsc_window(run_dir, t_clear)
    psspy, dyntools = _require_psspy()
    version = psspy.psseversion()
    print(f"PSS/E {version[1]}.{version[2]}.{version[3]}  {version[0].strip()}")
    if args.flat:
        folder = result_dir(tag, args.delt, flat=True, probe=False, fault_b=None)
        return run_flat(psspy, dyntools, tables, sav, args.delt, folder)
    if args.probe:
        if args.fault_b is None:
            raise SystemExit("--probe requires --fault-b.")
        folder = result_dir(tag, args.delt, flat=False, probe=True, fault_b=args.fault_b)
        run_fault(psspy, dyntools, tables, sav, args.delt, t_clear, args.fault_b, folder, probe=True)
        return 0
    sign = args.fault_b if args.fault_b is not None else _probe_sign(
        psspy, dyntools, tables, sav, args.delt, t_clear, tag
    )
    folder = result_dir(tag, args.delt, flat=False, probe=False, fault_b=None)
    # A new initialisation. The probe, if it ran, has left PSS/E in dynamics.
    run_fault(psspy, dyntools, tables, sav, args.delt, t_clear, sign, folder, probe=False)
    note = cd.CASE_DIR / "model_adjustments.txt"
    if not note.exists():
        note.write_text(
            "No GENROU data-check change. Xl = 0 and X\" = X' were accepted.\n"
            "ZX on the synchronous power-flow record is X'd so CONG has a source "
            "impedance. It is not in the solved power-flow admittance.\n"
            "Converter ZX = 0.15 on InvBase, WMOD = 1, so that impedance is the "
            "REGFMA1 coupling.\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
