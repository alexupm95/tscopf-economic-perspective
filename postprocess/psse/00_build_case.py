"""
Build the PSS/E 36 RAW/DYR and solve one TSC-OPF dispatch.

Writes PSSE/case/<tag>/ieee39.raw, ieee39.sav and the shared DYR.
Stops if voltage, angle, or generator P/Q miss the dispatch.
Does not run a dynamic simulation and does not modify the Results folders.

Usage:
  python 00_build_case.py --run "Results - 2026-10-06 174237"
  python 00_build_case.py --run pqon
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import pandas as pd

import case_data as cd

# PSS/E 36.6 install; edit both paths for another install.
PSSE_BIN = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSBIN")
PSSE_PY = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSPY313")
MACH_ID = "1 "


def _require_psspy():
    if (sys.version_info.major, sys.version_info.minor) != (3, 13):
        raise SystemExit(
            f"PSS/E 36.6 bindings need Python 3.13. This interpreter is {sys.version.split()[0]}."
        )
    if not PSSE_PY.is_dir():
        raise SystemExit(f"PSS/E Python path is missing: {PSSE_PY}")
    os.environ["PATH"] = str(PSSE_BIN) + os.pathsep + os.environ.get("PATH", "")
    sys.path.insert(0, str(PSSE_PY))
    try:
        import psse36  # noqa: F401
        import psspy
    except Exception as exc:
        raise SystemExit(f"psspy could not be imported: {exc}") from exc
    return psspy


def _check(ierr, what: str) -> None:
    code = ierr[0] if isinstance(ierr, (tuple, list)) else ierr
    if code != 0:
        raise RuntimeError(f"{what} failed with ierr={ierr}")


def _silence(psspy, log: Path) -> None:
    _check(psspy.psseinit(200), "psseinit")
    log.parent.mkdir(parents=True, exist_ok=True)
    if log.exists():
        log.unlink()
    psspy.progress_output(2, str(log), [0, 0])
    psspy.alert_output(2, str(log), [0, 0])
    psspy.prompt_output(6, "", [0, 0])
    psspy.report_output(6, "", [0, 0])


def _tighten(psspy) -> None:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    ints = [_i, 80, 40, _i, _i, _i, _i, _i, _i]
    reals = [_f] * 21
    reals[5] = 1e-10
    reals[7] = 1e-10
    reals[16] = 1e-10
    _check(psspy.solution_parameters_5(ints, reals), "solution_parameters_5")


def _snapshot(psspy, tables: dict) -> pd.DataFrame:
    rows = []
    for bus in sorted(tables["bus_v"]):
        ierr_v, voltage = psspy.busdat(bus, "PU")
        ierr_a, angle = psspy.busdat(bus, "ANGLED")
        _check(ierr_v, f"busdat PU {bus}")
        _check(ierr_a, f"busdat ANGLED {bus}")
        rows.append(
            {
                "quantity": "V_pu",
                "bus": bus,
                "unit": "",
                "psse": float(voltage),
                "tscopf": tables["bus_v"][bus],
                "error": float(voltage) - tables["bus_v"][bus],
                "unit_of_error": "pu",
            }
        )
        error_rad = math.radians(float(angle) - tables["bus_a_shift"][bus])
        rows.append(
            {
                "quantity": "theta_rad",
                "bus": bus,
                "unit": "",
                "psse": math.radians(float(angle)),
                "tscopf": math.radians(tables["bus_a_shift"][bus]),
                "error": error_rad,
                "unit_of_error": "rad",
            }
        )
    for gid in cd.UNIT_IDS:
        bus = cd.bus_of(gid)
        ierr_p, power = psspy.macdat(bus, MACH_ID, "P")
        ierr_q, reactive = psspy.macdat(bus, MACH_ID, "Q")
        _check(ierr_p, f"macdat P {bus}")
        _check(ierr_q, f"macdat Q {bus}")
        rows.append(
            {
                "quantity": "P_MW",
                "bus": bus,
                "unit": f"G{gid}",
                "psse": float(power),
                "tscopf": tables["gen_p"][bus],
                "error": float(power) - tables["gen_p"][bus],
                "unit_of_error": "MW",
            }
        )
        rows.append(
            {
                "quantity": "Q_MVAr",
                "bus": bus,
                "unit": f"G{gid}",
                "psse": float(reactive),
                "tscopf": tables["gen_q"][bus],
                "error": float(reactive) - tables["gen_q"][bus],
                "unit_of_error": "MVAr",
            }
        )
    return pd.DataFrame(rows)


def _worst(df: pd.DataFrame) -> tuple[float, str]:
    limits = {
        "V_pu": cd.V_TOL,
        "theta_rad": cd.ANG_TOL_RAD,
        "P_MW": cd.PQ_TOL_MW,
        "Q_MVAr": cd.PQ_TOL_MW,
    }
    messages = []
    failed = False
    for quantity, limit in limits.items():
        part = df.loc[df["quantity"] == quantity]
        idx = part["error"].abs().idxmax()
        row = part.loc[idx]
        err = abs(float(row.error))
        flag = "FAIL" if err > limit else "ok"
        if err > limit:
            failed = True
        where = f"bus {int(row.bus)}"
        if row.unit:
            where = f"{row.unit} {where}"
        messages.append(
            f"  {quantity:10s}  max |err|={err:.3e} {row.unit_of_error}  "
            f"limit {limit:.0e}  {flag}  at {where}  "
            f"(PSS/E {float(row.psse):.8g}, TSC-OPF {float(row.tscopf):.8g})"
        )
    return (1.0 if failed else 0.0), "\n".join(messages)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and solve one 39-bus GFM PSS/E case.")
    parser.add_argument("--run", required=True, help="Run tag: pqon")
    args = parser.parse_args()
    run_dir, tag, _t_clear = cd.resolve_run(args.run)
    tables = cd.load_tables(run_dir)
    case_dir = cd.CASE_DIR / tag
    case_dir.mkdir(parents=True, exist_ok=True)
    result_dir = cd.RESULTS / tag
    result_dir.mkdir(parents=True, exist_ok=True)
    raw = case_dir / "ieee39.raw"
    sav = case_dir / "ieee39.sav"
    log = result_dir / "pf_progress.txt"

    cd.write_raw(raw, tables, tag)
    xl, gap = cd.genrou_adjustment()
    dyr = case_dir / "ieee39.dyr"
    cd.write_dyr(dyr, tables, xl=xl, xpp_gap=gap)
    mapping = cd.mapping_text(tables)
    if tables["pq_overload"]:
        (case_dir / "mapping.txt").write_text(mapping, encoding="utf-8")
    else:
        cd.DYR_PATH.write_text(dyr.read_text(encoding="ascii"), encoding="ascii")
        cd.MAPPING_PATH.write_text(mapping, encoding="utf-8")
    print(f"Slack bus {cd.SWING} angle in the TSC-OPF dispatch is {tables['swing_deg']:.6f} deg.")
    print("PSS/E holds that bus at 0 deg. Relative angles are unchanged.")

    psspy = _require_psspy()
    version = psspy.psseversion()
    print(f"PSS/E {version[1]}.{version[2]}.{version[3]}  {version[0].strip()}")
    _silence(psspy, log)
    _check(psspy.read(0, str(raw)), "read raw")
    for gid in cd.GFM_IDS:
        bus = cd.bus_of(gid)
        ierr, wmod = psspy.macint(bus, MACH_ID, "WMOD")
        _check(ierr, f"macint WMOD {bus}")
        print(f"bus {bus} WMOD = {int(wmod)}")
        if int(wmod) != 1:
            raise SystemExit(
                f"Bus {bus} WMOD = {int(wmod)}. REGFMA1 requires a renewable machine (WMOD = 1)."
            )
    raw_log = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    if "Input error" in raw_log:
        excerpt = "\n".join(line for line in raw_log.splitlines() if "error" in line.lower())
        raise SystemExit("RAW was not accepted:\n" + excerpt[:4000])
    _tighten(psspy)
    # Same flags as the working 39-bus replay: taps, area interchange and
    # phase shifters locked; no flat start; VAR limits applied.
    _check(psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0]), "fnsl")
    snap = _snapshot(psspy, tables)
    failed, report = _worst(snap)
    print("=== power flow vs TSC-OPF dispatch ===")
    print(report)
    snap.to_csv(result_dir / "pf_mismatch.csv", index=False)
    if failed:
        raise SystemExit(
            "Power-flow mismatch exceeds the tolerance. Dynamics will not be run. "
            f"See {result_dir / 'pf_mismatch.csv'}"
        )
    _check(psspy.save(str(sav)), "save")
    print(f"Wrote {raw}")
    print(f"Wrote {dyr}")
    print(f"Wrote {sav}")
    print(f"Wrote {result_dir / 'pf_mismatch.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
