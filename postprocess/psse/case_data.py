"""Shared paths, run table, and RAW/DYR writer for the 39-bus GFM PSS/E check.

Does not import psspy. 00_build_case.py solves the power flow; 01_run_psse.py
runs the dynamics; 02_compare_psse_tscopf.py scores the trajectories.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

SBASE = 100.0
FREQ = 50.0
SWING = 39
FAULT_BUS = 23
TRIP = (23, 24)
T_FAULT_TABLE = 0.010
T_END = 5.0
TDO2 = 0.03
TQO2 = 0.04
XL_FALLBACK = 1e-5
XPP_GAP = 1e-6
# gfm_pq_overload = false: zero overload gains and wide P/Q limits.
PQ_OPEN = 10.0
SYNC_IDS = (1, 2, 7, 8, 9, 10)
GFM_IDS = (3, 4, 5, 6)
UNIT_IDS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
EXPECTED_GFM_H = {3: 4.2185, 4: 5.874, 5: 5.401, 6: 5.4285}

# Tag -> (TSC-OPF run folder, first post-fault sample). The run of the paper.
EXTERNAL: dict[str, tuple[Path, float]] = {
    "pqon": (HERE.parents[1] / "results" / "39bus" / "Results - 2026-10-07 124340", 0.110),
}

V_TOL = 1e-4
ANG_TOL_RAD = 1e-4
PQ_TOL_MW = 0.1

CASE_DIR = HERE / "case"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
DYR_PATH = CASE_DIR / "ieee39_gfm.dyr"
MAPPING_PATH = CASE_DIR / "mapping.txt"
FAULT_SIGN_PATH = CASE_DIR / "fault_b_mva.txt"
GENROU_PATH = CASE_DIR / "genrou_adjustment.txt"
VMAX_PATH = CASE_DIR / "tgov1_vmax.txt"


def bus_of(unit_id: int) -> int:
    return 29 + int(unit_id)


def unit_of(bus: int) -> str:
    return f"G{int(bus) - 29}"


def read_semi(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=";", encoding="utf-8-sig")


def resolve_run(name: str) -> tuple[Path, str, float]:
    """Accept a run tag (``pqon``, the run of the paper)."""
    key = name.strip().rstrip("\\/")
    if key not in EXTERNAL:
        raise SystemExit(f"Unknown run {name!r}. Known runs: {', '.join(EXTERNAL)}")
    path, t_clear = EXTERNAL[key]
    if not path.is_dir():
        raise SystemExit(f"Run folder is missing: {path}")
    return path, key, t_clear


def load_tables(run_dir: Path) -> dict:
    inputs = run_dir / "Inputs"
    buses = read_semi(inputs / "bus_data.csv")
    lines = read_semi(inputs / "line_data.csv")
    gens = read_semi(inputs / "generators_data.csv")
    dyn = read_semi(inputs / "gen_dynamic_data_full.csv")
    gfm = read_semi(inputs / "gfm_dynamic_data.csv")
    contingencies = read_semi(inputs / "contingencies.csv")
    bus_rep = read_semi(run_dir / "Dispatch" / "CSV" / "buses_report.csv")
    gen_rep = read_semi(run_dir / "Dispatch" / "CSV" / "generators_report.csv")
    with (run_dir / "run_manifest.toml").open("rb") as handle:
        manifest = tomllib.load(handle)

    _assert_manifest(manifest)
    contingency = contingencies.loc[contingencies["contingency_nr"] == 57]
    if len(contingency) != 1:
        raise SystemExit("contingencies.csv has no unique row for contingency 57.")
    row = contingency.iloc[0]
    if int(row.faulted_bus) != FAULT_BUS or int(row.from_bus) != 23 or int(row.to_bus) != 24:
        raise SystemExit(
            f"Contingency 57 is bus {int(row.faulted_bus)}, "
            f"line {int(row.from_bus)}-{int(row.to_bus)}, expected bus 23, line 23-24."
        )

    pairs = [(min(int(r.fbus), int(r.tbus)), max(int(r.fbus), int(r.tbus))) for _, r in lines.iterrows()]
    if len(pairs) != len(set(pairs)):
        raise SystemExit("line_data.csv has a parallel branch. Circuit ids are not assigned.")

    dyn_bus = set(int(b) for b in dyn["bus"])
    gfm_bus = set(int(b) for b in gfm["bus"])
    if dyn_bus != {bus_of(i) for i in SYNC_IDS}:
        raise SystemExit(f"Synchronous buses {sorted(dyn_bus)} are not { [bus_of(i) for i in SYNC_IDS] }.")
    if gfm_bus != {bus_of(i) for i in GFM_IDS}:
        raise SystemExit(f"Converter buses {sorted(gfm_bus)} are not buses 32-35.")
    if dyn_bus & gfm_bus:
        raise SystemExit("A bus is both a synchronous machine and a converter.")

    for _, machine in dyn.iterrows():
        if abs(float(machine.Xd_tr) - float(machine.Xq_tr)) > 1e-9:
            raise SystemExit(
                f"Bus {int(machine.bus)}: Xd_tr != Xq_tr. GENROU has one subtransient reactance."
            )
        if abs(float(machine.D)) > 1e-12:
            raise SystemExit(f"Bus {int(machine.bus)}: D = {float(machine.D)}, expected 0.")

    gfm_ids = [int(i) for i in manifest["resolved"]["gfm_ids"]]
    gfm_h = [float(h) for h in manifest["resolved"]["gfm_H_eff"]]
    if gfm_ids != list(GFM_IDS):
        raise SystemExit(f"Manifest gfm_ids {gfm_ids} are not {list(GFM_IDS)}.")
    for gid, h_eff in zip(gfm_ids, gfm_h):
        if abs(h_eff - EXPECTED_GFM_H[gid]) > 1e-6:
            raise SystemExit(
                f"G{gid} gfm_H_eff is {h_eff}, expected {EXPECTED_GFM_H[gid]}."
            )

    bus_v = {int(r.BUS): float(r.V_pu) for _, r in bus_rep.iterrows()}
    bus_a = {int(r.BUS): float(r.Theta_deg) for _, r in bus_rep.iterrows()}
    bus_pd = {int(r.BUS): float(r.PD_MW) for _, r in bus_rep.iterrows()}
    bus_qd = {int(r.BUS): float(r.QD_MVAr) for _, r in bus_rep.iterrows()}
    gen_p = {int(r.BUS): float(r.P_MW) for _, r in gen_rep.iterrows()}
    gen_q = {int(r.BUS): float(r.Q_MVAr) for _, r in gen_rep.iterrows()}
    if set(bus_v) != set(int(b) for b in buses["bus_i"]):
        raise SystemExit("buses_report.csv does not cover the same buses as bus_data.csv.")
    for _, bus in buses.iterrows():
        number = int(bus.bus_i)
        if abs(bus_pd[number] - float(bus.Pd)) > 1e-6 or abs(bus_qd[number] - float(bus.Qd)) > 1e-6:
            raise SystemExit(
                f"Bus {number}: dispatch load ({bus_pd[number]}, {bus_qd[number]}) "
                f"differs from bus_data ({float(bus.Pd)}, {float(bus.Qd)})."
            )
    pd_sum = sum(bus_pd.values())
    if abs(pd_sum - 5856.4) > 1e-3:
        raise SystemExit(f"Total load is {pd_sum} MW, expected 5856.4 MW.")

    swing = bus_a[SWING]
    bus_a_shift = {b: (0.0 if b == SWING else a - swing) for b, a in bus_a.items()}
    for gid in GFM_IDS:
        grow = gfm.loc[gfm["id"] == gid].iloc[0]
        for name, expected in (
            ("Xl", 0.15),
            ("mp", 0.01),
            ("mq", 0.05),
            ("Kpv", 0.0),
            ("Kiv", 5.86),
            ("Emax", 1.2),
            ("Emin", 0.0),
            ("Imax", 1.2),
            ("TPf", 0.01),
            ("TQf", 0.01),
            ("TVf", 0.01),
        ):
            if abs(float(grow[name]) - expected) > 1e-9:
                raise SystemExit(f"G{gid} {name} = {float(grow[name])}, expected {expected}.")

    return {
        "buses": buses,
        "lines": lines,
        "gens": gens.set_index("bus", drop=False),
        "dyn": dyn.set_index("bus", drop=False),
        "gfm": gfm.set_index("bus", drop=False),
        "bus_v": bus_v,
        "bus_a": bus_a,
        "bus_a_shift": bus_a_shift,
        "swing_deg": swing,
        "bus_pd": bus_pd,
        "bus_qd": bus_qd,
        "gen_p": gen_p,
        "gen_q": gen_q,
        "gfm_h": {gid: h for gid, h in zip(gfm_ids, gfm_h)},
        "e_max": float(manifest["transient"]["builder"]["limits"]["E_max_pu"]),
        "e_min": float(manifest["transient"]["builder"]["limits"]["E_min_pu"]),
        "pq_overload": bool(manifest["transient"]["dyn_model"]["gfm_pq_overload"]),
        "manifest": manifest,
    }


def _assert_manifest(manifest: dict) -> None:
    if abs(float(manifest["resolved"]["f_syn"]) - FREQ) > 1e-9:
        raise SystemExit(f"f_syn is {manifest['resolved']['f_syn']}, expected {FREQ}.")
    if abs(float(manifest["run"]["base_MVA"]) - SBASE) > 1e-9:
        raise SystemExit(f"base_MVA is {manifest['run']['base_MVA']}, expected {SBASE}.")
    if manifest["resolved"]["gfm_in_coi"] is not True:
        raise SystemExit("gfm_in_coi is not true. The centre of inertia would be wrong.")
    if manifest["transient"]["dyn_model"]["gfm_pq_overload"] not in (True, False):
        raise SystemExit("gfm_pq_overload is neither true nor false.")
    y = [float(v) for v in manifest["transient"]["dyn_model"]["fault"]["sc_fault_admittance_pu"]]
    if abs(y[0]) > 1e-9 or abs(y[1] + 1.0e6) > 1.0:
        raise SystemExit(f"Fault admittance {y} is not [0, -1e6].")
    e_max = float(manifest["transient"]["builder"]["limits"]["E_max_pu"])
    e_min = float(manifest["transient"]["builder"]["limits"]["E_min_pu"])
    if abs(e_max - 6.0) > 1e-9 or abs(e_min) > 1e-9:
        raise SystemExit(f"Exciter limits are ({e_min}, {e_max}), expected (0, 6).")


def genrou_adjustment() -> tuple[float, float]:
    """Xl and X'' gap. Absent file means the CSV mapping Xl = 0, X'' = X'."""
    if not GENROU_PATH.exists():
        return 0.0, 0.0
    xl_text, gap_text = GENROU_PATH.read_text(encoding="ascii").split()
    return float(xl_text), float(gap_text)


def vmax_override() -> dict[int, float]:
    """TGOV1 VMAX in pu on MBASE, only where PSS/E rejected Pmax/MBASE."""
    if not VMAX_PATH.exists():
        return {}
    override = {}
    for line in VMAX_PATH.read_text(encoding="ascii").splitlines():
        if not line.strip():
            continue
        bus_text, value_text = line.split(",")
        override[int(bus_text)] = float(value_text)
    return override


def inertia_weights(tables: dict) -> dict[str, float]:
    weights = {f"G{gid}": float(tables["dyn"].loc[bus_of(gid), "H"]) for gid in SYNC_IDS}
    for gid, h_eff in tables["gfm_h"].items():
        weights[f"G{gid}"] = float(h_eff)
    return weights


def mapping_text(tables: dict) -> str:
    lines = [
        "Parameter mapping from the TSC-OPF CSVs to PSS/E. Bases are stated per row.",
        "System base 100 MVA, system frequency 50 Hz.",
        "Synchronous reactances and H are already on the system base, and MBASE = 100 MVA,",
        "so they are copied unchanged. Converter CONs are on InvBase, and MBASE = InvBase.",
        "",
        "GENROU (system base, seconds and pu):",
        "  Td -> T'd0,  T''d0 = 0.03 s (not in the CSV; GENROU requires it)",
        "  Tq -> T'q0,  T''q0 = 0.04 s (not in the CSV; GENROU requires it)",
        "  H -> H,  D -> D = 0",
        "  Xd -> Xd,  Xq -> Xq,  Xd_tr -> X'd,  Xq_tr -> X'q",
        "  X''d = Xd_tr. Xd_tr = Xq_tr on every machine, and GENROU has one X''.",
        "  Xl = 0,  S(1.0) = 0,  S(1.2) = 0",
        "  Ra -> power-flow ZR (GENROU has no Ra CON). Used by CONG, not by the power-flow network.",
        "  Xd_tr -> power-flow ZX for the same reason. WMOD = 0.",
        "",
        "SEXS (pu, seconds):",
        "  Ta_exc = 0 so TA/TB = 0,  Tb_exc -> TB,  K_exc -> K,  T_exc -> TE",
        f"  EMIN = {tables['e_min']:.1f},  EMAX = {tables['e_max']:.1f}",
        "",
        "TGOV1 (pu on MBASE = 100 MVA, seconds):",
        "  R -> R,  T1 -> T1,  T2 -> T2,  T3 -> T3,  Dt = 0",
        "  VMAX = Pmax/100,  VMIN = Pmin/100. Every Pmin in generators_data.csv is 0.",
        "",
        "REGFMA1, one record per converter, CONs on that converter's InvBase:",
        "  Icons VFlag = 1, QVFlag = 1",
        "  TPf, TQf, TVf, Imax, Emax, Emin, mp, mq, Kpv, Kiv from gfm_dynamic_data.csv",
        (
            "  Pmax, Pmin, Qmax, Qmin, Kppmax, Kipmax, Kpqmax, Kiqmax from gfm_dynamic_data.csv"
            if tables["pq_overload"]
            else "  Pmax, Pmin, Qmax, Qmin = +10, -10, +10, -10 pu"
        ),
        (
            "  Those eight values are the gfm_pq_overload = true record."
            if tables["pq_overload"]
            else "  Kppmax = Kipmax = Kpqmax = Kiqmax = 0"
        ),
        (
            "  Imax = 1.2 stays on."
            if tables["pq_overload"]
            else "  Those eight values are the gfm_pq_overload = false record, not the CSV limits.\n"
            "  The CSV overload gains stay unused. Imax = 1.2 stays on."
        ),
        "  Power flow: MBASE = InvBase, ZR = 0, ZX = Xl = 0.15, WMOD = 1, WPF = 1.",
        "",
        "InvBase (MVA) and centre-of-inertia weight gfm_H_eff (s on the system base):",
    ]
    for gid in GFM_IDS:
        bus = bus_of(gid)
        base = float(tables["gfm"].loc[bus, "InvBase"])
        lines.append(f"  G{gid} bus {bus}: InvBase = {base:.1f}, H_eff = {tables['gfm_h'][gid]:.4f}")
    lines.extend(
        [
            "",
            "Loads: constant MVA in the power flow, then conl to 100% admittance on P and Q.",
            "Fault: dist_bus_fault units=1, [G, B] = [0, -1e8] MVA, from Y = [0, -1e6] pu.",
            "Bus 39 is the slack, so its angle is shifted to 0 before the solve.",
        ]
    )
    return "\n".join(lines) + "\n"


def _winding(tap: float, nomkv: float, rate: float) -> str:
    rates = [0.0] * 12
    rates[0] = rate
    rate_txt = ", ".join(f"{value:7.2f}" for value in rates)
    return (
        f"{tap:8.5f}, {nomkv:7.3f},   0.000, {rate_txt}, 0,     0, 0, "
        "1.10000, 0.90000, 1.10000, 0.90000, 33, 0, 0.00000, 0.00000, 0.00000"
    )


def write_raw(path: Path, tables: dict, tag: str) -> None:
    buses = tables["buses"]
    kv = {int(r.bus_i): float(r.baseKV) for _, r in buses.iterrows()}
    bus_lines = []
    for _, bus in buses.iterrows():
        number = int(bus.bus_i)
        name = f"BUS{number}"
        bus_lines.append(
            f"{number:6d},'{name:<12}', {float(bus.baseKV):8.4f},{int(bus.type)},"
            f"{int(bus.area):4d},{int(bus.zone):4d},   1,"
            f"{tables['bus_v'][number]:12.8f},{tables['bus_a_shift'][number]:12.6f},"
            f"{float(bus.Vmax):7.5f},{float(bus.Vmin):7.5f},"
            f"{float(bus.Vmax):7.5f},{float(bus.Vmin):7.5f}"
        )
    load_lines = []
    shunt_lines = []
    for _, bus in buses.iterrows():
        number = int(bus.bus_i)
        pd_mw = tables["bus_pd"][number]
        qd = tables["bus_qd"][number]
        if abs(pd_mw) > 1e-9 or abs(qd) > 1e-9:
            load_lines.append(
                f"{number:6d},'1 ',1,{int(bus.area):4d},{int(bus.zone):4d},"
                f"{pd_mw:12.5f},{qd:12.5f},     0.000,     0.000,     0.000,     0.000,"
                "   1,    1,0,      0.000,     0.000,     0,''"
            )
        if abs(float(bus.Gs)) > 1e-9 or abs(float(bus.Bs)) > 1e-9:
            shunt_lines.append(
                f"{number:6d},'1 ',1,{float(bus.Gs):10.3f},{float(bus.Bs):10.3f},'SHUNT'"
            )
    gen_lines = []
    for bus in [bus_of(i) for i in UNIT_IDS]:
        static = tables["gens"].loc[bus]
        gfm = bus in set(tables["gfm"].index)
        if gfm:
            mbase = float(tables["gfm"].loc[bus, "InvBase"])
            zr = 0.0
            zx = float(tables["gfm"].loc[bus, "Xl"])
            wmod = 1
        else:
            mbase = SBASE
            zr = float(tables["dyn"].loc[bus, "Ra"])
            zx = float(tables["dyn"].loc[bus, "Xd_tr"])
            wmod = 0
        gen_lines.append(
            f"{bus:6d},'1 ',{tables['gen_p'][bus]:12.6f},{tables['gen_q'][bus]:12.6f},"
            f"{float(static.Qmax):10.4f},{float(static.Qmin):10.4f},{tables['bus_v'][bus]:10.7f},"
            f"{bus:6d},   0,{mbase:10.3f},{zr:11.5E},{zx:11.5E}, 0.00000E+0, 0.00000E+0,"
            f"1.00000,1,  100.0,{float(static.Pmax):10.3f},{float(static.Pmin):10.3f},"
            f" 0,   1,1.0000,   0,   0.0,   0,   0.0,   0,   0.0, {wmod}, 1.0000"
        )
    branch_lines = []
    transformer_blocks = []
    for _, line in tables["lines"].iterrows():
        fbus, tbus = int(line.fbus), int(line.tbus)
        if abs(float(line.ratio) - 1.0) > 1e-8:
            transformer_blocks.append(
                "\n".join(
                    [
                        (
                            f"{fbus:6d},{tbus:6d},     0,'1 ',1,1,1,  0.00000,  0.00000,1,"
                            "'                                        ',1,   1,1.0000,   0,0.0000,"
                            "   0,0.0000,   0,0.0000,'            ',0"
                        ),
                        f"  {float(line.r):11.5E},{float(line.x):11.5E},   {SBASE:6.2f}",
                        _winding(float(line.ratio), kv[fbus], float(line.rateA)),
                        _winding(1.0, kv[tbus], 0.0),
                    ]
                )
            )
            continue
        branch_lines.append(
            f"{fbus:6d},{tbus:6d},'1 ',{float(line.r):11.5E},{float(line.x):11.5E},"
            f"{float(line.b):11.5E},'                                        ',{float(line.rateA):7.2f},"
            "   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,"
            " 0.00000, 0.00000, 0.00000, 0.00000, 1, 0, 1,  0.00,   1,1.0000"
        )
    areas = sorted({int(a) for a in buses["area"]})
    zones = sorted({int(z) for z in buses["zone"]})
    area_lines = []
    for area in areas:
        members = [int(r.bus_i) for _, r in buses.iterrows() if int(r.area) == area]
        isw = SWING if SWING in members else members[0]
        area_lines.append(f"  {area:4d},{isw:6d},     0.000,   10.000,'AREA{area}'")
    zone_lines = [f"  {zone:4d},'ZONE{zone}'" for zone in zones]
    shunt_body = "\n".join(shunt_lines) if shunt_lines else ""
    text = f"""@!IC,SBASE,REV,XFRRAT,NXFRAT,BASFRQ
0,  {SBASE:.2f}, 36,     0,     1, {FREQ:.2f}     / 39-bus GFM {tag}, 50 Hz, TSC-OPF dispatch
39-BUS SIX SYNCHRONOUS MACHINES AND FOUR REGFMA1 CONVERTERS
39 BUSES
GENERAL, THRSHZ=0.0001, PQBRAK=0.7, BLOWUP=5.0, MaxIsolLvls=4, CAMaxReptSln=20, ChkDupCntLbl=0
GAUSS, ITMX=100, ACCP=1.6, ACCQ=1.6, ACCM=1.0, TOL=0.00000001
NEWTON, ITMXN=80, ACCN=1.0, TOLN=0.00000001, VCTOLQ=0.00000001, VCTOLV=0.00000001, DVLIM=0.99, NDVFCT=0.99
ADJUST, ADJTHR=0.005, ACCTAP=1.0, TAPLIM=0.05, SWVBND=100.0, MXTPSS=99, MXSWIM=10
TYSL, ITMXTY=40, ACCTY=1.0, TOLTY=0.00000001
SOLVER, FDNS, ACTAPS=0, AREAIN=0, PHSHFT=0, DCTAPS=1, SWSHNT=1, FLATST=0, VARLIM=99, NONDIV=0
RATING, 1, "RATE1 ", "RATING SET 1                    "
0 / END OF SYSTEM-WIDE DATA, BEGIN BUS DATA
@!   I,'NAME        ', BASKV, IDE,AREA,ZONE,OWNER, VM,        VA,    NVHI,   NVLO,   EVHI,   EVLO
{chr(10).join(bus_lines)}
0 / END OF BUS DATA, BEGIN LOAD DATA
@!   I,'ID',STAT,AREA,ZONE,      PL,        QL,        IP,        IQ,        YP,        YQ, OWNER,SCALE,INTRPT
{chr(10).join(load_lines)}
0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA
@!   I,'ID',STATUS,  GL,        BL
{shunt_body}
0 / END OF FIXED SHUNT DATA, BEGIN VOLTAGE DROOP CONTROL DATA
0 / END OF VOLTAGE DROOP CONTROL DATA, BEGIN GENERATOR DATA
@!   I,'ID',      PG,        QG,        QT,        QB,     VS,    IREG,NREG,     MBASE,     ZR,         ZX,         RT,         XT,     GTAP,STAT, RMPCT,      PT,        PB,BASLOD,O1,    F1,WMOD, WPF
{chr(10).join(gen_lines)}
0 / END OF GENERATOR DATA, BEGIN SWITCHING DEVICE RATING SET DATA
0 / END OF SWITCHING DEVICE RATING SET DATA, BEGIN BRANCH DATA
@!   I,     J,'CKT',      R,           X,       B, 'N A M E', RATE1
{chr(10).join(branch_lines)}
0 / END OF BRANCH DATA, BEGIN SYSTEM SWITCHING DEVICE DATA
0 / END OF SYSTEM SWITCHING DEVICE DATA, BEGIN TRANSFORMER DATA
@!   I,     J,     K,'CKT',CW,CZ,CM,     MAG1,        MAG2,NMETR,'N A M E',STAT,O1,  F1
@!   R1-2,       X1-2, SBASE1-2
@!WINDV1, NOMV1,   ANG1, RATA1, COD1
@!WINDV2, NOMV2
{chr(10).join(transformer_blocks)}
0 / END OF TRANSFORMER DATA, BEGIN AREA DATA
@! I,   ISW,    PDES,     PTOL,    'ARNAME'
{chr(10).join(area_lines)}
0 / END OF AREA DATA, BEGIN TWO-TERMINAL DC DATA
0 / END OF TWO-TERMINAL DC DATA, BEGIN VSC DC LINE DATA
0 / END OF VSC DC LINE DATA, BEGIN IMPEDANCE CORRECTION DATA
0 / END OF IMPEDANCE CORRECTION DATA, BEGIN MULTI-TERMINAL DC DATA
0 / END OF MULTI-TERMINAL DC DATA, BEGIN MULTI-SECTION LINE DATA
0 / END OF MULTI-SECTION LINE DATA, BEGIN ZONE DATA
@! I,   'ZONAME'
{chr(10).join(zone_lines)}
0 / END OF ZONE DATA, BEGIN INTER-AREA TRANSFER DATA
0 / END OF INTER-AREA TRANSFER DATA, BEGIN OWNER DATA
@! I,   'OWNAME'
   1,'OWNER1'
0 / END OF OWNER DATA, BEGIN FACTS DEVICE DATA
0 / END OF FACTS DEVICE DATA, BEGIN SWITCHED SHUNT DATA
0 / END OF SWITCHED SHUNT DATA, BEGIN GNE DATA
0 / END OF GNE DATA, BEGIN INDUCTION MACHINE DATA
0 / END OF INDUCTION MACHINE DATA, BEGIN LOAD TYPE DATA
0 / END OF LOAD TYPE DATA, BEGIN INTERFACE DATA
0 / END OF INTERFACE DATA, BEGIN INTERFACE ELEMENT DATA
0 / END OF INTERFACE ELEMENT DATA, BEGIN SUBSTATION DATA
0 / END OF SUBSTATION DATA
Q
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="ascii")


def write_dyr(path: Path, tables: dict, *, xl: float, xpp_gap: float) -> None:
    blocks: list[str] = ["@!  BUS, 'MODEL', ID,  CON... /"]
    raised_vmax = vmax_override()
    for gid in SYNC_IDS:
        bus = bus_of(gid)
        machine = tables["dyn"].loc[bus]
        static = tables["gens"].loc[bus]
        xd1 = float(machine.Xd_tr)
        xq1 = float(machine.Xq_tr)
        xpp = xd1 - xpp_gap
        if not xpp > xl:
            raise SystemExit(f"Bus {bus}: X'' = {xpp} is not above Xl = {xl}.")
        pmax = raised_vmax.get(bus, float(static.Pmax) / SBASE)
        pmin = float(static.Pmin) / SBASE
        blocks.append(
            f"  {bus} 'GENROU' 1  {float(machine.Td):.8g} {TDO2:.8g} {float(machine.Tq):.8g} "
            f"{TQO2:.8g} {float(machine.H):.8g} {float(machine.D):.8g} {float(machine.Xd):.8g}\n"
            f"  {float(machine.Xq):.8g} {xd1:.8g} {xq1:.8g} {xpp:.8g} {xl:.8g} 0 0  /"
        )
        blocks.append(
            f"  {bus} 'SEXS' 1  0 {float(machine.Tb_exc):.8g} {float(machine.K_exc):.8g} "
            f"{float(machine.T_exc):.8g} {tables['e_min']:.8g} {tables['e_max']:.8g}  /"
        )
        blocks.append(
            f"  {bus} 'TGOV1' 1  {float(machine.R):.8g} {float(machine.T1):.8g} {pmax:.8g} "
            f"{pmin:.8g} {float(machine.T2):.8g} {float(machine.T3):.8g} 0  /"
        )
    for gid in GFM_IDS:
        bus = bus_of(gid)
        converter = tables["gfm"].loc[bus]
        if tables["pq_overload"]:
            limits = [
                float(converter.Pmax),
                float(converter.Pmin),
                float(converter.Qmax),
                float(converter.Qmin),
            ]
            gains = [
                float(converter.Kppmax),
                float(converter.Kipmax),
                float(converter.Kpqmax),
                float(converter.Kiqmax),
            ]
        else:
            limits = [PQ_OPEN, -PQ_OPEN, PQ_OPEN, -PQ_OPEN]
            gains = [0.0, 0.0, 0.0, 0.0]
        cons = [
            float(converter.TPf),
            float(converter.TQf),
            float(converter.TVf),
            float(converter.Imax),
            float(converter.Emax),
            float(converter.Emin),
            *limits,
            float(converter.mp),
            float(converter.mq),
            *gains,
            float(converter.Kpv),
            float(converter.Kiv),
        ]
        body = " ".join(f"{value:.8g}" for value in cons)
        blocks.append(f"  {bus} 'REGFMA1' 1 1 1 {body}/")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(blocks) + "\n", encoding="ascii")
