# PSS/E replay of the 39-bus run

PSS/E Xplore 36.6 replays the dispatch of the 39-bus run of the paper,
`results/39bus/Results - 2026-10-07 124340` (tag `pqon`: converter P/Q overload loops on). The
network, the machine data and the pre-fault point are read from that folder; TSC-OPF is not solved
again. `aligned.csv` puts both trajectories on the TSC-OPF time grid and is the file
`postprocess/tables_figure_39bus.py` reads.

PSS/E is commercial software and is not included. To rerun the PSS/E steps you need your own
licensed PSS/E 36.6 and a Python that its `psspy` supports, plus `pandas` and `matplotlib`
(Python 3.11 or later; `02_compare_psse_tscopf.py` needs no PSS/E).

`00_build_case.py` and `01_run_psse.py` load `psspy` from a fixed install folder, set in
`PSSE_BIN` and `PSSE_PY` at the top of each script:

```python
PSSE_BIN = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSBIN")
PSSE_PY = Path(r"C:\Program Files\PTI\PSSE36\36.6\PSSPY313")
```

PSS/E is not a Python package, so the scripts cannot find it on their own: they add these two
folders to `PATH` and `sys.path` before `import psspy`. Set both to your installation. `PSSPY313`
is the `psspy` build for Python 3.13; use the `PSSPY3xx` folder that matches your Python, or the
import fails. Every other path is relative to this folder and to the repository.

`results/pqon/psse.log` and `pf_progress.txt` are the PSS/E logs of the run, with the local folder
of the original machine replaced by `postprocess\psse\`.

From this folder:

```powershell
python 00_build_case.py --run pqon             # writes case/pqon/ieee39.{raw,dyr,sav}, solves the power flow
python 01_run_psse.py --run pqon --flat        # 5 s with no disturbance: initialisation check
python 01_run_psse.py --run pqon               # the fault
python 02_compare_psse_tscopf.py --run pqon    # aligned.csv, metrics and figures; needs no PSS/E
```

## Files

| Path | Contents |
|---|---|
| `case_data.py` | Run table, RAW and DYR writer, shared constants |
| `case/pqon/ieee39.raw`, `ieee39.dyr`, `ieee39.sav` | The PSS/E case of the run |
| `case/pqon/mapping.txt` | TSC-OPF column to PSS/E CON, with bases |
| `case/tgov1_vmax.txt`, `case/model_adjustments.txt` | The two raised TGOV1 ceilings |
| `case/fault_b_mva.txt` | Fault susceptance that shorts bus 23 |
| `results/pqon/psse_channels.csv`, `psse.out` | Every PSS/E channel, as CSV and as the binary output file |
| `results/pqon/psse.log`, `init_notes.txt` | Dynamics log and initialisation check |
| `results/pqon/pf_mismatch.csv`, `pf_progress.txt` | Power flow against the TSC-OPF dispatch |
| `results/pqon/bus23_fault_voltage.txt` | Bus 23 voltage during the fault |
| `results/pqon/gfm_max_diff.csv` | Largest converter differences |
| `results/comparison_metrics.csv` | RMSE and largest absolute error, every unit and quantity |
| `results/g7_peak.csv` | First-swing peak of G7 |
| `results/valve_margin.csv` | TSC-OPF valve command against Pmin and Pmax |
| `figures/` | Angle and converter figures, PDF and SVG |
| `aligned.csv` | Both trajectories on the TSC-OPF time grid |

## The case

System base 100 MVA, 50 Hz. G1, G2 and G7–G10 are synchronous machines on buses 30, 31 and 36–39
(GENROU, SEXS, TGOV1). G3–G6 are converters on buses 32–35 (REGFMA1), with `InvBase` 843.7,
1174.8, 1080.2 and 1085.7 MVA. Bus 39 is the slack. Loads are the constant-MVA dispatch in the
power flow and constant admittance in the dynamics.

`case/pqon/mapping.txt` is the full parameter map. In short: the TSC-OPF machine data are copied
unchanged; GENROU also needs `T''d0 = 0.03` s and `T''q0 = 0.04` s, which TSC-OPF does not have,
with `X'' = X'd`; SEXS has `EMIN = 0`, `EMAX = 6`; TGOV1 has `Dt = 0` and `VMAX = Pmax/100`. The
REGFMA1 records carry every converter column of `gfm_dynamic_data.csv`, overload limits and gains
included, because the run has `gfm_pq_overload = true`. The centre of inertia weights the machines
by `H` and the converters by `gfm_H_eff` (4.2185, 5.874, 5.401, 5.4285 s).

The fault is applied at bus 23 at 0.009 s and cleared, with line 23–24 opened, at 0.109 s. PSS/E
steps at 1 ms and applies each event one step before the first TSC-OPF sample that shows it
(0.010 s faulted, 0.110 s open), so the sample at that instant is the post-event state.

Two TGOV1 ceilings were raised once, by 0.001 pu above the mechanical power PSS/E initialises
(the scheduled electrical power plus stator `I²Ra`): G2 to 6.670 pu and G9 to 8.514 pu.

## Results

The power flow matches the dispatch to `1.5×10⁻⁷` pu, `1.3×10⁻⁷` rad, `6.6×10⁻⁵` MW and
`8.2×10⁻⁴` MVAr. The dynamics report `Initial conditions check O.K.` Bus 23 falls to
`5.3×10⁻⁵` pu during the fault.

G7 peaks at 100.00 deg at 0.225 s in TSC-OPF and at 99.07 deg at 0.227 s in PSS/E. Angle RMSE
over the five seconds runs from 0.038 deg (G10) to 1.38 deg (G2); the largest angle difference
is 1.96 deg (G7). Converter current loading differs by at most 0.012 of `Imax`.

## Differences that remain

- GENROU needs `T''d0` and `T''q0`; TSC-OPF has no subtransient dynamics.
- REGFMA1 and the TSC-OPF converter do not share current-limit logic. The internal voltage
  during the fault differs most (up to 0.36 pu); the current loading agrees.
- TSC-OPF uses a smooth governor limit and TGOV1 a hard one. The TSC-OPF valve command stays
  inside `[Pmin, Pmax]` on every machine (`results/valve_margin.csv`), so the hard limit does not
  act. G2's mechanical power in PSS/E stays about 17 MW above TSC-OPF, which is the `I²Ra`
  initialisation.
- TSC-OPF integrates at 5 ms, PSS/E at 1 ms.

No parameter was tuned to reduce these differences.
