# Dynamic controls: AVR and turbine governor

**How to configure** the optional **AVR** (DQ_4TH path) and **TGOV1 turbine governor** (FULL_BUS, classical or DQ_4TH), and what each writes to disk.

> **Equations, block diagrams, discretisation, and the limiter trade-offs live in Part I, [8. Machine controls: AVR and turbine governor](model/08_controls_avr_governor.md).** This page is the operational counterpart: presets, output files, and dual exports.

Implementation: [`functions_4_TS_avr.jl`](https://github.com/alexupm95/tscopf-economic-perspective/blob/main/TSCOPF/src/_transient_stability/functions_4_TS_avr.jl), [`functions_4_TS_governor.jl`](https://github.com/alexupm95/tscopf-economic-perspective/blob/main/TSCOPF/src/_transient_stability/functions_4_TS_governor.jl).

---

## AVR (automatic voltage regulator) — `include_avr=true`, `gen_order=DQ_4TH`

SEXS exciter driving field voltage `E_fd` from terminal-voltage feedback: an optional
lead-lag stage `(1+s·Ta_exc)/(1+s·Tb_exc)` feeding a gain-lag `K_exc/(1+s·T_exc)`,
saturated by a smooth sqrt clamp. Parameters `T_exc`, `K_exc`, `Ta_exc`, `Tb_exc`
come from `gen_dynamic_data_full.csv`; the clamp limits come from
`TsBoundLimitsConfig.E_min_pu` / `E_max_pu` — the **same pair** that bounds the classical
internal emf, so raising the field ceiling also widens that box.

`Ta_exc = Tb_exc = 0` bypasses the lead-lag (no variables, no rows) — every shipped
fixture uses this. `Tb_exc = 0` with `Ta_exc > 0` is rejected at validation, as are a
non-positive `T_exc` or `K_exc` (both are divided by).

Pre-fault explicit `bound_E` on `E_fd` still applies when `TsBuilderConfig.bound_E=true` (same `[E_min, E_max]` box as the transient AVR clamp limits). Optional research boxes on `E_fd_unlim_*`, `V_ref`, `Pv_*`, and `Pm_*` use separate `bound_*` toggles (default off).

### Config preset

```julia
DynModelConfig(
    gen_order = DQ_4TH,
    network_form = FULL_BUS,
    mech_power_mode = USE_PM,
    bound_style_δ = :coi_box,
    include_avr = true,
    # transient.gen_dynamic_filename = "gen_dynamic_data_full.csv"
)
```

### Outputs

| Artifact | Description |
|---|---|
| `Transient_Stability/CSV/avr_V_ref.csv` | Constant set-point per generator (p.u.) |
| `Transient_Stability/CSV/dq_E_fd_pu.csv` | Time-varying **saturated** `E_fd` trajectory |
| `Transient_Stability/CSV/dq_E_fd_unlim_pu.csv` | **Pre-saturation** `E_fd_unlim` trajectory; the gap to `E_fd` is the active clamp |
| `Transient_Stability/CSV/dq_E_LL_pu.csv` | Lead-lag output `E_LL` (absent when every generator sets `Ta_exc=Tb_exc=0`) |
| `Transient_Stability/Figures/dq_E_fd_traj.svg`, `dq_E_fd_unlim.svg`, `dq_E_LL.svg` | Trajectories when `save_ts_plots=true` |
| `Transient_Stability/CSV_duals/dual_Vref_init.csv` | Dual of pre-fault exciter link |
| `Transient_Stability/CSV_duals/dual_avr_leadlag.csv` | Dual of lead-lag stage (when present) |
| `Transient_Stability/CSV_duals/dual_avr_E_fd.csv` | Dual of exciter ODE on `E_fd_unlim` (fault + post-fault) |
| `Transient_Stability/CSV_duals/dual_avr_E_fd_sat.csv` | Dual of smooth field-voltage clamp |

---

## TGOV1 turbine governor — `include_governor=true`, FULL_BUS (classical or DQ_4TH)

IEEE TGOV1 primary-frequency control. Replaces constant `P_m` in the swing equation with time-varying `P_mech(t)`. Parameters `R, T1, T2, T3` from `gen_dynamic_data_full.csv`.

Available on **classical FULL_BUS** and **DQ_4TH FULL_BUS** (Milestone 3). Combine with `include_avr=true` on the DQ path for reference-style exciter + governor.

Note that the set-point is `R`-scaled — `P_ref = R·P_m` at the equilibrium — so the values in `governor_P_ref.csv` are not powers. See [Part I §8](model/08_controls_avr_governor.md) for why.

**Valve saturation** (`governor_limiter`): `GOV_NO_LIMIT` (default), `GOV_SMOOTH`, or `GOV_HARD_BOUND`. The choice trades physical fidelity against dual quality and, for `GOV_HARD_BOUND`, against feasibility — the comparison table is in [Part I §8](model/08_controls_avr_governor.md).

### Outputs

| Artifact | Description |
|---|---|
| `Transient_Stability/CSV/governor_P_mech.csv` | Mechanical power trajectory |
| `Transient_Stability/CSV/governor_P_valve.csv` | Valve output trajectory (post-limiter) |
| `Transient_Stability/CSV/governor_P_valve_raw.csv` | Raw integrator state ahead of the limiter; equals `P_valve` unless `GOV_SMOOTH` is clamping |
| `Transient_Stability/CSV/governor_P_ref.csv` | Set-point (p.u.) |
| `Transient_Stability/Figures/governor_P_mech.svg`, `governor_P_valve.svg`, `governor_P_valve_raw.svg` | The three trajectories when `save_ts_plots=true` |
| `Transient_Stability/CSV_duals/dual_Pref_init.csv`, `dual_gov_valve.csv`, `dual_gov_mech.csv` | Constraint duals |
| `Transient_Stability/CSV_duals/dual_gov_valve_sat.csv` | Dual of the `GOV_SMOOTH` clamp equality (this mode only) |
| `Transient_Stability/CSV_duals/dual_LB_gov_valve.csv`, `dual_UB_gov_valve.csv` | Duals of the `GOV_HARD_BOUND` valve bounds (this mode only) |

> `mechanical_power.csv` is not written when a governor is present: with `USE_PM` the
> mechanical-power trajectory *is* `governor_P_mech.csv`, and the two files were byte-identical.
> Runs without a governor still get `mechanical_power.csv`.

---

## First-step discretisation

TSCOPF anchors the `t = 1` row of the exciter and governor windows **trapezoidally**, where other implementations take a backward-Euler first step. The physics is identical either way; only the first row of each window differs. Set `DynModelConfig.ode_first_step = :backward_euler` to reproduce that first step exactly. The switch reaches the swing, EMF, exciter and governor rows on both FULL_BUS paths (classical and DQ); on `KRON_REDUCED` it is rejected, because those swing rows have no backward-Euler form.

## See also

- [Part I §8 — Machine controls](model/08_controls_avr_governor.md): equations, block diagrams, discretisation, limiter comparison
- [Complete example](complete_example.md): a run with AVR + `GOV_SMOOTH` and every field written out
- [Parameter reference §4](parameter_reference.md): the `DynModelConfig` field tables
