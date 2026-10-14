# Running a case

How to configure **RunConfig**, call **run_case!**, and interpret the results folder — from a single scenario through d_tol sweeps, MATPOWER import, and dual exports.

For every struct field and CSV column, see [Parameter reference](parameter_reference.md). For steady-state OPF equations and builder toggles in context, see [Steady-state OPF](model/02_steady_state_opf.md).

---

## How runs are configured

```
run script (runs/run_9bus.jl, runs/run_39bus.jl)
    │
    ▼
RunConfig(...)                    ← orchestration + two nested avenues
    ├── dispatch::DispatchConfig     steady-state OPF (type_model, builder toggles)
    └── transient::TransientConfig   simulation, TS builder, dyn_model (when trans_stab=true)
    │
    ├── load_system(cfg, path_main)                    reads INPUT_FILES/<case>/
    └── run_case!(cfg, sys, path_main, path_results)   build → solve → save RESULTS/
```

Both path arguments of `run_case!` are **positional and required** — `run_case!(cfg, sys)` is a `MethodError`. `load_system` is the lenient one: `path_main` defaults to `project_root()`.

**Copy helpers:**

- `reconfigure(cfg; field=value)` — shallow override on `RunConfig` (deep-copies `dispatch` / `transient` by default)
- `reconfigure_transient(tc; simulation=…, builder=…, dyn_model=…)` — override nested transient fields (δ_tol sweeps)

---


## 1. Example configurations

### 1.1 Plain ACOPF (no TS)

```julia
cfg = RunConfig(
    trans_stab = false,
    case = "9bus",
    load_factor = 1.5,
    solver_name = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF"),
)
```

### 1.2 TSC-ACOPF — Kron, SC contingency 2

```julia
cfg = RunConfig(
    trans_stab = true,
    case = "9bus",
    load_factor = 1.5,
    solver_name = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF", use_matrix = true),
    transient = TransientConfig(
        simulation = TsSimulationConfig(δ_tol_deg = 100.0),
        dyn_model = DynModelConfig(
            fault = FaultConfig(contingency_id = 2),
        ),
    ),
)
```

### 1.3 TSC-ACOPF — `FULL_BUS`, SC

```julia
cfg = RunConfig(
    trans_stab = true,
    case = "9bus",
    load_factor = 1.5,
    solver_name = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF"),
    transient = TransientConfig(
        simulation = TsSimulationConfig(δ_tol_deg = 100.0),
        dyn_model = DynModelConfig(
            network_form = FULL_BUS,
            mech_power_mode = USE_PM,
            bound_style_δ = :coi_box,
            zip_load_p = (1.0, 0.0, 0.0),  # (Z, I, P) — active demand, constant impedance
            zip_load_q = (1.0, 0.0, 0.0),  # (Z, I, P) — reactive demand, constant impedance
            fault = FaultConfig(contingency_id = 2),
        ),
    ),
)
```

### 1.4 TSC-ACOPF — mixed SG + GFM (`39bus_gfm`)

Requires `DQ_4TH` + `FULL_BUS` and split CSVs (`gen_dynamic_data_full.csv` for SGs, `gfm_dynamic_data.csv` for GFMs):

```julia
cfg = RunConfig(
    trans_stab = true,
    case = "39bus_gfm",
    load_factor = 1.0,
    solver_name = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF"),
    transient = TransientConfig(
        gen_dynamic_filename = "gen_dynamic_data_full.csv",
        simulation = TsSimulationConfig(δ_tol_deg = 100.0, clearing_time = 0.095, t_step = 0.005),
        dyn_model = DynModelConfig(
            allow_gfm = true,
            gen_order = DQ_4TH,
            network_form = FULL_BUS,
            mech_power_mode = USE_PM,
            bound_style_δ = :coi_box,
            zip_load_p = (1.0, 0.0, 0.0),  # (Z, I, P) — active demand, constant impedance
            zip_load_q = (1.0, 0.0, 0.0),  # (Z, I, P) — reactive demand, constant impedance
            fault = FaultConfig(contingency_id = 57),
        ),
    ),
)
```

The case of the paper has six synchronous machines and four converters (ids 3, 4, 5 and 6),
all with `Imax = 1.2` pu, so the current limiter is
built rather than bypassed. Contingency 57 is a fault at bus 23 cleared by opening branch
23–24; `clearing_time = 0.095` on a 5 ms grid is a 100 ms fault. `runs/run_39bus.jl` holds
the full configuration of the paper's run.

The converter's variable boxes are guard-rails derived from `gfm_dynamic_data.csv` and the
`gfm_*` margins on `TsBoundLimitsConfig`; they follow `bound_encoding` like every other
family, so binding GFM limits appear in the dual export (`dual_gfm_limiter_Id`,
`dual_gfm_droop`, `dual_LB/UB_gfm_Id_tf`, …). To widen or tighten one, set the margin
rather than editing the builder:

```julia
builder = TsBuilderConfig(
    limits = TsBoundLimitsConfig(gfm_I_floor_pu = 2.0, gfm_V_meas_max_pu = 1.8),
)
```

### 1.5 `FULL_BUS` — generator 3 trip

```julia
cfg = RunConfig(
    trans_stab = true,
    case = "9bus",
    solver_name = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF"),
    transient = TransientConfig(
        simulation = TsSimulationConfig(δ_tol_deg = 100.0),
        dyn_model = DynModelConfig(
            network_form = FULL_BUS,
            mech_power_mode = USE_PM,
            bound_style_δ = :coi_box,
            zip_load_p = (1.0, 0.0, 0.0),  # (Z, I, P) — active demand, constant impedance
            zip_load_q = (1.0, 0.0, 0.0),  # (Z, I, P) — reactive demand, constant impedance
            fault = FaultConfig(fault_type = GL, gl_gen_ids = [3]),
        ),
    ),
)
```

### 1.6 Custom TS builder toggles

```julia
transient = TransientConfig(
    builder = TsBuilderConfig(
        bound_δ_tf = false,
        ineq_δ_COI_tf_lower = true,
    ),
    dyn_model = DynModelConfig(…),
)
```

Pre-fault `P`/`Q`/`P_m` init equalities are always built by the model path
(no `TsBuilderConfig` fields). Use `mech_power_mode = USE_PG` when the swing
equation should use `P_g` instead of an independent `P_m` (then `Pm` init is omitted).

### 1.7 δ_tol sweep

```julia
base_cfg = RunConfig(
    trans_stab = true,
    dispatch = DispatchConfig(type_model = "ACOPF"),
    transient = TransientConfig(
        dyn_model = DynModelConfig(fault = FaultConfig(contingency_id = 2)),
    ),
)
sys = load_system(base_cfg, path_main)

for δ_val in 90:1:100
    tc = reconfigure_transient(base_cfg.transient;
        simulation = TsSimulationConfig(δ_tol_deg = Float64(δ_val)))
    cfg = reconfigure(base_cfg; transient = tc)
    run_case!(cfg, sys, path_main, path_folder_results)
end
```

---

## 2. Outputs

Each `run_case!` writes under `RESULTS/Results - <timestamp>/` (unless `overwrite_results=true`,
which uses the flat `RESULTS/` tree with the same conditional rules).

**Folder creation** happens at the start of `run_case!` via `build_results_paths` (`_common/auxiliar_functions.jl`).
**File writes** happen after a successful solve (`OPTIMAL`, `LOCALLY_SOLVED`, `ALMOST_LOCALLY_SOLVED` or `ITERATION_LIMIT`), except
`Inputs/` which is populated immediately from `INPUT_FILES/<case>/`.

| Subfolder | Created when (`mkpath` at run start) | Populated when |
|---|---|---|
| `Inputs/` | always | run start (`Copy_Input_CSVs_To_Results!`) |
| `<script>.jl` (run-folder root, not a subfolder) | `run_script` set | run start — byte-for-byte copy of the `.jl` file that built the `RunConfig`, written before the solve so it survives a failure |
| `run_manifest.toml` (run-folder root) | always | after the solve, before the model is released — the machine-readable description of the run (see below) |
| `Dispatch/CSV/` | always | successful solve (primal dispatch CSV/XLSX/TXT) |
| `Dispatch/model_details.txt`, `Dispatch/model_summary.txt` | always | at model build (objective, variables, constraints). Box limits appear as `ineq_const` rows when `bound_encoding = CONSTRAINT`, or as equivalent ≤-form lines from `bound_manifest` when `bound_encoding = VARIABLE`. |
| `Dispatch/CSV_duals/` | `save_duals=true` | successful solve (JuMP primal dual CSV + `Dispatch_Duals.xlsx`) |
| `Dispatch_WarmStart/`, `Dispatch_WarmStart/CSV/` | `save_warmstart_dispatch=true` and either `FULL_BUS` TSC-ACOPF or Kron TSC-DCOPF | successful pre-solve (before TS assembly); same reports as the standalone formulation |
| `Dispatch_WarmStart/CSV/delta_ref.csv` | `save_warmstart_dispatch=true` and Kron TSC-DCOPF | the Taylor anchor $\delta_{ref}$ the linearised $P_e$ rows expand around |
| `Dispatch_WarmStart/CSV_duals/` | `save_warmstart_dispatch=true` and `save_duals=true` | pre-solve JuMP duals |
| `Dispatch_WarmStart/prefault_coupling_starts.{txt,csv}` | `FULL_BUS` TSC-ACOPF (independent of `save_warmstart_dispatch`) | after TS assembly, before the joint solve: the `start=` values the coupling variables were seeded with |
| `Dispatch_Dual/` | `solve_explicit_dual=true` and `type_model` in `{"ED","DCOPF"}` | successful primal + explicit dual LP solve |
| `Dispatch_Dual/CSV/`, `Dispatch_Dual/CSV_duals/` | same as `Dispatch_Dual/` | explicit dual LP outputs |
| `Bus_Matrices/` | `type_model` in `{"ACOPF","DCOPF"}` and (`save_matrices=true` or `trans_stab=true`) | steady-state: Ybus/Bbus if `save_matrices`; TSC: fault/post-fault Y always when `trans_stab` |
| `Transient_Stability/CSV/` | `trans_stab=true` | successful TSC solve (classical + `FULL_BUS` trajectories) |
| `Transient_Stability/CSV/electrical_reactive_power.csv` | Kron TSC-ACOPF (`expressions[:Qe_tf]`) or `FULL_BUS` (`vars[:Qe_tf]`) | per-generator $Q_e(t)$ [MVAr] |
| `Transient_Stability/CSV/generator_*.csv` | `network_form=FULL_BUS` | per-generator $P_e, Q_e$ |
| `Transient_Stability/CSV/generator_internal_voltage.csv` | classical models (`gen_order=CLASSICAL_2ND`), Kron **and** `FULL_BUS` | pre-fault internal EMF magnitude $E$ per generator |
| `Transient_Stability/CSV/bus_*.csv` | `network_form=FULL_BUS` | nodal $P, Q, V, \theta$ trajectories |
| `Transient_Stability/CSV/total_energy.csv` | `trans_stab=true` | $V_e = V_{ke} + V_{pe}$ per SG, alongside `kinetic_energy.csv` / `potential_energy.csv` |
| `Transient_Stability/CSV/dq_*.csv` | `gen_order=DQ_4TH` | $E_d, E_q, I_d, I_q, T_e, E_{fd}$ and the terminal $V_d, V_q$ trajectories |
| `Transient_Stability/CSV/gfm_*.csv` | `allow_gfm=true` | converter states, `gfm_current_loading.csv`, and `gfm_setpoints.csv` (pre-fault $V_{set}$ and the $t=0$ filter states) |
| `Transient_Stability/CSV_duals/` | `trans_stab=true` and `save_duals=true` | successful TSC solve **with a dual certificate** (dynamic dual CSV/XLSX) |
| `Transient_Stability/Figures/` | `trans_stab=true` and `save_ts_plots=true` | successful TSC solve (trajectory SVG plots) |
| `Transient_Stability/Figures_Duals/` | `trans_stab=true`, `save_duals=true` **and** `save_ts_plots=true` | stability-corridor dual SVGs, one file per generator and bound side, for whichever corridors the run built: `Duals Trans. Stab. Const G<id> <side>.svg` for the δ corridor (COI- or machine-referenced) and `Duals Trans. Stab. Const Domega G<id> <side>.svg` for the Δω corridor |
| `Transient_Stability/CSV/Debug/` | `trans_stab=true` and `save_ts_debug_csv=true` | per-(window, gen, step) diagnostics: value, margin and dual on one row |
| `Optim_Matrices/` | `save_optim_matrices=true` and steady-state dispatch | successful solve (Jacobian/Hessian/gradient COO CSV). Row count depends on `bound_encoding`: `VARIABLE` mode has fewer inequality rows because box limits are not explicit `@constraint`s. |

`Transient_Stability/CSV/Debug/` and `Optim_Matrices/` are **not** in `path_names`: they are
created on demand by their writers (`Save_TS_Debug_CSV!` and
`Compute_and_Save_Optimal_Matrices`), so they appear only when their flag is on and the
solve reached the point of writing them.

**Not created:** `Dispatch/CSV_duals/` when `save_duals=false`; `Bus_Matrices/` for ED/UC or ACOPF/DCOPF with
`save_matrices=false` and `trans_stab=false`; any `Transient_Stability/` subtree when `trans_stab=false`; `Dispatch_Dual/` unless
`solve_explicit_dual=true` on a DC-OPF case.

**Duals and non-converged solves.** `run_case!` accepts `ITERATION_LIMIT` and
`ALMOST_LOCALLY_SOLVED` as "solved enough" to export the primal trajectories, but a solve
that stopped early carries no dual certificate. Both dual writers therefore require
`dual_status ∈ {FEASIBLE_POINT, NEARLY_FEASIBLE_POINT}` and skip the export with a warning
otherwise. Note this is stricter than `JuMP.has_duals`, which only asks for
`dual_status != NO_SOLUTION` and so accepts the `UNKNOWN_RESULT_STATUS` that Ipopt reports
after hitting `max_iter` — multipliers that exist as iterates and certify nothing. When
the export does run, `duals.txt` and `dynamic_model_duals.txt` open with the termination /
primal / dual status of the solve that produced them.

**Known limit.** The Kron-linear path (TSC-DCOPF) has no reactive-power variable or
expression, so `electrical_reactive_power.csv` is absent there by construction — not a
missing export.

### 2.1 `run_manifest.toml` — reading a run back with a script

`input_parameters.txt` is written for people, and it does not record `load_factor`;
`Inputs/bus_data.csv` archives the demand *before* that scaling is applied
(`load_system` scales after the copy). A `run_script` copy only contains the fields
someone typed. So none of the three lets a script recover what was actually solved.

`run_manifest.toml` does, and it is written on **every** run — failed ones included,
where `[status]` records why the folder is otherwise empty. Read it with any TOML
parser (`TOML.parsefile` in Julia, stdlib `tomllib` in Python):

| Table | Contents |
|---|---|
| `[run]` | every `RunConfig` scalar, `load_factor` included, plus `[run.ipopt]` / `[run.highs]` / … solver subtables |
| `[dispatch]` | `DispatchConfig`, enums as strings |
| `[transient.simulation]`, `[transient.builder]`, `[transient.dyn_model]`, `[transient.dyn_model.fault]` | the transient config tree |
| `[resolved]` | what only the built model knows: `delta_ref_gen` (the machine `:highest_H` picked), `delta_pair_count` / `delta_pairs` (the unordered pairs a `:pairwise` corridor covered), `sg_ids` / `gfm_ids`, `H` / `D` / `Xd_tr`, `delta_tol` / `delta_pair_tol`, `t_step`, `n_steps_fault` / `n_steps_postfault` / `n_steps_total`, and the cost curve `c0` / `c1` / `c2` per generator |
| `[exports]` | the dual CSV basenames this run actually wrote, so a reader does not have to probe the filesystem |
| `[status]` | `termination_status`, `objective_MVA`, `t_build`, `t_solve` |

The config tables are reflected off the structs with `fieldnames`, so a field added
to `RunConfig` or `DynModelConfig` appears in the manifest without touching the
writer. Greek field names are transliterated on the way out (`δ_tol` → `delta_tol`,
`Δω_tol` → `Delta_omega_tol`) so attribute access on the reading side stays sane.
Fields with no serialisable value are omitted rather than written as a placeholder:
`nothing`, and callbacks such as `post_solve_hook`.

### 2.2 Export flags

There is no master "save everything" switch. What a run writes is decided by the flags
below plus the model-shape gates (`network_form`, `gen_order`, `include_avr`,
`include_governor`, `allow_gfm`), which decide which quantities exist at all.

| Flag | Struct | Default | Controls |
|---|---|---|---|
| `save_duals` | `RunConfig` | `true` | `Dispatch/CSV_duals/`, `duals.txt`, `Dispatch_Duals.xlsx`, and the whole `Transient_Stability/CSV_duals/` + `OPF_Duals_Results.xlsx` layer |
| `save_matrices` | `RunConfig` | `true` | `Bus_Matrices/` (Ybus / Bbus dumps) |
| `save_ts_plots` | `RunConfig` | `false` | `Transient_Stability/Figures/` and `Figures_Duals/`; requires the Plots extension |
| `save_ts_debug_csv` | `RunConfig` | `false` | `Transient_Stability/CSV/Debug/` |
| `save_optim_matrices` | `RunConfig` | `false` | `Optim_Matrices/` (steady-state only; forced off for TSC runs) |
| `save_warmstart_dispatch` | `RunConfig` | `false` | `Dispatch_WarmStart/` (`FULL_BUS` TSC-ACOPF only) |
| `run_script` | `RunConfig` | `nothing` | Copy of the configuring `.jl` file in the run-folder root; set it to `@__FILE__` in the run script |
| `solve_explicit_dual` | `DispatchConfig` | `false` | `Dispatch_Dual/` (ED / DC-OPF only) |

`TransientConfig`, `TsBuilderConfig` and `DynModelConfig` carry **no** export flags —
their toggles change the model, and therefore which files have content to hold, but never
whether an existing quantity is written.

**MATPOWER import** (`Import_Matpower_Case`) creates only `Results - <ts>/Inputs/` — see [Parameter reference §6.1](parameter_reference.md).

### 2.3 Solver logs

Each run writes optimiser text logs at the **timestamped run root** (alongside `input_parameters.txt`):

| File | When written |
|---|---|
| `solver_log.txt` | Main solve (steady-state OPF or full TSC problem) |
| `solver_log_warmstart.txt` | Optional pre-solves only: DCOPF steady-state before TSC-DCOPF, or `FULL_BUS` ACOPF warm start before TSC assembly |

Logs are configured via `set_solver_log_path!` (`_manage_inputs/functions_2_setup_optim.jl`):
Ipopt/MadNLP/PARDISO use `output_file`; Gurobi uses `LogFile`; HiGHS uses `log_file`.

### 2.3.1 PARDISO on Windows (optional)

For large ACOPF / TSC-ACOPF KKT systems, you can use Intel MKL PARDISO as Ipopt's
linear solver instead of the default MUMPS:

```julia
cfg = RunConfig(
    solver_name = "Ipopt-pardiso",
    ipopt = IpoptSolverConfig(
        pardiso_lib_path = raw"C:\Libraries\Pardiso\lib\libpardiso.dll",
    ),
    dispatch = DispatchConfig(type_model = "ACOPF"),
    # ...
)
```

Requirements:

- A licensed `libpardiso.dll` (or platform equivalent).
- Set `IpoptSolverConfig.pardiso_lib_path` **or** `ENV["JULIA_PARDISO_LIB"]` before the run.
- By default `pardiso_license_message = true` sets `ENV["PARDISOLICMESSAGE"] = "1"`.

`"Ipopt-pardiso"` is an Ipopt backend (`is_ipopt_backend_solver`); all other `IpoptSolverConfig`
fields (tolerances, `max_iter`, Hessian mode) apply as for `solver_name = "Ipopt"`.

### 2.3.2 MadNLP + HSL (optional)

MadNLP can use HSL via `MadNLPHSL`. `solver_name` is namespaced so NLP solver and
linear engine stay explicit (`"Ipopt-ma57"` vs `"MadNLP-ma57"`). Pardiso remains
Ipopt-only (`"Ipopt-pardiso"`).

| `solver_name` | Linear solver | Extra packages |
|---|---|---|
| `"MadNLP"` | MUMPS (default) | `MadNLP` |
| `"MadNLP-ma57"` | HSL MA57 | `MadNLP`, `MadNLPHSL`, `HSL_jll` |
| `"MadNLP-ma97"` | HSL MA97 | `MadNLP`, `MadNLPHSL`, `HSL_jll` |

```julia
using TSCOPF, MadNLP, MadNLPHSL

cfg = RunConfig(
    solver_name = "MadNLP-ma57",
    madnlp = MadNLPSolverConfig(tol = 1e-8),
    dispatch = DispatchConfig(type_model = "ACOPF"),
)
```

Notes:

- MadNLP does **not** use Ipopt's `hsllib` attribute.
- HSL: install/configure `HSL_jll` as for Ipopt HSL, then `Pkg.add("MadNLPHSL")`.
  If you downloaded a licensed `HSL_jll` tree (e.g. from the HSL portal),
  **develop that path in every Julia environment** that should use MA57/MA97 —
  including the simulations project. Example:

  ```julia
  using Pkg
  Pkg.develop(path = raw"C:\path\to\HSL_jll.jl.v2025.7.21")   # wherever you unpacked it
  Pkg.add("MadNLPHSL")
  ```

  A stock registry `HSL_jll` (e.g. v4.0.6) may provide a `libhsl.dll` that loads
  but lacks full MA57 symbols; MadNLP then fails at factorize with missing
  `MA57_version`. Ipopt-only workflows often already have the licensed package
  developed in `~/.julia/environments/v1.x`, which is why Ipopt+ma57 can work
  without an explicit license environment variable.
- Optional MA97 threads: `MadNLPSolverConfig(ma97_num_threads = 8)`.
- Further MA57/MA97 knobs go in `madnlp.raw_options` (e.g. `"ma57_pivtol" => 1e-8`).

!!! note "Windows file locks"
    Ipopt keeps `output_file` open until the native `IpoptProblem` is freed.
    `run_case!` calls `release_solver_backend!(model)` **after** all result exports:
    finalize any Ipopt inner problem, `empty!(model)`, then `GC.gc()`. You can move
    or archive `solver_log*.txt` without closing Julia. Do **not** call
    `release_solver_backend!` yourself while you still need `JuMP.value` / `JuMP.dual`
    on that model. The returned `dyn_model_dict` variable references are stale after
    release — use the saved CSV/TXT outputs. In batch loops, discard the result tuple
    between iterations so nothing pins the emptied model.

### 2.4 Dual variable exports

!!! warning "Sign convention — the easiest place to get an economic reading backwards"
    Inequalities are coded as `(LHS − RHS) ≤ 0`, so `JuMP.dual()` on an *active*
    inequality constraint is **non-positive** by construction — a negative dual on
    a `P_g ≤ P_g^max` constraint is the expected sign for a binding capacity
    limit, not an error. Power-balance LMPs use **π_k = +λ_k** on the primal JuMP
    duals: the balance is coded `P_g − P_d − Σflows == 0` and the Lagrangian is
    `f − Σλ(LHS − RHS)`, so the two minus signs cancel. The **explicit dual LP**
    under `Dispatch_Dual/` is a different object with the opposite convention,
    `π_k = −λ_k`. When exporting a
    new dual family, cross-check the sign on a small case (case9, contingency 2)
    before trusting the CSV: flipping this once and propagating it through a
    stability-adjusted-LMP figure is a quiet way to publish a wrong-signed rent.

**Steady-state — `Dispatch/Dispatch_Duals.xlsx`**

Registry-driven by `_manage_outputs/DispatchDualRegistry.jl` (`STEADY_STATE_DUAL_SPECS`):
one generic exporter writes whichever constraint families exist in the solved model, so
ED, DC-OPF and AC-OPF share the same sheet/column naming. Sheets present depend on the
model and on the `DispatchConfig` toggles (`ineq_sbranch_upper`, `ineq_ang_diff_branch`,
`bound_*`, `ineq_sg_upper`).

| Sheet | Constraint | Notation / id column |
|---|---|---|
| `P_Balance` / `Q_Balance` | Active / reactive power balance | λ_k (LMP = +λ_k), Bus_ID |
| `Pd_Pin` | Demand pin, `P_d[k] - p_d⁰[k] == 0` (`pd_as_variable`, ACOPF) | σ_k = ∂L/∂P_d(k), Bus_ID |
| `Sg_Upper` | Generator capability curve | Gen_ID |
| `Sik_Upper` / `Ski_Upper` | Branch thermal limits (i→k / k→i) | Branch_ID |
| `Ang_Diff_Lo` / `Ang_Diff_Up` | Angle difference limits | ρ_km⁻ / ρ_km⁺, Branch_ID |
| `V_Mag_Lo` / `V_Mag_Up` | Bus voltage magnitude limits | Bus_ID |
| `V_Ang_Lo` / `V_Ang_Up` | Bus angle limits | μ_k⁻ / μ_k⁺, Bus_ID |
| `Pg_Lo` / `Pg_Up` | Generator P limits | η_g⁻ / η_g⁺, Gen_ID |
| `Qg_Lo` / `Qg_Up` | Generator Q limits | Gen_ID |
| `Eq_Pik` / `Eq_Qik` / `Eq_Pki` / `Eq_Qki` | Branch flow definitions (explicit-flow models) | Branch_ID |
| `Pik_Lo`/`Up`, `Qik_Lo`/`Up`, `Pki_Lo`/`Up`, `Qki_Lo`/`Up` | Branch flow bounds (when `bound_P_ik` … enabled) | Branch_ID |
| `GFM_Imax` / `GFM_Emax` | Grid-forming converter capability at the dispatch point (`allow_gfm`, ACOPF) | Gen_ID |

CSV mirrors live in `Dispatch/CSV_duals/` (one file per family, e.g. `dual_P_balance.csv`,
`dual_Ski_Upper.csv`, `dual_diff_ang_Upper.csv`). Column headers match the table above; the
generator-bound CSVs use `Gen_ID` and the `k→i` / angle-upper families carry their own
`Dual_Ski_Upper` / `Dual_diff_ang_Upper` columns.

**Transient — `Transient_Stability/OPF_Duals_Results.xlsx`**

Each sheet is named after its registry export with the `dual_` prefix dropped and `_Duals`
appended (`dual_δ_COI_upper` → `δ_COI_upper_Duals`). The names keep the Greek letters.

| Sheet | Constraint family |
|---|---|
| `Pe_init_Duals` | Initial P–δ–E link |
| `δCOI_Duals` | COI angle definition |
| `Pe_Duals` | Electrical power during transient |
| `δ_Duals` | Rotor angle (trapezoidal) |
| `Δω_Duals` | Speed deviation (swing) |
| `δ_COI_lower_Duals` / `δ_COI_upper_Duals` | δ w.r.t. COI stability bounds |
| `δ_ref_lower_Duals` / `δ_ref_upper_Duals` | δ w.r.t. a reference machine (`:highest_H`, `:ref_gen`) |
| `δ_pair_lower_Duals` / `δ_pair_upper_Duals` | pairwise δ separation (`:pairwise`) |

The CSV files under `Transient_Stability/CSV_duals/` use ASCII names instead
(`dual_delta_COI_upper.csv`, `dual_Delta_Omega.csv`, …), except the box families below.

Optional box families (`TsBuilderConfig.bound_*`) export under their own names on both
windows — `dual_LB_δ_tf` / `dual_UB_δ_tf`, and the same pattern for `Δω`, `Pe`, `Qe`,
`δCOI`, and the dq states `Ed`, `Eq`, `Id`, `Iq`, `Te`. With `allow_gfm` the converter
families follow: the measurement filters (`filter_P`, `filter_Q`, `filter_V`), the δ
integration step, the raw and clipped Q–V PI states, the current limiter, `Pe`/`Qe`, and
the pre-fault `δ` / `P_set` boxes.

Every one of these is available on **GL and OB** disturbances too. Those build a single
window, and the exported series covers that one window.

Authoritative mapping: `_transient_stability/DynDualRegistry.jl` (the three `*_DUAL_SPECS`
constants).

**Id columns in the TS dual CSVs.** Every file under `Transient_Stability/CSV_duals/`
names what its rows and columns are:

- per-generator trajectories use `Gen_<id>` headers, per-**bus** ones use `Bus_<id>`.
  The nodal families are `dual_Pbalance.csv`, `dual_Qbalance.csv`, `dual_V_lower.csv`, ….
- the one-value-per-generator files (`dual_Pe_init.csv`, `dual_Pm_init.csv`, the
  `dual_LB_*` / `dual_UB_*` boxes) carry a leading `Gen_ID` column, and the merged
  time series (`dual_delta_COI.csv`, …) carry a running `Index`. Select rows by these id
  columns, not by position: positions do not line up between families. On a mixed fleet
  `dual_Pe_init.csv` spans SG **and** converters while `dual_Pe.csv` spans the machines alone.

---

## 3. Validation before solve

`run_case!` calls, in order:

1. `Check_Coherence_Input_Data` — model/solver pairing  
2. `validate_run_config!` — dispatch + transient consistency; calls `validate_dispatch_config!` and, when `trans_stab=true`, `validate_dyn_config!` (`DynModelConfig` rules)  

When `trans_stab=true`, then:

3. `validate_fault_config!` — disturbance vs network, for the main fault and each of `transient.extra_contingency_ids`  
4. `validate_dyn_data!` — the machine CSV has every column the selected model needs  
5. `validate_δ_reference!` — the reference machine of a `:highest_H` / `:ref_gen` δ corridor against the actual generator data, and the two-survivor requirement shared with `:pairwise`  
6. `validate_all_gfm_fleet!` — the extra rules for a converter-only fleet (no-op otherwise)  
7. `validate_gfm_dispatch_band!` — each converter's `DGEN` P/Q interval sits inside its `DGFM` band (mixed fleets only)

---

## 4. Quick lookup — where do I change X?

| I want to change… | Edit… |
|---|---|
| Case, load level, solver, I/O flags | `RunConfig` in the run script |
| ACOPF vs DCOPF, fuel cost, OPF toggles | `RunConfig.dispatch` (`DispatchConfig`) |
| DC branch susceptance (`1/x` vs PowerModels) | `RunConfig.dispatch.susceptance_model` |
| δ_tol, simulation timing, f_syn | `RunConfig.transient.simulation` |
| Which TS constraints are built | `RunConfig.transient.builder` (`TsBuilderConfig`) |
| Kron vs `FULL_BUS`, ZIP loads, P_m mode | `RunConfig.transient.dyn_model` |
| AVR / governor on, valve limiter | `RunConfig.transient.dyn_model` — physics in [Part I §8](model/08_controls_avr_governor.md) |
| Grid-forming fleet on, converter data file | `dyn_model.allow_gfm`, `transient.gfm_dynamic_filename` — physics in [Part I §9](model/09_grid_forming.md) |
| SC vs GL vs OB, which gen/load/branch trips | `RunConfig.transient.dyn_model.fault` |
| Which bus/line faults (SC) | `contingencies.csv` + `fault.contingency_id` |
| Network topology, costs, limits | CSV files in `INPUT_FILES/<case>/` |
| Several loads on one bus, each with its own ZIP split | `load_data.csv` in the case folder ([parameter reference §6.2](parameter_reference.md#62-load_datacsv--several-loads-on-one-bus-each-with-its-own-zip-split)) |
| Use a MATPOWER `.m` directly (no CSV conversion) | `RunConfig.matpower_file = "case.m"` ([Parameter reference §6.1.1](parameter_reference.md)) |
| Convert a MATPOWER `.m` to in-house CSVs | `Import_Matpower_Case(path_m)` ([Parameter reference §6.1.2](parameter_reference.md)) |
| Everything at once, explicitly | [Complete example](complete_example.md) — every field of every struct |
