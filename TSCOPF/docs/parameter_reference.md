# Parameter reference

Field-by-field reference for **RunConfig**, **DispatchConfig**, **TransientConfig**, **DynModelConfig**, **FaultConfig**, and network CSV inputs: what each knob does, where to set it, and valid combinations.

- [Running a case](running_a_case.md) — workflow, examples, outputs, and validation.
- [Configuration map](configuration_map.md) — developer-oriented struct field map.

---

## 1. `RunConfig` — orchestration

Edit in **`main.jl`** (single scenario) or **`main_loop.jl`** (δ_tol sweep).  
Defaults: `engine.jl` (`Base.@kwdef struct RunConfig`).

These fields apply to **both** plain OPF and TSC runs.

| Field | Type | Default | Description |
|---|---|---|---|
| `trans_stab` | `Bool` | `false` | `true` → attach TS constraints; requires `transient` |
| `case` | `String` | `"9bus"` | Subfolder under `INPUT_FILES/` |
| `base_MVA` | `Float64` | `100.0` | System power base [MVA] |
| `load_factor` | `Float64` | `1.5` | Multiplier on `p_d` and `q_d` after read |
| `solver_name` | `String` | `"Ipopt"` | `"Ipopt"`, `"Ipopt-ma57"`, `"Ipopt-ma97"`, `"Ipopt-pardiso"`, `"MadNLP"`, `"MadNLP-ma57"`, `"MadNLP-ma97"`, `"Gurobi"`, `"HiGHS"` |
| `silent_solver` | `Bool` | `false` | Suppress optimiser console output (`print_level = 0` for Ipopt) |
| `time_limit_sec` | `Float64` | `600.0` | JuMP solver time limit [s] |
| `ipopt` | `IpoptSolverConfig` | defaults below | Ipopt / HSL NLP options — see §1.1 |
| `highs` | `HiGHSSolverConfig` | defaults below | HiGHS LP/QP options — see §1.2 |
| `gurobi` | `GurobiSolverConfig` | defaults below | Gurobi LP/QP/MILP options — see §1.3 |
| `madnlp` | `MadNLPSolverConfig` | defaults below | MadNLP NLP options — see §1.4 |
| `overwrite_results` | `Bool` | `false` | `false` → timestamped `RESULTS/` folder |
| `save_duals` | `Bool` | `true` | Export JuMP duals from the **primal** solve after success. Creates `Dispatch/CSV_duals/` at run start (all formulations: ED, DCOPF, ACOPF). When `trans_stab=true`, also creates `Transient_Stability/CSV_duals/`. Independent of `solve_explicit_dual`. A solve with no dual certificate (`dual_status` outside `{FEASIBLE_POINT, NEARLY_FEASIBLE_POINT}` — e.g. the `UNKNOWN_RESULT_STATUS` Ipopt reports on `ITERATION_LIMIT`) is skipped with a warning instead of exporting numbers that are not prices; when the export runs, the solve status is stamped at the top of `duals.txt` / `dynamic_model_duals.txt`. |
| `save_optim_matrices` | `Bool` | `false` | Export Jacobian, Hessian, Lagrangian gradient as sparse COO CSV. **Steady-state only** — auto-forced `false` (with a warning) for TSC runs, where these matrices are huge. |
| `save_matrices` | `Bool` | `true` | Steady-state only: write Ybus/Bbus XLSX to `Bus_Matrices/` when `type_model` is ACOPF or DCOPF. Ignored for ED and UC. TSC fault/post-fault matrices are always exported when `trans_stab=true` (folder created when `trans_stab` or steady-state `save_matrices` on ACOPF/DCOPF). |
| `save_ts_plots` | `Bool` | `false` | TSC only: trajectory SVG figures under `Transient_Stability/Figures/`, plus the stability-corridor dual figures under `Figures_Duals/` (δ and Δω, whichever styles the run built) when `save_duals` is also on (opt-in; uses Plots.jl). CSV/TXT transient results are still saved when `trans_stab=true`. On a headless machine `load_plots_extension!` sets `GKSwstype=100` so GR writes to file instead of failing to open a window. |
| `save_ts_debug_csv` | `Bool` | `false` | TSC only: per-(window, generator, step) diagnostic dumps under `Transient_Stability/CSV/Debug/` — value, distance to the limit and dual on one row. Writes `gfm_filter_debug.csv`, `gfm_limiter_debug.csv`, `gfm_voltage_debug.csv` (mixed-fleet runs) and `swing_debug.csv`. Opt-in: the filter file alone is three rows per converter per step. **`swing_debug.csv` uses the Greek variable names of the model (`δ_curr`, `Δω_step`, …) in its header.** The file is UTF-8 without a BOM, so opening it by double-click in Excel decodes it as ANSI and mangles those names; import it instead (Data → From Text/CSV → 65001 Unicode UTF-8). Every other exported CSV uses ASCII headers. |
| `save_warmstart_dispatch` | `Bool` | `false` | TSC runs that pre-solve a steady-state OPF before assembling the dynamics: writes that solution to `Dispatch_WarmStart/` (same reports as `Dispatch/`; dual export follows `save_duals`). Two such paths — **FULL_BUS TSC-ACOPF** (the mandatory ACOPF warm start) and **Kron TSC-DCOPF** (the DC solve that fixes the Taylor anchor, which also writes `CSV/delta_ref.csv`). Throws on Kron TSC-ACOPF, which solves once jointly and has no pre-solve, and on steady-state runs. |
| `run_script` | `String` or `nothing` | `nothing` | Archive the script that configured the run. Set it to `@__FILE__` inside the run script and `run_case!` copies that `.jl` file, byte for byte, into the run-folder root next to `input_parameters.txt` and `solver_log.txt` — so a `RESULTS/` folder carries the exact configuration that produced it, comments included. Any readable path is accepted (a sweep driver can archive itself); a path that is not an existing file throws in `validate_run_config!`, before any model is built. The copy is written right after the folder is created, so it survives a failed or iteration-limited solve. Sweeps that call `run_case!` per `reconfigure`d iteration get one copy per timestamped folder. |
| `post_solve_hook` | `Function` or `nothing` | `nothing` | Called once, with the model still assembled, immediately before `release_solver_backend!` invalidates every `VariableRef`/`ConstraintRef` — the last point at which a caller can read duals or variable values off the live model. Receives **keyword arguments only** (`cfg`, `sys`, `path_names`, `model`, `opf_dict`, `dyn_model_dict`, `dyn_parameters_dict`, `status`, `status_solved`, `obj_MVA`, `t_build`, `t_solve`), so declare it as `f(; kwargs...)` and pick out what you need; more data may be handed over in future versions without breaking a hook written today. Runs on failed solves too — gate on `status_solved` if that matters. Exceptions propagate: an export that failed silently is worse than a run that stopped and said so. |
| `dispatch` | `DispatchConfig` | `DispatchConfig()` | Avenue 1 — see §2 |
| `transient` | `TransientConfig` or `nothing` | `nothing` | Avenue 2 — see §3; **required** when `trans_stab=true` |

### Solver compatibility (enforced automatically)

Checked via `Check_Coherence_Input_Data(trans_stab, dispatch.type_model, solver_name)`.

| `dispatch.type_model` | Allowed `solver_name` | Notes |
|---|---|---|
| `"ACOPF"` | Ipopt, MadNLP (+ HSL aliases); Pardiso via `"Ipopt-pardiso"` only | Nonlinear AC — **not** Gurobi/HiGHS |
| `"DCOPF"` | Gurobi, HiGHS | Linear DC OPF |
| TSC-ACOPF | Ipopt / MadNLP (+ `Ipopt-*` / `MadNLP-*` HSL); `"Ipopt-pardiso"` | `trans_stab=true` + `"ACOPF"` |
| TSC-DCOPF | Gurobi, HiGHS | `trans_stab=true` + `"DCOPF"` |

### 1.1 `IpoptSolverConfig` — NLP solver knobs

Set on `RunConfig.ipopt` (or pass to `Setup_Optim_Model`). Used when
`solver_name` is `"Ipopt"`, an Ipopt+HSL linear solver (`"Ipopt-ma57"`, `"Ipopt-ma97"`), or
`"Ipopt-pardiso"`. Replaces the legacy `config_file.jl` Ipopt block. Log paths are **not** here:
`run_case!` sets `output_file` to `RESULTS/.../solver_log.txt` at solve time.

| Field | Type | Default | Description |
|---|---|---|---|
| `tol` | `Float64` | `1e-8` | Ipopt overall convergence tolerance |
| `print_level` | `Int` | `5` | Verbosity 0–12; overridden to `0` when `silent_solver=true` |
| `max_iter` | `Int` | `5000` | Maximum Ipopt iterations |
| `constr_viol_tol` | `Float64` | `1e-8` | Constraint violation tolerance |
| `dual_inf_tol` | `Float64` | `1e-8` | Dual infeasibility tolerance |
| `compl_inf_tol` | `Float64` | `1e-8` | Complementarity tolerance |
| `hessian_approximation` | `String` | `"exact"` | `"exact"` or `"limited-memory"` (L-BFGS) |
| `mu_strategy` | `String` or `nothing` | `nothing` | `"monotone"`, `"adaptive"`, or `nothing`. `nothing` leaves Ipopt's monotone barrier |
| `mu_init` | `Float64` or `nothing` | `nothing` | Positive initial barrier parameter, or `nothing` to leave Ipopt's 0.1 — except when `gfm_in_coi` and `constrain_Δω` are both on, where `run_case!` uses `1e-3`. The converter reactive nameplate makes the default barrier cycle in restoration on that model even though the optimum sits inside the nameplate |
| `limited_memory_max_history` | `Int` | `50` | L-BFGS history (only when `hessian_approximation = "limited-memory"`) |
| `pardiso_lib_path` | `String` or `nothing` | `nothing` | Path to `libpardiso.dll` (or `.so`); falls back to `ENV["JULIA_PARDISO_LIB"]` when `solver_name = "Ipopt-pardiso"` |
| `pardiso_license_message` | `Bool` | `true` | When `solver_name = "Ipopt-pardiso"`, sets `ENV["PARDISOLICMESSAGE"] = "1"` before model build |

Example — PARDISO linear solver (Windows):

```julia
cfg = RunConfig(
    solver_name = "Ipopt-pardiso",
    ipopt = IpoptSolverConfig(
        pardiso_lib_path = raw"C:\Libraries\Pardiso\lib\libpardiso.dll",
    ),
    dispatch = DispatchConfig(type_model = "ACOPF"),
)
```

Example — limited-memory Hessian:

```julia
cfg = RunConfig(
    solver_name = "Ipopt",
    ipopt = IpoptSolverConfig(
        hessian_approximation = "limited-memory",
        limited_memory_max_history = 50,
    ),
)
```

### 1.2 `HiGHSSolverConfig` — LP/QP solver knobs

Set on `RunConfig.highs`. Used when `solver_name = "HiGHS"` (DCOPF, TSC-DCOPF, explicit dual LP when Gurobi is unavailable). Log path via `log_file` at solve time.

| Field | Type | Default | HiGHS option | Description |
|---|---|---|---|---|
| `output_flag` | `Bool` | `true` | `output_flag` | Console logging; overridden by `silent_solver` |
| `primal_feasibility_tolerance` | `Float64` | `1e-8` | `primal_feasibility_tolerance` | Primal feasibility tolerance |
| `dual_feasibility_tolerance` | `Float64` | `1e-8` | `dual_feasibility_tolerance` | Dual feasibility tolerance |
| `ipm_optimality_tolerance` | `Float64` | `1e-8` | `ipm_optimality_tolerance` | IPM optimality tolerance |
| `simplex_iteration_limit` | `Int` | `5000` | `simplex_iteration_limit` | Simplex iteration cap |
| `ipm_iteration_limit` | `Int` | `5000` | `ipm_iteration_limit` | IPM iteration cap |
| `solver` | `String` | `"choose"` | `solver` | `"choose"`, `"simplex"`, `"ipm"`, or `"pdlp"` |
| `raw_options` | `Dict{String,Any}` | `{}` | — | Extra HiGHS MOI attributes (e.g. `"time_limit" => 120.0`) |

Example:

```julia
cfg = RunConfig(
    solver_name = "HiGHS",
    highs = HiGHSSolverConfig(
        solver = "ipm",
        ipm_iteration_limit = 10_000,
    ),
)
```

### 1.3 `GurobiSolverConfig` — LP/QP/MILP/NLP solver knobs

Set on `RunConfig.gurobi`. Used when `solver_name = "Gurobi"` (DCOPF, UC, explicit dual LP; nonlinear when licensed). Log path via `LogFile` at solve time.

| Field | Type | Default | Gurobi parameter | Description |
|---|---|---|---|---|
| `output_flag` | `Int` | `1` | `OutputFlag` | Console output (`0`/`1`) |
| `feasibility_tol` | `Float64` | `1e-8` | `FeasibilityTol` | Primal feasibility tolerance |
| `optimality_tol` | `Float64` | `1e-8` | `OptimalityTol` | Dual feasibility tolerance |
| `bar_iter_limit` | `Int` | `5000` | `BarIterLimit` | Barrier iteration limit (LP/QP) |
| `mip_gap` | `Float64` | `1e-8` | `MIPGap` | MIP relative gap (UC MILP) |
| `nl_bar_iter_limit` | `Int` | `5000` | `NLBarIterLimit` | NL barrier iteration limit (Gurobi 13+) |
| `nl_bar_p_feas_tol` | `Float64` | `1e-8` | `NLBarPFeasTol` | NL barrier primal residual tolerance |
| `nl_bar_d_feas_tol` | `Float64` | `1e-8` | `NLBarDFeasTol` | NL barrier dual residual tolerance |
| `nl_bar_c_feas_tol` | `Float64` | `1e-8` | `NLBarCFeasTol` | NL barrier complementarity tolerance |
| `optimality_target` | `Int` | `1` | `OptimalityTarget` | `1` selects NL barrier for continuous nonconvex (Gurobi 13+) |
| `raw_options` | `Dict{String,Any}` | `{}` | — | Extra Gurobi parameters (e.g. `"Threads" => 4`) |

Example:

```julia
cfg = RunConfig(
    solver_name = "Gurobi",
    gurobi = GurobiSolverConfig(
        mip_gap = 1e-6,
        nl_bar_iter_limit = 10_000,
        optimality_target = 1,
    ),
)
```

### 1.4 `MadNLPSolverConfig` — NLP solver knobs

Set on `RunConfig.madnlp`. Used when `solver_name` is `"MadNLP"` or a MadNLP
HSL alias (`"MadNLP-ma57"`, `"MadNLP-ma97"`). Requires `using MadNLP` after
`using TSCOPF` (package extension). HSL aliases also need `MadNLPHSL` + `HSL_jll`
in the environment — see [Running a case](running_a_case.md). Pardiso is
Ipopt-only (`"Ipopt-pardiso"`).

| Field | Type | Default | Description |
|---|---|---|---|
| `tol` | `Float64` | `1e-8` | KKT residual tolerance |
| `max_iter` | `Int` | `3000` | Maximum IPM iterations |
| `print_level` | `Int` | `3` | MadNLP `LogLevels` (3 = INFO); quiet when `silent_solver=true` |
| `hessian_approximation` | `String` | `"exact"` | `"exact"`, `"compact-lbfgs"`, `"bfgs"`, or `"damped-bfgs"` (see note below) |
| `ma97_num_threads` | `Int` or `nothing` | `nothing` | When set, MadNLP `ma97_num_threads` (useful with `"MadNLP-ma97"`) |
| `raw_options` | `Dict{String,Any}` | `{}` | Advanced MadNLP options (e.g. `"acceptable_tol" => 1e-5`, `"ma57_pivtol" => 1e-8`) |

Unlike Ipopt, MadNLP does **not** take an `hsllib` attribute. Select the KKT
linear solver with `solver_name`; library loading is handled by `MadNLPHSL`.

Example — MadNLP with HSL MA57:

```julia
using TSCOPF, MadNLP, MadNLPHSL

cfg = RunConfig(
    solver_name = "MadNLP-ma57",
    madnlp = MadNLPSolverConfig(
        tol = 1e-8,
        raw_options = Dict("ma57_pivtol" => 1e-8),
    ),
)
```

`hessian_approximation = "exact"` uses analytic second derivatives. Of the
quasi-Newton modes, only `"compact-lbfgs"` works with the sparse KKT system
MadNLP builds for JuMP models; `"bfgs"` and `"damped-bfgs"` are dense and
additionally require a dense KKT system (practical on small cases only).
Quasi-Newton tuning travels through `raw_options`.

Example — MadNLP with MA97 threads and a limited-memory Hessian:

```julia
using TSCOPF, MadNLP

cfg = RunConfig(
    solver_name = "MadNLP-ma97",
    madnlp = MadNLPSolverConfig(
        hessian_approximation = "compact-lbfgs",
        ma97_num_threads = 8,
        raw_options = Dict(
            "acceptable_tol" => 1e-5,
            "quasi_newton_options" => MadNLP.QuasiNewtonOptions(max_history = 20),
        ),
    ),
)
```

The dense variants need the KKT system switched over as well:

```julia
madnlp = MadNLPSolverConfig(
    hessian_approximation = "bfgs",          # or "damped-bfgs"
    raw_options = Dict(
        "callback"      => MadNLP.DenseCallback,
        "kkt_system"    => MadNLP.DenseCondensedKKTSystem,
        "linear_solver" => MadNLP.LapackCPUSolver,
    ),
)
```

---

## 2. `DispatchConfig` — steady-state OPF (avenue 1)

Defined in `_common/DispatchConfig.jl`. Set on `RunConfig.dispatch`.

### 2.1 Core selection

| Field | Type | Default | Description |
|---|---|---|---|
| `type_model` | `String` | `"ACOPF"` | `"ACOPF"`, `"DCOPF"`, `"ED"`, `"UC"` |
| `use_matrix` | `Bool` | `true` | Ybus/Bbus matrix form for power balance |
| `cost_type` | `String` | `"quadratic"` | `"quadratic"` (MATPOWER) or `"linear"` |
| `susceptance_model` | `SusceptanceModel` | `SIMPLE` | **DC-OPF only.** Branch susceptance convention: `SIMPLE` = `1/x` (textbook DC); `POWERMODELS` = `imag(inv(r+jx))` = `x/(r²+x²)` (matches PowerModels). Identical when `r = 0`; diverge otherwise. Ignored by AC/ED. |
| `pd_as_variable` | `Bool` | `false` | **ACOPF only, active power only.** Build the bus demand as a JuMP variable pinned by `P_d[k] - p_d⁰[k] == 0` instead of inlining the data value. The feasible set, the optimum and every other dual are unchanged; what is gained is that the pin's dual σ_k **is** `∂L/∂P_d(k)`, accumulated by the solver over every row that contains the demand, with no census of constraint families. Exported as `dual_Pd_pin.csv`. On a dispatch-only ACOPF it reproduces `dual_P_balance.csv`. `q_d` stays a constant. |
| `solve_explicit_dual` | `Bool` | `false` | **ED or DC-OPF (linear).** After a successful primal solve, build and solve the in-house dual LP; write to `Dispatch_Dual/`. See §2.3. |

### 2.2 Optional builder toggles

Mandatory equalities (slack-bus angle, P/Q balance, and branch flow equations where
required by the formulation) are **always built** — they are not configurable.

When a `bound_*` flag is `true`, limits are applied according to `bound_encoding`:

| `bound_encoding` | Primal | Dual export | `model_details.txt` box limits |
|---|---|---|---|
| `CONSTRAINT` (default) | explicit `@constraint` inequalities | `JuMP.dual` on `ineq_const_*` containers | constraint rows from `ineq_const` |
| `VARIABLE` | `lower_bound` / `upper_bound` on `@variable` | `bound_manifest` adapter → same CSV column names | equivalent ≤-form lines from `bound_manifest` |

| Group | Fields | Default highlights |
|---|---|---|
| Variable bounds | `bound_V`, `bound_θ`, `bound_P_g`, `bound_Q_g`, `bound_P_ik`, … | V/θ/P_g/Q_g on; branch flows off |
| Box encoding | `bound_encoding` | `CONSTRAINT` |
| Inequalities | `ineq_sbranch_upper`, `ineq_ang_diff_branch`, `ineq_sg_upper` | branch thermal + angle-diff on; gen capability off |
| Limit values | `limits::DispatchLimitsConfig` | see below |

See `DispatchConfig.jl` for the full field list.

**Limit sources.** `V`, `P_g`, `Q_g`, branch thermal ratings, and per-branch angle-difference limits always come from the case CSVs (`DBUS`, `DGEN`, `DCIR`). The remaining values live on `DispatchConfig.limits` (`DispatchLimitsConfig`):

| Field | Default | Meaning |
|---|---|---|
| `θ_min_rad`, `θ_max_rad` | ±π | Bus voltage-angle box when `bound_θ=true` (AC and DC) |
| `ang_diff_clamp_ac_deg` | `60.0` | Safety clamp on `DCIR.ang_min/ang_max` in the AC angle-diff constraint (values beyond ±clamp are tightened, with a printed warning) |
| `ang_diff_clamp_dc_deg` | `30.0` | Same clamp for the DC angle-diff constraint and the DC dual LP |

**Branch thermal / flow limits** are **not** fields of `DispatchLimitsConfig`. They use case data `DCIR.l_cap_1` under two layers:

1. Master toggle `DispatchConfig.ineq_sbranch_upper` — when `false`, no thermal inequalities for any branch (even if every `l_cap_1` is nonzero). Optional flow variable boxes use `bound_P_ik` / `bound_Q_ik` / … (default off).
2. Per-branch rating — when the toggle is on, only branches with `l_cap_1 != 0` get `±l_cap_1 / base_MVA`; unrated (`l_cap_1 == 0`) are skipped (no finite sentinel). The DC dual LP mirrors the same rated-branch set.

The explicit DC dual LP (§2.3) reads the same `limits` for θ / angle-diff clamps, so primal RHS and dual objective coefficients cannot drift apart.

### 2.3 Explicit dispatch dual (`solve_explicit_dual`)

Two ways to obtain dual information after a dispatch solve:

| Mechanism | Flag | Output folder | When to use |
|---|---|---|---|
| JuMP constraint duals from the **primal** | `RunConfig.save_duals = true` | `Dispatch/CSV_duals/` (+ `Dispatch/Dispatch_Duals.xlsx`) | Default; ED / DCOPF / ACOPF; folder created at run start, files written on successful solve |
| In-house **dual LP** (standalone formulation) | `dispatch.solve_explicit_dual = true` | `Dispatch_Dual/` (`CSV/`, `CSV_duals/`) | Strong-duality check; prices as decision variables; independent of `save_duals` |

**Requirements** (enforced in `validate_dispatch_config!`):

| Setting | ED | DC-OPF |
|---|---|---|
| `type_model` | `"ED"` | `"DCOPF"` |
| `cost_type` | `"linear"` | `"linear"` |
| `use_matrix` | (ignored) | `true` (Bbus primal) |

**ED dual:** scalar balance multiplier `λ`; system marginal price `π_SMP = -λ`; generator bound multipliers `α_lower`, `α_upper`. Outputs include `dual_ed_SMP.csv`, `dual_ed_alpha.csv`, `ED_Dual_Results.xlsx`.

**DC-OPF dual:** nodal `λ_k`, LMPs `π_k = -λ_k`, plus network multipliers when constraints are active. Outputs include `dual_dcopf_LMP.csv`, `DC_OPF_Dual_Results.xlsx`.

Quadratic primal cost → use `save_duals` on the primal NLP/QP instead.

**Safe defaults:** `solve_explicit_dual` defaults to `false`, so new users running the stock `main.jl` / smoke tests are unaffected.

**Example (ED, dispatch-only, 9-bus):**

```julia
cfg = RunConfig(
    trans_stab = false,
    dispatch = DispatchConfig(
        type_model = "ED",
        cost_type = "linear",
        solve_explicit_dual = true,
    ),
    solver_name = "Gurobi",
)
```

**Example (DC-OPF, dispatch-only, 9-bus):**

```julia
dispatch = DispatchConfig(
    type_model = "DCOPF",
    use_matrix = true,
    cost_type = "linear",
    solve_explicit_dual = true,
)
```

Results: primal under `RESULTS/.../Dispatch/`, explicit dual under `RESULTS/.../Dispatch_Dual/` (model dump, LMP CSV, duality gap report).

---

## 3. `TransientConfig` — transient analysis (avenue 2)

Defined in `_transient_stability/TsConfig.jl`. Set on `RunConfig.transient`.

```julia
transient = TransientConfig(
    simulation = TsSimulationConfig(δ_tol_deg = 100.0),
    builder    = TsBuilderConfig(),          # TS eq/ineq toggles
    dyn_model  = DynModelConfig(…),          # physics + fault
    gen_dynamic_filename = "gen_dynamic_data.csv",
)
```

### 3.1 `TsSimulationConfig` — timing and stability limits

| Field | Type | Default | Description |
|---|---|---|---|
| `δ_tol_deg` | `Float64` | `90.0` | Symmetric max \|δ − δ_COI\| [degrees] when lower/upper unset |
| `δ_tol_deg_lower` | `Float64` or `nothing` | `nothing` | Below-COI limit [degrees]; default → `δ_tol_deg` |
| `δ_tol_deg_upper` | `Float64` or `nothing` | `nothing` | Above-COI limit [degrees]; default → `δ_tol_deg` |
| `δ_pair_tol_deg` | `Float64` | `180.0` | Symmetric max \|δᵢ − δⱼ\| [degrees] for `bound_style_δ = :pairwise`; ignored by every other style |
| `t_start_sim` | `Float64` | `0.0` | Simulation origin [s] |
| `t_end_sim` | `Float64` | `5.0` | Simulation horizon [s] |
| `t_step` | `Float64` | `0.01` | Trapezoidal time step [s] |
| `t_start_fault` | `Float64` | `0.01` | Disturbance application [s] |
| `clearing_time` | `Float64` | `0.3` | SC fault duration before clearing [s] |
| `f_syn` | `Float64` | `50.0` | Synchronous frequency [Hz] |
| `Δω_0` | `Float64` | `0.0` | Initial speed deviation [p.u.] |

**SC timing (defaults):** fault at 0.01 s, clearing 0.3 s, post-fault to 5 s, Δt = 0.01 s.  
**GL timing:** single window from `t_start_fault` to `t_end_sim`.

### 3.2 `TsBuilderConfig` — which TS constraints are built

User-accessible toggles for **δ-COI stability inequalities** and **optional explicit
variable bounds**. Pre-fault init equalities that couple the OPF point to the
initial dynamic state are **always built by the path** (physics, not toggles):

| Constraint | When built |
|---|---|
| `eq_const_P_init` | Always (AC classical / DQ and Kron linear) |
| `eq_const_Q_init` | AC classical / DQ only — **omitted** on TSC-DCOPF / Kron linear |
| `eq_const_Pm_init` | Only when `mech_power_mode = USE_PM` — **omitted** under `USE_PG` |

Fault- and post-fault-window equalities (swing, Pe, COI, nodal balance) are likewise always built.

| Field | Default | Meaning |
|---|---|---|
| `ineq_δ_COI_tf_lower`, `ineq_δ_COI_tf_upper` | `true` | δ w.r.t. COI during fault (each side optional) |
| `ineq_δ_COI_tpf_lower`, `ineq_δ_COI_tpf_upper` | `true` | Same for post-fault window |
| `bound_E`, `bound_δ`, `bound_P_m` | `false` | Pre-fault explicit bound inequalities |
| `bound_δ_tf`, `bound_Δω_tf`, `bound_Pe_tf`, `bound_Qe_tf`, `bound_δCOI_tf` | `false` | Optional fault-window explicit bounds |
| `bound_δ_tpf`, `bound_Δω_tpf`, `bound_Pe_tpf`, `bound_Qe_tpf`, `bound_δCOI_tpf` | `false` | Optional post-fault explicit bounds |

Each family above exports one dual CSV per side and per window pair —
`dual_LB_δ_tf.csv` / `dual_UB_δ_tf.csv`, and likewise for `Δω`, `Pe`, `Qe` (per generator)
and `δCOI` (a COI scalar per time step). The `_tf` file name covers both windows: the
series is the fault window followed by the post-fault window, exactly as `dual_Pe.csv` is.
The dq boxes `Ed`, `Eq`, `Id`, `Iq`, `Te` follow the same rule.
| `bound_Ed`, `bound_Eq`, `bound_Id`, `bound_Iq` | `false` | DQ pre-fault algebraic boxes |
| `bound_V_ref`, `bound_P_ref` | `false` | AVR / governor set-point boxes |
| `bound_Ed_tf` … `bound_Te_tf` (+ `_tpf`) | `false` | DQ machine-state tf/tpf boxes |
| `bound_E_fd_unlim_tf`, `bound_E_fd_unlim_tpf` | `false` | AVR research boxes on `E_fd_unlim` |
| `bound_Pv_tf`, `bound_Pv_tpf`, `bound_Pm_tf`, `bound_Pm_tpf` | `false` | Governor research boxes (independent of `governor_limiter`) |
| `limits` | `TsBoundLimitsConfig()` | Limit values for the bounds above (see below) |
| `bound_encoding` | `CONSTRAINT` | How box limits attach (`CONSTRAINT` vs `VARIABLE`) |

**`TsBoundLimitsConfig`** (`builder.limits`): physical defaults — `E` [0, 2] pu, `δ` [-π, π], `P_m` from `DGEN` (`P_m_source = :dgen_pg_limits`), FullBus `V_bus_min_pu = 0` (forced on `V_tf`/`V_tpf`), governor valve `[gov_valve_min_pu, pg_max/base_MVA]` when `include_governor` and limiter ≠ `GOV_NO_LIMIT`. Optional tf/tpf `TsBoundLimitPair` fields default to `(-Inf, Inf)` = inactive (not `±9999`). `VARIABLE` encoding omits JuMP bounds; `CONSTRAINT` omits the ≤-row — Inf never enters an explicit inequality. Toggle on + Inf is a no-op until you set finite limits. DQ pre-fault `Ed/Eq/Id/Iq` and AVR `V_ref` have finite defaults when their toggles are on.

**GFM box knobs** (`builder.limits`, read only when `allow_gfm`): `gfm_δ_min_rad`/`gfm_δ_max_rad` (±π), `gfm_V_meas_min_pu`/`gfm_V_meas_max_pu` (0 / 2.5), `gfm_E_raw_extra_pu` (0.25), `gfm_E_raw_slack_pu` (0.5), `gfm_E_clip_slack_pu` (0.05), `gfm_PQ_bound_scale` (1.5), `gfm_PQ_bound_offset_pu` (0.05), `gfm_P_meas_floor_pu`/`gfm_Q_meas_floor_pu` (3.0), `gfm_I_floor_pu`/`gfm_I_ceiling_pu` (5.0 / 20.0), `gfm_pq_raw_span_pu` (20), `gfm_pq_clip_slack_pu` (0.05). These are margins and floors, not final numbers: `resolve_gfm_bound_limits` derives the per-unit box from them plus the `DGFM` row. The P and Q *measurement* boxes are solver guard-rails, not device ratings. Each box has a `TsBuilderConfig` toggle — `bound_gfm_δ` plus `bound_gfm_{P_meas,Q_meas,V_meas,E_int_raw,E_int,E_droop_raw,E_droop,Id,Iq}_{tf,tpf}`, and `bound_gfm_pq_{raw,clip}_{tf,tpf}` for the overload states — which unlike the SG `bound_*` families default to **`true`**, since the nonconvex limiter needs the guard-rails. They follow `bound_encoding`. The dispatch range on a converter is the `DGEN` box (`pg_min`/`pg_max`, `qg_min`/`qg_max`). `DGFM.Pmin`/`Pmax`/`Qmin`/`Qmax` are comparator inputs to the overload PIs and do not box `P_set`. `validate_gfm_dispatch_band!` requires the `DGEN` interval to sit inside the `DGFM` band. Defaults reproduce the values the reference port hard-coded. The smoothing constants `_GFM_EPS_E`, `_GFM_PQLIM_EPS_STEP_FRAC` and `_GFM_PQLIM_OUT_ERR_TOL` (overload state-clip ε as a fraction of one integrator step, output-clip ε as `Kp` × an error tolerance; both resolved per converter and axis into `run_manifest.toml`), `_GFM_EPS_NORM_V`, `_GFM_EPS_LIM_I` and the bypass threshold `_GFM_IMAX_NO_LIMIT` are module constants, not configuration. They live in `src/_transient_stability/numerics.jl` with the AVR and governor smoothing constants. They are `const`, so restart Julia after editing one: a live session, Revise included, keeps the compiled value.

Full toggle inventory and “wired vs forced” notes: `OUTPUTS/reference/constraint_toggle_map.md`.

Flip these in `main.jl` under `transient.builder` (the `TsBuilderConfig` struct).

### 3.3 `gen_dynamic_filename` / `gfm_dynamic_filename`

| Field | Default | Description |
|---|---|---|
| `gen_dynamic_filename` | `"gen_dynamic_data.csv"` | **SG** machine-parameter CSV in the case folder |
| `gfm_dynamic_filename` | `"gfm_dynamic_data.csv"` | **GFM** machine-parameter CSV (only when `allow_gfm=true`) |

**SG-only cases** (`allow_gfm=false`, default): row `i` in `gen_dynamic_data` must correspond to generator `id = i` and the same `bus` as row `i` of `generators_data.csv`. `load_system` calls `validate_dgen_dyn_row_ids!`.

**Mixed SG+GFM** (`allow_gfm=true`): do **not** put GFM columns in `gen_dynamic_data`. Split instead:

- `generators_data.csv` — all units (SG and GFM)
- `gen_dynamic_data*.csv` — **SG rows only** (optional `id` column; required when SG ids are not `1:n_sg`)
- `gfm_dynamic_data.csv` — **GFM rows only**; see the column table below

Every `DGEN.id` must appear in exactly one of the two dynamic files. Example case: `INPUT_FILES/9bus_gfm/`. Phases G0–G2 implement ingest, init/ACOPF limits, and transient GFM dynamics. The physics behind each column is in [Part I §9 — Grid-forming inverters](model/09_grid_forming.md).

| Column | Required | Meaning |
|---|---|---|
| `id` | yes | Generator id; must match a `DGEN` row and must **not** appear in `gen_dynamic_data*` |
| `bus` | yes | Terminal bus |
| `Xl` | yes | Coupling (filter) reactance behind the terminal |
| `mq` | yes | Q–V droop gain |
| `mp` | yes | P–f droop gain; $\Delta\omega = m_p(P_{\mathrm{set}} - P_{\mathrm{meas}})$ |
| `Kpv`, `Kiv` | yes | Q–V PI proportional and integral gains |
| `Emax`, `Emin` | yes | Internal-voltage clip, and the ACOPF-side internal-voltage limit |
| `Pmax`, `Pmin` | yes | Active-power device rating. Comparator inputs to the P overload PI when `gfm_pq_overload` is on. Not a box on `P_set` — `DGEN` owns the dispatch range, and it must sit inside this band |
| `TPf`, `TQf`, `TVf` | yes | Measurement-filter time constants for the $P$, $Q$ and $V$ channels (REGFM_A1 eq. 4–6), independent of one another; `0` means unfiltered (pass-through). A legacy single `Tf` column is still accepted in their place and is broadcast into all three, so case folders and archived `RESULTS/*/Inputs/` copies written before the split keep loading. Only `TPf` feeds the COI weight $H^{\mathrm{eff}} = T_{Pf}/(2 m_p)$ used by `gfm_in_coi` |
| `Imax` | yes | Current-limiter threshold **for the transient only** — brief overcurrent up to the junction thermal limit (e.g. 1.2 pu). The steady-state ACOPF circle uses the nameplate rating instead (decision D4, [Part I §9](model/09_grid_forming.md)). **A value ≥ 20 disables the limiter entirely** (`_GFM_IMAX_NO_LIMIT`), removing those rows from the model rather than merely relaxing them; bypassed units are logged and recorded on `meta[:gfm_limiter_bypassed]` |
| `Qmax`, `Qmin` | yes | Reactive-power device rating (REGFM_A1 Table 1, `VFlag = 1`: `±0.44` pu on inverter base). Comparator inputs to the Q overload PI. Divided by `base_ratio` like `Pmax`/`Pmin`. Not a second box on `Q_g` |
| `Kppmax`, `Kipmax` | yes | P-overload PI gains, shared by the upper and lower branch. Map a power error to per-unit `Δω`, so they are multiplied by `base_ratio` like `mp`. Table 1: `0.01` and `0.1` |
| `Kpqmax`, `Kiqmax` | yes | Q-overload PI gains, shared by the upper and lower branch. `Kiqmax` is also accepted as `kiqmax`. Map a power error to a voltage command, so they are multiplied by `base_ratio` like `mq`. Table 1, `VFlag = 1`: `3` and `20` |
| `InvBase` | optional | Converter MVA base (aliases `mach_base_MVA`, `inv_base`). Defaults to the system base; a `0` entry also falls back to it. `apply_gfm_base_conversion!` rescales `Xl`, `mq`, `mp`, `Kppmax`, `Kipmax`, `Kpqmax`, `Kiqmax`, `Pmax`, `Pmin`, `Qmax`, `Qmin`, `Imax` to system base at read time |

Note also that a legacy `Eg` column in `gen_dynamic_data*.csv` is accepted and **ignored** — the internal emf is solved from the init equalities, not read from CSV.

---

## 4. `DynModelConfig` — physics and network form

Nested in `RunConfig.transient.dyn_model`. Defined in `_transient_stability/DynModelConfig.jl`.

### 4.1 Core selection

| Field | Options | Default | Description |
|---|---|---|---|
| `gen_order` | `CLASSICAL_2ND`, `DQ_4TH` | `CLASSICAL_2ND` | `DQ_4TH` = 4th-order dq machine (FULL_BUS only; see preset below) |
| `network_form` | `KRON_REDUCED`, `FULL_BUS` | `KRON_REDUCED` | Kron vs full-network nodal |
| `mech_power_mode` | `USE_PG`, `USE_PM` | `USE_PG` | Swing uses `P_g` or explicit `P_m` |
| `include_governor` | `Bool` | `false` | TGOV1 turbine governor → time-varying `P_mech(t)`; **FULL_BUS + USE_PM only** |
| `include_avr` | `Bool` | `false` | SEXS AVR → time-varying `E_fd(t)`; **DQ_4TH + FULL_BUS + USE_PM only** |
| `governor_limiter` | `GOV_NO_LIMIT`, `GOV_SMOOTH`, `GOV_HARD_BOUND` | `GOV_NO_LIMIT` | Valve saturation treatment (only read when `include_governor`) |
| `allow_gfm` | `Bool` | `false` | Mixed fleet with separate `gfm_dynamic_data.csv`; **requires `DQ_4TH` + `FULL_BUS`** (Phase G2: init + ACOPF limits + transient filters/droop/PI/limiter) |
| `gfm_in_coi` | `Bool` | `false` | Converters take part in the Centre of Inertia, weighted by `H_eff = TPf/(2·mp)`, and gain `:coi_box` corridor rows under the same tolerances as the machines. Governs `bound_style_δ = :coi_box` and `bound_style_Δω = :coi_box` together; warns when inert. **Required** on a converter-only fleet whenever either corridor references the COI |
| `gfm_pq_overload` | `Bool` | `false` | REGFM_A1 Fig. 3 P and Q overload PIs. Off, the droop row and the Q–V voltage error are unchanged. On, eight states per converter-step. The integrators follow `gfm_integrator`. Warns when `allow_gfm` is false |
| `allow_all_gfm` | `Bool` | `false` | Permit a **converter-only** fleet (every active generator in `gfm_dynamic_data.csv`). Default `false` throws, so an all-converter run can never happen by accident. When true, `validate_all_gfm_fleet!` still rejects `bound_style_δ ∈ (:swing_propagated, :highest_H)` (both need machine inertia), a COI-referenced corridor without `gfm_in_coi`, and `include_avr`/`include_governor` (SG-only families, so they would be silent no-ops). See `INPUT_FILES/9bus_gfm_all/` |
| `dq_speed_dev_in_algebra` | `Bool` | `true` | **DQ_4TH only:** include `(1+Δω)` in `Pe` and stator `Vd`/`Vq` algebra (RMS-style); `false` neglects speed in those equations |

**Common presets:**

| Scenario | Suggested settings |
|---|---|
| Classic TSC-ACOPF (Kron) | defaults: `KRON_REDUCED`, `USE_PG`, `:swing_propagated` |
| FULL_BUS TSC-ACOPF | `network_form=FULL_BUS`, `mech_power_mode=USE_PM`, `bound_style_δ=:coi_box` |
| FULL_BUS + turbine governor | above **plus** `include_governor=true`, `gen_dynamic_filename="gen_dynamic_data_full.csv"` (needs `R,T1,T2,T3`) |
| DQ_4TH FULL_BUS TSC-ACOPF | `gen_order=DQ_4TH`, `FULL_BUS`, `USE_PM`, `bound_style_δ=:coi_box`, `gen_dynamic_filename="gen_dynamic_data_full.csv"` (needs `Xq_tr,Xd,Xq,Td,Tq,Ra`); optional `dq_speed_dev_in_algebra=false` |
| DQ_4TH + AVR | DQ preset **plus** `include_avr=true` (needs `T_exc,K_exc,Ta_exc,Tb_exc` in full CSV; `Ta_exc=Tb_exc=0` bypasses the lead-lag) |
| DQ_4TH + governor | DQ preset **plus** `include_governor=true` (needs `R,T1,T2,T3` in full CSV) |
| DQ_4TH + AVR + governor | Both flags true (reference-style full control stack) |
| Mixed SG + GFM (Phase G2) | DQ preset **plus** `allow_gfm=true`, case `9bus_gfm`, split dyn CSVs (see §3.3) |
| TSC-DCOPF (linearised Pe) | `dispatch.type_model="DCOPF"`, Kron defaults |

#### Turbine governor (TGOV1)

When `include_governor=true`, mechanical power becomes a time-varying state driven by speed
deviation through a droop/lead-lag governor + turbine (`T1·dPv/dt=(P_ref−Δω)/R−Pv`,
`T3·dPm/dt+Pm=T2·dPv/dt+Pv`), replacing the constant `P_m` in the swing
equation. Requires machine columns `R, T1, T2, T3` (supplied by `gen_dynamic_data_full.csv`).
The set-point is pinned to the dispatch (`P_ref=R·P_m`, `P_m=P_g`). Valve saturation:

- `GOV_NO_LIMIT` — linear governor, valve unbounded (smoothest, cleanest duals).
- `GOV_SMOOTH` — sqrt smooth min/max anti-windup clamp (physically saturating, nonconvex).
- `GOV_HARD_BOUND` — explicit `0 ≤ P_valve ≤ pg_max/base_MVA` (clean bound duals; can be
  infeasible on a strong disturbance — that infeasibility is a modelling signal, not a bug).

Governor trajectories are written to `Transient_Stability/CSV/governor_*.csv`
(`P_mech`, `P_valve`, `P_valve_raw`, `P_ref`), with SVG figures for the three time series when
`save_ts_plots=true`. `P_valve_raw` is the integrator state ahead of the limiter, so
`P_valve_raw − P_valve` is the clamp actually applied under `GOV_SMOOTH`. The integrator is
anti-windup — each step advances from the *limited* `P_valve` of the previous step — so that
gap stays within one step's valve travel instead of growing while the valve sits on its limit.

Duals `dual_Pref_init`, `dual_gov_valve`, `dual_gov_mech` join the standard TS dual export, plus
one limiter family depending on the mode: `dual_gov_valve_sat` (`GOV_SMOOTH` clamp equality) or
`dual_LB_gov_valve` / `dual_UB_gov_valve` (`GOV_HARD_BOUND` ≤-form pair). `GOV_NO_LIMIT` adds neither.

With a governor on an SG-only fleet, `P_m(t)` is the governor mechanical-power state, so
`mechanical_power.csv` is **not** written — it would duplicate `governor_P_mech.csv`. Runs
without a governor still get it, and so do mixed SG + GFM runs: `governor_P_mech.csv` covers
the governed machines only, while `mechanical_power.csv` also carries each converter's
constant `P_set`.

Under `bound_encoding = VARIABLE`, `GOV_HARD_BOUND` stamps the limits on `P_valve_raw`
itself and exports `dual_LB_gov_valve` / `dual_UB_gov_valve` through the bound manifest —
same column names as the `CONSTRAINT` ≤-rows.

Block diagram and equations: [`docs/dynamic_controls_avr_governor.md`](dynamic_controls_avr_governor.md) (governor section).

#### AVR (SEXS exciter) — DQ_4TH only

When `include_avr=true`, field voltage `E_fd` becomes a time-varying state driven by
terminal-voltage feedback through the SEXS chain:

```text
u = V_ref − V_terminal
E_LL = (1 + Ta_exc·s) / (1 + Tb_exc·s) · u          # lead-lag (optional)
T_exc · d(E_fd_unlim)/dt = K_exc · E_LL − E_fd_unlim  # gain-lag
E_fd = smoothclamp(E_fd_unlim; E_min, E_max)
```

Requires `T_exc, K_exc, Ta_exc, Tb_exc` in `gen_dynamic_data_full.csv`. Saturation
limits come from `TsBoundLimitsConfig` (default `[0, 2]` pu). The saturated `E_fd(t)`
couples into the Eq EMF ODE. Pre-fault: `E_fd = K_exc·(V_ref − V_bus)` (unity DC gain
through the lead-lag, so the steady-state link is unchanged). Can be combined with
`include_governor=true` on DQ_4TH.

**Lead-lag bypass.** `Ta_exc = Tb_exc = 0` skips the lead-lag stage entirely for that
generator — no `E_LL` variables, no rows — matching ANDES `zero_out=True`. Every
shipped fixture uses this, so default AVR runs are identical to the old
first-order-only model. `Tb_exc = 0` with `Ta_exc > 0` is rejected
(`validate_dyn_data!`): that is a bare differentiator with no state-space realization.
Both time constants must be ≥ 0. When `Tb_exc > 0` and `t_step > 2·Tb_exc`, the
builder warns that the trapezoidal update will alternate in sign (stable ringing).

**Gain-lag values.** `T_exc` and `K_exc` must both be **> 0**, also enforced by
`validate_dyn_data!`. The exciter row divides by `T_exc` (`c = Δt/(2·T_exc)`), and
`K_exc = 0` would pin the pre-fault link to `E_fd = 0` and strip the exciter of its
voltage feedback — a run that converges to something meaningless rather than failing.
Both checks apply only under `include_avr = true`, so a DQ-only run may still leave
these columns blank.

Trajectories: `avr_V_ref.csv`, `dq_E_fd_pu.csv` (saturated), `dq_E_fd_unlim_pu.csv`
(pre-saturation), and `dq_E_LL_pu.csv` when the lead-lag is active. The gap between
`E_fd` and `E_fd_unlim` shows when the exciter is clamping.
Duals: `dual_Vref_init`, `dual_avr_leadlag` (when present), `dual_avr_E_fd` (ODE on
`E_fd_unlim`), `dual_avr_E_fd_sat` (clamp).

See [`docs/dynamic_controls_avr_governor.md`](dynamic_controls_avr_governor.md) for the AVR block diagram.

#### DQ_4TH machine (machine core + optional AVR / governor)

`gen_order=DQ_4TH` adds subtransient emfs `Ed`/`Eq`, dq currents `Id`/`Iq`, and field voltage
`E_fd` (constant by default; dynamic when `include_avr=true`). Mechanical power is constant
`P_m` by default; time-varying `P_mech(t)` when `include_governor=true`. Requires
`gen_dynamic_data_full.csv` with columns `Xq_tr, Xd, Xq, Td, Tq, Ra` in addition to the usual
`Xd_tr, H, D` (plus `T_exc,K_exc,Ta_exc,Tb_exc` for AVR and `R,T1,T2,T3` for governor). SC and GL faults are
supported on the FULL_BUS network path.
`TSC-DCOPF` + `DQ_4TH` is not implemented.

### 4.2 Loads and stability corridors

| Field | Type | Default | Description |
|---|---|---|---|
| `zip_load_p` | `(Float64,Float64,Float64)` | `(1,0,0)` | **Active**-demand `(Z,I,P)` = impedance / current / power fractions; **must sum to 1**; default = constant impedance |
| `zip_load_q` | `(Float64,Float64,Float64)` | `(1,0,0)` | **Reactive**-demand `(Z,I,P)` split; independent of `zip_load_p`; **must sum to 1** |
| `constrain_δ` | `Bool` | `true` | Build the rotor-angle corridor |
| `bound_style_δ` | `Symbol` | `:swing_propagated` | Reference and form of that corridor — see the table below |
| `δ_ref_gen_id` | `Int` or `nothing` | `nothing` | Reference machine id; **required** by `bound_style_δ = :ref_gen`, ignored (with a warning) otherwise |
| `constrain_Δω` | `Bool` | `false` | Build the speed corridor |
| `bound_style_Δω` | `Symbol` | `:coi_box` | `:coi_box` = box on `Δωᵢ − Δω_COI`; `:abs` = box on the raw `Δωᵢ`, no COI formed |
| `Δω_tol_pu` | `Float64` | `0.5` | Symmetric half-width [p.u.] when lower/upper unset |
| `Δω_tol_pu_lower` | `Float64` or `nothing` | `nothing` | Below-reference limit [p.u.]; default → `Δω_tol_pu` |
| `Δω_tol_pu_upper` | `Float64` or `nothing` | `nothing` | Above-reference limit [p.u.]; default → `Δω_tol_pu` |

#### The rotor-angle corridor

The two corridors are independent: each has its own toggle and its own reference, and a
transient run must enable at least one of them. Angle half-widths come from
`TsSimulationConfig.δ_tol_deg*`; the reference is chosen by `bound_style_δ`:

| `bound_style_δ` | Reference | Rows built |
|---|---|---|
| `:swing_propagated` | δ_COI | `δ_tol[1] ≤ δ_g − δ_COI ≤ δ_tol[2]` with one swing step substituted for `δ_g[t]`, so the rows also carry `P_mech`, `Pe`, `Δω`, `ω_syn`, `Δt` |
| `:coi_box` | δ_COI | the same band, boxed directly |
| `:highest_H` | largest-`H` **surviving** machine | `δ_tol[1] ≤ δ_g − δ_ref ≤ δ_tol[2]` |
| `:ref_gen` | `δ_ref_gen_id` | same |
| `:pairwise` | none | `−Δ ≤ δᵢ − δⱼ ≤ +Δ` for every unordered pair `i < j`, with `Δ = δ_pair_tol_deg` |

The machine-referenced styles follow the convention used in much of the
transient-stability literature: bound the angular separation of the n−1 remaining machines
against one live reference unit. The reference carries **no row of its own** — `δ_ref − δ_ref ≡ 0`
can never bind and its dual would mean nothing — so a fleet of `n` machines yields `n−1`
corridors per window.

`:highest_H` picks over the machines that **survive the disturbance**, so a unit the run is
about to trip can never become the reference. The resolved id is recorded in
`dynamic_model_details.txt` and in the `δ_ref_gen` column of `swing_debug.csv`.

On a mixed SG + GFM fleet the two machine-referenced styles **also bound the converters**,
against the same reference and the same `δ_tol`: `δ_g − δ_ref` needs no inertia from `g`,
only an angle in the same frame. `:swing_propagated` never does, because it propagates a
per-machine `H` *and* `D`. `:coi_box` does so when `gfm_in_coi = true`, which also weights
the converters into the COI mean itself with `H_eff = Tf/(2·mp)`, the equivalent inertia of
the LPF-droop frequency loop; left `false` (the default) it stays SG-only. The same
split applies on the speed side — `bound_style_Δω = :abs` always spans the converters,
`:coi_box` only under `gfm_in_coi`. A converter can be the reference under `:ref_gen`, but never under `:highest_H`,
which ranks by `H`. Converter duals arrive as extra `Gen_*` columns in the existing
`dual_delta_ref_*.csv` / `dual_Delta_Omega_abs_*.csv`; `swing_debug.csv` remains SG-only.

Every style is honoured on **every** network form (Kron, Kron-linear, FULL_BUS), and they
price differently — a Kron-vs-FULL_BUS dual comparison must fix the same style on both
sides. Note that `USE_PM` rejects `:swing_propagated`, and `FULL_BUS` requires `USE_PM`, so
in practice the propagated form is a Kron + `USE_PG` option.

##### `:pairwise` — the separation corridor

`:pairwise` drops the reference entirely and states the loss-of-synchronism criterion in its
direct form: every unordered pair of synchronised units is held within `δ_pair_tol_deg` of
each other. That knob is **separate from `δ_tol_deg*`** and symmetric, because `|δᵢ − δⱼ|`
has no "above" and "below" reading — the two sides of the box already cover both signs of
the difference.

Three consequences are worth knowing before running one:

- **It spans SG + GFM.** A separation bound needs two angles in the same frame and no
  inertia from either, so converters are corridor members exactly as under `:highest_H` /
  `:ref_gen`.
- **It is quadratic in fleet size.** `n(n−1)/2` pairs per time step against `n−1` for the
  reference-framed styles — 3 pairs on the 9-bus, 45 on a 10-machine 39-bus.
- **Only the widest-spread pair can bind.** The other rows sit inactive with a zero dual,
  which is exactly what makes the export readable: the non-zero column *names* the two units
  the corridor is holding apart. When two pairs bind at the same instant with the same
  spread (symmetric fleets) the split between their multipliers is not unique — read the sum
  over pairs as the corridor's rent at that step.

A pairwise threshold `Δ` is **not** comparable to a reference-framed threshold `Δ`:
`|δᵢ − δ_ref| ≤ Δ` for all `i` only implies `|δᵢ − δⱼ| ≤ 2Δ`. Comparisons across styles must
be made on the realised spread, not on the knob.

Duals export under `dual_δ_COI_lower` / `dual_δ_COI_upper` for the COI-referenced styles,
under `dual_δ_ref_lower` / `dual_δ_ref_upper` for the machine-referenced ones, and under
`dual_δ_pair_lower` / `dual_δ_pair_upper` for `:pairwise`. Only one family exists in a given
run; the names are kept apart because the shadow price of "stay within δ_tol of the COI", of
"stay within δ_tol of machine k" and of "stay within Δ of every other unit" are different
quantities. The pairwise CSV and Excel columns are labelled `Pair_G<i>_G<j>` rather than
`Gen_<id>`, since each row belongs to two units and to neither alone.

When the corridor does not reference the COI, no `δ_COI` variable and no defining equality
are built — the COI survives as an unconstrained expression, so `angle_rel_COI.csv` and the
`δ_rel_COI` column of `swing_debug.csv` are still written, but `dual_δCOI` is absent and
`TsBuilderConfig.bound_δCOI_tf` / `bound_δCOI_tpf` are rejected at validation. A
machine-referenced run additionally writes `angle_rel_ref.csv`, which is the trajectory the
corridor actually bounds. With `save_ts_plots = true` the same trajectory is drawn in
`Figures/Delta (ref G<k>) vs time.svg` (reference machine flat at zero, ±δ_tol as horizontal
lines), and `Delta vs time.svg` anchors its ±δ_tol band on machine `k` rather than on the COI —
`Delta (ref COI) vs time.svg` is still written, but it is not the picture of the constraint here.

#### The speed corridor

`constrain_Δω = true` boxes the speed deviation, and `bound_style_Δω` picks against what.
`:coi_box` builds the `Δω_COI` variable and bounds `Δωᵢ − Δω_COI`; `:abs` bounds the raw
`Δωᵢ` of the swing equation and forms no COI at all.

!!! warning "`:abs` also prices the common-mode excursion"
    These are not the same constraint with the COI held at zero. `Δω_COI` is a free
    variable the optimiser can shift, so the COI box prices only the *spread* between
    machines, while the absolute box additionally pins how far the whole fleet drifts from
    synchronous speed. During a bus short circuit that drift is largely redispatch-invariant
    — mechanical power is fixed while electrical power collapses — so an absolute band
    tighter than it makes the problem **infeasible**, not merely expensive. On the 9-bus
    case at `load_factor = 1.5`, contingency 2, `Δω_tol_pu = 0.02` is infeasible while the
    `0.5` default solves. Start loose and tighten.

Like `δ_tol_deg_lower` / `δ_tol_deg_upper`, both Δω overrides are **positive magnitudes**, so
`Δω_tol_pu_lower = 0.05` with `Δω_tol_pu_upper = 0.03` gives the signed corridor
`−0.05 ≤ Δωᵢ − Δω_ref ≤ 0.03`. Setting only one leaves the other at `Δω_tol_pu`. Useful when
under- and over-frequency excursions carry different consequences. Both styles read the same
three fields, and duals export as `dual_Δω_COI_*` or `dual_Δω_abs_*` respectively.

#### ZIP load splits

The two ZIP vectors are **independent** because TSOs do not generally apply the same load
model to both channels. The Spanish TSO, for instance, recommends representing active demand
as constant current and reactive demand as constant admittance during RMS stability studies:

```julia
DynModelConfig(
    network_form = FULL_BUS,
    zip_load_p   = (0.0, 1.0, 0.0),   # (Z, I, P) — constant current
    zip_load_q   = (1.0, 0.0, 0.0),   # (Z, I, P) — constant admittance
)
```

!!! note "Order is (Z, I, P) in both vectors"
    Index 1 is **impedance**, index 2 is **current**, index 3 is **power**. The reference
    implementation stores the vector as `(P, I, Z)`, so its `[0,0,1]` (constant-Z) is our `(1,0,0)`.

Both splits apply only on `network_form = FULL_BUS` (classical and `DQ_4TH`, SG-only or mixed
SG+GFM). On `KRON_REDUCED` loads are folded into `Y_red` as constant admittance, so the two
fields are ignored and a non-default value raises a warning.

These two fields set the load model for the **whole system**. To give individual loads at
the same bus different splits — an industrial constant-impedance load beside a
constant-power one — add a `load_data.csv` to the case folder (§6.2); each bus then enters
its balance with the demand-weighted mean of the shares of the loads on it, and these
fields become the default for every bus the table does not mention.

Because the balance is scaled by the steady-state voltage `V`, the load term reads
`p_d·(Z·V_t² + I·V_t·V + P·V²)`, which collapses to `p_d·V²` at the pre-fault point
`V_t = V` for any split summing to 1. A mixed P/Q pair therefore stays consistent with the
constant-power ACOPF warm start at `t = 0`; the two channels only diverge off equilibrium.

### 4.3 Validation rules (important)

!!! warning "Enforced in `validate_dyn_config!`, not just convention"
    These are hard `ArgumentError`s raised before any model is built, not style
    guidance — a run that violates one fails fast rather than building a model
    that silently means something different than you intended.

- `USE_PM` **rejects** `bound_style_δ = :swing_propagated` (that row substitutes the dispatch power into the angle band, which is inconsistent with an explicit `P_m`); the three box styles are all valid
- a transient run **requires at least one corridor**: `constrain_δ = true` or `constrain_Δω = true`
- `bound_style_δ = :ref_gen` **requires** `δ_ref_gen_id`; the id is checked against the system data by `validate_δ_reference!` (in service, not tripped, at least two surviving units) before the warm-start solve. A GFM id is accepted; the `H > 0` and id-range checks then do not apply, since `gen_dynamic_data` holds machine rows only. `:highest_H` additionally requires one surviving machine with `H > 0` — it ranks by inertia, so its reference is never a converter
- `bound_style_δ = :pairwise` requires `δ_pair_tol_deg > 0` and, like the machine-referenced styles, **at least two surviving units** after the disturbance (`validate_δ_reference!`); it resolves no reference machine, so `δ_ref_gen_id` is inert on such a run
- `bound_style_δ ∈ (:highest_H, :ref_gen, :pairwise)` **rejects** `TsBuilderConfig.bound_δCOI_tf` / `bound_δCOI_tpf`, which need a `δ_COI` variable those styles do not build
- `FULL_BUS` **requires** `mech_power_mode = USE_PM`, and always pre-solves the steady-state ACOPF that seeds the coupling variables (there is no flat-start alternative; a failed pre-solve aborts the run)
- `ode_first_step = :backward_euler` **requires** `network_form = FULL_BUS` — the Kron swing rows are trapezoidal at every step, so the Kron paths reject it instead of silently ignoring it
- `TSC-DCOPF` + `FULL_BUS` → **not implemented**
- `zip_load_p` and `zip_load_q` coefficients must each sum to 1 (checked per vector)
- `DQ_4TH` **requires** `FULL_BUS`, `mech_power_mode=USE_PM`, a box-form `bound_style_δ`, and full machine columns in `gen_dynamic_data` (use `gen_dynamic_data_full.csv`)
- `include_avr=true` **requires** `gen_order=DQ_4TH` (and `T_exc,K_exc,Ta_exc,Tb_exc` in the dynamic CSV; `Ta_exc=Tb_exc=0` is the lead-lag bypass)
- `include_governor=true` on DQ_4TH uses the same TGOV1 layer as classical FULL_BUS (`R,T1,T2,T3` in CSV)
- `TSC-DCOPF` + `DQ_4TH` → **not implemented**
- `include_avr` **requires** `gen_order = DQ_4TH` when enabled (reserved)
- `include_governor` **requires** `mech_power_mode = USE_PM` **and** `network_form = FULL_BUS` (not supported on the Kron path)

---

## 5. `FaultConfig` — disturbances

Nested in `RunConfig.transient.dyn_model.fault`. Defined in `_transient_stability/FaultConfig.jl`.

### 5.1 Short-circuit (SC) — default

```julia
dyn_model = DynModelConfig(
    fault = FaultConfig(
        fault_type = SC,
        contingency_id = 2,
    ),
)
```

**Contingency table** — `INPUT_FILES/<case>/contingencies.csv`:

| Column | Meaning |
|---|---|
| `contingency_nr` | Row id (matches `contingency_id`) |
| `faulted_bus` | Bus where the fault is applied |
| `circuit_id` | Branch index to disconnect after clearing |
| `from_bus`, `to_bus` | Branch endpoints (informational) |

Timing is controlled by `transient.simulation` (§3.1), not `contingencies.csv`.

| Field | Default | Meaning |
|---|---|---|
| `sc_fault_admittance_pu` | `0 - im*FAULT_BUS_SHUNT` (`FAULT_BUS_SHUNT = 1e6`) | Faulted-bus shunt, pu on `base_MVA`. Read only for `fault_type = SC` |
| `sc_fault_stamp` | `:add` | `:add` adds the shunt to the network `Y_kk`; `:overwrite` replaces the diagonal |

The default is an inductive solid short: the residual voltage stays near the pre-fault phasor and the fault draws no active power. The resistive solid short remains available as `FAULT_BUS_SHUNT + 0im` with `:overwrite`: the residual voltage sits near −90° and the fault is an artificial active-power sink during the fault, with a small shadow price on that bus's fault-window P balance. A positive imaginary part is capacitive (`+jB` in the Ybus) and produces a warning.

### 5.2 Generator / load disconnection (GL)

**Generator trip** (9-bus — trip gen 3):

```julia
dyn_model = DynModelConfig(
    network_form = FULL_BUS,
    mech_power_mode = USE_PM,
    bound_style_δ = :coi_box,
    fault = FaultConfig(fault_type = GL, gl_gen_ids = [3]),
)
```

**Load trip** (9-bus — remove 100% load at bus 5):

```julia
fault = FaultConfig(
    fault_type = GL,
    gl_load_bus_ids = [5],
    gl_percent_power = [-1.0],   # alpha = -1  =>  full trip (zero demand)
)
```

**Load scaling convention:** for each disturbed bus $k$,

```math
p_{d,k}^{\mathrm{new}} = p_{d,k}\,(1 + \alpha_k), \qquad
q_{d,k}^{\mathrm{new}} = q_{d,k}\,(1 + \alpha_k)
```

where $\alpha_k$ is the parallel entry in `gl_percent_power`.

| `gl_percent_power` ($\alpha$) | Effect |
|---|---|
| `-1.0` | Disconnect 100% of load ($p_d, q_d \to 0$) |
| `-0.5` | Shed 50% (half of pre-fault demand retained) |
| `0.0` | Unchanged |
| `+1.0` | Double demand |

Values $\alpha < -1$ are rejected (would imply negative demand).

| Field | Description |
|---|---|
| `gl_gen_ids` | Generator row indices in `generators_data.csv` |
| `gl_load_bus_ids` | Bus indices whose demand is scaled (parallel to `gl_percent_power`) |
| `gl_percent_power` | Per-bus scaling factor $\alpha$ in the formula above |

GL modifies **copies** of network DataFrames inside the TS builder; CSV files on disk are unchanged.

### 5.3 Open branch (OB)

Open one or more in-service lines/transformers **without** a short-circuit stage.
Same single-window timeline as GL. Branch indices are **DCIR row indices**
(1-based; same numbering as `circuit_id` in `contingencies.csv` / `line_data.csv` order).

```julia
fault = FaultConfig(fault_type = OB, ob_branch_ids = [6])       # one circuit (9-bus 5–7)
# Multi-branch is allowed when the remaining network stays connected from the slack
# (on the classic 9-bus, most N-2 opens island — validation hard-fails those).
```

| Field | Description |
|---|---|
| `ob_branch_ids` | DCIR row indices to open (`l_status → 0` in the TS copy) |

Validation requires every listed branch to be in service, rejects mixing with `gl_*`
fields, and **hard-fails** if the open set islands the network or isolates the slack.

FULL_BUS classical, Kron, and DQ FULL_BUS builders implement OB (same single
dynamic window as GL).

---

## 6. Network input files (CSV)

Place case data in `INPUT_FILES/<case>/`. Selected by `RunConfig.case`.

| File | Contents |
|---|---|
| `bus_data.csv` | Bus id, type (3=slack), `Pd`/`Qd`, voltage limits |
| `generators_data.csv` | Bus, P/Q limits, cost coeffs, `status` |
| `line_data.csv` | From/to bus, `r,x,b`, ratings, `status` |
| `gen_dynamic_data.csv` | `bus`, `Xd_tr`, `H`, `D` (minimal) |
| `contingencies.csv` | SC contingency definitions |
| `load_data.csv` | *Optional* — splits a bus's demand into individual loads, each with its own ZIP shares (§6.2) |

Use `load_factor` for parametric demand studies rather than editing CSVs.

### 6.1 MATPOWER `.m` input — two workflows

There are two ways to use a MATPOWER `.m` file as input.

#### 6.1.1 Direct `RunConfig.matpower_file` (recommended for TSC runs)

Set `RunConfig.matpower_file` to the `.m` filename (relative to the case folder) or an absolute
path. `load_system` will parse the file in-memory — no CSV round-trip — remap buses to 1..nBUS,
and scale demand via `load_factor` exactly as in CSV mode. Dynamic data (`gen_dynamic_data.csv`)
must still be present in the case folder.

```julia
cfg = RunConfig(
    trans_stab    = true,
    case          = "9bus_test_mfile",   # INPUT_FILES/9bus_test_mfile/
    matpower_file = "case9.m",           # relative to INPUT_FILES/9bus_test_mfile/
    load_factor   = 1.5,
    solver_name   = "Ipopt",
    dispatch = DispatchConfig(type_model = "ACOPF"),
    transient = TransientConfig(
        simulation = TsSimulationConfig(δ_tol_deg = 100.0),
        dyn_model = DynModelConfig(
            network_form    = FULL_BUS,
            mech_power_mode = USE_PM,
            bound_style_δ     = :coi_box,
            zip_load_p        = (1.0, 0.0, 0.0),  # (Z, I, P) — active demand, constant impedance
            zip_load_q        = (1.0, 0.0, 0.0),  # (Z, I, P) — reactive demand, constant impedance
            fault = FaultConfig(contingency_id = 2),
        ),
    ),
)
sys    = load_system(cfg, path_main)
result = run_case!(cfg, sys, path_main, path_folder_results)
```

The case folder layout for this mode:

```
INPUT_FILES/9bus_test_mfile/
  case9.m                 ← steady-state (bus, gen, branch, gencost)
  gen_dynamic_data.csv    ← machine data (mandatory for TSC)
  contingencies.csv       ← SC contingency table (required for SC fault)
```

When results are archived, the `.m` file is copied to `RESULTS/.../Inputs/` instead of the three
CSVs. `gen_dynamic_data.csv` and `contingencies.csv` are always archived for TSC runs.

> **Note on `baseMVA`:** if the `.m` file's `baseMVA` differs from `RunConfig.base_MVA`, a warning
> is printed but `RunConfig.base_MVA` wins. Case9 uses 100 MVA throughout, so this is usually a
> non-issue.

> **Load scaling:** MATPOWER mirror files may carry operating-level loads already. If so, set
> `load_factor = 1.0` to avoid double-scaling.

#### 6.1.2 `Import_Matpower_Case` — write CSVs to disk then run

To start from a MATPOWER case and generate the in-house CSV files explicitly:

```julia
imp = Import_Matpower_Case("INPUT_FILES/PowerModels/case39.m")
imp.DBUS, imp.DGEN, imp.DCIR     # in-house DataFrames (1..nBUS remapped)
imp.path_names[:pf_inputs]       # RESULTS/Results - <ts>/Inputs/ with the 3 CSVs
```

`Import_Matpower_Case` (`_manage_inputs/functions_2_parse_matpower.jl`) is a standalone parser (no
PowerModels dependency). It maps `mpc.bus`/`mpc.gen`+`mpc.gencost`/`mpc.branch` to the in-house schema
and writes `bus_data.csv` / `generators_data.csv` / `line_data.csv` into a timestamped
`RESULTS/Results - <ts>/Inputs/` folder only — **no** `Dispatch/`, `Bus_Matrices/`, or
`Transient_Stability/` subfolders are created at import time.
Only **polynomial** gen costs (model 2, order ≤ 2) are supported; **dynamic data is not produced** —
add `gen_dynamic_data.csv` (and `contingencies.csv` for SC) yourself to run a TS case. To run a parsed
case through the pipeline, copy the three CSVs into an `INPUT_FILES/<case>/` folder and set
`RunConfig.case`, or use `RunConfig.matpower_file` directly (§6.1.1).

---

### 6.2 `load_data.csv` — several loads on one bus, each with its own ZIP split

Optional, per case. Without it a bus has one aggregate demand carrying the system-wide
`zip_load_p` / `zip_load_q`. With it, that demand is split into named loads and each one
gets its own `(Z, I, P)` shares — the input a sensitivity study on load composition needs.

```
id;bus;p_d;q_d;aZ_p;aI_p;aP_p;aZ_q;aI_q;aP_q
1;5;75;30;1.0;0.0;0.0;1.0;0.0;0.0
2;5;50;20;0.0;0.0;1.0;0.0;0.0;1.0
3;6;90;30;0.5;0.0;0.5;1.0;0.0;0.0
```

| Column | Meaning |
|---|---|
| `id` | Unique load id (any integers; used to report per-load results) |
| `bus` | Bus label as written in `bus_data.csv` — several rows may name the same bus |
| `p_d`, `q_d` | MW / MVAr, **unscaled**, same convention as `bus_data.csv` |
| `aZ_p`, `aI_p`, `aP_p` | Active-demand `(Z, I, P)` shares of *this load*; must sum to 1 |
| `aZ_q`, `aI_q`, `aP_q` | Reactive-demand shares of this load; must sum to 1, independent of the active ones |

Rules enforced at `load_system` (all throw, none silently correct the data):

- the rows of a bus must sum to that bus's `Pd` and `Qd` — the table says how a bus's demand
  splits, it never redefines it;
- each triple sums to 1, the same rule as `zip_load_p`;
- `id` is unique, and `bus` exists in `bus_data.csv`.

`load_factor` scales the table and the bus demand together, so the two stay consistent.

**What honours it.** The FULL_BUS nodal balances, classical and DQ_4TH. Each bus enters its
balance with the demand-weighted mean of the shares of the loads sitting on it, which is exact
rather than an approximation: every load at a bus sees the same trajectory voltage and the same
nominal voltage, so the sum of their ZIP polynomials is the mean polynomial times the bus total.

**What ignores it.** `KRON_REDUCED`, which folds demand into `Y_red` as a constant admittance —
only the Z share could ever survive the Schur elimination, so the run warns and continues, the
same way it warns about a non-default `zip_load_p`. The steady-state ACOPF / DCOPF / ED see only
the bus aggregate, which the table leaves unchanged by construction; the split is a refinement of
the *dynamics*, not of the dispatch.

**GL load trips** still scale a whole bus (`gl_load_bus_ids`, `gl_percent_power`). That scales
every load at the bus by the same factor, so the mix is unchanged and stays consistent with the
table; shedding one individual load is not supported.

The run archives the file it read into `Inputs/load_data.csv`, unscaled, and
`Transient_Stability/dynamic_model_details.txt` records how many loads it used.

