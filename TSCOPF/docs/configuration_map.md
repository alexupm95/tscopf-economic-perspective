# Configuration map — structs

This document tracks **where each input parameter lives**: steady-state dispatch and transient analysis are split into nested `@kwdef` structs on `RunConfig`. (The legacy `config_file.jl` shim layer was **removed** — all configuration now lives in the structs below.)

!!! tip "This page vs. the User guide"
    This page maps *where* a field lives (struct → file). The
    [Parameter reference](parameter_reference.md) explains *what each field does*,
    valid combinations, and worked examples. Start there; come back here when you
    need the exact struct/file for an edit.

## Architecture (current)

```
RunConfig                          ← edit in main.jl / main_loop.jl / tests
├── orchestration (case, solver, I/O, scaling)
├── dispatch::DispatchConfig       ← avenue 1: steady-state OPF builder
└── transient::TransientConfig     ← avenue 2 (required when trans_stab=true)
    ├── simulation::TsSimulationConfig   (δ_tol, timing, f_syn, Δω_0)
    ├── builder::TsBuilderConfig         (TS eq/ineq toggles)
    ├── dyn_model::DynModelConfig
    │   └── fault::FaultConfig
    │       ├── SC  → contingencies.csv (via fault.contingency_id)
    │       ├── GL  → gl_gen_ids / gl_load_bus_ids / gl_percent_power
    │       └── OB  → ob_branch_ids (open branch; no short-circuit)
    └── gen_dynamic_filename
```

**Entry point:** build a `RunConfig` and call `run_case!(cfg, sys, …)`.

---

## 1. `RunConfig` — orchestration (`engine.jl`)

Shared knobs for dispatch-only and TSC runs. Does **not** carry `type_model`, `δ_tol_deg`, or `contingency_id` directly.

| Parameter | Field | Default | Notes |
|---|---|---|---|
| Enable transient stability | `trans_stab` | `false` | `true` requires `transient::TransientConfig` |
| Case folder name | `case` | `"9bus"` | `INPUT_FILES/<case>/` |
| Power base [MVA] | `base_MVA` | `100.0` | |
| Demand scale factor | `load_factor` | `1.5` | Scales `p_d`, `q_d` in memory |
| Optimiser | `solver_name` | `"Ipopt"` | `"Ipopt"`, `"Ipopt-ma57"`, `"Ipopt-ma97"`, `"Ipopt-pardiso"`, `"MadNLP"`, `"MadNLP-ma57"`, `"MadNLP-ma97"`, `"Gurobi"`, `"HiGHS"` — coherence checked vs `dispatch.type_model` |
| Suppress solver console | `silent_solver` | `false` | Forces Ipopt `print_level = 0` when true |
| Solver time limit [s] | `time_limit_sec` | `600.0` | |
| Ipopt NLP options | `ipopt` | `IpoptSolverConfig()` | §1.1 |
| HiGHS LP/QP options | `highs` | `HiGHSSolverConfig()` | §1.2 |
| Gurobi LP/MILP/NL options | `gurobi` | `GurobiSolverConfig()` | §1.3 |
| MadNLP NLP options | `madnlp` | `MadNLPSolverConfig()` | §1.4 |
| Overwrite vs timestamped results | `overwrite_results` | `false` | |
| Save dual variables (primal JuMP) | `save_duals` | `true` | Creates `Dispatch/CSV_duals/`; with `trans_stab`, also `Transient_Stability/CSV_duals/`. Files written on successful solve. |
| Save Jacobian / Hessian / gradient | `save_optim_matrices` | `false` | Sparse COO CSV in `RESULTS/.../Optim_Matrices/`. Steady-state only; auto-forced `false` for TSC runs via `resolve_save_optim_matrices`. |
| Admittance XLSX dumps | `save_matrices` | `true` | Steady-state Ybus/Bbus (ACOPF/DCOPF only; not ED/UC). TSC fault/post-fault matrices always when `trans_stab=true`. |
| TSC trajectory figures | `save_ts_plots` | `false` | `Transient_Stability/Figures/` (Plots.jl); requires `trans_stab=true` |
| TS diagnostic dumps | `save_ts_debug_csv` | `false` | `Transient_Stability/CSV/Debug/`: `gfm_{filter,limiter,voltage}_debug.csv` + `swing_debug.csv`, one row per (window, gen, step) with value, margin and dual |
| Pre-solve dispatch archive | `save_warmstart_dispatch` | `false` | `Dispatch_WarmStart/` on the two TSC paths that pre-solve a steady-state OPF: FULL_BUS TSC-ACOPF (warm start) and Kron TSC-DCOPF (δ_ref anchor, plus `CSV/delta_ref.csv`). Throws elsewhere. |
| Run-script archive | `run_script` | `nothing` | `nothing` = off. Set to `@__FILE__` in the run script to copy that `.jl` file, byte for byte, into the run-folder root next to `input_parameters.txt`. Any readable path is accepted; a non-existent one throws in `validate_run_config!`. |
| Steady-state builder config | `dispatch` | `DispatchConfig()` | Avenue 1 |
| Transient bundle | `transient` | `nothing` | Avenue 2; required when `trans_stab=true` |
| MATPOWER input file | `matpower_file` | `nothing` | `nothing` = CSV mode (default). Set to `"case9.m"` (relative to case folder) or an absolute path to load steady-state data from a MATPOWER `.m` file instead of the three CSVs. Dynamic data (`gen_dynamic_data.csv`) is always read from the case folder regardless. See §6.1 of the user guide. |

**Helpers:** `reconfigure(cfg; kwargs…)`, `validate_run_config!`.

---

## 2. `DispatchConfig` — steady-state OPF (`_common/DispatchConfig.jl`)

User-facing dispatch configuration.

| Parameter | Field | Default | Notes |
|---|---|---|---|
| Steady-state formulation | `type_model` | `"ACOPF"` | `"ACOPF"`, `"DCOPF"`, `"ED"`, `"UC"` |
| Matrix-form power balance | `use_matrix` | `true` | Ybus/Bbus |
| Fuel cost form | `cost_type` | `"quadratic"` | `"quadratic"` or `"linear"` |
| Branch susceptance model (DC) | `susceptance_model` | `SIMPLE` | `@enum SusceptanceModel`: `SIMPLE` = `1/x`; `POWERMODELS` = `imag(inv(r+jx))` = `x/(r²+x²)`. **DC-OPF only** (ignored by AC/ED); the two agree when `r=0` |
| Demand as a pinned variable | `pd_as_variable` | `false` | **ACOPF only.** `P_d[k]` becomes a variable held by `P_d[k] - p_d⁰[k] == 0`, whose dual σ_k is `∂L/∂P_d(k)` directly. Solution-invariant; adds one variable and one equality per bus; exports `dual_Pd_pin.csv`. Active power only |
| Explicit ≤ bounds on V, θ, P_g, Q_g, branch flows | `bound_*` | see source | Dual-friendly when `true` |
| Box limit encoding | `bound_encoding` | `CONSTRAINT` | `CONSTRAINT` = explicit `@constraint` inequalities; `VARIABLE` = JuMP `lower_bound`/`upper_bound` on `@variable` (fewer rows; dual export via `bound_manifest`). UC ignores this field. |
| Optional inequality families | `ineq_*` | see source | branch thermal, angle-diff, gen capability |
| Explicit dispatch dual LP | `solve_explicit_dual` | `false` | Requires `ED` or `DCOPF` + `linear` cost; DCOPF also needs `use_matrix=true`; results in `Dispatch_Dual/`. Quadratic primal → use `save_duals` on primal instead. |
| Non-case-data limit values | `limits::DispatchLimitsConfig` | `θ_min/max_rad = ∓π`, `ang_diff_clamp_ac_deg = 60`, `ang_diff_clamp_dc_deg = 30` | θ box and angle-diff safety clamps; branch thermal uses `DCIR.l_cap_1` only when `ineq_sbranch_upper` (unrated `l_cap_1 == 0` skipped). DC dual LP reads the same struct (user guide §2.2). Full default-state inventory: `OUTPUTS/reference/constraint_bound_defaults.md`. |

Mandatory equalities (not in `DispatchConfig`): slack angle, P/Q balance, and branch-flow
eq when `use_matrix=false` — see `Make_*OPF_Model*!` in `_acopf/` / `_dcopf/`.

**Helpers:** `build_opf_input_param(dispatch)`, `copy_dispatch_config`, `validate_dispatch_config!`.

---

## 3. `TransientConfig` — transient bundle (`_transient_stability/TsConfig.jl`)

Required when `RunConfig.trans_stab = true`.

| Parameter | Field | Default | Notes |
|---|---|---|---|
| Simulation timing & limits | `simulation` | `TsSimulationConfig()` | See §3.1 |
| TS builder toggles | `builder` | `TsBuilderConfig()` | See §3.2 |
| Physics / network / disturbance | `dyn_model` | `DynModelConfig()` | See §4 |
| Dynamic generator CSV (SG) | `gen_dynamic_filename` | `"gen_dynamic_data.csv"` | In case folder |
| GFM dynamic CSV | `gfm_dynamic_filename` | `"gfm_dynamic_data.csv"` | Used when `allow_gfm=true` |

**Helpers:** `default_transient_config()`, `copy_transient_config`, `reconfigure_transient`, `validate_transient_config!`.

### 3.1 `TsSimulationConfig`

Replaces `Time_Parameters_4_TS_*` and the `δ_tol_deg` argument to `Common_Parameters_4_TS`.

| Parameter | Field | Default | Notes |
|---|---|---|---|
| Max \|δ − δ_COI\| [deg] | `δ_tol_deg` | `90.0` | Symmetric half-width when `δ_tol_deg_lower` / `δ_tol_deg_upper` unset |
| Below δ_COI [deg] | `δ_tol_deg_lower` | `nothing` | Optional; defaults to `δ_tol_deg` |
| Above δ_COI [deg] | `δ_tol_deg_upper` | `nothing` | Optional; defaults to `δ_tol_deg` |
| Max \|δᵢ − δⱼ\| [deg] | `δ_pair_tol_deg` | `180.0` | Symmetric; read only by `bound_style_δ = :pairwise` |
| Simulation start [s] | `t_start_sim` | `0.0` | |
| Simulation end [s] | `t_end_sim` | `5.0` | |
| Time step [s] | `t_step` | `0.01` | |
| Fault application time [s] | `t_start_fault` | `0.01` | |
| SC clearing duration [s] | `clearing_time` | `0.3` | SC only |
| Synchronous frequency [Hz] | `f_syn` | `50.0` | |
| Initial speed deviation [p.u.] | `Δω_0` | `0.0` | |

**Helpers:** `common_ts_parameters(simulation)`, `time_windows_sc(simulation)`, `time_windows_gld(simulation)`.

**Legacy shims:**

- `Time_Parameters_4_TS_SC()` → `time_windows_sc(TsSimulationConfig())`
- `Time_Parameters_4_TS_GLD()` → `time_windows_gld(TsSimulationConfig())`
- `Common_Parameters_4_TS(; δ_tol_deg)` → `common_ts_parameters(TsSimulationConfig(δ_tol_deg=…))`

### 3.2 `TsBuilderConfig`

Replaces `Define_Inputs_4TS_Optim()` as the user-facing TS builder configuration.

| Group | Fields | Default pattern | Meaning |
|---|---|---|---|
| Variable bounds | `bound_E`, `bound_δ`, `bound_P_m`, `bound_*_tf`, `bound_*_tpf` | pre-fault `false`; tf/tpf `false` | Explicit ≤ inequalities when `true` (all wired) |
| Box limit encoding | `bound_encoding` | `CONSTRAINT` | Same semantics as `DispatchConfig.bound_encoding`; δ-COI stability and Kron `Qe` expressions always use `CONSTRAINT`. |
| Bound limit values | `limits::TsBoundLimitsConfig` | physical E/δ/P_m; tf/tpf default `(-Inf, Inf)` | Numbers used by bound inequalities; see below |
| DQ / AVR / governor optional bounds | `bound_Ed` … `bound_Pm_tpf` | all `false` | Wired on DQ / AVR / governor paths only |
| δ-COI stability | `ineq_δ_COI_tf_lower`, `ineq_δ_COI_tf_upper`, `ineq_δ_COI_tpf_*` | all `true` | Optional box constraints w.r.t. COI (wired) |

**`TsBoundLimitsConfig`:** `E_min_pu`/`E_max_pu`, `δ_min_rad`/`δ_max_rad`, `P_m_source` (`:dgen_pg_limits`), `gov_valve_min_pu`/`gov_valve_max_source` (`:dgen_pg_limits`), `V_bus_min_pu` (FullBus forced lower on bus voltages), DQ `Ed/Eq/Id/Iq` and AVR `V_ref` scalars, and `δ_tf`, `Δω_tf`, … as `TsBoundLimitPair` (default `-Inf`/`Inf` = inactive; CONSTRAINT skips ≤-rows, VARIABLE omits JuMP bounds). Plus the GFM box knobs `gfm_δ_*`, `gfm_V_meas_*`, `gfm_E_raw_extra_pu`, `gfm_E_raw_slack_pu`, `gfm_E_clip_slack_pu`, `gfm_PQ_bound_scale`, `gfm_PQ_bound_offset_pu`, `gfm_P_meas_floor_pu`, `gfm_Q_meas_floor_pu`, `gfm_I_floor_pu`, `gfm_I_ceiling_pu` — margins and floors that `resolve_gfm_bound_limits` combines with the `DGFM` row (read only when `allow_gfm`). Their toggles are `TsBuilderConfig.bound_gfm_δ` and `bound_gfm_{P_meas,Q_meas,V_meas,E_int_raw,E_int,E_droop_raw,E_droop,Id,Iq}_{tf,tpf}`, all defaulting `true` (the SG `bound_*` families default `false`). Full default-state inventory: `OUTPUTS/reference/constraint_bound_defaults.md`.

Pre-fault init equalities (`eq_const_P_init` / `Q_init` / `Pm_init`) and fault/post-fault equalities (swing, Pe, COI, KCL) are **mandatory physics** — not in `TsBuilderConfig`. Path rules: `Q_init` is omitted on Kron linear / TSC-DCOPF; `Pm_init` only when `mech_power_mode = USE_PM`. See `OUTPUTS/reference/constraint_toggle_map.md`.

**Helper:** `build_ts_input_param(builder)`.

**Legacy shim:** `Define_Inputs_4TS_Optim()` → `build_ts_input_param(TsBuilderConfig())`.

---

## 4. `DynModelConfig` — physics and network form

(`_transient_stability/DynModelConfig.jl`; nested in `TransientConfig.dyn_model`)

| Parameter | Field | Default | Notes |
|---|---|---|---|
| Generator order | `gen_order` | `CLASSICAL_2ND` | `DQ_4TH` = 4th-order dq (FULL_BUS) |
| Network form | `network_form` | `KRON_REDUCED` | `FULL_BUS` for full Ybus TSC-ACOPF |
| Mechanical power in swing | `mech_power_mode` | `USE_PG` | `USE_PM` required for FULL_BUS |
| Include AVR | `include_avr` | `false` | SEXS exciter (lead-lag + gain-lag); requires `DQ_4TH` + full CSV (`T_exc,K_exc,Ta_exc,Tb_exc`; `0;0` bypasses lead-lag) |
| Include governor | `include_governor` | `false` | TGOV1 turbine governor; requires `USE_PM` + `FULL_BUS` (classical or DQ_4TH) |
| Governor valve limiter | `governor_limiter` | `GOV_NO_LIMIT` | `GOV_NO_LIMIT` / `GOV_SMOOTH` / `GOV_HARD_BOUND` (only when `include_governor`) |
| Allow GFM fleet | `allow_gfm` | `false` | Separate `gfm_dynamic_data.csv`; requires `DQ_4TH` + `FULL_BUS` (G2: full transient GFM) |
| GFM converters in the COI | `gfm_in_coi` | `false` | Weights converters into both COI means with `H_eff = Tf/(2*mp)` and gives them `:coi_box` corridor rows under the same tolerances. Governs `bound_style_δ = :coi_box` and `bound_style_Δω = :coi_box` only; inert (warns) without `allow_gfm` or without a COI-referenced corridor |
| GFM P/Q overload PIs | `gfm_pq_overload` | `false` | REGFM_A1 Fig. 3. Off, droop and the Q–V voltage error are unchanged. On, both branch pairs (eight states per converter-step). Integrators follow `gfm_integrator`. Warns when `allow_gfm` is false |
| ZIP load split `(Z,I,P)` — active | `zip_load_p` | `(1,0,0)` | impedance / current / power fractions; must sum to 1; default = constant impedance |
| ZIP load split `(Z,I,P)` — reactive | `zip_load_q` | `(1,0,0)` | independent of `zip_load_p`; must sum to 1 (e.g. REE: `zip_load_p=(0,1,0)`, `zip_load_q=(1,0,0)`) |
| Build the δ corridor | `constrain_δ` | `true` | At least one of `constrain_δ` / `constrain_Δω` must be true |
| δ corridor reference | `bound_style_δ` | `:swing_propagated` | `:coi_box` \| `:highest_H` \| `:ref_gen` \| `:pairwise`; USE_PM / DQ_4TH reject `:swing_propagated`; all build on every network form. The machine-referenced and pairwise styles bound GFM converters too; `:coi_box` does so under `gfm_in_coi`, and `:swing_propagated` never does |
| δ reference machine id | `δ_ref_gen_id` | `nothing` | Required by `:ref_gen`; validated against the system data before the warm start. May name a GFM unit — `:highest_H` may not, since it ranks by `H` |
| Build the Δω corridor | `constrain_Δω` | `false` | |
| Δω corridor reference | `bound_style_Δω` | `:coi_box` | `:abs` bounds the raw Δω, forms no COI, and spans GFM converters; `:coi_box` is SG-only unless `gfm_in_coi=true` |
| Δω tolerance [p.u.] | `Δω_tol_pu` | `0.5` | Symmetric half-width when lower/upper unset; requires `constrain_Δω=true` |
| Below reference [p.u.] | `Δω_tol_pu_lower` | `nothing` | Optional; defaults to `Δω_tol_pu` |
| Above reference [p.u.] | `Δω_tol_pu_upper` | `nothing` | Optional; defaults to `Δω_tol_pu` |
| Disturbance spec | `fault` | `FaultConfig()` | See §5 |

**Validation:** `validate_dyn_config!(cfg::RunConfig)` (reads `cfg.transient.dyn_model` when `trans_stab=true`), plus `validate_δ_reference!(cfg, DGEN, DGEN_DYN, DGFM)` for the machine-referenced and pairwise δ styles, which needs the system data and so runs from `run_case!`.

---

## 5. `FaultConfig` — disturbances

(`_transient_stability/FaultConfig.jl`; nested in `TransientConfig.dyn_model.fault`)

| Parameter | Field | Default | Applies when |
|---|---|---|---|
| Disturbance class | `fault_type` | `SC` | `SC`, `GL`, or `OB` |
| Contingency table row | `contingency_id` | `2` | `fault_type == SC` |
| SC fault admittance [pu] | `sc_fault_admittance_pu` | `0 - im*FAULT_BUS_SHUNT` | `fault_type == SC`. Default inductive solid short (`FAULT_BUS_SHUNT = 1e6`). Resistive option: `FAULT_BUS_SHUNT + 0im`. `real(y) ≥ 0`, nonzero, finite. `imag(y) > 0` warns (capacitive) |
| SC diagonal stamp | `sc_fault_stamp` | `:add` | `fault_type == SC`. `:add` adds the shunt to `Y_kk`; `:overwrite` replaces it (the resistive option) |
| Generators to trip | `gl_gen_ids` | `Int[]` | `fault_type == GL` (gen trip) |
| Load buses to modify | `gl_load_bus_ids` | `Int[]` | `fault_type == GL` (load trip) |
| Load scaling factor α | `gl_percent_power` | `Float64[]` | parallel to `gl_load_bus_ids`; $p_d^{\mathrm{new}} = p_d(1+\alpha)$, same for $q_d$; $\alpha \ge -1$ |
| Branches to open | `ob_branch_ids` | `Int[]` | `fault_type == OB` (DCIR row indices; lines/transformers) |

**SC details** are read from `INPUT_FILES/<case>/contingencies.csv` via `build_fault_details(fault, path)` (also used by `Print_Input_Parameters` for run archival).

**OB (Open Branch)** uses the same single-window timeline as GL (no fault-on short-circuit). Validation rejects opens that island the network or isolate the slack. Implemented on classical FULL_BUS, Kron, and DQ FULL_BUS.

---

## 6. `config_file.jl` — removed

The legacy shim layer (`Define_Initial_Parameters`, `Define_Inputs_4OPF`, `Define_Inputs_4TS_Optim`, `Time_Parameters_4_TS_*`, `Common_Parameters_4_TS`, `Define_Fault_Details`, …) has been **deleted**. Use the structs and their builders directly:

| Old function | Replacement |
|---|---|
| `Define_Inputs_4OPF()` | `build_opf_input_param(default_dispatch_config())` |
| `Define_Inputs_4TS_Optim()` | `build_ts_input_param(TsBuilderConfig())` |
| `Time_Parameters_4_TS_SC()` | `time_windows_sc(TsSimulationConfig())` |
| `Time_Parameters_4_TS_GLD()` | `time_windows_gld(TsSimulationConfig())` |
| `Common_Parameters_4_TS(; δ_tol_deg)` | `common_ts_parameters(TsSimulationConfig(δ_tol_deg=…))` |
| `Define_Fault_Details` | `build_fault_details(fault, path)` |

**Add all new user parameters to the structs** (`DispatchConfig`, `TsSimulationConfig`, `TsBuilderConfig`, `DynModelConfig`).

---

## 7. Network data (CSV — not in structs)

Read from `INPUT_FILES/<case>/` by `load_system`. Scaled in memory by `RunConfig.load_factor`.

| File | Role |
|---|---|
| `bus_data.csv` | Buses: demand, limits, type (slack=3) |
| `generators_data.csv` | Generators: costs, limits, status |
| `line_data.csv` | Branches: impedance, ratings, status |
| `gen_dynamic_data.csv` | Machine parameters: `Xd_tr`, `H`, `D`, … |
| `contingencies.csv` | SC fault location + branch to disconnect |

**MATPOWER import:** a standalone parser ingests a MATPOWER `.m` case into this schema —
`Import_Matpower_Case(path_m)` (`_manage_inputs/functions_2_parse_matpower.jl`) returns the in-house
`DBUS`/`DGEN`/`DCIR` DataFrames and writes the three steady-state CSVs into
`RESULTS/Results - <ts>/Inputs/` only (no dispatch/TSC subfolders). No PowerModels dependency;
`mpc.gencost` (polynomial, model 2) maps to `c2,c1,c0`. Dynamic data (`gen_dynamic_data.csv`) is
**not** produced — author it separately for TS runs.

**Results folder layout:** conditional subfolders are chosen by `build_results_paths(cfg)` from
`RunConfig` flags (`save_duals`, `save_matrices`, `save_ts_plots`, `save_warmstart_dispatch`, `trans_stab`,
`dispatch.solve_explicit_dual`, `dispatch.type_model`). See [Running a case](running_a_case.md) §8.

---

## 8. File cross-reference

| Concern | Primary file |
|---|---|
| Run orchestration | `engine.jl` |
| Steady-state builder config | `_common/DispatchConfig.jl` |
| Branch susceptance variants (DC) | `_common/functions_4_admittance_matrices.jl` (`dc_branch_susceptance`, `Calculate_Matrix_B`, `Calculate_Matrix_B_PowerModels`) |
| MATPOWER `.m` import | `_manage_inputs/functions_2_parse_matpower.jl` |
| Transient simulation + builder | `_transient_stability/TsConfig.jl` |
| Dynamic physics selection | `_transient_stability/DynModelConfig.jl`, `AbstractDynamicGenModel.jl` |
| Disturbances | `_transient_stability/FaultConfig.jl` |
| Include order | `src/TSCOPF.jl` |
| Single run | `main.jl` |
| δ_tol sweep | `main_loop.jl` |
| Solver / model coherence | `_common/functions_4_sanity_checks.jl` |
| Optimiser factory + log paths + backend release | `_manage_inputs/functions_2_setup_optim.jl` (`Setup_Optim_Model`, `set_solver_log_path!`, `release_solver_backend!`) |
| Ipopt NLP / HSL / PARDISO options | `_manage_inputs/IpoptSolverConfig.jl`, `_common/ipopt_hsl.jl` |
| MadNLP HSL linear solvers | `_common/madnlp_linear_solvers.jl` |
| HiGHS / Gurobi / MadNLP options | `_manage_inputs/HiGHSSolverConfig.jl`, `GurobiSolverConfig.jl`, `MadNLPSolverConfig.jl` |
| Numerical constants (`FAULT_BUS_SHUNT`) | `_common/constants.jl` |
| Results path layout (`build_results_paths`, `results_folder_keys`) | `_common/auxiliar_functions.jl` |
| Archived / not-loaded experimental code | `_archive/` |
