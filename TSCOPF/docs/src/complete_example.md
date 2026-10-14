# Complete example — every parameter

The examples in [Running a case](running_a_case.md) are deliberately partial: each sets the handful of fields that matter for one scenario and leaves everything else at its default. This page is the opposite. It writes out **one run with every user-facing field of every configuration struct set explicitly**, so you can copy it, delete what you do not need, and know that nothing is hiding behind a default you never saw.

The scenario is the most demanding path the package supports today: the 39-bus case of the paper (`39bus_gfm`) with a mixed synchronous / grid-forming fleet, 4th-order machines on the full nodal network, AVR and TGOV1 governor active, smooth valve saturation, and an ACOPF warm start.

!!! warning "Struct-verified, not run-verified"
    Every field name, type, and default on this page was checked against the struct definitions in `src/`. The configuration below was **not executed**, so no objective value, solve time, or dual is quoted. Treat it as a syntactically and semantically valid configuration, not as a benchmark result. If you run it, expect a large NLP: 500 time steps × two windows × (4 machine states + 2 control states + GFM states) per unit.

!!! note "One import gotcha before you start"
    `BoundEncoding` is exported as a *type*, but its two values are not. Write `TSCOPF.CONSTRAINT` / `TSCOPF.VARIABLE`, not bare `CONSTRAINT` / `VARIABLE`. Every other enum used below (`DQ_4TH`, `FULL_BUS`, `USE_PM`, `GOV_SMOOTH`, `SC`, `SIMPLE`, …) *is* exported and can be written bare.

---

## The configuration

```julia
using TSCOPF
using JuMP

const PATH_MAIN    = @__DIR__
const PATH_RESULTS = joinpath(PATH_MAIN, "RESULTS")
mkpath(PATH_RESULTS)

# Optional boxes that should stay inactive: (-Inf, Inf) is the documented "off"
# value, not a placeholder. Under either encoding no row and no JuMP bound is
# created for an infinite side.
const OPEN = TsBoundLimitPair(min = -Inf, max = Inf)

# ── 1. Steady-state limits that are not case data ────────────────────────────
const dispatch_limits = DispatchLimitsConfig(
    θ_min_rad             = -π,      # bus angle box, only used when bound_θ = true
    θ_max_rad             =  π,
    ang_diff_clamp_ac_deg = 60.0,    # safety clamp on DCIR.ang_min/ang_max (AC)
    ang_diff_clamp_dc_deg = 30.0,    # same for the DC angle-difference rows
)

# ── 2. Steady-state OPF builder ──────────────────────────────────────────────
const dispatch_cfg = DispatchConfig(
    type_model           = "ACOPF",           # ACOPF | DCOPF | ED | UC
    use_matrix           = true,              # Ybus/Bbus form of the power balance
    cost_type            = "quadratic",       # quadratic (MATPOWER) | linear
    susceptance_model    = SIMPLE,            # DC-OPF only; ignored by ACOPF
    bound_V              = true,
    bound_θ              = true,
    bound_P_g            = true,
    bound_Q_g            = true,
    bound_P_ik           = false,             # per-direction branch flow boxes
    bound_Q_ik           = false,
    bound_P_ki           = false,
    bound_Q_ki           = false,
    ineq_sg_upper        = false,             # generator capability circle
    ineq_sbranch_upper   = true,              # branch thermal limits (DCIR.l_cap_1)
    ineq_ang_diff_branch = true,              # per-branch angle-difference limits
    pd_as_variable       = false,             # ACOPF only: demand as a pinned variable (dual = ∂L/∂P_d)
    solve_explicit_dual  = false,             # ED/DCOPF + linear cost only
    limits               = dispatch_limits,
    bound_encoding       = TSCOPF.CONSTRAINT, # CONSTRAINT (dual-friendly) | VARIABLE
)

# ── 3. Transient timing and the stability tolerance ──────────────────────────
const ts_simulation = TsSimulationConfig(
    δ_tol_deg       = 100.0,   # symmetric |δ − δ_COI| corridor [deg]
    δ_tol_deg_lower = nothing, # positive magnitude below COI; nothing → δ_tol_deg
    δ_tol_deg_upper = nothing, # positive magnitude above COI; nothing → δ_tol_deg
    δ_pair_tol_deg  = 180.0,   # |δᵢ − δⱼ| limit, read only by bound_style_δ = :pairwise
    t_start_sim     = 0.0,
    t_end_sim       = 5.0,     # end of the post-fault window [s]
    t_step          = 0.01,    # trapezoidal Δt [s]
    t_start_fault   = 0.01,    # disturbance instant [s]
    clearing_time   = 0.15,    # SC only: fault-on duration [s]; the fault lasts clearing_time + t_step
    f_syn           = 50.0,    # synchronous frequency [Hz] → ω_s = 2π f
    Δω_0            = 0.0,     # initial speed deviation [p.u.]
)

# ── 4. Numeric values behind the optional transient boxes ────────────────────
const ts_limits = TsBoundLimitsConfig(
    # Classical / shared physical limits
    E_min_pu             = 0.0,
    E_max_pu             = 4.0,    # also the AVR field-voltage clamp — see note below
    δ_min_rad            = -π,
    δ_max_rad            =  π,
    P_m_source           = :dgen_pg_limits,   # only accepted value
    V_bus_min_pu         = 0.0,    # forced lower bound on V_tf / V_tpf (FULL_BUS)
    # Governor valve limits (read when include_governor and limiter ≠ GOV_NO_LIMIT)
    gov_valve_min_pu     = 0.0,
    gov_valve_max_source = :dgen_pg_limits,   # only accepted value
    # Optional classical / shared trajectory boxes — inactive
    δ_tf = OPEN, Δω_tf = OPEN, Pe_tf = OPEN, Qe_tf = OPEN, δCOI_tf = OPEN,
    δ_tpf = OPEN, Δω_tpf = OPEN, Pe_tpf = OPEN, Qe_tpf = OPEN, δCOI_tpf = OPEN,
    # DQ pre-fault algebraic states (finite defaults; only used if the toggle is on)
    Ed_min_pu = -2.0, Ed_max_pu = 2.0,
    Eq_min_pu = -2.0, Eq_max_pu = 2.0,
    Id_min_pu = -5.0, Id_max_pu = 5.0,
    Iq_min_pu = -5.0, Iq_max_pu = 5.0,
    # AVR / governor pre-fault set-points
    V_ref_min_pu = 0.0,
    V_ref_max_pu = Inf,            # Inf ⇒ upper side never materialised
    P_ref_source = :dgen_pg_limits,
    # DQ machine-state boxes per window — inactive
    Ed_tf = OPEN, Eq_tf = OPEN, Id_tf = OPEN, Iq_tf = OPEN, Te_tf = OPEN,
    Ed_tpf = OPEN, Eq_tpf = OPEN, Id_tpf = OPEN, Iq_tpf = OPEN, Te_tpf = OPEN,
    # AVR / governor research boxes — inactive (saturation physics is unaffected)
    E_fd_unlim_tf = OPEN, E_fd_unlim_tpf = OPEN,
    Pv_tf = OPEN, Pv_tpf = OPEN,
    Pm_tf = OPEN, Pm_tpf = OPEN,
    # GFM guard-rail margins and floors (read only when allow_gfm)
    gfm_δ_min_rad = -π, gfm_δ_max_rad = π,
    gfm_V_meas_min_pu = 0.0, gfm_V_meas_max_pu = 2.5,
    gfm_E_raw_extra_pu = 0.25, gfm_E_raw_slack_pu = 0.5, gfm_E_clip_slack_pu = 0.05,
    gfm_PQ_bound_scale = 1.5, gfm_PQ_bound_offset_pu = 0.05,
    gfm_P_meas_floor_pu = 3.0, gfm_Q_meas_floor_pu = 3.0,
    gfm_I_floor_pu = 5.0, gfm_I_ceiling_pu = 20.0,
    gfm_pq_raw_span_pu = 20.0, gfm_pq_clip_slack_pu = 0.05,
)

# ── 5. Which transient constraints and boxes get built ───────────────────────
# The pre-fault init equalities (P/Q/Pm links) and the window equalities (swing,
# Pe, COI, nodal balance) are always built by the path — they are physics, not
# toggles. Only optional boxes and the δ–COI corridor are switchable here.
const ts_builder = TsBuilderConfig(
    # Pre-fault classical boxes
    bound_E   = true,
    bound_δ   = true,
    bound_P_m = true,
    # Optional classical trajectory boxes
    bound_δ_tf = false, bound_Δω_tf = false, bound_Pe_tf = false,
    bound_Qe_tf = false, bound_δCOI_tf = false,
    bound_δ_tpf = false, bound_Δω_tpf = false, bound_Pe_tpf = false,
    bound_Qe_tpf = false, bound_δCOI_tpf = false,
    # DQ pre-fault algebraic boxes
    bound_Ed = false, bound_Eq = false, bound_Id = false, bound_Iq = false,
    # AVR / governor pre-fault set-point boxes
    bound_V_ref = true,
    bound_P_ref = false,
    # DQ machine-state boxes per window
    bound_Ed_tf = false, bound_Eq_tf = false, bound_Id_tf = false,
    bound_Iq_tf = false, bound_Te_tf = false,
    bound_Ed_tpf = false, bound_Eq_tpf = false, bound_Id_tpf = false,
    bound_Iq_tpf = false, bound_Te_tpf = false,
    # AVR / governor research boxes
    bound_E_fd_unlim_tf = false, bound_E_fd_unlim_tpf = false,
    bound_Pv_tf = false, bound_Pv_tpf = false,
    bound_Pm_tf = false, bound_Pm_tpf = false,
    # GFM guard-rail boxes (default true; read only when allow_gfm)
    bound_gfm_δ = true,
    bound_gfm_P_meas_tf = true, bound_gfm_Q_meas_tf = true, bound_gfm_V_meas_tf = true,
    bound_gfm_E_int_raw_tf = true, bound_gfm_E_int_tf = true,
    bound_gfm_E_droop_raw_tf = true, bound_gfm_E_droop_tf = true,
    bound_gfm_Id_tf = true, bound_gfm_Iq_tf = true,
    bound_gfm_P_meas_tpf = true, bound_gfm_Q_meas_tpf = true, bound_gfm_V_meas_tpf = true,
    bound_gfm_E_int_raw_tpf = true, bound_gfm_E_int_tpf = true,
    bound_gfm_E_droop_raw_tpf = true, bound_gfm_E_droop_tpf = true,
    bound_gfm_Id_tpf = true, bound_gfm_Iq_tpf = true,
    bound_gfm_pq_raw_tf = true, bound_gfm_pq_raw_tpf = true,
    bound_gfm_pq_clip_tf = true, bound_gfm_pq_clip_tpf = true,
    # The δ–COI stability corridor itself (this is the TSC in TSC-OPF)
    ineq_δ_COI_tf_lower  = true,
    ineq_δ_COI_tf_upper  = true,
    ineq_δ_COI_tpf_lower = true,
    ineq_δ_COI_tpf_upper = true,
    limits         = ts_limits,
    bound_encoding = TSCOPF.CONSTRAINT,
)

# ── 6. The disturbance ───────────────────────────────────────────────────────
# contingencies.csv row 57 of 39bus_gfm: three-phase fault at bus 23, trip circuit 29 (23–24).
const fault_cfg = FaultConfig(
    fault_type       = SC,          # SC | GL | OB
    contingency_id   = 57,          # SC only: row of contingencies.csv
    sc_fault_admittance_pu = -im * TSCOPF.FAULT_BUS_SHUNT,  # SC only: inductive solid short
    sc_fault_stamp   = :add,        # SC only: :add to Y_kk, or :overwrite
    gl_gen_ids       = Int[],       # GL only: DGEN row indices to trip
    gl_load_bus_ids  = Int[],       # GL only: buses whose demand is scaled
    gl_percent_power = Float64[],   # GL only: α per bus, p_d ← p_d(1+α)
    ob_branch_ids    = Int[],       # OB only: DCIR row indices to open
)

# ── 7. Physics and network form ──────────────────────────────────────────────
const dyn_model_cfg = DynModelConfig(
    gen_order               = DQ_4TH,          # CLASSICAL_2ND | DQ_4TH
    network_form            = FULL_BUS,        # KRON_REDUCED | FULL_BUS
    mech_power_mode         = USE_PM,          # USE_PG | USE_PM
    dq_speed_dev_in_algebra = true,            # keep (1+Δω) in the stator/Pe algebra
    include_avr             = true,            # DQ_4TH only
    include_governor        = true,            # FULL_BUS only, needs USE_PM
    governor_limiter        = GOV_SMOOTH,      # GOV_NO_LIMIT | GOV_SMOOTH | GOV_HARD_BOUND
    allow_gfm               = true,            # DQ_4TH + FULL_BUS only
    gfm_in_coi              = false,           # converters in the COI + :coi_box corridors
    gfm_pq_overload         = false,           # REGFM_A1 P/Q overload PIs
    allow_all_gfm           = false,           # opt-in for a converter-only fleet
    ode_first_step          = :trapezoidal,    # :trapezoidal | :backward_euler
    gfm_integrator          = :backward_euler, # GFM filters + Q–V PI only
    zip_load_p              = (1.0, 0.0, 0.0), # ACTIVE demand (Z, I, P), sums to 1
    zip_load_q              = (1.0, 0.0, 0.0), # REACTIVE demand (Z, I, P), sums to 1
    constrain_δ             = true,            # build the rotor-angle corridor
    bound_style_δ           = :coi_box,        # :coi_box | :swing_propagated | :highest_H | :ref_gen | :pairwise
    δ_ref_gen_id            = nothing,         # required by :ref_gen only
    constrain_Δω            = false,           # optional speed corridor
    bound_style_Δω          = :coi_box,        # :coi_box | :abs
    Δω_tol_pu               = 0.5,             # symmetric half-width [p.u.]
    Δω_tol_pu_lower         = nothing,         # positive magnitude below COI
    Δω_tol_pu_upper         = nothing,         # positive magnitude above COI
    fault                   = fault_cfg,
)

# ── 8. The transient bundle ──────────────────────────────────────────────────
const transient_cfg = TransientConfig(
    simulation           = ts_simulation,
    builder              = ts_builder,
    dyn_model            = dyn_model_cfg,
    gen_dynamic_filename = "gen_dynamic_data_full.csv",  # DQ needs the extended file
    gfm_dynamic_filename = "gfm_dynamic_data.csv",       # read only when allow_gfm
    extra_contingency_ids = Int[],                       # further SC contingencies in the same solve
)

# ── 9. Solver option blocks (all four; only the selected one is applied) ─────
const ipopt_cfg = IpoptSolverConfig(
    tol                        = 1e-8,
    print_level                = 5,       # 0…12; silent_solver = true forces 0
    max_iter                   = 3_000,
    constr_viol_tol            = 1e-8,
    dual_inf_tol               = 1e-8,
    compl_inf_tol              = 1e-8,
    hessian_approximation      = "exact", # "exact" | "limited-memory"
    limited_memory_max_history = 50,      # only read for "limited-memory"
    acceptable_tol             = nothing, # nothing ⇒ leave the Ipopt default
    obj_scaling_factor         = nothing,
    nlp_scaling_method         = nothing, # e.g. "gradient-based"
    mu_strategy                = nothing, # "monotone" | "adaptive"
    mu_init                    = nothing, # nothing ⇒ Ipopt's 0.1 (1e-3 with gfm_in_coi + constrain_Δω)
    pardiso_lib_path           = nothing, # required for solver_name "Ipopt-pardiso"
    pardiso_license_message    = true,
)

const highs_cfg = HiGHSSolverConfig(
    output_flag                  = true,
    primal_feasibility_tolerance = 1e-8,
    dual_feasibility_tolerance   = 1e-8,
    ipm_optimality_tolerance     = 1e-8,
    simplex_iteration_limit      = 5_000,
    ipm_iteration_limit          = 5_000,
    solver                       = "choose",  # choose | simplex | ipm | pdlp
    raw_options                  = Dict{String, Any}(),
)

const gurobi_cfg = GurobiSolverConfig(
    output_flag       = 1,
    feasibility_tol   = 1e-8,
    optimality_tol    = 1e-8,
    bar_iter_limit    = 5_000,
    mip_gap           = 1e-8,
    nl_bar_iter_limit = 5_000,
    nl_bar_p_feas_tol = 1e-8,
    nl_bar_d_feas_tol = 1e-8,
    nl_bar_c_feas_tol = 1e-8,
    optimality_target = 1,
    raw_options       = Dict{String, Any}(),
)

const madnlp_cfg = MadNLPSolverConfig(
    tol                   = 1e-8,
    max_iter              = 3_000,
    print_level           = MADNLP_LOG_INFO,  # 3; MADNLP_LOG_ERROR is 6
    hessian_approximation = "exact",
    ma97_num_threads      = nothing,          # e.g. 8 with solver_name "MadNLP-ma97"
    raw_options           = Dict{String, Any}(),
)

# ── 10. The run ──────────────────────────────────────────────────────────────
const cfg = RunConfig(
    # Avenue selector
    trans_stab = true,

    # Orchestration
    case           = "39bus_gfm",
    base_MVA       = 100.0,
    load_factor    = 1.0,
    solver_name    = "Ipopt",
    silent_solver  = false,
    time_limit_sec = 3600.0,

    # Solver option blocks
    ipopt  = ipopt_cfg,
    highs  = highs_cfg,
    gurobi = gurobi_cfg,
    madnlp = madnlp_cfg,

    # Results / I/O
    overwrite_results       = false,
    save_duals              = true,
    save_optim_matrices     = false,  # forced false on TSC paths regardless
    save_matrices           = true,
    save_ts_plots           = false,  # true also needs Plots + load_plots_extension!
    save_ts_debug_csv       = false,  # per-step diagnostics under CSV/Debug/
    save_warmstart_dispatch = true,
    run_script              = nothing,  # @__FILE__ archives this script in the run folder
    post_solve_hook         = nothing,  # f(; kwargs...) called before the model is released

    # Nested avenues
    dispatch  = dispatch_cfg,
    transient = transient_cfg,

    # Optional MATPOWER override (nothing ⇒ read the case CSVs)
    matpower_file = nothing,
)

validate_run_config!(cfg)

sys    = load_system(cfg, PATH_MAIN)
result = run_case!(cfg, sys, PATH_MAIN, PATH_RESULTS)

println("Termination: ", result.status)
println("Objective:   ", result.obj_MVA)          # a plain number, or nothing
println("Results in:  ", result.path_names[:pf_results_date])
```

!!! warning "`run_case!` takes four positional arguments"
    `run_case!(cfg, sys)` is a `MethodError` — the only method is `run_case!(cfg, sys, path_main, path_folder_results)`. `load_system` is the lenient one: `load_system(cfg)` works, defaulting `path_main` to `project_root()`.

!!! warning "The returned model is dead"
    `run_case!` calls `release_solver_backend!`, which empties the JuMP model on Ipopt paths. Every `VariableRef` and `ConstraintRef` inside the returned `dyn_model_dict` is invalid once the call returns, and `result.obj_MVA` is already a plain `Float64` (or `nothing` if the solve did not reach a saved status) — never pass it to `JuMP.value`. Read numbers from the exported CSVs under `result.path_names`, or from `result.RGEN` / `result.RBUS` / `result.RCIR`, which are `DataFrame`s and unaffected.

---

## Field reference

Every table below lists **all** fields of the struct. "Default" is the value in the `@kwdef` definition; the example above overrides some of them.

### `RunConfig` — `engine.jl`

| Field | Type | Default | Notes |
|---|---|---|---|
| `trans_stab` | `Bool` | `false` | Selects the avenue. `true` requires a non-`nothing` `transient`; `false` requires `transient === nothing` |
| `case` | `String` | `"9bus"` | Folder name under `INPUT_FILES/` |
| `base_MVA` | `Float64` | `100.0` | System base. Wins over a MATPOWER file's `baseMVA` (with a warning) |
| `load_factor` | `Float64` | `1.5` | Multiplies `p_d`/`q_d` at load time. Changing it needs a fresh `load_system` |
| `solver_name` | `String` | `"Ipopt"` | `Ipopt`, `Ipopt-ma57`, `Ipopt-ma97`, `Ipopt-pardiso`, `HiGHS`, `Gurobi`, `MadNLP`, `MadNLP-ma57`, `MadNLP-ma97` |
| `silent_solver` | `Bool` | `false` | Forces `print_level = 0` / `set_silent`, overriding the per-solver setting |
| `time_limit_sec` | `Float64` | `600.0` | Wall-clock cap passed to the optimiser |
| `ipopt` | `IpoptSolverConfig` | `IpoptSolverConfig()` | Applied only when `solver_name` is an Ipopt backend |
| `highs` | `HiGHSSolverConfig` | `HiGHSSolverConfig()` | Applied only for `"HiGHS"` |
| `gurobi` | `GurobiSolverConfig` | `GurobiSolverConfig()` | Applied only for `"Gurobi"` |
| `madnlp` | `MadNLPSolverConfig` | `MadNLPSolverConfig()` | Applied only for a MadNLP backend |
| `overwrite_results` | `Bool` | `false` | `true` reuses one folder instead of a new timestamped one |
| `save_duals` | `Bool` | `true` | Writes the registry-driven dual exports. The reason this package exists |
| `save_optim_matrices` | `Bool` | `false` | Steady-state only; `resolve_save_optim_matrices` forces `false` on TSC paths |
| `save_matrices` | `Bool` | `true` | Ybus / Bbus / Kron matrices to `Bus_Matrices/` |
| `save_ts_plots` | `Bool` | `false` | SVG trajectories. `true` without Plots loaded raises at validation |
| `save_ts_debug_csv` | `Bool` | `false` | Per-(window, generator, step) diagnostics under `Transient_Stability/CSV/Debug/` |
| `save_warmstart_dispatch` | `Bool` | `false` | Archives the steady-state pre-solve. Requires TSC + (ACOPF + `FULL_BUS`) or (DCOPF + Kron) |
| `run_script` | `String` or `nothing` | `nothing` | Path of a file to copy into the run folder; usually `@__FILE__` |
| `post_solve_hook` | `Function` or `nothing` | `nothing` | Called once with the model still assembled; keyword arguments only |
| `dispatch` | `DispatchConfig` | `DispatchConfig()` | Avenue 1 |
| `transient` | `TransientConfig` or `nothing` | `nothing` | Avenue 2 |
| `matpower_file` | `String` or `nothing` | `nothing` | Filename relative to the case folder, or an absolute path. Dynamic and contingency CSVs are still read from the case folder |

### `DispatchConfig` — `_common/DispatchConfig.jl`

| Field | Type | Default | Notes |
|---|---|---|---|
| `type_model` | `String` | `"ACOPF"` | `"ACOPF"`, `"DCOPF"`, `"ED"`, `"UC"` |
| `use_matrix` | `Bool` | `true` | Matrix power balance. Required by the explicit DC dual |
| `cost_type` | `String` | `"quadratic"` | `"quadratic"` or `"linear"`. UC requires linear |
| `susceptance_model` | `SusceptanceModel` | `SIMPLE` | `SIMPLE` = $1/x$; `POWERMODELS` = $\Im(1/(r+jx))$. DC-OPF only |
| `bound_V`, `bound_θ`, `bound_P_g`, `bound_Q_g` | `Bool` | `true` | Explicit boxes on the primary variables |
| `bound_P_ik`, `bound_Q_ik`, `bound_P_ki`, `bound_Q_ki` | `Bool` | `false` | Per-direction branch-flow boxes |
| `ineq_sg_upper` | `Bool` | `false` | Generator apparent-power capability |
| `ineq_sbranch_upper` | `Bool` | `true` | Branch thermal limits from `DCIR.l_cap_1`; unrated branches (`0`) are skipped |
| `ineq_ang_diff_branch` | `Bool` | `true` | Angle-difference rows, clamped by `DispatchLimitsConfig` |
| `pd_as_variable` | `Bool` | `false` | ACOPF only: demand as a variable pinned to the data; its dual is `∂L/∂P_d` |
| `solve_explicit_dual` | `Bool` | `false` | Solves the hand-written dual LP. ED or DCOPF, linear cost, and DCOPF also needs `use_matrix = true` |
| `limits` | `DispatchLimitsConfig` | `DispatchLimitsConfig()` | Non-case-data limits |
| `bound_encoding` | `BoundEncoding` | `CONSTRAINT` | `TSCOPF.CONSTRAINT` or `TSCOPF.VARIABLE`. UC ignores it |

### `DispatchLimitsConfig`

| Field | Type | Default | Notes |
|---|---|---|---|
| `θ_min_rad` | `Float64` | `-π` | Only used when `bound_θ = true`; must be `<` the max |
| `θ_max_rad` | `Float64` | `π` | |
| `ang_diff_clamp_ac_deg` | `Float64` | `60.0` | Per-branch angle limits beyond this are tightened to it, with a printed warning |
| `ang_diff_clamp_dc_deg` | `Float64` | `30.0` | Same for the DC path |

Voltage, generator, thermal, and per-branch angle limits are **case data** and come from the CSVs — they deliberately have no fields here.

### `TsSimulationConfig` — `_transient_stability/TsConfig.jl`

| Field | Type | Default | Notes |
|---|---|---|---|
| `δ_tol_deg` | `Float64` | `90.0` | Symmetric corridor half-width, degrees. Must be positive. **The knob to sweep** |
| `δ_tol_deg_lower` | `Float64` or `nothing` | `nothing` | Positive magnitude below COI; `nothing` falls back to `δ_tol_deg` |
| `δ_tol_deg_upper` | `Float64` or `nothing` | `nothing` | Positive magnitude above COI |
| `δ_pair_tol_deg` | `Float64` | `180.0` | Symmetric \|δᵢ − δⱼ\| limit, degrees; read only by `:pairwise` |
| `t_start_sim` | `Float64` | `0.0` | Recorded for output; the grids start at `t_start_fault` |
| `t_end_sim` | `Float64` | `5.0` | Must exceed the clearing instant (SC) or `t_start_fault` (GL/OB) |
| `t_step` | `Float64` | `0.01` | Trapezoidal $\Delta t$. Halving it roughly doubles the NLP |
| `t_start_fault` | `Float64` | `0.01` | Disturbance instant |
| `clearing_time` | `Float64` | `0.3` | SC only. Fault-on duration; ignored for GL and OB |
| `f_syn` | `Float64` | `50.0` | $\omega_s = 2\pi f_{\mathrm{syn}}$ in the angle update |
| `Δω_0` | `Float64` | `0.0` | Initial speed deviation |

### `TsBuilderConfig`

Sixty-two `Bool` toggles plus `limits` and `bound_encoding`. The synchronous-machine boxes default to `false`; the four corridor flags and the GFM guard-rail boxes default to `true`.

| Group | Fields | Default | Effect |
|---|---|---|---|
| Pre-fault classical boxes | `bound_E`, `bound_δ`, `bound_P_m` | `false` | Boxes on the initial internal emf, rotor angle, and mechanical power |
| Classical trajectory boxes | `bound_δ_tf`, `bound_Δω_tf`, `bound_Pe_tf`, `bound_Qe_tf`, `bound_δCOI_tf` and the five `_tpf` twins | `false` | Per-step boxes on the classical states in each window |
| DQ pre-fault algebraic | `bound_Ed`, `bound_Eq`, `bound_Id`, `bound_Iq` | `false` | Boxes on the initial two-axis emfs and stator currents |
| DQ trajectory boxes | `bound_Ed_tf`, `bound_Eq_tf`, `bound_Id_tf`, `bound_Iq_tf`, `bound_Te_tf` and the five `_tpf` twins | `false` | Per-step DQ machine-state boxes |
| Control set-points | `bound_V_ref`, `bound_P_ref` | `false` | Boxes on the AVR and governor set-points |
| Control research boxes | `bound_E_fd_unlim_tf/_tpf`, `bound_Pv_tf/_tpf`, `bound_Pm_tf/_tpf` | `false` | Boxes on the pre-saturation exciter state and the governor trajectories. Saturation physics is unchanged either way |
| δ–COI corridor | `ineq_δ_COI_tf_lower`, `ineq_δ_COI_tf_upper`, `ineq_δ_COI_tpf_lower`, `ineq_δ_COI_tpf_upper` | `true` | The transient-stability constraint itself. Turning all four off gives a plain OPF with dynamics attached but nothing binding |
| GFM guard-rails | `bound_gfm_δ`, `bound_gfm_{P_meas,Q_meas,V_meas,E_int_raw,E_int,E_droop_raw,E_droop,Id,Iq}_{tf,tpf}`, `bound_gfm_pq_{raw,clip}_{tf,tpf}` | `true` | Boxes on the converter states; read only when `allow_gfm` |
| Values | `limits::TsBoundLimitsConfig` | `TsBoundLimitsConfig()` | The numbers behind every box above |
| Encoding | `bound_encoding::BoundEncoding` | `CONSTRAINT` | δ–COI rows and expression targets always use `CONSTRAINT` regardless |

The pre-fault init equalities that link the OPF point to the initial dynamic state (`eq_const_P_init`, `eq_const_Q_init`, `eq_const_Pm_init`) are **not** toggles — they are always built by whichever path applies.

### `TsBoundLimitsConfig`

`TsBoundLimitPair(min, max)` defaults to `(-Inf, Inf)`, which means **inactive**, not "very large". Under `VARIABLE` encoding JuMP simply omits that side; under `CONSTRAINT` no ≤-row is written. A toggle set to `true` with an infinite limit is a no-op.

| Field | Type | Default | Notes |
|---|---|---|---|
| `E_min_pu`, `E_max_pu` | `Float64` | `0.0`, `2.0` | Classical internal emf box **and** the AVR field-voltage clamp. See the warning below |
| `δ_min_rad`, `δ_max_rad` | `Float64` | `-π`, `π` | Pre-fault rotor-angle box |
| `P_m_source` | `Symbol` | `:dgen_pg_limits` | Only accepted value; takes `DGEN.pg_min/pg_max ÷ base_MVA` |
| `V_bus_min_pu` | `Float64` | `0.0` | Forced lower bound on `V_tf` / `V_tpf` on `FULL_BUS` |
| `gov_valve_min_pu` | `Float64` | `0.0` | Valve floor when a limiter is active |
| `gov_valve_max_source` | `Symbol` | `:dgen_pg_limits` | Only accepted value; valve ceiling per generator |
| `δ_tf`, `Δω_tf`, `Pe_tf`, `Qe_tf`, `δCOI_tf` | `TsBoundLimitPair` | `(-Inf, Inf)` | Classical fault-window boxes |
| `δ_tpf`, `Δω_tpf`, `Pe_tpf`, `Qe_tpf`, `δCOI_tpf` | `TsBoundLimitPair` | `(-Inf, Inf)` | Post-fault twins |
| `Ed_min_pu`/`Ed_max_pu`, `Eq_min_pu`/`Eq_max_pu` | `Float64` | `∓2.0` | DQ pre-fault emf boxes |
| `Id_min_pu`/`Id_max_pu`, `Iq_min_pu`/`Iq_max_pu` | `Float64` | `∓5.0` | DQ pre-fault current boxes |
| `V_ref_min_pu`, `V_ref_max_pu` | `Float64` | `0.8`, `1.2` | AVR set-point box |
| `P_ref_source` | `Symbol` | `:dgen_pg_limits` | Governor set-point box source. Note $P_{\mathrm{ref}} = R P_m$ is $R$-scaled, so a power-unit box rarely binds |
| `Ed_tf … Te_tf`, `Ed_tpf … Te_tpf` | `TsBoundLimitPair` | `(-Inf, Inf)` | DQ machine states per window (5 + 5) |
| `E_fd_unlim_tf`, `E_fd_unlim_tpf` | `TsBoundLimitPair` | `(-Inf, Inf)` | Pre-saturation exciter state |
| `Pv_tf`, `Pv_tpf`, `Pm_tf`, `Pm_tpf` | `TsBoundLimitPair` | `(-Inf, Inf)` | Governor valve and mechanical trajectories |
| `gfm_*` (14 fields) | `Float64` | see the example | GFM guard-rail margins and floors; [Parameter reference §3.2](parameter_reference.md) lists each default |

!!! warning "`E_min_pu` / `E_max_pu` do double duty"
    On the classical path they bound the internal emf $E'$. On the DQ path with `include_avr = true` the **same pair** supplies the field-voltage clamp $[E^{\min}, E^{\max}]$ of Model 8.4. The default `[0, 2]` is tight for a field voltage, which is why the example uses `4.0`. There is no separate knob today; see [chapter 8](model/08_controls_avr_governor.md).

### `DynModelConfig` — `_transient_stability/DynModelConfig.jl`

| Field | Type | Default | Notes |
|---|---|---|---|
| `gen_order` | `GenOrder` | `CLASSICAL_2ND` | `CLASSICAL_2ND` or `DQ_4TH`. `DQ_4TH` needs `gen_dynamic_data_full.csv` |
| `network_form` | `NetworkForm` | `KRON_REDUCED` | `KRON_REDUCED` or `FULL_BUS` |
| `mech_power_mode` | `MechPowerMode` | `USE_PG` | `USE_PG` puts the dispatch $P_g$ in the swing; `USE_PM` introduces an explicit $P_m$ |
| `dq_speed_dev_in_algebra` | `Bool` | `true` | Keeps the $(1+\Delta\omega)$ factor in the stator and $P_e$ algebra (RMS convention). `false` neglects it |
| `include_avr` | `Bool` | `false` | Chapter 8. Requires `DQ_4TH` and `T_exc`, `K_exc`, `Ta_exc`, `Tb_exc` columns |
| `include_governor` | `Bool` | `false` | Chapter 8. Requires `USE_PM`, `FULL_BUS`, and `R`, `T1`, `T2`, `T3` columns |
| `governor_limiter` | `GovernorLimiter` | `GOV_NO_LIMIT` | `GOV_NO_LIMIT`, `GOV_SMOOTH`, `GOV_HARD_BOUND`. Only read when the governor is on |
| `allow_gfm` | `Bool` | `false` | Chapter 9. Requires `DQ_4TH` + `FULL_BUS`. At least one active SG must remain unless `allow_all_gfm = true` |
| `allow_all_gfm` | `Bool` | `false` | Chapter 9. Opt-in for a converter-only fleet |
| `gfm_in_coi` | `Bool` | `false` | Chapter 9. Weights converters into the COI with `H_eff = TPf/(2·mp)` and gives them `:coi_box` corridor rows. Governs `bound_style_δ`/`bound_style_Δω` `= :coi_box` only |
| `gfm_pq_overload` | `Bool` | `false` | Chapter 9. `REGFM_A1` Fig. 3 P and Q overload PIs. Off, the droop and the Q–V voltage error are unchanged |
| `ode_first_step` | `Symbol` | `:trapezoidal` | `:trapezoidal` or `:backward_euler` for the first row of each window (swing, EMF, AVR, governor, and the GFM angle). Every later row stays trapezoidal |
| `gfm_integrator` | `Symbol` | `:backward_euler` | GFM measurement filters and Q–V PI **only**. `:backward_euler` is BE at every step (the default — the trapezoidal variants cost ~70× solve time on a 3 s horizon); `:follow_ode_first_step` applies the `ode_first_step` rule above; `:trapezoidal` is trapezoidal at every step. Trapezoidal rings once `t_step > 2·Tf` |
| `zip_load_p` | `NTuple{3,Float64}` | `(1.0, 0.0, 0.0)` | **Active** demand split `(Z, I, P)`; must sum to 1 |
| `zip_load_q` | `NTuple{3,Float64}` | `(1.0, 0.0, 0.0)` | **Reactive** demand split `(Z, I, P)`; independent of `zip_load_p`; must sum to 1 |
| `constrain_δ` | `Bool` | `true` | Build the rotor-angle corridor. At least one of `constrain_δ` / `constrain_Δω` must be on |
| `δ_ref_gen_id` | `Int` or `nothing` | `nothing` | Reference machine; required by `:ref_gen` |
| `bound_style_δ` | `Symbol` | `:swing_propagated` | `:coi_box` writes the corridor directly on $\delta - \delta^{\mathrm{COI}}$; `:swing_propagated` substitutes the trapezoidal update first; `:highest_H` / `:ref_gen` reference one live machine instead of the COI; `:pairwise` bounds every unordered pair of units against `δ_pair_tol_deg` and references nothing |
| `constrain_Δω` | `Bool` | `false` | Adds the speed corridor |
| `bound_style_Δω` | `Symbol` | `:coi_box` | `:coi_box` bounds Δω − Δω_COI; `:abs` bounds the raw Δω and forms no COI |
| `Δω_tol_pu` | `Float64` | `0.5` | Symmetric half-width, p.u. Must be positive when the corridor is on |
| `Δω_tol_pu_lower` | `Float64` or `nothing` | `nothing` | Positive magnitude below COI |
| `Δω_tol_pu_upper` | `Float64` or `nothing` | `nothing` | Positive magnitude above COI |
| `fault` | `FaultConfig` | `FaultConfig()` | The disturbance |

!!! warning "The ZIP ordering trap"
    Both tuples are `(Z, I, P)`: index 1 is constant **impedance**, index 2 constant **current**, index 3 constant **power**. A code that stores the same vector as `(P, I, Z)` writes constant impedance as `[0,0,1]`, which is this package's `(1.0, 0.0, 0.0)`, **not** `(0.0, 0.0, 1.0)`. Copying that ordering silently gives constant-power loads, and the resulting TSC solve typically runs to the iteration limit rather than failing cleanly. The two splits are independent: set **both** fields.

    The two are separate because TSOs do not use one model for both channels — Red Eléctrica, for example, recommends constant-current active demand with constant-admittance reactive demand, which is `zip_load_p = (0.0, 1.0, 0.0)` with `zip_load_q = (1.0, 0.0, 0.0)`.

    On `KRON_REDUCED` both fields are ignored, because loads are folded into `Y_red` as constant admittance. A non-default value there produces a warning, not an error.

### `FaultConfig`

| Field | Type | Default | Notes |
|---|---|---|---|
| `fault_type` | `FaultType` | `SC` | `SC` bus short circuit, `GL` generator or load trip, `OB` open branch |
| `contingency_id` | `Int` | `2` | SC only: row of `contingencies.csv`, giving faulted bus and branch to trip |
| `sc_fault_admittance_pu` | `ComplexF64` | `-im * FAULT_BUS_SHUNT` | SC only: faulted-bus shunt (inductive solid short; `FAULT_BUS_SHUNT = 1e6`) |
| `sc_fault_stamp` | `Symbol` | `:add` | SC only: `:add` adds the shunt to `Y_kk`; `:overwrite` replaces it |
| `gl_gen_ids` | `Vector{Int}` | `Int[]` | GL only: `DGEN` row indices. Cannot trip the slack; at least one unit must survive |
| `gl_load_bus_ids` | `Vector{Int}` | `Int[]` | GL only: buses whose demand is scaled. Mutually exclusive with `gl_gen_ids` in one run |
| `gl_percent_power` | `Vector{Float64}` | `Float64[]` | GL only: $\alpha$ per bus, $p_d \leftarrow p_d(1+\alpha)$. $\alpha = -1$ is a full trip, $-0.5$ sheds half. Values below $-1$ are rejected |
| `ob_branch_ids` | `Vector{Int}` | `Int[]` | OB only: `DCIR` row indices of in-service branches. Duplicates and islanding opens are rejected |

SC uses the two-window timeline (`time_windows_sc`: fault-on for `clearing_time`, then post-fault). GL and OB use a single window from `t_start_fault` to `t_end_sim` with no fault-on stage, and nothing reconnects.

To switch disturbance, replace step 6:

```julia
# Trip generator 3 at t = t_start_fault
FaultConfig(fault_type = GL, gl_gen_ids = [3])

# Shed 40 % of the demand at buses 5 and 7
FaultConfig(fault_type = GL, gl_load_bus_ids = [5, 7], gl_percent_power = [-0.4, -0.4])

# Open circuits 4 and 5 (DCIR row indices)
FaultConfig(fault_type = OB, ob_branch_ids = [4, 5])
```

### `TransientConfig`

| Field | Type | Default | Notes |
|---|---|---|---|
| `simulation` | `TsSimulationConfig` | default | Timing and the δ tolerance |
| `builder` | `TsBuilderConfig` | default | Which constraints get built |
| `dyn_model` | `DynModelConfig` | default | Physics, network form, disturbance |
| `gen_dynamic_filename` | `String` | `"gen_dynamic_data.csv"` | SG machine data. `DQ_4TH` needs the `_full` variant |
| `gfm_dynamic_filename` | `String` | `"gfm_dynamic_data.csv"` | Read only when `allow_gfm = true` |
| `extra_contingency_ids` | `Vector{Int}` | `Int[]` | Further SC contingencies in the same solve; results in `Transient_Stability_c<id>/` |

### Solver configs

| Struct | Fields | Notable defaults |
|---|---|---|
| `IpoptSolverConfig` | `tol`, `print_level`, `max_iter`, `constr_viol_tol`, `dual_inf_tol`, `compl_inf_tol`, `hessian_approximation`, `limited_memory_max_history`, `acceptable_tol`, `obj_scaling_factor`, `nlp_scaling_method`, `mu_strategy`, `mu_init`, `pardiso_lib_path`, `pardiso_license_message` | `tol = 1e-8`, `print_level = 5` (range 0–12), `max_iter = 5000`, `hessian_approximation = "exact"`. The six `nothing`-able fields leave Ipopt's own default in place when unset |
| `HiGHSSolverConfig` | `output_flag`, `primal_feasibility_tolerance`, `dual_feasibility_tolerance`, `ipm_optimality_tolerance`, `simplex_iteration_limit`, `ipm_iteration_limit`, `solver`, `raw_options` | `solver = "choose"` (`simplex`, `ipm`, `pdlp` also valid); iteration limits `5000` |
| `GurobiSolverConfig` | `output_flag`, `feasibility_tol`, `optimality_tol`, `bar_iter_limit`, `mip_gap`, `nl_bar_iter_limit`, `nl_bar_p_feas_tol`, `nl_bar_d_feas_tol`, `nl_bar_c_feas_tol`, `optimality_target`, `raw_options` | `output_flag = 1` (an `Int`, not a `Bool`); the four `nl_bar_*` fields need Gurobi 13+ |
| `MadNLPSolverConfig` | `tol`, `max_iter`, `print_level`, `hessian_approximation`, `ma97_num_threads`, `raw_options` | `max_iter = 3000`, `print_level = MADNLP_LOG_INFO` (`3`; `MADNLP_LOG_ERROR` is `6`) |

`raw_options` is the escape hatch: any solver attribute not surfaced as a field can go in there as a `String => value` pair. `IpoptSolverConfig` has no `raw_options` — and no `output_file`, because `run_case!` redirects logs to `RESULTS/.../solver_log.txt` and `solver_log_warmstart.txt` itself.

Solver choice is partly enforced: `Check_Coherence_Input_Data` rejects ACOPF or TSC-ACOPF with an LP solver, UC with any solver but Gurobi, and `trans_stab = true` on ED or UC. A DCOPF or ED run with an NLP solver is accepted. Gurobi and MadNLP live in package extensions — `using Gurobi` or `using MadNLP` must happen before `run_case!` so the extension loads.

---

## Combinations the validators reject

Each of these raises an `ArgumentError` before any variable is created. The check runs in `validate_run_config!`, which fans out to the others.

| Configuration | Rejected by | Message (abridged) |
|---|---|---|
| `trans_stab = true`, `transient = nothing` | `validate_run_config!` | `trans_stab=true requires transient::TransientConfig.` |
| `trans_stab = false`, `transient` set | `validate_run_config!` | `trans_stab=false requires transient=nothing.` |
| `save_ts_plots = true` without Plots loaded | `validate_run_config!` | `save_ts_plots=true requires Plots.jl in the active environment …` |
| `save_warmstart_dispatch = true` on a run with no steady-state pre-solve | `validate_run_config!` | `save_warmstart_dispatch=true requires a TSC run with a steady-state pre-solve …` |
| `FULL_BUS` with `mech_power_mode = USE_PG` | `validate_dyn_config!` | `network_form=FULL_BUS requires mech_power_mode=USE_PM …` |
| `KRON_REDUCED` with `ode_first_step = :backward_euler` | `validate_dyn_config!` | `ode_first_step=:backward_euler is not implemented on network_form=KRON_REDUCED …` |
| Unknown `type_model` or `cost_type` | `validate_dispatch_config!` | `Unknown type_model "…"` |
| `type_model = "UC"` with quadratic cost | `validate_dispatch_config!` | `UC requires cost_type="linear" (MILP).` |
| `solve_explicit_dual = true` outside ED/DCOPF, or with quadratic cost, or DCOPF without `use_matrix` | `validate_dispatch_config!` | `solve_explicit_dual=true requires type_model="DCOPF" or "ED".` |
| `θ_min_rad ≥ θ_max_rad`, or a non-positive angle clamp | `validate_dispatch_limits!` | `require θ_min_rad < θ_max_rad.` |
| A ZIP tuple not summing to 1 | `validate_dyn_config!` | `zip_load_p (Z, I, P) coefficients must sum to 1 …` |
| `DQ_4TH` with `KRON_REDUCED` | `validate_dyn_config!` | `DQ_4TH requires FULL_BUS network_form …` |
| `DQ_4TH` without `USE_PM` | `validate_dyn_config!` | `DQ_4TH requires mech_power_mode=USE_PM.` |
| `DQ_4TH` with `:swing_propagated` | `validate_dyn_config!` | `DQ_4TH requires a box-form bound_style_δ (:coi_box), not :swing_propagated.` |
| `include_avr` without `DQ_4TH` | `validate_dyn_config!` | `include_avr=true requires gen_order=DQ_4TH.` |
| `include_governor` without `USE_PM` | `validate_dyn_config!` | `include_governor=true requires mech_power_mode=USE_PM …` |
| `include_governor` on `KRON_REDUCED` | `validate_dyn_config!` | `include_governor=true is currently supported only on network_form=FULL_BUS …` |
| `USE_PM` with `:swing_propagated` | `validate_dyn_config!` | `USE_PM requires a box-form bound_style_δ (:coi_box) …` |
| `bound_style_δ` outside the five symbols | `validate_dyn_config!` | `bound_style_δ must be :swing_propagated, :coi_box, :highest_H, :ref_gen or :pairwise …` |
| `:ref_gen` without an id | `validate_dyn_config!` | `bound_style_δ=:ref_gen requires δ_ref_gen_id …` |
| `:ref_gen` id that is tripped or out of service (a GFM id is legal) | `validate_δ_reference!` | `δ_ref_gen_id=… is tripped by the GL disturbance …` |
| Neither corridor enabled | `validate_dyn_config!` | `A transient-stability run needs at least one corridor …` |
| `ode_first_step` outside the two symbols | `validate_dyn_config!` | `ode_first_step must be :trapezoidal or :backward_euler …` |
| `gfm_integrator` outside the three symbols | `validate_dyn_config!` | `gfm_integrator must be :follow_ode_first_step, :backward_euler or :trapezoidal …` |
| TSC-DCOPF with `FULL_BUS` | `validate_dyn_config!` | `TSC-DCOPF with FULL_BUS is not implemented yet …` |
| TSC-DCOPF with `DQ_4TH` | `validate_dyn_config!` | `TSC-DCOPF with DQ_4TH is not implemented …` |
| `allow_gfm` outside `DQ_4TH` + `FULL_BUS`, or under DCOPF | `validate_dyn_config!` | `allow_gfm=true requires network_form=FULL_BUS.` |
| `constrain_Δω` with a non-positive tolerance | `validate_dyn_config!` | `constrain_Δω=true requires Δω_tol_pu > 0.` |
| Non-positive `δ_tol_deg` or `t_step` | `validate_transient_config!` | `δ_tol_deg must be positive.` |
| Any limit pair with `min ≥ max`, or a `*_source` other than `:dgen_pg_limits` | `validate_ts_bound_limits!` | `require E_min_pu < E_max_pu.` |
| SC with a non-empty `ob_branch_ids`, GL specifying both gens and loads, OB with an empty list, an OB set that islands the network, tripping the slack | `validate_fault_config!` | `OB open of branches […] would island the network …` |
| An all-GFM fleet (no active SG) without `allow_all_gfm = true` | `resolve_sg_gfm_gens!` | `DQ FULL_BUS requires at least one active SG …` |
| An all-GFM fleet with `:swing_propagated` / `:highest_H`, a COI corridor without `gfm_in_coi`, or AVR / governor on | `validate_all_gfm_fleet!` | raised before the warm-start solve |
| ACOPF with an LP solver, UC without Gurobi, `trans_stab = true` on ED or UC | `Check_Coherence_Input_Data` | raised at the top of `run_case!` |

---

## What is *not* configurable here

Three groups of numbers look like parameters but are not reachable from `RunConfig`:

- **Smoothing constants** — `_AVR_SMOOTH_RHO` and `_GOV_SMOOTH_RHO`, both $10^{-4}$, fix the corner rounding of the AVR and governor clamps; `_GFM_EPS_E`, `_GFM_EPS_NORM_V` and `_GFM_EPS_LIM_I` do the same for the converter clip and current limiter.
- **`_GFM_IMAX_NO_LIMIT`** — the 20.0 pu threshold above which the GFM current limiter is replaced by unsaturated stator algebra. A formulation switch driven by input data rather than a tuning knob; bypassed units are logged. See [chapter 9](model/09_grid_forming.md). The GFM *box* limits, by contrast, **are** configurable — `TsBoundLimitsConfig` fields `gfm_*`.
- **Machine and network data** — everything in `INPUT_FILES/<case>/*.csv`. Use `load_factor` for demand studies rather than editing the CSVs in place; the builders `deepcopy` the DataFrames before applying GL or OB disturbances precisely so the originals stay clean.

## Varying one field

For sweeps, do not rebuild the config by hand — `reconfigure` copies it with overrides, deep-copying the nested `dispatch` and `transient` so the original is untouched:

```julia
sys = load_system(cfg, PATH_MAIN)          # load the network once

for δ_val in 80.0:10.0:110.0
    tc = reconfigure_transient(cfg.transient;
        simulation = TsSimulationConfig(
            (f => f === :δ_tol_deg ? δ_val : getfield(cfg.transient.simulation, f)
             for f in fieldnames(TsSimulationConfig))...))
    run_case!(reconfigure(cfg; transient = tc), sys, PATH_MAIN, PATH_RESULTS)
end
```

!!! warning "`TsSimulationConfig(δ_tol_deg = …)` resets the other nine fields"
    The shorter idiom `simulation = TsSimulationConfig(δ_tol_deg = δ_val)` — used in [Running a case](running_a_case.md) §1.7 — constructs a **fresh** struct, so `t_step`, `clearing_time`, `t_end_sim` and everything else silently revert to their `@kwdef` defaults. That is harmless when the base config already uses defaults, and a silent corruption of the sweep when it does not. The splat above carries the base values through; use it whenever the base `TsSimulationConfig` is non-default. The same caveat applies to `TsBuilderConfig` and `DynModelConfig` overrides.

`load_factor` cannot be swept this way: it is applied to `p_d`/`q_d` inside `load_system`, so a demand sweep needs a fresh `SystemData` per point.

---

## See also

- [Running a case](running_a_case.md) — output folder layout, solver logs, dual exports
- [Parameter reference](parameter_reference.md) — the same fields grouped by topic, with the CSV schemas
- [Configuration map](configuration_map.md) — struct nesting and where each field is consumed
- [Chapter 8](model/08_controls_avr_governor.md), [chapter 9](model/09_grid_forming.md) — the physics behind the control and converter fields
