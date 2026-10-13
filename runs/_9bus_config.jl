#=
================================================================================
 SHARED MODEL CONFIGURATION — 9-BUS RUNS (SECTION IV-A)
================================================================================
  Included by `run_9bus.jl`, which must define `CASE` (a folder under
  INPUT_FILES/) first. Defines `cfg` (TSCOPF, every saver off) and `ipopt_acopf`.

    * three-phase short circuit at bus 7 (contingency 2), applied at t = 0.01 s
      and cleared by opening branch 5–7;
    * 0–5 s, δ − δ_COI within ±100°;
    * DQ_4TH + AVR (SEXS) + governor (TGOV1), FULL_BUS network, ZIP (1,0,0) default;
    * lossless case, no branch-flow limits (ineq_sbranch_upper = false);
    * Ipopt tol 1e-8, limited-memory Hessian with history 50.

  `t_step` and `clearing_time` below are placeholders: `run_9bus.jl` overrides
  both, and the fault model, with the values of the paper.
================================================================================
=#

OPEN = TsBoundLimitPair(min = -Inf, max = Inf)

dispatch_cfg = DispatchConfig(
    type_model = "ACOPF", use_matrix = true, cost_type = "linear",
    bound_V = true, bound_θ = true, bound_P_g = true, bound_Q_g = true,
    bound_P_ik = false, bound_Q_ik = false, bound_P_ki = false, bound_Q_ki = false,
    ineq_sg_upper = false, ineq_sbranch_upper = false, ineq_ang_diff_branch = false,
    solve_explicit_dual = false,
    limits = DispatchLimitsConfig(θ_min_rad = -π, θ_max_rad = π, ang_diff_clamp_ac_deg = 60.0),
    bound_encoding = TSCOPF.CONSTRAINT,
)

ts_simulation = TsSimulationConfig(
    δ_tol_deg = 100.0, δ_tol_deg_lower = nothing, δ_tol_deg_upper = nothing,
    t_start_sim = 0.0, t_end_sim = 5.0, t_step = 0.01,
    t_start_fault = 0.01, clearing_time = 0.10, f_syn = 50.0, Δω_0 = 0.0,
)

ts_limits = TsBoundLimitsConfig(
    E_min_pu = 0.0, E_max_pu = 6.0, δ_min_rad = -π, δ_max_rad = π,
    P_m_source = :dgen_pg_limits, V_bus_min_pu = 0.0,
    δ_tf = OPEN, Δω_tf = OPEN, Pe_tf = OPEN, Qe_tf = OPEN, δCOI_tf = OPEN,
    δ_tpf = OPEN, Δω_tpf = OPEN, Pe_tpf = OPEN, Qe_tpf = OPEN, δCOI_tpf = OPEN,
    Ed_min_pu = -2.0, Ed_max_pu = 2.0, Eq_min_pu = -2.0, Eq_max_pu = 2.0,
    Id_min_pu = -5.0, Id_max_pu = 5.0, Iq_min_pu = -5.0, Iq_max_pu = 5.0,
    Ed_tf = OPEN, Eq_tf = OPEN, Id_tf = OPEN, Iq_tf = OPEN, Te_tf = OPEN,
    Ed_tpf = OPEN, Eq_tpf = OPEN, Id_tpf = OPEN, Iq_tpf = OPEN, Te_tpf = OPEN,
    gov_valve_min_pu = 0.0, gov_valve_max_source = :dgen_pg_limits,
    P_ref_source = :dgen_pg_limits,
    Pv_tf = OPEN, Pv_tpf = OPEN, Pm_tf = OPEN, Pm_tpf = OPEN,
    V_ref_min_pu = 0.0, V_ref_max_pu = 6.0,
    E_fd_unlim_tf = OPEN, E_fd_unlim_tpf = OPEN,
)

ts_builder = TsBuilderConfig(
    bound_E = false, bound_δ = true, bound_P_m = false,
    bound_δ_tf = false, bound_Δω_tf = false, bound_Pe_tf = false,
    bound_Qe_tf = false, bound_δCOI_tf = false,
    bound_δ_tpf = false, bound_Δω_tpf = false, bound_Pe_tpf = false,
    bound_Qe_tpf = false, bound_δCOI_tpf = false,
    bound_Ed = false, bound_Eq = false, bound_Id = false, bound_Iq = false,
    bound_Ed_tf = false, bound_Eq_tf = false, bound_Id_tf = false,
    bound_Iq_tf = false, bound_Te_tf = false,
    bound_Ed_tpf = false, bound_Eq_tpf = false, bound_Id_tpf = false,
    bound_Iq_tpf = false, bound_Te_tpf = false,
    bound_V_ref = false, bound_E_fd_unlim_tf = false, bound_E_fd_unlim_tpf = false,
    ineq_δ_COI_tf_lower = true, ineq_δ_COI_tf_upper = true,
    ineq_δ_COI_tpf_lower = true, ineq_δ_COI_tpf_upper = true,
    limits = ts_limits, bound_encoding = TSCOPF.CONSTRAINT,
)

fault_cfg = FaultConfig(
    fault_type = SC, contingency_id = 2,
    gl_gen_ids = Int[], gl_load_bus_ids = Int[], gl_percent_power = Float64[],
    ob_branch_ids = Int[],
)

dyn_model_cfg = DynModelConfig(
    gen_order = DQ_4TH, network_form = FULL_BUS, mech_power_mode = USE_PM,
    dq_speed_dev_in_algebra = false,
    include_avr = true, include_governor = true, governor_limiter = GOV_SMOOTH,
    ode_first_step = :trapezoidal,
    zip_load_p = (1.0, 0.0, 0.0), zip_load_q = (1.0, 0.0, 0.0),
    constrain_δ = true, bound_style_δ = :coi_box, δ_ref_gen_id = 1,
    constrain_Δω = false, bound_style_Δω = :coi_box, Δω_tol_pu = 0.03,
    Δω_tol_pu_lower = nothing, Δω_tol_pu_upper = nothing,
    fault = fault_cfg,
)

# TSCOPF solver
ipopt_cfg = IpoptSolverConfig(
    tol = 1e-8, print_level = 5, max_iter = 1_000,
    constr_viol_tol = 1e-8, dual_inf_tol = 1e-8, compl_inf_tol = 1e-8,
    hessian_approximation = "limited-memory", limited_memory_max_history = 50,
    pardiso_lib_path = get(ENV, "JULIA_PARDISO_LIB", nothing),
    pardiso_license_message = false,
)

# ACOPF benchmark solver
ipopt_acopf = IpoptSolverConfig(
    tol = 1e-8, print_level = 5, max_iter = 1_000,
    constr_viol_tol = 1e-8, dual_inf_tol = 1e-8, compl_inf_tol = 1e-8,
    hessian_approximation = "exact", limited_memory_max_history = 50,
)

cfg = RunConfig(
    trans_stab = true, case = CASE, base_MVA = 100.0, load_factor = 1.0,
    solver_name = "Ipopt", silent_solver = true, time_limit_sec = 3600.0,
    ipopt = ipopt_cfg, overwrite_results = true,
    save_duals = false, save_optim_matrices = false, save_matrices = false,
    save_ts_plots = false, save_ts_debug_csv = false, save_warmstart_dispatch = false,
    dispatch = dispatch_cfg,
    transient = TransientConfig(
        simulation = ts_simulation, builder = ts_builder, dyn_model = dyn_model_cfg,
        gen_dynamic_filename = "gen_dynamic_data_full.csv", gfm_dynamic_filename = "",
    ),
)
