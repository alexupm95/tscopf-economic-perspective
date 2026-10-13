#=
================================================================================
 39-BUS RUN OF THE PAPER (SECTION IV-B)
================================================================================
  Usage, from the root of the repository:

      julia --project=TSCOPF runs/run_39bus.jl

  With no environment variable set, the script solves the run of the paper,
  results/39bus/Results - 2026-10-07 124340. A new `Results - <timestamp>`
  folder is written under new_runs/, and one row per solve is appended to
  new_runs/summary.csv.

  Settings (environment variables; lists are comma-separated):

      J_CASE      case folder(s) under INPUT_FILES                    (39bus_gfm)
      J_CONT      contingency row(s) of contingencies.csv             (57)
      J_CLEAR     clearing time(s) [s]                                (0.095)
      J_DTOL      COI corridor half-width [deg]                       (100)
      J_TSTEP     time step [s]                                       (0.005)
      J_TEND      end of the window [s]                               (5.0)
      J_FAULT     inductive (−j1e6, :add) | resistive (1e10, :overwrite)
      J_GFM_INT   gfm_integrator                                      (trapezoidal)
      J_GFM_PQ    1 = converter P/Q overload loops on                 (1)
      J_KKT       0 none | 1 kkt_Pg_terms | 2 kkt_Pg_terms + sensitivity (1)
      J_LIMITS    branch thermal limits, 1 | 0                        (1)
      J_MAXITER   Ipopt iteration cap                                 (300)
      J_SOLVER    Ipopt | Ipopt-ma57 (needs HSL_jll in the environment)
      J_PLOTS     1 = save the trajectory plots                       (1)

  The fault window keeps the node at t_start_fault + clearing_time, so the fault
  lasts clearing_time + t_step: 0.095 s on a 5 ms grid is a 100 ms fault.
================================================================================
=#

env(k, d) = get(ENV, k, d)
const PLOTS = env("J_PLOTS", "1") == "1"
PLOTS && (ENV["GKSwstype"] = "100")          # headless GR, before Plots loads

using TSCOPF
using Dates
PLOTS && @eval using Plots

const PATH_MAIN = normpath(joinpath(@__DIR__, ".."))
const PATH_RES  = joinpath(PATH_MAIN, "new_runs")
mkpath(PATH_RES)
const SUMMARY = joinpath(PATH_RES, "summary.csv")

const CASES     = String.(split(env("J_CASE", "39bus_gfm"), ','))
const CONTS     = parse.(Int, split(env("J_CONT", "57"), ','))
const CLEARINGS = parse.(Float64, split(env("J_CLEAR", "0.095"), ','))
const DTOL      = parse(Float64, env("J_DTOL", "100"))
const TSTEP     = parse(Float64, env("J_TSTEP", "0.005"))
const TEND      = parse(Float64, env("J_TEND", "5.0"))
const FAULT     = env("J_FAULT", "inductive")
const GFM_INT   = Symbol(env("J_GFM_INT", "trapezoidal"))
const GFM_PQ    = env("J_GFM_PQ", "1") == "1"
const KKT       = parse(Int, env("J_KKT", "1"))
const LIMITS    = env("J_LIMITS", "1") == "1"
const MAXITER   = parse(Int, env("J_MAXITER", "300"))
const SOLVER    = env("J_SOLVER", "Ipopt")

FAULT in ("inductive", "resistive") || error("J_FAULT must be inductive or resistive (got $FAULT)")

# ── KKT exports, attached through the post-solve hook ─────────────────────────────
const KKT_SOURCE = joinpath(PATH_MAIN, "kkt")
KKT >= 1 && include(joinpath(KKT_SOURCE, "kkt_pg_terms.jl"))
KKT >= 2 && include(joinpath(KKT_SOURCE, "kkt_sensitivity.jl"))

function kkt_post_solve_hook(; cfg, sys, path_names, opf_dict, dyn_model_dict, model = nothing,
                               status_solved, kwargs...)
    (cfg.save_duals && status_solved) || return nothing
    KKT >= 1 && Save_KKT_Pg_Terms!(path_names, cfg, sys, opf_dict, dyn_model_dict)
    KKT >= 2 && Save_KKT_Sensitivity!(path_names, cfg, sys, opf_dict, dyn_model_dict, model)
    return nothing
end

const OPEN = TsBoundLimitPair(min = -Inf, max = Inf)

# E_max_pu is also the AVR field ceiling:
# at the default 2.0 the units near the fault sit on it before the disturbance.
const TS_LIMITS = TsBoundLimitsConfig(
    E_min_pu = 0.0, E_max_pu = 6.0, δ_min_rad = -π, δ_max_rad = π,
    P_m_source = :dgen_pg_limits, V_bus_min_pu = 0.0,
    Ed_min_pu = -2.0, Ed_max_pu = 2.0, Eq_min_pu = -2.0, Eq_max_pu = 2.0,
    Id_min_pu = -5.0, Id_max_pu = 5.0, Iq_min_pu = -5.0, Iq_max_pu = 5.0,
    gov_valve_min_pu = 0.0, gov_valve_max_source = :dgen_pg_limits,
    P_ref_source = :dgen_pg_limits, V_ref_min_pu = 0.0, V_ref_max_pu = 6.0,
)

# Every optional box off except the pre-fault δ box. That covers the 19
# `bound_gfm_*` guard-rails too, which default to true and stall Ipopt on a GFM fleet.
const BUILDER = TsBuilderConfig(;
    (f => (f == :bound_δ) for f in fieldnames(TsBuilderConfig)
     if startswith(string(f), "bound_") && f != :bound_encoding)...,
    ineq_δ_COI_tf_lower = true, ineq_δ_COI_tf_upper = true,
    ineq_δ_COI_tpf_lower = true, ineq_δ_COI_tpf_upper = true,
    limits = TS_LIMITS, bound_encoding = TSCOPF.CONSTRAINT,
)

"""Run configuration for one (case, contingency, clearing time) cell."""
function journal_cfg(case::String, cid::Int, clearing::Float64)
    gfm = isfile(joinpath(PATH_MAIN, "INPUT_FILES", case, "gfm_dynamic_data.csv"))

    # The resistive 1e10 short is an active-power sink
    # with a small shadow price on the faulted bus. The inductive short is the clean one
    # for pricing and is the package default.
    fault = FAULT == "resistive" ?
        FaultConfig(fault_type = SC, contingency_id = cid,
                    sc_fault_admittance_pu = 1e10 + 0im, sc_fault_stamp = :overwrite) :
        FaultConfig(fault_type = SC, contingency_id = cid)

    dyn = DynModelConfig(
        gen_order = DQ_4TH, network_form = FULL_BUS, mech_power_mode = USE_PM,
        dq_speed_dev_in_algebra = false,
        include_avr = true, include_governor = true, governor_limiter = GOV_SMOOTH,
        ode_first_step = :trapezoidal,
        allow_gfm = gfm, gfm_in_coi = gfm, gfm_pq_overload = gfm && GFM_PQ,
        gfm_integrator = gfm ? GFM_INT : :backward_euler,
        zip_load_p = (1.0, 0.0, 0.0), zip_load_q = (1.0, 0.0, 0.0),   # constant impedance
        constrain_δ = true, bound_style_δ = :coi_box,                  # the one frame
        constrain_Δω = false,                                          # never both
        fault = fault,
    )
    sim = TsSimulationConfig(
        δ_tol_deg = DTOL, t_start_sim = 0.0, t_end_sim = TEND, t_step = TSTEP,
        t_start_fault = 0.01, clearing_time = clearing, f_syn = 50.0, Δω_0 = 0.0,
    )
    return RunConfig(
        trans_stab = true, case = case, base_MVA = 100.0, load_factor = 1.0,
        solver_name = SOLVER, silent_solver = false, time_limit_sec = 3600.0,
        # Quasi-Newton Hessian, history 50.
        ipopt = IpoptSolverConfig(
            tol = 1e-8, print_level = 5, max_iter = MAXITER,
            constr_viol_tol = 1e-8, dual_inf_tol = 1e-8, compl_inf_tol = 1e-8,
            hessian_approximation = "limited-memory", limited_memory_max_history = 50),
        overwrite_results = false,
        save_duals = true, save_optim_matrices = false, save_matrices = false,
        save_ts_plots = PLOTS, save_ts_debug_csv = gfm, save_warmstart_dispatch = true,
        run_script = @__FILE__,
        post_solve_hook = KKT >= 1 ? kkt_post_solve_hook : nothing,
        dispatch = DispatchConfig(
            type_model = "ACOPF", use_matrix = true, cost_type = "linear",
            bound_V = true, bound_θ = true, bound_P_g = true, bound_Q_g = true,
            ineq_sg_upper = false, ineq_sbranch_upper = LIMITS, ineq_ang_diff_branch = false,
            limits = DispatchLimitsConfig(θ_min_rad = -π, θ_max_rad = π,
                                          ang_diff_clamp_ac_deg = 60.0),
            bound_encoding = TSCOPF.CONSTRAINT),
        transient = TransientConfig(
            simulation = sim, builder = BUILDER, dyn_model = dyn,
            gen_dynamic_filename = "gen_dynamic_data_full.csv",
            gfm_dynamic_filename = gfm ? "gfm_dynamic_data.csv" : ""),
    )
end

"""Ipopt iteration count of the TSC solve, read from the run's `solver_log.txt`."""
function ipopt_iterations(run_dir::String)
    log = joinpath(run_dir, "solver_log.txt")
    isfile(log) || return ""
    m = match(r"Number of Iterations\.+:\s*(\d+)", read(log, String))
    return m === nothing ? "" : m.captures[1]
end

isfile(SUMMARY) || write(SUMMARY,
    "time;case;contingency;clearing_s;dtol_deg;t_step;fault;gfm_int;status;obj;t_build_s;t_solve_s;iterations;run_dir\n")

function summary_row(case, cid, ct, cols...)
    open(SUMMARY, "a") do io
        println(io, join([string(now()), case, cid, ct, DTOL, TSTEP, FAULT, GFM_INT, cols...], ';'))
    end
end

for case in CASES, cid in CONTS, ct in CLEARINGS
    summary_row(case, cid, ct, "START", "", "", "", "", "")
    row = try
        cfg = journal_cfg(case, cid, ct)
        validate_run_config!(cfg)
        validate_dyn_config!(cfg)
        r = run_case!(cfg, load_system(cfg, PATH_MAIN), PATH_MAIN, PATH_RES)
        dir = r.path_names[:pf_results_date]
        (string(r.status), string(r.obj_MVA), round(r.t_build; digits = 1),
         round(r.t_solve; digits = 1), ipopt_iterations(dir), dir)
    catch err
        ("ERROR: " * first(replace(sprint(showerror, err), r"[\r\n;]+" => " "), 300),
         "", "", "", "", "")
    end
    summary_row(case, cid, ct, row...)
    println("<<< ", case, "  cont ", cid, "  clearing ", ct, " s : ", row[1], "  obj ", row[2],
            "  iterations ", row[5])
    flush(stdout)
    GC.gc()                     # release the solved model before the next cell is built
end
