#=
================================================================================
 engine.jl  —  single orchestration core for the TSC-OPF pipeline
================================================================================
 Phase 1 refactor (MF-1 + MF-3 of the audit blueprint).

 Purpose
 -------
 Replace the ~140-line orchestration body that was copy-pasted across main.jl,
 main_loop.jl, main_mcp.jl and main_loop_mcp.jl with ONE function, `run_case!`,
 driven by ONE configuration object, `RunConfig`.

 This file defines:
   • RunConfig    — every run knob in one immutable struct (was scattered across
                    Define_Initial_Parameters(), the swept δ_tol global, and the
                    _contingency_id global).
   • reconfigure  — copy a RunConfig overriding a few fields (used by sweeps).
   • SystemData   — the loaded network data, built once and reused across a sweep.
   • load_system  — read input files + organise the network (was inlined per run).
   • run_case!    — build → solve → (δ_ref) → dynamics → solve → save, for ONE
                    configuration. Returns a NamedTuple with the outcome.

 What it deliberately does NOT change (behaviour preservation)
 ------------------------------------------------------------
   • The model builders, constraint code, and save functions are called exactly
     as before. Only the *plumbing* around them is consolidated.
   • The pre-solve that produces δ_ref runs only for trans_stab + "DCOPF", which
     is what main.jl did. (Old main_loop.jl also pre-solved for ACOPF, a
     redundant extra solve that, for the nonconvex ACOPF, could warm-start the
     final solve into a different local optimum. Dropping it here is a fix, and
     it is flagged in the Phase 1 notes. It does not affect the DCOPF δ_tol
     sweep, where the pre-solve is required and happens identically.)

 Author: Alex Junior da Cunha Coelho  ·  refactor scaffolding, June 2026
================================================================================
=#

# ==============================================================================
#  RunConfig — the single source of truth for one run
# ==============================================================================
# Defaults reproduce the values returned by Define_Initial_Parameters() in
# config_file.jl, so `run_case!(default RunConfig)` == the old default run.
# Override only what a given scenario needs (see main.jl / main_loop.jl wrappers).
"""
    RunConfig

Top-level configuration for one TSC-OPF run: every knob shared by plain
dispatch and transient-stability runs, plus the two nested avenues
`dispatch::DispatchConfig` (steady-state OPF) and
`transient::Union{Nothing,TransientConfig}` (required when `trans_stab = true`).

Build one with keyword arguments (defaults reproduce the legacy integrated
9-bus TSC-ACOPF run), pass it to [`load_system`](@ref) and [`run_case!`](@ref),
and use [`reconfigure`](@ref) to vary a handful of fields (for example a
δ_tol sweep) without mutating the original or touching global state.

See the parameter reference, section 1, and the configuration map, section 1, for the
full field reference.
"""
Base.@kwdef struct RunConfig
    # --- avenue selector -----------------------------------------------------
    trans_stab::Bool = false

    # --- orchestration (shared by dispatch and TSC) --------------------------
    case::String = "9bus"
    base_MVA::Float64 = 100.0
    load_factor::Float64 = 1.5
    solver_name::String = "Ipopt"
    silent_solver::Bool = false
    time_limit_sec::Float64 = 600.0
    ipopt::IpoptSolverConfig = IpoptSolverConfig()
    highs::HiGHSSolverConfig = HiGHSSolverConfig()
    gurobi::GurobiSolverConfig = GurobiSolverConfig()
    madnlp::MadNLPSolverConfig = MadNLPSolverConfig()

    # --- results / I/O -------------------------------------------------------
    overwrite_results::Bool = false
    save_duals::Bool = true
    save_optim_matrices::Bool = false   # steady-state only; resolve_save_optim_matrices forces false for TSC
    save_matrices::Bool = true
    save_ts_plots::Bool = false   # TSC trajectory figures (Plots.jl); opt-in for sweeps/CI
    # Per-(window, gen, step) diagnostic dumps under Transient_Stability/CSV/Debug/:
    # value, distance to the limit and dual on one row. Off by default — the filter file
    # alone is three rows per generator per step, which adds up fast over a sweep.
    save_ts_debug_csv::Bool = false
    # TSC runs that pre-solve a steady-state OPF before TS assembly (FULL_BUS ACOPF
    # warm start, TSC-DCOPF δ_ref anchor): archive that solution to Dispatch_WarmStart/.
    save_warmstart_dispatch::Bool = false
    # Archive the script that configured this run. `nothing` = off; set it to
    # `@__FILE__` in the run script to drop a byte-for-byte copy of that file into
    # the run folder, next to input_parameters.txt. Any readable path works (a
    # sweep driver can archive its own driver file).
    run_script::Union{Nothing, String} = nothing
    # Called once per run with the model still assembled, immediately before
    # `release_solver_backend!` invalidates every VariableRef/ConstraintRef. Receives
    # keyword arguments only (`cfg`, `sys`, `path_names`, `model`, `opf_dict`,
    # `dyn_model_dict`, `dyn_parameters_dict`, `status`, `status_solved`, `obj_MVA`,
    # `t_build`, `t_solve`), so a hook written today keeps working when more data is
    # handed over later — declare it as `f(; kwargs...)` and pick out what it needs.
    # `nothing` = off. Exceptions are not caught: an export that silently failed is
    # worse than a run that stopped and said so.
    post_solve_hook::Union{Nothing, Function} = nothing

    # --- avenue 1: steady-state dispatch -------------------------------------
    dispatch::DispatchConfig = DispatchConfig()

    # --- avenue 2: transient analysis (required when trans_stab=true) --------
    transient::Union{Nothing, TransientConfig} = nothing

    # --- input-data source (optional MATPOWER override) ----------------------
    # When `nothing` (default) steady-state data is read from the in-house CSVs
    # in INPUT_FILES/<case>/ (bus_data.csv, generators_data.csv, line_data.csv).
    # When set to a filename (relative to INPUT_FILES/<case>/) or an absolute
    # path, DBUS/DGEN/DCIR are parsed from that MATPOWER .m file instead.
    # Dynamic data (gen_dynamic_data.csv) and contingency tables are always
    # read from INPUT_FILES/<case>/ regardless of this setting.
    matpower_file::Union{Nothing, String} = nothing
end

"""Deep-copy a `DynModelConfig` (avoids @kwdef shared-default mutation)."""
function copy_dyn_model_config(dm::DynModelConfig)::DynModelConfig
    return DynModelConfig(;
        (f => f == :fault ? copy_fault_config(dm.fault) : getfield(dm, f)
         for f in fieldnames(DynModelConfig))...)
end

"""
    reconfigure(cfg::RunConfig; kwargs...) -> RunConfig

Copy `cfg`, overriding the fields given as keyword arguments. `dispatch` and
`transient` are deep-copied by default (via `copy_dispatch_config` /
`copy_transient_config`) unless explicitly passed in `kwargs`, so sweeps can
vary one knob (for example `load_factor` or `δ_tol_deg` through
[`reconfigure_transient`](@ref)) without mutating the original `cfg` or
touching global state.
"""
function reconfigure(cfg::RunConfig; kwargs...)
    base = (; (f => getfield(cfg, f) for f in fieldnames(RunConfig))...)
    merged = pairs((; base..., kwargs...))
    args = Dict{Symbol, Any}(merged)
    if !haskey(args, :dispatch)
        args[:dispatch] = copy_dispatch_config(cfg.dispatch)
    end
    if !haskey(args, :transient) && cfg.transient !== nothing
        args[:transient] = copy_transient_config(cfg.transient)
    end
    return RunConfig(; args...)
end

"""
    validate_δCOI_bound_request!(tc::TransientConfig)

Reject `TsBuilderConfig.bound_δCOI_tf` / `bound_δCOI_tpf` on a run whose δ corridor does
not reference the centre of inertia.

Those toggles put explicit bounds on the `δCOI_*` **variable**, and a machine-referenced
run (or one with `constrain_δ=false`) builds the COI as an expression instead — there is
no variable to bound. The attach helpers skip absent variables silently, so without this
check the request would simply evaporate.
"""
function validate_δCOI_bound_request!(tc::TransientConfig)
    dyn = tc.dyn_model
    coi_referenced = dyn.constrain_δ &&
        dyn.bound_style_δ ∈ (:coi_box, :swing_propagated)
    coi_referenced && return nothing
    for field in (:bound_δCOI_tf, :bound_δCOI_tpf)
        getfield(tc.builder, field) && throw(ArgumentError(
            "$field=true needs a δ_COI variable, which only exists when the δ corridor " *
            "references the COI (constrain_δ=true with bound_style_δ ∈ (:coi_box, " *
            ":swing_propagated)). This run has constrain_δ=$(dyn.constrain_δ), " *
            "bound_style_δ=:$(dyn.bound_style_δ)."))
    end
    return nothing
end

"""
    validate_run_config!(cfg::RunConfig)

Fail fast on incoherent run configuration (avenue split, dispatch, transient).
"""
function validate_run_config!(cfg::RunConfig)
    validate_dispatch_config!(cfg.dispatch)
    if cfg.trans_stab
        cfg.transient === nothing &&
            throw(ArgumentError("trans_stab=true requires transient::TransientConfig."))
        validate_transient_config!(cfg.transient)
        validate_dyn_config!(cfg)
        validate_δCOI_bound_request!(cfg.transient)
    elseif cfg.transient !== nothing
        throw(ArgumentError("trans_stab=false requires transient=nothing."))
    end
    if cfg.trans_stab && cfg.save_ts_plots && !load_plots_extension!()
        throw(ArgumentError(
            "save_ts_plots=true requires Plots.jl in the active environment. " *
            "Add it with `Pkg.add(\"Plots\")` and call `load_plots_extension!()` " *
            "before `run_case!` (loads Plots + Measures for the extension)."))
    end
    if cfg.save_warmstart_dispatch && !warmstart_dispatch_applicable(cfg)
        throw(ArgumentError(
            "save_warmstart_dispatch=true requires a TSC run with a steady-state pre-solve: " *
            "either dispatch.type_model=\"ACOPF\" with network_form=FULL_BUS (ACOPF warm start) " *
            "or dispatch.type_model=\"DCOPF\" with network_form=KRON_REDUCED (δ_ref anchor)."))
    end
    # Caught here rather than at export time: a typo in the path should not surface
    # after a long solve, when the run folder already exists.
    if cfg.run_script !== nothing && !isfile(cfg.run_script)
        throw(ArgumentError(
            "run_script points at no readable file: $(cfg.run_script). " *
            "Use `run_script = @__FILE__` inside the run script, or `nothing` to disable."))
    end
    if is_ipopt_backend_solver(cfg.solver_name)
        validate_ipopt_solver_config!(cfg.ipopt, cfg.solver_name)
    elseif cfg.solver_name == "HiGHS"
        validate_highs_solver_config!(cfg.highs)
    elseif cfg.solver_name == "Gurobi"
        validate_gurobi_solver_config!(cfg.gurobi)
    elseif cfg.solver_name == "MadNLP"
        validate_madnlp_solver_config!(cfg.madnlp)
    end
    return nothing
end

"""
    validate_dyn_config!(cfg::RunConfig)

Fail fast on incoherent model-selection combinations. Default `DynModelConfig()`
always passes so existing runs are unchanged.
"""
function validate_dyn_config!(cfg::RunConfig)
    cfg.transient === nothing && return nothing
    dyn = cfg.transient.dyn_model
    type_model = cfg.dispatch.type_model

    # Active and reactive demand carry independent (Z, I, P) splits; each must sum to 1
    # on its own — a shared check would let a mismatched pair through.
    for (name, split) in ((:zip_load_p, dyn.zip_load_p), (:zip_load_q, dyn.zip_load_q))
        z, i, p = split   # (Z, I, P): impedance, current, power
        if !isapprox(z + i + p, 1.0; atol=1e-9)
            throw(ArgumentError(
                "$name (Z, I, P) coefficients must sum to 1 (got Z=$z, I=$i, P=$p)."))
        end
    end

    # Kron folds loads into Y_red as constant admittance, so the splits are dead knobs
    # there. Warn rather than throw — a non-default split is a misconfiguration, not a
    # reason to block an otherwise valid Kron run.
    if dyn.network_form == KRON_REDUCED &&
       (dyn.zip_load_p != (1.0, 0.0, 0.0) || dyn.zip_load_q != (1.0, 0.0, 0.0))
        @warn "zip_load_p / zip_load_q are ignored on KRON_REDUCED (loads are folded " *
              "into Y_red as constant admittance). Use network_form=FULL_BUS to apply them." *
              " Got zip_load_p=$(dyn.zip_load_p), zip_load_q=$(dyn.zip_load_q)."
    end

    if dyn.gen_order == DQ_4TH && dyn.network_form != FULL_BUS
        throw(ArgumentError(
            "DQ_4TH requires FULL_BUS network_form (Kron reduction is inconsistent with dq dynamics)."))
    end

    if dyn.include_avr && dyn.gen_order != DQ_4TH
        throw(ArgumentError("include_avr=true requires gen_order=DQ_4TH."))
    end

    # The FULL_BUS builders need mechanical power as its own state (the swing is driven
    # by an explicit nodal Pe, so pinning P_mech ≡ P_g would over-constrain the coupling).
    # Both builders already throw, but only after the warm-start ACOPF has been solved —
    # check here so the run dies before paying for that solve.
    if dyn.network_form == FULL_BUS && dyn.mech_power_mode != USE_PM
        throw(ArgumentError(
            "network_form=FULL_BUS requires mech_power_mode=USE_PM " *
            "(P_m must be an independent variable on the full-network paths)."))
    end

    if dyn.gen_order == DQ_4TH
        # Milestone 1/2/3: machine core; optional AVR and governor.
        dyn.mech_power_mode != USE_PM && throw(ArgumentError(
            "DQ_4TH requires mech_power_mode=USE_PM."))
        dyn.bound_style_δ == :swing_propagated && throw(ArgumentError(
            "DQ_4TH requires a box-form bound_style_δ (:coi_box), not :swing_propagated."))
    end

    if dyn.include_governor && dyn.mech_power_mode != USE_PM
        throw(ArgumentError(
            "include_governor=true requires mech_power_mode=USE_PM (governor acts on P_m / P_mech)."))
    end

    if dyn.include_governor && dyn.network_form != FULL_BUS
        throw(ArgumentError(
            "include_governor=true is currently supported only on network_form=FULL_BUS " *
            "(the Kron path does not run the FULL_BUS coupling initialisation pipeline)."))
    end

    if cfg.trans_stab && type_model == "DCOPF" && dyn.gen_order == DQ_4TH
        throw(ArgumentError(
            "TSC-DCOPF with DQ_4TH is not implemented (no linearised dq formulation yet)."))
    end

    # --- stability corridors: two independent knobs, at least one of them on -------
    # The rotor-angle and speed corridors are what makes a run transient-stability
    # *constrained*; with both off the TS block is pure simulation bolted onto an OPF,
    # which is never what the caller meant.
    if dyn.bound_style_δ ∉ (:swing_propagated, :coi_box, :highest_H, :ref_gen, :pairwise)
        throw(ArgumentError(
            "bound_style_δ must be :swing_propagated, :coi_box, :highest_H, :ref_gen or " *
            ":pairwise (got $(dyn.bound_style_δ))."))
    end

    # The id itself is checked against the system data in `validate_δ_reference!` — here we
    # only catch the two errors that need no DataFrames.
    if dyn.bound_style_δ == :ref_gen
        dyn.δ_ref_gen_id === nothing && throw(ArgumentError(
            "bound_style_δ=:ref_gen requires δ_ref_gen_id (a generator id from gen_dynamic_data)."))
        dyn.δ_ref_gen_id > 0 || throw(ArgumentError(
            "δ_ref_gen_id must be a positive generator id (got $(dyn.δ_ref_gen_id))."))
    elseif dyn.δ_ref_gen_id !== nothing
        @warn "δ_ref_gen_id=$(dyn.δ_ref_gen_id) is ignored when " *
              "bound_style_δ=:$(dyn.bound_style_δ) (it only selects the reference machine " *
              "for :ref_gen)."
    end

    if dyn.bound_style_Δω ∉ (:coi_box, :abs)
        throw(ArgumentError(
            "bound_style_Δω must be :coi_box or :abs (got $(dyn.bound_style_Δω))."))
    end

    if cfg.trans_stab && !dyn.constrain_δ && !dyn.constrain_Δω
        throw(ArgumentError(
            "A transient-stability run needs at least one corridor: set constrain_δ=true " *
            "(rotor angle) or constrain_Δω=true (speed deviation)."))
    end

    if dyn.ode_first_step ∉ (:trapezoidal, :backward_euler)
        throw(ArgumentError(
            "ode_first_step must be :trapezoidal or :backward_euler (got $(dyn.ode_first_step))."))
    end

    # The Kron swing rows hard-code the trapezoidal average at every step, including
    # t=1, so backward Euler is unimplemented there rather than merely unwired. Throw
    # instead of silently integrating with a scheme the user did not ask for.
    if dyn.ode_first_step === :backward_euler && dyn.network_form == KRON_REDUCED
        throw(ArgumentError(
            "ode_first_step=:backward_euler is not implemented on network_form=KRON_REDUCED " *
            "(the Kron swing equations are trapezoidal at every step). " *
            "Use network_form=FULL_BUS for the backward-Euler first step."))
    end

    if dyn.gfm_integrator ∉ (:follow_ode_first_step, :backward_euler, :trapezoidal)
        throw(ArgumentError(
            "gfm_integrator must be :follow_ode_first_step, :backward_euler or " *
            ":trapezoidal (got $(dyn.gfm_integrator))."))
    end

    # The GFM knob only reaches the DQ FULL_BUS converter path. Warn rather than throw:
    # a stale non-default on a Kron or classical run is a misconfiguration, not a reason
    # to block an otherwise valid case.
    if dyn.gfm_integrator != :backward_euler && !dyn.allow_gfm
        @warn "gfm_integrator=$(dyn.gfm_integrator) is ignored when allow_gfm=false " *
              "(it only governs the GFM measurement filters and the Q–V PI integrator)."
    end

    if dyn.gfm_pq_overload && !dyn.allow_gfm
        @warn "gfm_pq_overload=true is ignored when allow_gfm=false " *
              "(it only stamps the REGFM_A1 P and Q overload PIs)."
    end

    # Same idiom for the COI membership knob, which is inert on two counts: no converters
    # in the fleet, or no COI-referenced corridor switched on for them to join.
    if dyn.gfm_in_coi
        if !dyn.allow_gfm
            @warn "gfm_in_coi=true is ignored when allow_gfm=false (there are no " *
                  "converters to weight into the COI)."
        else
            δ_coi = dyn.constrain_δ && dyn.bound_style_δ == :coi_box
            Δω_coi = dyn.constrain_Δω && dyn.bound_style_Δω == :coi_box
            if !δ_coi && !Δω_coi
                @warn "gfm_in_coi=true has no effect: no COI-referenced corridor is " *
                      "active (bound_style_δ=$(dyn.bound_style_δ), constrain_δ=" *
                      "$(dyn.constrain_δ); bound_style_Δω=$(dyn.bound_style_Δω), " *
                      "constrain_Δω=$(dyn.constrain_Δω))."
            end
        end
    end

    if dyn.mech_power_mode == USE_PM && dyn.bound_style_δ == :swing_propagated
        throw(ArgumentError(
            "USE_PM requires a box-form bound_style_δ (:coi_box); :swing_propagated " *
            "substitutes the dispatch power into the angle band, which is invalid with " *
            "an explicit P_m."))
    end

    if cfg.trans_stab && type_model == "DCOPF" && dyn.network_form == FULL_BUS
        throw(ArgumentError(
            "TSC-DCOPF with FULL_BUS is not implemented yet (Phase 3 is TSC-ACOPF only)."))
    end

    if dyn.constrain_Δω
        dyn.Δω_tol_pu > 0.0 ||
            throw(ArgumentError("constrain_Δω=true requires Δω_tol_pu > 0."))
        # Overrides are positive magnitudes below/above the reference (see Δω_tol_tuple).
        for (name, v) in ((:Δω_tol_pu_lower, dyn.Δω_tol_pu_lower),
                          (:Δω_tol_pu_upper, dyn.Δω_tol_pu_upper))
            v === nothing || v > 0.0 ||
                throw(ArgumentError("$name must be positive when set."))
        end
    end

    if dyn.allow_gfm
        dyn.network_form == FULL_BUS || throw(ArgumentError(
            "allow_gfm=true requires network_form=FULL_BUS."))
        dyn.gen_order == DQ_4TH || throw(ArgumentError(
            "allow_gfm=true requires gen_order=DQ_4TH."))
        type_model == "DCOPF" && throw(ArgumentError(
            "allow_gfm=true is not supported with TSC-DCOPF (linearize path)."))
    end

    return nothing
end

"""
    validate_δ_reference!(cfg, DGEN, DGEN_DYN, DGFM=nothing)

Check the machine-referenced δ corridor (`bound_style_δ = :highest_H | :ref_gen`)
against the actual system data.

`validate_dyn_config!` cannot do this: whether a generator id is in service, tripped by
the disturbance, a GFM unit, or carries inertia is a property of the CSVs, not of the
config. Run it from `run_case!` **before** the mandatory FULL_BUS ACOPF warm start —
otherwise a typo in `δ_ref_gen_id` costs a full NLP solve before it is caught.

The corridor is built over the units that survive the disturbance, so the reference is
selected from that set: `:highest_H` can never land on a unit the run is about to trip.

Grid-forming converters are corridor *members* here — the machine-referenced styles bound
their angle like any other — but they are not `:highest_H` *candidates*, because that style
ranks by inertia and `DGFM` carries none. A converter becomes the reference only through an
explicit `:ref_gen`, and the `H > 0` check then does not apply to it: `DGEN_DYN` holds
machine rows only, which is also why `H` is looked up by id rather than row.
"""
function validate_δ_reference!(
    cfg::RunConfig,
    DGEN::DataFrame,
    DGEN_DYN::Union{DataFrame, Nothing},
    DGFM::Union{DataFrame, Nothing}=nothing,
)
    cfg.trans_stab || return nothing
    cfg.transient === nothing && return nothing
    dyn = cfg.transient.dyn_model
    dyn.constrain_δ || return nothing
    # `:pairwise` resolves no reference, but it shares the "two units or the corridor is
    # empty" failure mode and must be caught here rather than at the solve.
    dyn.bound_style_δ ∈ (:highest_H, :ref_gen, :pairwise) || return nothing
    DGEN_DYN === nothing && throw(ArgumentError(
        "bound_style_δ=:$(dyn.bound_style_δ) needs machine dynamic data (DGEN_DYN)."))

    # Units the corridor will actually span: in service and not tripped, converters included.
    gfm_ids = DGFM === nothing ? Set{Int}() : gfm_id_set(DGFM)
    tripped = dyn.fault.fault_type == GL ? Set{Int}(dyn.fault.gl_gen_ids) : Set{Int}()
    survivors = [g for g in findall(x -> x == 1, DGEN.g_status) if g ∉ tripped]
    # The narrower set `:highest_H` ranks over, and the only one `DGEN_DYN.H` may be indexed by.
    sg_survivors = [g for g in survivors if g ∉ gfm_ids]

    length(survivors) ≥ 2 || throw(ArgumentError(
        "bound_style_δ=:$(dyn.bound_style_δ) needs at least two units " *
        "left after the disturbance (found $(length(survivors)): $(survivors)). The " *
        "corridor is built between units, so a single survivor yields no rows at all."))

    # Pairwise has no reference machine to validate: every surviving unit carries rows, and
    # the count check above is the whole contract.
    dyn.bound_style_δ == :pairwise && return nothing

    if dyn.bound_style_δ == :ref_gen
        ref = dyn.δ_ref_gen_id
        # A converter is a legal reference, but it has no row in the machine CSV, so the
        # id-range and inertia checks below are for synchronous machines only.
        1 ≤ ref ≤ nrow(DGEN) || throw(ArgumentError(
            "δ_ref_gen_id=$ref is not a generator id (generators_data has " *
            "$(nrow(DGEN)) rows)."))
        ref ∈ tripped && throw(ArgumentError(
            "δ_ref_gen_id=$ref is tripped by the GL disturbance (gl_gen_ids=" *
            "$(dyn.fault.gl_gen_ids)); the reference must stay synchronised."))
        ref ∈ survivors || throw(ArgumentError(
            "δ_ref_gen_id=$ref is out of service (DGEN.g_status = 0)."))
        if ref ∉ gfm_ids
            # By id: under GFM `DGEN_DYN` holds SG rows only, so row ≠ id past a converter.
            _dgen_dyn_H(DGEN_DYN, ref) > 0.0 || throw(ArgumentError(
                "δ_ref_gen_id=$ref has H = 0 in gen_dynamic_data; pick a machine with inertia."))
        end
    else
        any(g -> _dgen_dyn_H(DGEN_DYN, g) > 0.0, sg_survivors) || throw(ArgumentError(
            "bound_style_δ=:highest_H found no surviving machine with H > 0 " *
            "(candidates: $(sg_survivors)). The reference is ranked by inertia and so is " *
            "always a synchronous machine, even when converters are bounded by the corridor."))
    end

    return nothing
end

"""
    validate_all_gfm_fleet!(cfg, DGEN, DGEN_DYN, DGFM)

Screen a **converter-only** fleet — every in-service generator listed in
`gfm_dynamic_data.csv` — against the combinations that cannot work without at least one
synchronous machine.

No-ops unless the active fleet really is all grid-forming, so nothing here can affect a
mixed or SG-only run. Fleet membership is a property of the CSVs rather than the config,
so this runs from `run_case!` with `SystemData` in hand, next to
[`validate_δ_reference!`](@ref) and **before** the mandatory FULL_BUS ACOPF warm start.

Why each rejection exists:

  * `allow_all_gfm = false` — the default. A case that happens to list every unit in the
    GFM file would otherwise silently become a converter-only run whose "COI" is a defined
    weighting (`H_eff = TPf/(2·mp)`) rather than a physical centre of inertia. Opting in is
    a statement that this was intended.
  * `bound_style_δ = :swing_propagated` — substitutes a swing step into the angle band and
    needs a per-machine `H` *and* `D`; a converter has neither.
  * `bound_style_δ = :highest_H` — ranks its reference by machine inertia, and there are no
    machines to rank.
  * a COI-referenced corridor with `gfm_in_coi = false` — `_coi_gens` would return an empty
    set and `eq_const_kron_COI_generic!` would die reducing over an empty collection. A
    clear message here beats an obscure one three layers down.
  * `include_avr` / `include_governor` — both attach to `sg_gens` only, so on this fleet
    they are silent no-ops. Rejecting is better than letting a run *look* like it has an
    exciter when it does not.
  * `δ_ref_gen_id` naming a non-converter — nothing else exists to reference.
"""
function validate_all_gfm_fleet!(
    cfg::RunConfig,
    DGEN::DataFrame,
    DGEN_DYN::Union{DataFrame, Nothing},
    DGFM::Union{DataFrame, Nothing}=nothing,
)
    cfg.trans_stab || return nothing
    cfg.transient === nothing && return nothing
    dyn = cfg.transient.dyn_model

    gfm_ids = DGFM === nothing ? Set{Int}() : gfm_id_set(DGFM)
    isempty(gfm_ids) && return nothing
    active = findall(x -> x == 1, DGEN.g_status)
    isempty(active) && return nothing
    # Only a fleet with no synchronous machine at all reaches the checks below.
    all(g -> g ∈ gfm_ids, active) || return nothing

    n = length(active)
    dyn.allow_all_gfm || throw(ArgumentError(
        "every in-service generator ($n of $n: $(active)) is grid-forming, and " *
        "converter-only fleets are opt-in. Set DynModelConfig.allow_all_gfm=true if that " *
        "is intended. Note that the COI of such a fleet is a defined weighting " *
        "(H_eff = TPf/(2·mp)), not a physical centre of inertia."))

    if dyn.constrain_δ && dyn.bound_style_δ === :swing_propagated
        throw(ArgumentError(
            "bound_style_δ=:swing_propagated is unavailable on a converter-only fleet: " *
            "it substitutes a swing step into the angle band and needs a per-machine H " *
            "and D, which grid-forming units do not have. Use :coi_box (with " *
            "gfm_in_coi=true), :ref_gen or :pairwise."))
    end
    if dyn.constrain_δ && dyn.bound_style_δ === :highest_H
        throw(ArgumentError(
            "bound_style_δ=:highest_H is unavailable on a converter-only fleet: it ranks " *
            "its reference machine by inertia and there are no machines to rank. Use " *
            ":coi_box (with gfm_in_coi=true), :ref_gen or :pairwise."))
    end

    # `:coi_box` on either corridor makes the COI the reference; with no SGs the COI set is
    # empty unless the converters are weighted into it.
    coi_referenced = (dyn.constrain_δ && dyn.bound_style_δ === :coi_box) ||
                     (dyn.constrain_Δω && dyn.bound_style_Δω === :coi_box)
    if coi_referenced && !dyn.gfm_in_coi
        throw(ArgumentError(
            "a COI-referenced corridor on a converter-only fleet requires " *
            "gfm_in_coi=true; otherwise the COI spans synchronous machines only and the " *
            "set is empty. Either set gfm_in_coi=true or pick a corridor that does not " *
            "reference the COI (:ref_gen, :pairwise, or bound_style_Δω=:abs)."))
    end

    for (flag, name) in ((dyn.include_avr, "include_avr"),
                         (dyn.include_governor, "include_governor"))
        flag && throw(ArgumentError(
            "$name=true has no effect on a converter-only fleet — the AVR and governor " *
            "families are built over synchronous machines only. Set it false rather than " *
            "running a model that appears to carry a control it does not."))
    end

    if dyn.constrain_δ && dyn.bound_style_δ === :ref_gen && dyn.δ_ref_gen_id !== nothing
        dyn.δ_ref_gen_id ∈ gfm_ids || throw(ArgumentError(
            "δ_ref_gen_id=$(dyn.δ_ref_gen_id) is not a grid-forming unit, and this fleet " *
            "has no synchronous machines. Reference one of $(sort(collect(gfm_ids)))."))
    end
    return nothing
end


"""
    validate_dyn_data!(cfg::RunConfig, DGEN_DYN)

Fail fast when the machine CSV named by `TransientConfig.gen_dynamic_filename` does
not carry the columns the selected `DynModelConfig` needs — `R, T1, T2, T3` under
`include_governor`, `T_exc, K_exc, Ta_exc, Tb_exc` under `include_avr`, the dq set
under `DQ_4TH` — or carries AVR parameters that have no realization: a non-positive
`T_exc`/`K_exc`, or lead-lag time constants that do not describe a proper block
(see [`validate_avr_data!`](@ref)).

Separate from [`validate_dyn_config!`](@ref) because that one is config-only and
runs before any data is read; this needs `DGEN_DYN`, so it runs from `run_case!`
where `SystemData` is in hand — the same shape as `validate_fault_config!`.

Without this the run would not fail at all on the classical path: the parser
pre-fills absent optional columns with `NaN`, so `DGEN_DYN.R` exists whatever the
file contains, and those `NaN`s get stamped straight into governor constraint
coefficients.

Every row is checked, out-of-service generators included — the same rule
`has_required_dyn_columns` has always applied. Mixed SG + GFM fleets cannot
false-positive here: grid-forming units live in `gfm_dynamic_data.csv` and are not
rows of `DGEN_DYN`.
"""
function validate_dyn_data!(cfg::RunConfig, DGEN_DYN::Union{DataFrame, Nothing})
    cfg.trans_stab || return nothing
    dyn = cfg.transient.dyn_model
    fname = cfg.transient.gen_dynamic_filename

    DGEN_DYN === nothing && throw(ArgumentError(
        "trans_stab=true requires generator dynamic data, but none was loaded " *
        "(expected $fname). Build the SystemData with load_system on a TSC RunConfig."))

    missing_cols = missing_dyn_columns(dyn, DGEN_DYN)
    if !isempty(missing_cols)
        cols = join(string.(missing_cols), ", ")
        why = join(dyn_column_reasons(missing_cols), ", ")
        throw(ArgumentError(
            "$fname is missing machine data required by the selected dynamic model: " *
            "$cols (required by $why). Columns must be present and finite for every " *
            "generator row."))
    end

    dyn.include_avr && validate_avr_data!(DGEN_DYN, fname)
    return nothing
end

"""
    validate_avr_data!(DGEN_DYN, fname)

Reject SEXS exciter parameters that the AVR builders cannot represent.

Gain-lag stage, per generator row:

  * `T_exc > 0` and `K_exc > 0`. Both are divided by — `T_exc` in the exciter row's
    `c = Δt/(2·T_exc)` and its backward-Euler `Δt/T_exc`, `K_exc` in the lead-lag
    warm start — so a zero reaches the model as an infinite coefficient. `K_exc = 0`
    is the quieter of the two and the worse: the pre-fault link `E_fd = K_exc·(V_ref − V)`
    degenerates to `E_fd = 0` and the exciter loses its input, which converges to
    something meaningless rather than failing.

Lead-lag stage, per generator row:

  * `Ta_exc ≥ 0` and `Tb_exc ≥ 0` — negative time constants are not a model.
  * `Tb_exc == 0` requires `Ta_exc == 0`. With both zero the block is the exact
    pass-through `E_LL ≡ V_ref − V` and is skipped outright; with only `Tb_exc` zero
    it is a bare `1 + Ta_exc·s` differentiator, improper, whose discretization has
    coefficients that blow up as `Δt → 0` and whose alternating mode is genuinely
    excited. That second case used to be caught by the old `Tb_exc > 0` build gate;
    now that the bypass keys on *both* being zero, nothing else rejects it.

Called from [`validate_dyn_data!`](@ref) only when `include_avr = true`, after the
column-presence check has guaranteed all four columns exist and are finite — hence
no `isnan` guard on the comparisons below. Keeping it behind that gate is what lets
a non-AVR run carry `NaN` in these columns, which the DQ-only fixtures rely on.
"""
function validate_avr_data!(DGEN_DYN::DataFrame, fname::String)
    for i in 1:nrow(DGEN_DYN)
        T_exc, K_exc = DGEN_DYN.T_exc[i], DGEN_DYN.K_exc[i]
        (T_exc <= 0.0 || K_exc <= 0.0) && throw(ArgumentError(
            "$fname row $i (bus $(DGEN_DYN.bus[i])) has a non-positive AVR gain-lag " *
            "parameter (T_exc=$T_exc, K_exc=$K_exc). Both must be > 0: the exciter row " *
            "divides by T_exc, and K_exc = 0 leaves the field voltage pinned at zero " *
            "with no voltage feedback."))

        Ta, Tb = DGEN_DYN.Ta_exc[i], DGEN_DYN.Tb_exc[i]
        (Ta < 0.0 || Tb < 0.0) && throw(ArgumentError(
            "$fname row $i (bus $(DGEN_DYN.bus[i])) has a negative AVR lead-lag time " *
            "constant (Ta_exc=$Ta, Tb_exc=$Tb). Both must be ≥ 0."))
        if iszero(Tb) && !iszero(Ta)
            throw(ArgumentError(
                "$fname row $i (bus $(DGEN_DYN.bus[i])) sets Tb_exc=0 with Ta_exc=$Ta, " *
                "which is a pure differentiator (1 + Ta_exc·s) and has no state-space " *
                "realization. Use Ta_exc=Tb_exc=0 to bypass the lead-lag stage, or " *
                "Tb_exc > 0 to model it."))
        end
    end
    return nothing
end

# ==============================================================================
#  SystemData — network data loaded once, reused across a sweep
# ==============================================================================
# DGEN_DYN is only present when trans_stab is true; hence the Union with Nothing.
"""
    SystemData

Network data loaded once by [`load_system`](@ref) and reused across a run or
a parameter sweep: bus/generator/circuit `DataFrame`s (`DBUS`, `DGEN`,
`DCIR`), optional SG dynamic machine data (`DGEN_DYN`, `nothing` unless
`trans_stab = true`), optional GFM dynamic data (`DGFM`, `nothing` unless
`allow_gfm`), the optional per-load demand split (`DLOAD`, `nothing` unless
the case folder has a `load_data.csv`), the bus-index remap (`bus_mapping`,
`reverse_bus_mapping`), element counts (`nBUS`, `nGEN`, `nCIR`), adjacency
lookups (`bus_gen_circ_dict`, `bus_gen_circ_dict_ON`), and the resolved
input-file path (`path_input`).

Demand (`p_d`, `q_d`) is already scaled by `RunConfig.load_factor` at this
point — build a fresh `SystemData` if `load_factor` changes.
"""
struct SystemData
    DBUS::DataFrame
    DGEN::DataFrame
    DGEN_DYN::Union{DataFrame, Nothing}
    DGFM::Union{DataFrame, Nothing}
    DLOAD::Union{DataFrame, Nothing}
    DCIR::DataFrame
    bus_mapping::OrderedDict
    reverse_bus_mapping::OrderedDict
    nBUS::Int
    nGEN::Int
    nCIR::Int
    bus_gen_circ_dict::OrderedDict
    bus_gen_circ_dict_ON::OrderedDict
    path_input::String
end

"""
    load_system(cfg::RunConfig, path_main::String=project_root()) -> SystemData

Read the input files, scale demand by `cfg.load_factor`, and organise the
bus/generator/circuit adjacency into a [`SystemData`](@ref). `path_main`
defaults to the repository root ([`project_root`](@ref)).

Two input modes, governed by `cfg.matpower_file`:

- CSV mode (`matpower_file = nothing`, default): reads `bus_data.csv`,
  `generators_data.csv`, `line_data.csv` from `<path_main>/INPUT_FILES/<case>/`.
- MATPOWER mode (`matpower_file = "case9.m"` or an absolute path):
  steady-state data (`DBUS`, `DGEN`, `DCIR`) is parsed in-memory from the
  `.m` file; dynamic data (`DGEN_DYN`) and contingency tables are still read
  from `<path_main>/INPUT_FILES/<case>/`.

Every input path derives from that one join, so a case tree outside the
repository is used by passing its root — no other argument or environment
variable is involved.

Call once per case and reuse the returned `SystemData` across a sweep
(`load_factor` is only applied here, so build a fresh `SystemData` if it
changes between runs).
"""
function load_system(cfg::RunConfig, path_main::String=project_root())::SystemData
    path_input = joinpath(path_main, "INPUT_FILES", cfg.case)
    allow_gfm = cfg.trans_stab && cfg.transient !== nothing &&
                cfg.transient.dyn_model.allow_gfm
    gfm_dyn_file = cfg.trans_stab && cfg.transient !== nothing ?
        cfg.transient.gfm_dynamic_filename : "gfm_dynamic_data.csv"

    if isnothing(cfg.matpower_file)
        # ── CSV path (existing behaviour, unchanged) ───────────────────────────
        if cfg.trans_stab
            gen_dyn_file = cfg.transient.gen_dynamic_filename
            DBUS, DGEN, DGEN_DYN, DCIR, bus_mapping, reverse_bus_mapping =
                Read_Input_Data(path_input, cfg.trans_stab;
                    gen_dynamic_filename=gen_dyn_file)
            nBUS = length(DBUS.bus); nGEN = length(DGEN.id); nCIR = length(DCIR.from_bus)
            if allow_gfm
                DGFM = Read_GFM_Dynamic_Data(path_input;
                    filename=gfm_dyn_file, sys_base_MVA=cfg.base_MVA)
                # Prefer explicit `id` from SG CSV when present (partitioned fleets).
                DGEN_DYN = _ensure_gen_dyn_ids_from_raw_if_needed(path_input, gen_dyn_file, DGEN_DYN)
                validate_gen_dyn_partition!(DGEN, DGEN_DYN, DGFM)
            else
                DGFM = nothing
                Check_Gen_Sta_Dyn_Data(nGEN, nrow(DGEN_DYN))
                validate_dgen_dyn_row_ids!(DGEN, DGEN_DYN)
            end
        else
            DBUS, DGEN, DCIR, bus_mapping, reverse_bus_mapping =
                Read_Input_Data(path_input, cfg.trans_stab)
            nBUS = length(DBUS.bus); nGEN = length(DGEN.id); nCIR = length(DCIR.from_bus)
            DGEN_DYN = nothing
            DGFM = nothing
        end
    else
        # ── MATPOWER .m path ───────────────────────────────────────────────────
        path_m = isabspath(cfg.matpower_file) ?
                 cfg.matpower_file :
                 joinpath(path_input, cfg.matpower_file)
        isfile(path_m) || throw(ArgumentError(
            "matpower_file not found: \"$path_m\". " *
            "Verify RunConfig.matpower_file and RunConfig.case."))

        DBUS, DGEN, DCIR, parsed_baseMVA, bus_mapping, reverse_bus_mapping =
            matpower_to_inhouse_dataframes(path_m)

        if !isapprox(parsed_baseMVA, cfg.base_MVA; rtol=1e-6)
            @warn ("MATPOWER baseMVA=$(parsed_baseMVA) differs from " *
                   "RunConfig.base_MVA=$(cfg.base_MVA). Using RunConfig value. " *
                   "Set base_MVA=$(parsed_baseMVA) to suppress this warning.")
        end

        nBUS = length(DBUS.bus); nGEN = length(DGEN.id); nCIR = length(DCIR.from_bus)

        if cfg.trans_stab
            gen_dyn_file = cfg.transient.gen_dynamic_filename
            DGEN_DYN = Read_Gen_Dynamic_Data(path_input; filename=gen_dyn_file)
            if allow_gfm
                DGFM = Read_GFM_Dynamic_Data(path_input;
                    filename=gfm_dyn_file, sys_base_MVA=cfg.base_MVA)
                DGEN_DYN = _ensure_gen_dyn_ids_from_raw_if_needed(path_input, gen_dyn_file, DGEN_DYN)
                validate_gen_dyn_partition!(DGEN, DGEN_DYN, DGFM)
            else
                DGFM = nothing
                Check_Gen_Sta_Dyn_Data(nGEN, nrow(DGEN_DYN))
                validate_dgen_dyn_row_ids!(DGEN, DGEN_DYN)
            end
        else
            DGEN_DYN = nothing
            DGFM = nothing
        end
    end

    # Optional per-load table: splits a bus's demand into individual loads, each with
    # its own ZIP shares. Read after the bus remap (its `bus` column carries the original
    # CSV labels) and validated against the UNSCALED bus demand it claims to split.
    DLOAD = Read_Load_Data(path_input)
    if DLOAD !== nothing
        remap_load_buses!(DLOAD, bus_mapping)
        validate_load_table!(DBUS, DLOAD)
    end

    # Scale demand (vectorised, in place on the freshly-read DataFrame columns).
    # The load table moves with the bus aggregate, so their sums stay consistent.
    DBUS.p_d = cfg.load_factor .* DBUS.p_d
    DBUS.q_d = cfg.load_factor .* DBUS.q_d
    if DLOAD !== nothing
        DLOAD.p_d = cfg.load_factor .* DLOAD.p_d
        DLOAD.q_d = cfg.load_factor .* DLOAD.q_d
    end

    bus_gen_circ_dict, bus_gen_circ_dict_ON = Organize_Bus_Gen_Circ(DBUS, DGEN, DCIR)

    return SystemData(DBUS, DGEN, DGEN_DYN, DGFM, DLOAD, DCIR, bus_mapping, reverse_bus_mapping,
                      nBUS, nGEN, nCIR, bus_gen_circ_dict, bus_gen_circ_dict_ON,
                      path_input)
end

"""
If the SG dynamic CSV has an `id` column, replace auto `1:n` ids with those values.
Needed when `allow_gfm` partitions the fleet and SG rows are a subset of `DGEN`.
"""
function _ensure_gen_dyn_ids_from_raw_if_needed(
    folder_path::String,
    filename::String,
    DGEN_DYN::DataFrame,
)::DataFrame
    path = joinpath(folder_path, filename)
    isfile(path) || return DGEN_DYN
    df_raw = CSV.read(path, DataFrame; delim=';')
    headers = lowercase.(strip.(string.(names(df_raw))))
    id_j = findfirst(==("id"), headers)
    id_j === nothing && return DGEN_DYN
    nrow(df_raw) == nrow(DGEN_DYN) || return DGEN_DYN
    DGEN_DYN = copy(DGEN_DYN)
    DGEN_DYN.id = Int64.(df_raw[!, id_j])
    return DGEN_DYN
end

"""
    resolve_fullbus_coupling_hints!(cfg, model, opf_dict, obj_function_MVA,
                                    path_names, sys, solver_log_warmstart)

Solve the steady-state ACOPF that seeds every FULL_BUS coupling variable and return
its operating point as [`SteadyStateHints`](@ref).

The pre-solve is **mandatory** on FULL_BUS: the dynamic variables (E, δ, Ed/Eq,
Id/Iq, bus V/θ) are initialised from it, and a flat guess leaves the joint NLP far
enough from any equilibrium that Ipopt typically stalls at iteration 0. A failed
pre-solve therefore aborts the run rather than falling back.
"""
function resolve_fullbus_coupling_hints!(
    cfg::RunConfig,
    model::Model,
    opf_dict::OrderedDict{Symbol, Any},
    obj_function_MVA,
    path_names::OrderedDict{Symbol, String},
    sys::SystemData,
    solver_log_warmstart::String,
)::SteadyStateHints
    println("\n--- FULL_BUS: solving steady-state ACOPF for warm start ---")
    set_solver_log_path!(model, cfg.solver_name, solver_log_warmstart)
    optimize!(model)
    ws_status = termination_status(model)
    if ws_status != MOI.OPTIMAL && ws_status != MOI.LOCALLY_SOLVED
        throw(ArgumentError(
            "FULL_BUS requires a successful ACOPF warm start before TS assembly " *
            "(termination status: $ws_status). Check the case data and the steady-state " *
            "limits, or relax the Ipopt tolerances; see $solver_log_warmstart."))
    end
    hints = extract_opf_solved_hints(opf_dict)
    println("ACOPF warm start succeeded — injecting solved V/θ/P_g/Q_g into TS starts.")
    println("Warm-start objective: $(JuMP.value.(obj_function_MVA))\n")
    cfg.save_warmstart_dispatch && save_warmstart_dispatch!(
        cfg, path_names, model, obj_function_MVA, opf_dict, sys;
        type_model="ACOPF",
        model_label=cfg.dispatch.use_matrix ? "AC-OPF (Ybus, warm start)" : "AC-OPF (warm start)")
    return hints
end

"""
    save_warmstart_dispatch!(cfg, path_names, model, obj_function_MVA, opf_dict, sys;
                             type_model, model_label)

Persist the steady-state pre-solve that precedes TS assembly to `Dispatch_WarmStart/`.

Used by both pre-solving paths — the FULL_BUS ACOPF warm start and the TSC-DCOPF
solve that fixes the Taylor anchor `δ_ref` — because `Export_OPF_Model` and
`Save_Solution_Optimal_Dispatch` already select their content from `opf_dict`.
"""
function save_warmstart_dispatch!(
    cfg::RunConfig,
    path_names::OrderedDict{Symbol, String},
    model::Model,
    obj_function_MVA,
    opf_dict::OrderedDict{Symbol, Any},
    sys::SystemData;
    type_model::String,
    model_label::String,
)
    ws_paths = dispatch_warmstart_path_names(path_names)
    mkpath(ws_paths[:pf_dispatch])
    mkpath(ws_paths[:pf_dispatch_CSV])
    cfg.save_duals && mkpath(ws_paths[:pf_dispatch_CSV_duals])

    Export_OPF_Model(model, ws_paths, obj_function_MVA, opf_dict; model_label=model_label)
    Save_Solution_Optimal_Dispatch(
        ws_paths, model, type_model, cfg.dispatch.use_matrix,
        obj_function_MVA, opf_dict,
        sys.bus_gen_circ_dict, sys.DBUS, sys.DGEN, sys.DCIR,
        cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR,
        sys.bus_mapping, sys.reverse_bus_mapping;
        _save_duals=cfg.save_duals)
    println("Warm-start dispatch ($type_model) saved in: ", ws_paths[:pf_dispatch])
    return nothing
end

"""
    save_delta_ref_csv!(path_names, δ_ref)

Write the TSC-DCOPF linearisation anchor to `Dispatch_WarmStart/CSV/delta_ref.csv`.

`δ_ref[g] = P_g·X'_d + θ_bus(g)` from the DC pre-solve is the point the electrical
power is Taylor-expanded around, so it is a property of the *run*, not of the final
solution: nothing downstream re-derives it.
"""
function save_delta_ref_csv!(
    path_names::OrderedDict{Symbol, String},
    δ_ref::OrderedDict{Int64, Float64},
)
    isempty(δ_ref) && return nothing
    ws_csv = path_names[:pf_dispatch_warmstart_CSV]
    mkpath(ws_csv)
    gens = sort(collect(keys(δ_ref)))
    table = DataFrame(
        gen = gens,
        delta_ref_rad = [δ_ref[g] for g in gens],
        delta_ref_deg = [rad2deg(δ_ref[g]) for g in gens],
    )
    CSV.write(joinpath(ws_csv, "delta_ref.csv"), table; delim=';')
    println("TSC-DCOPF linearisation anchor saved to: ", joinpath(ws_csv, "delta_ref.csv"))
    return nothing
end

"""
    extra_sc_fault(primary::FaultConfig, id::Int)

Short-circuit `FaultConfig` for one row of `TransientConfig.extra_contingency_ids`. It keeps
the faulted-bus shunt of the primary fault (`sc_fault_admittance_pu`, `sc_fault_stamp`), so
every contingency of a run uses the same fault model instead of the package default.
"""
extra_sc_fault(primary::FaultConfig, id::Int) = FaultConfig(
    fault_type=SC, contingency_id=id,
    sc_fault_admittance_pu=primary.sc_fault_admittance_pu,
    sc_fault_stamp=primary.sc_fault_stamp)

"""
    run_case!(cfg::RunConfig, sys::SystemData, path_main::String, path_folder_results::String)

Build, solve, and save one TSC-OPF run: validate `cfg`, resolve the results
folder layout, build the steady-state (and, when `cfg.trans_stab = true`,
transient-stability) JuMP model, solve it, and write outputs under
`path_folder_results` (conditional subfolders — see Running a case, section 2).

Validation runs in order before any model is built: `Check_Coherence_Input_Data`
(model/solver pairing), `validate_run_config!`, and, for TSC runs,
`validate_fault_config!`.

Returns a `NamedTuple` with `status` (the `MOI` termination status),
`obj_MVA`, `RBUS`, `RGEN`, `RCIR` (result `DataFrame`s), `dyn_model_dict`,
`dyn_parameters_dict`, `t_build`, `t_solve`, and `path_names` (the resolved
results paths).

!!! warning "The JuMP model is released before this function returns"
    `release_solver_backend!` calls `empty!(model)` on Ipopt paths, so **every**
    `JuMP.VariableRef` created during the run is invalid once `run_case!` returns.
    `obj_MVA` is therefore a plain **number** (the already-evaluated objective, or
    `nothing` when the solve did not converge) — do not wrap it in `JuMP.value`.
    The same applies to `dyn_model_dict`: its `:vars` entries are dead references
    afterwards. Query variable values inside the builders, or read them back from
    the exported CSVs under `path_names`.

`path_folder_results` is typically `joinpath(path_main, "RESULTS")`.
"""
function run_case!(cfg::RunConfig, sys::SystemData,
                   path_main::String, path_folder_results::String)

    # --- coherence / sanity checks (fail fast before building anything) -------
    # Throws on an invalid (type_model, solver, trans_stab) combo, e.g. an
    # ACOPF / TSC-ACOPF run with Gurobi or HiGHS (NLP solver required).
    Check_Coherence_Input_Data(cfg.trans_stab, cfg.dispatch.type_model, cfg.solver_name)
    validate_run_config!(cfg)
    if cfg.trans_stab
        validate_fault_config!(
            cfg.transient.dyn_model.fault, sys.DGEN, sys.DBUS, sys.DCIR)
        for id in cfg.transient.extra_contingency_ids
            validate_fault_config!(extra_sc_fault(cfg.transient.dyn_model.fault, id),
                sys.DGEN, sys.DBUS, sys.DCIR)
        end
        # Machine data must satisfy the selected dyn model before anything is built:
        # a missing governor/AVR/dq column reaches the builders as NaN, not as an error.
        validate_dyn_data!(cfg, sys.DGEN_DYN)
        # A machine-referenced δ corridor names a generator that must survive the
        # disturbance — checked here so a bad id fails before the warm-start ACOPF solve.
        validate_δ_reference!(cfg, sys.DGEN, sys.DGEN_DYN, sys.DGFM)
        # A converter-only fleet needs an explicit opt-in and rules out a few corridor
        # and control combinations outright. No-op on every mixed or SG-only fleet.
        validate_all_gfm_fleet!(cfg, sys.DGEN, sys.DGEN_DYN, sys.DGFM)
        # The pre-fault dispatch has to sit inside the device band the overload PIs
        # initialise against, whether or not the knob is on. DGFM is nothing on an
        # SG-only run, so this is a no-op there.
        if sys.DGFM !== nothing
            validate_gfm_dispatch_band!(sys.DGEN, sys.DGFM, cfg.base_MVA)
        end
    end

    # Resolve the optimization-matrix export request. Forced false (with a warning)
    # for TSC runs, where the Jacobian/Hessian are huge. Non-interactive — never
    # blocks a sweep/test/batch job. The resolved value is what we both act on and
    # record in the input-parameters file, so the log reflects what actually ran.
    save_optim_matrices_eff =
        resolve_save_optim_matrices(cfg.save_optim_matrices, cfg.trans_stab)

    # --- results folder for this run (simulation-dependent subfolders only) ---
    path_names = build_results_paths(path_main, path_folder_results, cfg)
    path_names[:pf_input_files] = sys.path_input
    # Toy / in-memory systems may have no INPUT_FILES folder — skip archival.
    if !isempty(sys.path_input) && isdir(sys.path_input)
        Copy_Input_CSVs_To_Results!(
            path_names;
            trans_stab = cfg.trans_stab,
            gen_dynamic_filename = cfg.transient === nothing ?
                "gen_dynamic_data.csv" : cfg.transient.gen_dynamic_filename,
            gfm_dynamic_filename = (cfg.transient !== nothing &&
                cfg.transient.dyn_model.allow_gfm) ?
                cfg.transient.gfm_dynamic_filename : nothing,
            matpower_file = cfg.matpower_file,
        )
    end

    # Archive the script that configured this run, alongside input_parameters.txt.
    # Done before the solve on purpose: a run that fails or hits the iteration limit
    # is exactly when the source configuration matters most.
    if cfg.run_script !== nothing
        cp(cfg.run_script,
           joinpath(path_names[:pf_results_date], basename(cfg.run_script));
           force = true)   # force: overwrite_results=true reuses one folder across runs
    end

    # --- optimiser setup ------------------------------------------------------
    model = Setup_Optim_Model(
        cfg.solver_name;
        ipopt=cfg.ipopt,
        highs=cfg.highs,
        gurobi=cfg.gurobi,
        madnlp=cfg.madnlp,
        silent=cfg.silent_solver,
    )
    JuMP.set_time_limit_sec(model, cfg.time_limit_sec)
    # Converter reactive nameplates are far tighter than the machine Q boxes.
    # With the converter inside the COI and the Δω corridor on, Ipopt's default
    # barrier initialisation (mu_init = 0.1) cycles in restoration — the faulted-bus
    # balance stays ~17 pu off — even though the optimum is interior to that
    # nameplate. A smaller initial barrier parameter recovers the same point.
    # Example 12c (AVR + governor, same corridors) is unchanged at 1e-11.
    # An explicit `IpoptSolverConfig.mu_init` wins over this default.
    dyn = cfg.transient === nothing ? nothing : cfg.transient.dyn_model
    gfm_coi_speed = dyn !== nothing && dyn.allow_gfm && dyn.gfm_in_coi &&
        dyn.constrain_Δω && sys.DGFM !== nothing && nrow(sys.DGFM) > 0
    if gfm_coi_speed && cfg.ipopt.mu_init === nothing &&
            startswith(cfg.solver_name, "Ipopt")
        JuMP.set_optimizer_attribute(model, "mu_init", 1.0e-3)
    end
    log_dir = path_names[:pf_results_date]
    solver_log_main = joinpath(log_dir, "solver_log.txt")
    solver_log_warmstart = joinpath(log_dir, "solver_log_warmstart.txt")

    # --- build steady-state model --------------------------------------------
    t_build = time()
    model, obj_function, obj_function_MVA, opf_dict = Build_SteadyState_Model!(
        path_names, model, cfg.dispatch,
        sys.DBUS, sys.DGEN, sys.DCIR, sys.bus_gen_circ_dict_ON,
        cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR;
        save_matrices=cfg.save_matrices)

    type_model = cfg.dispatch.type_model

    # Phase G1: GFM steady-state ACOPF current / internal-E limits (system-base Imax).
    if type_model == "ACOPF" &&
       cfg.trans_stab &&
       cfg.transient !== nothing &&
       cfg.transient.dyn_model.allow_gfm &&
       sys.DGFM !== nothing &&
       nrow(sys.DGFM) > 0
        attach_gfm_acopf_limits!(model, opf_dict, sys.DGEN, sys.DGFM, cfg.base_MVA)
    end

    # --- DCOPF steady-state objective before TS (for explicit dual verification) ---
    primal_dcopf_f_star = nothing
    δ_ref = OrderedDict{Int64, Float64}()
    steady_state_hints = nothing
    if cfg.trans_stab && type_model == "DCOPF"
        set_solver_log_path!(model, cfg.solver_name, solver_log_warmstart)
        optimize!(model)
        status = termination_status(model)
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
            primal_dcopf_f_star = objective_value(model)
            Pg_sol = [JuMP.value(v) for (i, v) in opf_dict[:vars][:P_g]]
            θ_sol  = [JuMP.value(v) for (i, v) in opf_dict[:vars][:θ]]
            active_gen = findall(x -> x == 1, sys.DGEN.g_status)
            for gen in active_gen
                bus = sys.DGEN.bus[gen]
                Xd_prime = sys.DGEN_DYN.Xd_tr[gen]   # transient reactance
                δ_ref[gen] = Pg_sol[gen] * Xd_prime + θ_sol[bus]
            end
            println("Steady-state objective: $(JuMP.value.(obj_function_MVA)) \n")
            # This solve is a genuine pre-solve: its dispatch is what δ_ref is built
            # from, and the TS constraints below linearise around it.
            if cfg.save_warmstart_dispatch
                save_warmstart_dispatch!(
                    cfg, path_names, model, obj_function_MVA, opf_dict, sys;
                    type_model="DCOPF",
                    model_label=cfg.dispatch.use_matrix ? "DC-OPF (Bbus, warm start)" :
                                                          "DC-OPF (warm start)")
                save_delta_ref_csv!(path_names, δ_ref)
            end
        end
    elseif cfg.trans_stab && type_model == "ACOPF" &&
           cfg.transient.dyn_model.network_form == FULL_BUS
        steady_state_hints = resolve_fullbus_coupling_hints!(
            cfg, model, opf_dict, obj_function_MVA, path_names, sys, solver_log_warmstart)
    end

    # --- add transient-stability constraints ---------------------------------
    # δ_tol and contingency_id are threaded in explicitly (no globals, no runtime
    # redefinition of Common_Parameters_4_TS).
    dyn_model_dict = nothing
    dyn_parameters_dict = nothing
    extra_dyn_blocks = []   # (dyn_model_dict, dyn_parameters_dict, path_names) per extra contingency
    if cfg.trans_stab
        # Machine-data columns were checked in the sanity block above (validate_dyn_data!),
        # which covers every gen_order and control flag, not just DQ_4TH, and fires before
        # the warm-start ACOPF solve rather than after it.
        # Factory dispatch (Phase 1.2): DCOPF → linearized Pe; ACOPF → nonlinear Pe
        linearize = type_model == "DCOPF"
        model, dyn_model_dict, dyn_parameters_dict = Build_Dynamic_Model!(
            path_names, model, opf_dict,
            sys.DBUS, sys.DGEN, sys.DGEN_DYN, sys.DCIR, sys.bus_gen_circ_dict_ON,
            cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR, cfg.transient;
            linearize=linearize,
            δ_ref=linearize ? δ_ref : nothing,
            steady_state_hints=steady_state_hints,
            DGFM=sys.DGFM,
            DLOAD=sys.DLOAD,
        )
        # Further contingencies: one more transient block each, on the same pre-fault
        # variables, with its own results folder `Transient_Stability_c<id>/`.
        for id in cfg.transient.extra_contingency_ids
            dm = cfg.transient.dyn_model
            transient_c = reconfigure_transient(cfg.transient;
                dyn_model=DynModelConfig(;
                    (f => f == :fault ? extra_sc_fault(dm.fault, id) :
                          getfield(dm, f) for f in fieldnames(DynModelConfig))...),
                extra_contingency_ids=Int[])
            path_names_c = copy(path_names)
            ts_c = joinpath(path_names[:pf_results_date], "Transient_Stability_c$(id)")
            path_names_c[:pf_TS] = ts_c
            path_names_c[:pf_TS_figures] = joinpath(ts_c, "Figures")
            path_names_c[:pf_TS_figures_duals] = joinpath(ts_c, "Figures_Duals")
            path_names_c[:pf_TS_CSV] = joinpath(ts_c, "CSV")
            path_names_c[:pf_TS_CSV_duals] = joinpath(ts_c, "CSV_duals")
            path_names_c[:pf_bus_matrices] = joinpath(ts_c, "Bus_Matrices")
            for k in (:pf_TS_figures, :pf_TS_figures_duals, :pf_TS_CSV, :pf_TS_CSV_duals,
                      :pf_bus_matrices)
                mkpath(path_names_c[k])
            end
            model, dict_c, params_c = Build_Dynamic_Model!(
                path_names_c, model, opf_dict,
                sys.DBUS, sys.DGEN, sys.DGEN_DYN, sys.DCIR, sys.bus_gen_circ_dict_ON,
                cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR, transient_c;
                linearize=linearize,
                δ_ref=linearize ? δ_ref : nothing,
                steady_state_hints=steady_state_hints,
                DGFM=sys.DGFM,
                DLOAD=sys.DLOAD,
            )
            push!(extra_dyn_blocks, (dict_c, params_c, path_names_c))
        end
        # After TS assembly, Ipopt cold-starts from `start=` (not the prior ACOPF
        # primal). Re-stamp solved V/θ/P_g/Q_g so the joint NLP begins at the reference
        # warm-start operating point (avoids obj=0 / inf_pr~1e10 at iter 0).
        if steady_state_hints !== nothing
            apply_steady_state_hints_to_opf!(opf_dict, steady_state_hints)
        end
        # Snapshot pre-fault coupling JuMP start= values under Dispatch_WarmStart/
        # (before joint optimize!). Independent of save_warmstart_dispatch, but only
        # meaningful where the starts came from a solved ACOPF: off FULL_BUS the Kron
        # builders seed δ=0 / E=1 constants, so the file would record nothing.
        if steady_state_hints !== nothing && dyn_model_dict !== nothing
            Save_Prefault_Coupling_Starts!(dyn_model_dict, path_names)
        end
    end
    t_build = time() - t_build

    # --- solve the full problem ----------------------------------------------
    t_solve = time()
    set_solver_log_path!(model, cfg.solver_name, solver_log_main)
    if cfg.solver_name in ("Ipopt", "Ipopt-ma57", "Ipopt-ma97", "Ipopt-pardiso")
        # Prefer previous primal when available; starts from apply_* cover OPF vars.
        JuMP.set_optimizer_attribute(model, "warm_start_init_point", "yes")
    end
    JuMP.optimize!(model)
    t_solve = time() - t_solve
    status_model = JuMP.termination_status(model)
    println("\nTime to build model : $(round(t_build; digits=2)) sec")
    println("Time to solve model : $(round(t_solve; digits=2)) sec")
    println("Termination status  : $status_model \n")
    println("-" ^ 70)

    # --- save results ---------------------------------------------------------
    RBUS = nothing; RGEN = nothing; RCIR = nothing
    status_solved = status_model == OPTIMAL || status_model == LOCALLY_SOLVED ||
                    status_model == ITERATION_LIMIT || status_model == ALMOST_LOCALLY_SOLVED
    if status_solved
        println("Objective value: $(JuMP.value.(obj_function_MVA)) \n")

        RBUS, RGEN, RCIR = Save_Solution_Optimal_Dispatch(
            path_names, model, type_model, cfg.dispatch.use_matrix,
            obj_function_MVA, opf_dict,
            sys.bus_gen_circ_dict, sys.DBUS, sys.DGEN, sys.DCIR,
            cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR,
            sys.bus_mapping, sys.reverse_bus_mapping;
            _save_duals=cfg.save_duals)

        # Export the optimization matrices (Jacobian, Hessian, Lagrangian gradient)
        # as sparse COO CSV when requested. Already resolved to false for TSC runs.
        if save_optim_matrices_eff
            Compute_and_Save_Optimal_Matrices(model, path_names)
        end

        if cfg.trans_stab
            # Model-agnostic dynamic-results saver (handles both TSC-ACOPF and
            # TSC-DCOPF). Renamed from Save_Results_Dynamic_Model_tsred_linear (SF-2).
            Save_Results_Dynamic_Model(
                model, path_names, dyn_model_dict, dyn_parameters_dict,
                cfg.base_MVA, dyn_model_dict[:refs][:P_mech], sys.DGEN_DYN;
                save_ts_plots=cfg.save_ts_plots,
                save_ts_debug_csv=cfg.save_ts_debug_csv)
            if cfg.save_duals
                Save_Duals_Dynamic_Model_tsred(model, path_names, dyn_model_dict,
                    dyn_parameters_dict; save_ts_plots=cfg.save_ts_plots)
            end
            for (dict_c, params_c, path_names_c) in extra_dyn_blocks
                Save_Results_Dynamic_Model(
                    model, path_names_c, dict_c, params_c,
                    cfg.base_MVA, dict_c[:refs][:P_mech], sys.DGEN_DYN;
                    save_ts_plots=cfg.save_ts_plots,
                    save_ts_debug_csv=cfg.save_ts_debug_csv)
                if cfg.save_duals
                    Save_Duals_Dynamic_Model_tsred(model, path_names_c, dict_c,
                        params_c; save_ts_plots=cfg.save_ts_plots)
                end
            end
        end

        # --- explicit dual LP (in-house builder, separate Dispatch_Dual/) ---
        if cfg.dispatch.solve_explicit_dual && type_model == "DCOPF"
            f_star = primal_dcopf_f_star
            if f_star === nothing
                f_star = objective_value(model)
            end
            Run_Explicit_DC_OPF_Dual!(
                path_names, cfg.dispatch, model, opf_dict,
                sys.DBUS, sys.DGEN, sys.DCIR, sys.bus_gen_circ_dict_ON,
                cfg.base_MVA, sys.nBUS, sys.nGEN, sys.nCIR;
                primal_f_star = f_star,
                solver_name = cfg.solver_name,
                silent_solver = cfg.silent_solver,
                highs = cfg.highs,
                gurobi = cfg.gurobi,
                cross_check_primal_duals = !cfg.trans_stab,
            )
        elseif cfg.dispatch.solve_explicit_dual && type_model == "ED"
            Run_Explicit_ED_Dual!(
                path_names, cfg.dispatch, model, opf_dict,
                sys.DGEN, sys.bus_gen_circ_dict_ON,
                cfg.base_MVA, sys.nGEN;
                primal_f_star = objective_value(model),
                solver_name = cfg.solver_name,
                silent_solver = cfg.silent_solver,
                highs = cfg.highs,
                gurobi = cfg.gurobi,
            )
        end
    else
        @warn "Optimisation did not converge (status: $status_model). Results not saved."
    end

    # --- save the input parameters used in this run --------------------------
    fault = cfg.transient === nothing ? nothing : cfg.transient.dyn_model.fault
    if dyn_parameters_dict !== nothing
        Print_Input_Parameters(path_names, cfg.overwrite_results, cfg.trans_stab,
            type_model, cfg.solver_name, cfg.dispatch.use_matrix, cfg.silent_solver,
            cfg.case, cfg.base_MVA, cfg.load_factor, save_optim_matrices_eff;
            fault=fault, dyn_parameters_dict=dyn_parameters_dict)
    else
        Print_Input_Parameters(path_names, cfg.overwrite_results, cfg.trans_stab,
            type_model, cfg.solver_name, cfg.dispatch.use_matrix, cfg.silent_solver,
            cfg.case, cfg.base_MVA, cfg.load_factor, save_optim_matrices_eff;
            fault=fault)
    end

    # Snapshot the objective BEFORE the backend is released. `release_solver_backend!`
    # calls `empty!(model)` on Ipopt paths (to free the solver log handle), which
    # invalidates every `VariableRef` — including the ones inside `obj_function_MVA`
    # and inside `dyn_model_dict`. Returning the live expression would hand the caller
    # a reference that throws on `JuMP.value`. Broadcast keeps the shape for the
    # scalar and container forms alike.
    obj_MVA_value = status_solved ? JuMP.value.(obj_function_MVA) : nothing

    # --- machine-readable run description ------------------------------------
    # Written unconditionally, failed runs included: `[status]` then records why the
    # folder is empty. This is the only place `load_factor` and the resolved δ
    # reference / SG-GFM partition are archived in a form a script can read back.
    Save_Run_Manifest!(path_names, cfg, sys;
        opf_dict=opf_dict, dyn_model_dict=dyn_model_dict,
        dyn_parameters_dict=dyn_parameters_dict,
        status=status_model, obj_MVA=obj_MVA_value,
        t_build=t_build, t_solve=t_solve)

    # Last chance to read the assembled model: everything below invalidates it.
    if cfg.post_solve_hook !== nothing
        cfg.post_solve_hook(; cfg=cfg, sys=sys, path_names=path_names, model=model,
            opf_dict=opf_dict, dyn_model_dict=dyn_model_dict,
            dyn_parameters_dict=dyn_parameters_dict,
            status=status_model, status_solved=status_solved,
            obj_MVA=obj_MVA_value, t_build=t_build, t_solve=t_solve)
    end

    release_solver_backend!(model)
    model = nothing

    return (; status=status_model, obj_MVA=obj_MVA_value,
              RBUS, RGEN, RCIR, dyn_model_dict, dyn_parameters_dict,
              t_build, t_solve, path_names)
end
