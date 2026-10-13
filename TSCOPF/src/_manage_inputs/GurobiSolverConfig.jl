#=
================================================================================
 GurobiSolverConfig.jl  —  Gurobi LP/QP/MILP/NLP solver options
================================================================================
 Used when `solver_name` is `"Gurobi"`. Log paths via `set_solver_log_path!` (`LogFile`).
================================================================================
=#

"""
    GurobiSolverConfig

Gurobi options for LP/QP/MILP (and nonlinear when licensed) runs
(`solver_name = "Gurobi"`). Set on `RunConfig.gurobi` or pass to
[`Setup_Optim_Model`](@ref). Use `raw_options` for any additional Gurobi
parameter (e.g. `"Threads" => 4`).

Defaults match the typed knobs below (`MIPGap = 1e-8`, NL barrier tolerances,
`OptimalityTarget = -1`, i.e. left to Gurobi).

!!! warning "`OptimalityTarget` numbering is not stable across Gurobi versions"
    Gurobi 13 accepts only `-1` (automatic), `0` (global) and `1` (**local**),
    where older releases used `1` for global and `2`/`3` for local. Selecting
    local optimization makes Gurobi reject *any* discrete model with
    `Error 10016`, which breaks unit commitment. The default is therefore `-1`:
    do not pin a number unless you know what your Gurobi version means by it.
"""
Base.@kwdef struct GurobiSolverConfig
    # Logging
    output_flag::Int = 1
    # Continuous models
    feasibility_tol::Float64 = 1e-8
    optimality_tol::Float64 = 1e-8
    # Barrier (LP/QP)
    bar_iter_limit::Int = 5_000
    # Mixed-integer models
    mip_gap::Float64 = 1e-8
    # Nonlinear / NL barrier (Gurobi 13+)
    nl_bar_iter_limit::Int = 5_000
    nl_bar_p_feas_tol::Float64 = 1e-8
    nl_bar_d_feas_tol::Float64 = 1e-8
    nl_bar_c_feas_tol::Float64 = 1e-8
    # -1 = leave it to Gurobi. Do NOT default this to a number: the parameter was
    # renumbered in Gurobi 13, where 1 means LOCAL optimization and makes Gurobi
    # reject every discrete model with "Error 10016: Local optimization cannot be
    # used for discrete problems or SOS constraints" — which silently broke UC.
    optimality_target::Int = -1
    raw_options::Dict{String, Any} = Dict{String, Any}()
end

"""Fail fast on invalid `GurobiSolverConfig` field values."""
function validate_gurobi_solver_config!(cfg::GurobiSolverConfig)
    cfg.output_flag in (0, 1) ||
        throw(ArgumentError(
            "GurobiSolverConfig.output_flag must be 0 or 1 (got $(cfg.output_flag))."))
    cfg.optimality_tol > 0 ||
        throw(ArgumentError("GurobiSolverConfig.optimality_tol must be positive."))
    cfg.feasibility_tol > 0 ||
        throw(ArgumentError("GurobiSolverConfig.feasibility_tol must be positive."))
    cfg.bar_iter_limit > 0 ||
        throw(ArgumentError("GurobiSolverConfig.bar_iter_limit must be positive."))
    cfg.mip_gap >= 0 ||
        throw(ArgumentError("GurobiSolverConfig.mip_gap must be non-negative."))
    cfg.nl_bar_iter_limit > 0 ||
        throw(ArgumentError("GurobiSolverConfig.nl_bar_iter_limit must be positive."))
    cfg.nl_bar_p_feas_tol > 0 ||
        throw(ArgumentError("GurobiSolverConfig.nl_bar_p_feas_tol must be positive."))
    cfg.nl_bar_d_feas_tol > 0 ||
        throw(ArgumentError("GurobiSolverConfig.nl_bar_d_feas_tol must be positive."))
    cfg.nl_bar_c_feas_tol > 0 ||
        throw(ArgumentError("GurobiSolverConfig.nl_bar_c_feas_tol must be positive."))
    -1 <= cfg.optimality_target <= 3 ||
        throw(ArgumentError(
            "GurobiSolverConfig.optimality_target must be in -1:3 (got $(cfg.optimality_target))."))
    return nothing
end

"""
    apply_gurobi_options!(model::JuMP.Model, cfg::GurobiSolverConfig; silent::Bool=false)

Stamp Gurobi options from `cfg` onto `model`. When `silent=true`, calls `set_silent`.
"""
function apply_gurobi_options!(
    model::JuMP.Model,
    cfg::GurobiSolverConfig;
    silent::Bool = false,
)::JuMP.Model
    JuMP.set_optimizer_attribute(model, "OutputFlag", cfg.output_flag)
    JuMP.set_optimizer_attribute(model, "FeasibilityTol", cfg.feasibility_tol)
    JuMP.set_optimizer_attribute(model, "OptimalityTol", cfg.optimality_tol)
    JuMP.set_optimizer_attribute(model, "BarIterLimit", cfg.bar_iter_limit)
    JuMP.set_optimizer_attribute(model, "MIPGap", cfg.mip_gap)
    _try_set_gurobi_nl_attribute!(model, "NLBarIterLimit", cfg.nl_bar_iter_limit)
    _try_set_gurobi_nl_attribute!(model, "NLBarPFeasTol", cfg.nl_bar_p_feas_tol)
    _try_set_gurobi_nl_attribute!(model, "NLBarDFeasTol", cfg.nl_bar_d_feas_tol)
    _try_set_gurobi_nl_attribute!(model, "NLBarCFeasTol", cfg.nl_bar_c_feas_tol)
    _try_set_gurobi_nl_attribute!(model, "OptimalityTarget", cfg.optimality_target)
    silent && JuMP.set_silent(model)
    _apply_solver_raw_options!(model, cfg.raw_options)
    return model
end

"""
Best-effort stamp for Gurobi NL-only attributes.

Some Gurobi/MOI combinations reject these parameters on models that do not expose
the nonlinear barrier interface yet (for example MILP UC tests). In that case we
skip the attribute instead of failing the whole model build.

Gurobi.jl records the raw name in `Optimizer.params` *before* throwing
`UnsupportedAttribute`, and `MOI.empty!` later replays that dict. Scrub the
sticky entry so `optimize!` does not resurrect the rejected attribute.
"""
function _try_set_gurobi_nl_attribute!(
    model::JuMP.Model,
    attr::String,
    value,
)::Nothing
    try
        JuMP.set_optimizer_attribute(model, attr, value)
    catch err
        if err isa MOI.UnsupportedAttribute
            _scrub_gurobi_raw_param!(model, attr)
            return nothing
        end
        rethrow()
    end
    return nothing
end

"""Remove a sticky unsupported raw parameter from the Gurobi MOI wrapper."""
function _scrub_gurobi_raw_param!(model::JuMP.Model, attr::String)::Nothing
    try
        opt = JuMP.unsafe_backend(model)
        while hasproperty(opt, :inner) && !hasproperty(opt, :params)
            opt = getproperty(opt, :inner)
        end
        if hasproperty(opt, :params)
            delete!(getproperty(opt, :params), attr)
        end
    catch
        # Backend may not be attached yet; nothing sticky to scrub.
    end
    return nothing
end
