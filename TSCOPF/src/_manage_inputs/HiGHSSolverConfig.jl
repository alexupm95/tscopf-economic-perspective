#=
================================================================================
 HiGHSSolverConfig.jl  —  HiGHS LP/QP solver options
================================================================================
 Used when `solver_name` is `"HiGHS"`. Log paths via `set_solver_log_path!`.
================================================================================
=#

"""
    HiGHSSolverConfig

HiGHS options for LP/QP runs (`solver_name = "HiGHS"`). Set on `RunConfig.highs`
or pass to [`Setup_Optim_Model`](@ref). Use `raw_options` for any additional
HiGHS MOI attribute (e.g. `"time_limit" => 120.0`).
"""
Base.@kwdef struct HiGHSSolverConfig
    output_flag::Bool = true
    primal_feasibility_tolerance::Float64 = 1e-8
    dual_feasibility_tolerance::Float64 = 1e-8
    ipm_optimality_tolerance::Float64 = 1e-8
    simplex_iteration_limit::Int = 5_000
    ipm_iteration_limit::Int = 5_000
    solver::String = "choose"
    raw_options::Dict{String, Any} = Dict{String, Any}()
end

const _HIGHS_SOLVER_MODES = ("choose", "simplex", "ipm", "pdlp")

"""Fail fast on invalid `HiGHSSolverConfig` field values."""
function validate_highs_solver_config!(cfg::HiGHSSolverConfig)
    cfg.primal_feasibility_tolerance > 0 ||
        throw(ArgumentError("HiGHSSolverConfig.primal_feasibility_tolerance must be positive."))
    cfg.dual_feasibility_tolerance > 0 ||
        throw(ArgumentError("HiGHSSolverConfig.dual_feasibility_tolerance must be positive."))
    cfg.ipm_optimality_tolerance > 0 ||
        throw(ArgumentError("HiGHSSolverConfig.ipm_optimality_tolerance must be positive."))
    cfg.simplex_iteration_limit > 0 ||
        throw(ArgumentError("HiGHSSolverConfig.simplex_iteration_limit must be positive."))
    cfg.ipm_iteration_limit > 0 ||
        throw(ArgumentError("HiGHSSolverConfig.ipm_iteration_limit must be positive."))
    cfg.solver in _HIGHS_SOLVER_MODES ||
        throw(ArgumentError(
            "HiGHSSolverConfig.solver must be one of $(_HIGHS_SOLVER_MODES) (got \"$(cfg.solver)\")."))
    return nothing
end

"""
    apply_highs_options!(model::JuMP.Model, cfg::HiGHSSolverConfig; silent::Bool=false)

Stamp HiGHS options from `cfg` onto `model`. When `silent=true`, calls `set_silent`.
"""
function apply_highs_options!(
    model::JuMP.Model,
    cfg::HiGHSSolverConfig;
    silent::Bool = false,
)::JuMP.Model
    JuMP.set_optimizer_attribute(model, "output_flag", cfg.output_flag)
    JuMP.set_optimizer_attribute(model, "primal_feasibility_tolerance", cfg.primal_feasibility_tolerance)
    JuMP.set_optimizer_attribute(model, "dual_feasibility_tolerance", cfg.dual_feasibility_tolerance)
    JuMP.set_optimizer_attribute(model, "ipm_optimality_tolerance", cfg.ipm_optimality_tolerance)
    JuMP.set_optimizer_attribute(model, "simplex_iteration_limit", cfg.simplex_iteration_limit)
    JuMP.set_optimizer_attribute(model, "ipm_iteration_limit", cfg.ipm_iteration_limit)
    JuMP.set_optimizer_attribute(model, "solver", cfg.solver)
    silent && JuMP.set_silent(model)
    _apply_solver_raw_options!(model, cfg.raw_options)
    return model
end
