#=
================================================================================
 MadNLPSolverConfig.jl  —  MadNLP interior-point NLP solver options
================================================================================
 Used when `solver_name` is `"MadNLP"` or a MadNLP+HSL backend
 (`"MadNLP-ma57"`, `"MadNLP-ma97"`).
 Log paths via `set_solver_log_path!`.
================================================================================
=#

"""MadNLP `LogLevels` values (see MadNLP.jl `enums.jl`)."""
const MADNLP_LOG_INFO = 3
const MADNLP_LOG_ERROR = 6

# User-facing Hessian mode  =>  MadNLP type name (resolved at model-build time).
#   exact         : analytic second derivatives (ExactHessian), any KKT system
#   compact-lbfgs : limited-memory quasi-Newton, compatible with the *sparse*
#                   KKT system MadNLP uses by default for JuMP models
#   bfgs          : dense BFGS         — requires a DENSE KKT system
#   damped-bfgs   : dense damped BFGS  — requires a DENSE KKT system
# MadNLP's `check_option_sanity` rejects the two dense modes when the KKT system
# is sparse, so those need `raw_options` (see the docstring below).
const _MADNLP_HESSIAN_MODES = Dict(
    "exact" => "ExactHessian",
    "bfgs" => "BFGS",
    "damped-bfgs" => "DampedBFGS",
    "compact-lbfgs" => "CompactLBFGS",
)

"""
    MadNLPSolverConfig

MadNLP options for NLP runs (`solver_name = "MadNLP"` or MadNLP-HSL aliases).
Set on `RunConfig.madnlp` or pass to [`Setup_Optim_Model`](@ref).
`print_level` uses MadNLP `LogLevels` integers (default `3` = INFO). Use
`raw_options` for advanced MadNLP attributes (e.g. `"acceptable_tol" => 1e-5`,
`"ma57_pivtol" => 1e-8`). Linear-solver selection is via `solver_name`, not
`raw_options`. Pardiso is Ipopt-only (`"Ipopt-pardiso"`).

`hessian_approximation` accepts `"exact"` (default), `"compact-lbfgs"`,
`"bfgs"` and `"damped-bfgs"`. Only `"compact-lbfgs"` works with the sparse KKT
system MadNLP builds for JuMP models; the two dense BFGS variants additionally
need

```julia
raw_options = Dict(
    "callback"      => MadNLP.DenseCallback,
    "kkt_system"    => MadNLP.DenseCondensedKKTSystem,
    "linear_solver" => MadNLP.LapackCPUSolver,
)
```

which is only practical on small models. Quasi-Newton tuning (history depth,
initialisation) travels through `raw_options` as well, e.g.
`"quasi_newton_options" => MadNLP.QuasiNewtonOptions(max_history = 20)`.
"""
Base.@kwdef struct MadNLPSolverConfig
    tol::Float64 = 1e-8
    max_iter::Int = 3_000
    print_level::Int = MADNLP_LOG_INFO
    hessian_approximation::String = "exact"
    ma97_num_threads::Union{Nothing, Int} = nothing
    raw_options::Dict{String, Any} = Dict{String, Any}()
end

"""Map user-facing Hessian mode to MadNLP type name."""
function madnlp_hessian_type_name(mode::String)::String
    get(_MADNLP_HESSIAN_MODES, mode) do
        throw(ArgumentError(
            "MadNLPSolverConfig.hessian_approximation must be one of " *
            "$(join(sort(collect(keys(_MADNLP_HESSIAN_MODES))), ", ")) " *
            "(got \"$(mode)\")."))
    end
end

"""Fail fast on invalid `MadNLPSolverConfig` field values."""
function validate_madnlp_solver_config!(cfg::MadNLPSolverConfig)
    cfg.tol > 0 || throw(ArgumentError("MadNLPSolverConfig.tol must be positive."))
    cfg.max_iter > 0 ||
        throw(ArgumentError("MadNLPSolverConfig.max_iter must be positive."))
    1 <= cfg.print_level <= 6 ||
        throw(ArgumentError("MadNLPSolverConfig.print_level must be in 1:6 (MadNLP LogLevels)."))
    madnlp_hessian_type_name(cfg.hessian_approximation)
    if cfg.ma97_num_threads !== nothing
        cfg.ma97_num_threads > 0 || throw(ArgumentError(
            "MadNLPSolverConfig.ma97_num_threads must be positive when set."))
    end
    return nothing
end

"""
    apply_madnlp_options!(model::JuMP.Model, cfg::MadNLPSolverConfig; silent::Bool=false)

Stamp MadNLP options from `cfg` onto `model`. When `silent=true`, forces quiet
`print_level` (ERROR).
"""
function apply_madnlp_options!(
    model::JuMP.Model,
    cfg::MadNLPSolverConfig;
    silent::Bool = false,
)::JuMP.Model
    JuMP.set_optimizer_attribute(model, "tol", cfg.tol)
    JuMP.set_optimizer_attribute(model, "max_iter", cfg.max_iter)
    JuMP.set_optimizer_attribute(
        model, "print_level", silent ? MADNLP_LOG_ERROR : cfg.print_level)
    if cfg.hessian_approximation != "exact"
        # `MadNLPOptions.hessian_approximation` is a `::Type` field and MadNLP
        # has no String parser for it, so resolve the type object first.
        JuMP.set_optimizer_attribute(
            model,
            "hessian_approximation",
            resolve_madnlp_type(
                madnlp_hessian_type_name(cfg.hessian_approximation);
                hint = "Ensure the MadNLP version in the environment provides it.",
            ),
        )
    end
    _apply_solver_raw_options!(model, cfg.raw_options)
    return model
end
