#=
================================================================================
 IpoptSolverConfig.jl  —  Ipopt NLP solver options (replaces legacy config_file.jl)
================================================================================
 Used when `solver_name` is `"Ipopt"`, an Ipopt+HSL linear solver (`Ipopt-ma57`,
 `Ipopt-ma97`), or `"Ipopt-pardiso"`. Log paths are set at solve time via
 `set_solver_log_path!`.
================================================================================
=#

"""
    IpoptSolverConfig

Ipopt options for NLP runs (`solver_name = "Ipopt"`, `"Ipopt-ma57"`, `"Ipopt-ma97"`,
or `"Ipopt-pardiso"`). Set on `RunConfig.ipopt` or pass to [`Setup_Optim_Model`](@ref).

Defaults reproduce the legacy integrated `config_file.jl` tolerances and
iteration cap. `output_file` is **not** stored here: `run_case!` redirects
logs to `RESULTS/.../solver_log.txt` (and `solver_log_warmstart.txt` when a
pre-solve runs).

For PARDISO as the KKT linear solver:

```julia
RunConfig(
    solver_name = "Ipopt-pardiso",
    ipopt = IpoptSolverConfig(
        pardiso_lib_path = raw"C:\\Libraries\\Pardiso\\lib\\libpardiso.dll",
    ),
)
```

For a quasi-Newton Hessian instead of the exact KKT Hessian:

```julia
RunConfig(
    ipopt = IpoptSolverConfig(
        hessian_approximation = "limited-memory",
        limited_memory_max_history = 50,
    ),
)
```
"""
Base.@kwdef struct IpoptSolverConfig
    tol::Float64 = 1e-8
    print_level::Int = 5
    max_iter::Int = 5_000
    constr_viol_tol::Float64 = 1e-8
    dual_inf_tol::Float64 = 1e-8
    compl_inf_tol::Float64 = 1e-8
    hessian_approximation::String = "exact"
    limited_memory_max_history::Int = 50
    # Optional reference-style scaling / acceptability (nothing → leave Ipopt default)
    acceptable_tol::Union{Nothing, Float64} = nothing
    obj_scaling_factor::Union{Nothing, Float64} = nothing
    nlp_scaling_method::Union{Nothing, String} = nothing
    # `nothing` leaves Ipopt's monotone barrier. `"monotone"` or `"adaptive"` forces it.
    mu_strategy::Union{Nothing, String} = nothing
    # `nothing` leaves Ipopt's default 0.1, except when `gfm_in_coi` and
    # `constrain_Δω` are both on: `run_case!` then uses `1e-3`. The converter
    # reactive nameplate makes the default barrier stiff enough to cycle.
    mu_init::Union{Nothing, Float64} = nothing
    pardiso_lib_path::Union{Nothing, String} = nothing
    pardiso_license_message::Bool = true
end

const _IPOPT_HESSIAN_MODES = ("exact", "limited-memory")

"""
    validate_ipopt_solver_config!(cfg::IpoptSolverConfig, solver_name::String="Ipopt")

Fail fast on invalid `IpoptSolverConfig` field values. When `solver_name` is
`"Ipopt-pardiso"`, also checks that a PARDISO library path is configured and exists.
"""
function validate_ipopt_solver_config!(
    cfg::IpoptSolverConfig,
    solver_name::String = "Ipopt",
)
    cfg.tol > 0 || throw(ArgumentError("IpoptSolverConfig.tol must be positive (got $(cfg.tol))."))
    cfg.max_iter > 0 ||
        throw(ArgumentError("IpoptSolverConfig.max_iter must be positive (got $(cfg.max_iter))."))
    0 <= cfg.print_level <= 12 ||
        throw(ArgumentError("IpoptSolverConfig.print_level must be in 0:12 (got $(cfg.print_level))."))
    cfg.constr_viol_tol > 0 ||
        throw(ArgumentError("IpoptSolverConfig.constr_viol_tol must be positive."))
    cfg.dual_inf_tol > 0 ||
        throw(ArgumentError("IpoptSolverConfig.dual_inf_tol must be positive."))
    cfg.compl_inf_tol > 0 ||
        throw(ArgumentError("IpoptSolverConfig.compl_inf_tol must be positive."))
    if cfg.mu_init !== nothing && !(cfg.mu_init > 0)
        throw(ArgumentError(
            "IpoptSolverConfig.mu_init must be positive or nothing (got $(cfg.mu_init))."))
    end
    if cfg.mu_strategy !== nothing && !(cfg.mu_strategy in ("monotone", "adaptive"))
        throw(ArgumentError(
            "IpoptSolverConfig.mu_strategy must be \"monotone\", \"adaptive\", or nothing " *
            "(got $(repr(cfg.mu_strategy)))."))
    end
    cfg.hessian_approximation in _IPOPT_HESSIAN_MODES ||
        throw(ArgumentError(
            "IpoptSolverConfig.hessian_approximation must be \"exact\" or " *
            "\"limited-memory\" (got \"$(cfg.hessian_approximation)\")."))
    if cfg.hessian_approximation == "limited-memory"
        cfg.limited_memory_max_history > 0 ||
            throw(ArgumentError("IpoptSolverConfig.limited_memory_max_history must be positive."))
    end
    if ipopt_pardiso_requested(solver_name)
        path = cfg.pardiso_lib_path !== nothing ?
            cfg.pardiso_lib_path : get(ENV, "JULIA_PARDISO_LIB", "")
        isempty(path) && throw(ArgumentError(
            "solver_name \"Ipopt-pardiso\" requires IpoptSolverConfig.pardiso_lib_path or " *
            "ENV[\"JULIA_PARDISO_LIB\"] pointing to libpardiso (e.g. libpardiso.dll)."))
        isfile(path) || throw(ArgumentError(
            "PARDISO library not found at \"$(path)\". " *
            "Set IpoptSolverConfig.pardiso_lib_path or ENV[\"JULIA_PARDISO_LIB\"]."))
    end
    return nothing
end

"""
    apply_ipopt_options!(model::JuMP.Model, cfg::IpoptSolverConfig; silent::Bool=false)

Stamp Ipopt raw options from `cfg` onto `model`. When `silent=true`, forces
`print_level = 0` regardless of `cfg.print_level`.
"""
function apply_ipopt_options!(
    model::JuMP.Model,
    cfg::IpoptSolverConfig;
    silent::Bool = false,
)::JuMP.Model
    JuMP.set_optimizer_attribute(model, "tol", cfg.tol)
    JuMP.set_optimizer_attribute(model, "print_level", silent ? 0 : cfg.print_level)
    JuMP.set_optimizer_attribute(model, "max_iter", cfg.max_iter)
    JuMP.set_optimizer_attribute(model, "constr_viol_tol", cfg.constr_viol_tol)
    JuMP.set_optimizer_attribute(model, "dual_inf_tol", cfg.dual_inf_tol)
    JuMP.set_optimizer_attribute(model, "compl_inf_tol", cfg.compl_inf_tol)
    if cfg.hessian_approximation == "limited-memory"
        JuMP.set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
        JuMP.set_optimizer_attribute(
            model, "limited_memory_max_history", cfg.limited_memory_max_history)
    end
    if cfg.acceptable_tol !== nothing
        JuMP.set_optimizer_attribute(model, "acceptable_tol", cfg.acceptable_tol)
    end
    if cfg.obj_scaling_factor !== nothing
        JuMP.set_optimizer_attribute(model, "obj_scaling_factor", cfg.obj_scaling_factor)
    end
    if cfg.nlp_scaling_method !== nothing
        JuMP.set_optimizer_attribute(model, "nlp_scaling_method", cfg.nlp_scaling_method)
    end
    if cfg.mu_strategy !== nothing
        JuMP.set_optimizer_attribute(model, "mu_strategy", cfg.mu_strategy)
    end
    if cfg.mu_init !== nothing
        JuMP.set_optimizer_attribute(model, "mu_init", cfg.mu_init)
    end
    return model
end
