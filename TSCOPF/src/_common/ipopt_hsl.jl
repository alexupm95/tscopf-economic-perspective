# ==============================================================================
#  Optional Ipopt linear solvers (HSL ma57/ma97, PARDISO)
# ==============================================================================
# Ipopt's default sparse linear solver is MUMPS.  With a valid HSL license and
# `HSL_jll` in the active environment, `solver_name = "Ipopt-ma57"` or
# `"Ipopt-ma97"` selects Ipopt backed by the corresponding HSL routine.  With a
# licensed PARDISO library, `solver_name = "Ipopt-pardiso"` selects Intel MKL
# PARDISO via `pardisolib`.
#
# RunConfig `solver_name` values are namespaced (`Ipopt-*`) to mirror MadNLP
# aliases.  Ipopt's MOI `linear_solver` attribute still uses the native strings
# `"ma57"`, `"ma97"`, and `"pardiso"`.
# ==============================================================================

"""RunConfig `solver_name` values for Ipopt + HSL (not Ipopt's MOI attribute names)."""
const IPOPT_HSL_LINEAR_SOLVERS = ("Ipopt-ma57", "Ipopt-ma97")

"""RunConfig `solver_name` for Ipopt + Intel MKL PARDISO."""
const IPOPT_PARDISO_SOLVER = "Ipopt-pardiso"

"""Ipopt MOI `linear_solver` string for MKL PARDISO."""
const _IPOPT_PARDISO_LINEAR_SOLVER_ATTR = "pardiso"

const _IPOPT_HSL_ATTR = Dict{String, String}(
    "Ipopt-ma57" => "ma57",
    "Ipopt-ma97" => "ma97",
)

"""`true` when `solver_name` builds an `Ipopt.Optimizer` (default, HSL, or PARDISO)."""
function is_ipopt_backend_solver(solver_name::String)::Bool
    return solver_name == "Ipopt" ||
        solver_name in IPOPT_HSL_LINEAR_SOLVERS ||
        solver_name == IPOPT_PARDISO_SOLVER
end

"""Return Ipopt `linear_solver` `"ma57"` / `"ma97"` when requested; otherwise `nothing`."""
function ipopt_hsl_linear_solver(solver_name::String)::Union{Nothing, String}
    return get(_IPOPT_HSL_ATTR, solver_name, nothing)
end

"""`true` when `solver_name` requests PARDISO as Ipopt's linear solver."""
function ipopt_pardiso_requested(solver_name::String)::Bool
    return solver_name == IPOPT_PARDISO_SOLVER
end

"""`true` when the optional `HSL_jll` artifact is present in the active environment."""
function hsl_jll_available()::Bool
    return !isnothing(Base.identify_package("HSL_jll"))
end

"""Absolute path to `libhsl` for Ipopt's `hsllib` attribute."""
function hsl_jll_lib_path()::String
    pkgid = Base.identify_package("HSL_jll")
    pkgid === nothing && throw(ArgumentError(
        "solver_name \"Ipopt-ma57\" / \"Ipopt-ma97\" requires the optional HSL_jll package " *
        "(`using Pkg; Pkg.add(\"HSL_jll\")`) and a valid HSL license."))
    try
        mod = Base.require(pkgid)
        return mod.libhsl_path::String
    catch err
        throw(ArgumentError(
            "HSL_jll is listed in the environment but could not be loaded ($err). " *
            "Check that your HSL license is configured."))
    end
end

"""Resolve PARDISO shared-library path from config or `ENV["JULIA_PARDISO_LIB"]`."""
function resolve_pardiso_lib_path(ipopt::IpoptSolverConfig)::Union{Nothing, String}
    path = if ipopt.pardiso_lib_path !== nothing
        ipopt.pardiso_lib_path
    else
        get(ENV, "JULIA_PARDISO_LIB", "")
    end
    return isempty(path) ? nothing : path
end

"""Absolute path to `libpardiso` for Ipopt's `pardisolib` attribute."""
function pardiso_lib_path(ipopt::IpoptSolverConfig)::String
    path = resolve_pardiso_lib_path(ipopt)
    path === nothing && throw(ArgumentError(
        "solver_name \"Ipopt-pardiso\" requires IpoptSolverConfig.pardiso_lib_path or " *
        "ENV[\"JULIA_PARDISO_LIB\"] pointing to libpardiso (e.g. libpardiso.dll)."))
    isfile(path) || throw(ArgumentError(
        "PARDISO library not found at \"$(path)\". " *
        "Set IpoptSolverConfig.pardiso_lib_path or ENV[\"JULIA_PARDISO_LIB\"]."))
    return path
end

"""Set process-level PARDISO env vars before building the Ipopt model."""
function setup_pardiso_env!(ipopt::IpoptSolverConfig)::Nothing
    if ipopt.pardiso_license_message
        ENV["PARDISOLICMESSAGE"] = "1"
    end
    if ipopt.pardiso_lib_path !== nothing
        ENV["JULIA_PARDISO_LIB"] = ipopt.pardiso_lib_path
    end
    return nothing
end

"""`true` when a PARDISO library path is configured and the file exists."""
function pardiso_available(ipopt::IpoptSolverConfig = IpoptSolverConfig())::Bool
    path = resolve_pardiso_lib_path(ipopt)
    return path !== nothing && isfile(path)
end

"""Configure Ipopt HSL or PARDISO linear solver when requested."""
function configure_ipopt_linear_solver!(
    model::JuMP.Model,
    solver_name::String,
    ipopt::IpoptSolverConfig,
)::JuMP.Model
    hsl = ipopt_hsl_linear_solver(solver_name)
    if hsl !== nothing
        JuMP.set_optimizer_attribute(model, "hsllib", hsl_jll_lib_path())
        JuMP.set_optimizer_attribute(model, "linear_solver", hsl)
        # the reference implementation: more stable pivoting for difficult KKT systems
        if hsl == "ma57"
            JuMP.set_optimizer_attribute(model, "ma57_automatic_scaling", "yes")
        end
    elseif ipopt_pardiso_requested(solver_name)
        setup_pardiso_env!(ipopt)
        lib = pardiso_lib_path(ipopt)
        JuMP.set_optimizer_attribute(model, "pardisolib", lib)
        JuMP.set_optimizer_attribute(model, "linear_solver", _IPOPT_PARDISO_LINEAR_SOLVER_ATTR)
    end
    return model
end

"""Deprecated alias — use [`configure_ipopt_linear_solver!`](@ref)."""
function configure_ipopt_hsl_linear_solver!(model::JuMP.Model, solver_name::String)::JuMP.Model
    return configure_ipopt_linear_solver!(model, solver_name, IpoptSolverConfig())
end
