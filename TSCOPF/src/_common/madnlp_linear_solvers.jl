# ==============================================================================
#  Optional MadNLP linear solvers (HSL ma57/ma97)
# ==============================================================================
# MadNLP's default sparse linear solver is MUMPS.  With optional packages in the
# active environment:
#   - MadNLPHSL + HSL_jll  →  solver_name "MadNLP-ma57" / "MadNLP-ma97"
#
# Unlike Ipopt, MadNLP does not take an `hsllib` path attribute.  Library loading
# is delegated to MadNLPHSL / HSL_jll.  Pardiso remains Ipopt-only
# (`solver_name = "Ipopt-pardiso"`).
# ==============================================================================

"""MadNLP + HSL linear-solver `solver_name` values (no collision with Ipopt HSL)."""
const MADNLP_HSL_LINEAR_SOLVERS = ("MadNLP-ma57", "MadNLP-ma97")

"""All MadNLP-backend `solver_name` values including the default `"MadNLP"`."""
const MADNLP_BACKEND_SOLVERS = ("MadNLP", MADNLP_HSL_LINEAR_SOLVERS...)

const _MADNLP_LINEAR_SOLVER_TYPE = Dict{String, String}(
    "MadNLP-ma57" => "Ma57Solver",
    "MadNLP-ma97" => "Ma97Solver",
)

"""`true` when `solver_name` builds a `MadNLP.Optimizer` (default or HSL)."""
function is_madnlp_backend_solver(solver_name::String)::Bool
    return solver_name in MADNLP_BACKEND_SOLVERS
end

"""Return the MadNLP linear-solver type name for HSL backends; else `nothing`."""
function madnlp_linear_solver_type_name(solver_name::String)::Union{Nothing, String}
    return get(_MADNLP_LINEAR_SOLVER_TYPE, solver_name, nothing)
end

"""`true` when optional `MadNLPHSL` is present in the active environment."""
function madnlp_hsl_ext_available()::Bool
    return !isnothing(Base.identify_package("MadNLPHSL"))
end

"""Package name that must be loaded for a MadNLP HSL `solver_name`, if any."""
function madnlp_linear_solver_package_name(solver_name::String)::Union{Nothing, String}
    return solver_name in MADNLP_HSL_LINEAR_SOLVERS ? "MadNLPHSL" : nothing
end

"""Load `MadNLPHSL` when required by `solver_name`."""
function ensure_madnlp_linear_solver_package!(solver_name::String)
    pkg_name = madnlp_linear_solver_package_name(solver_name)
    pkg_name === nothing && return nothing
    pkgid = Base.identify_package(pkg_name)
    pkgid === nothing && throw(ArgumentError(
        "solver_name \"$solver_name\" requires the optional $pkg_name package " *
        "(`using Pkg; Pkg.add(\"$pkg_name\")`) plus HSL_jll and a valid HSL license."))
    try
        Base.require(pkgid)
    catch err
        throw(ArgumentError(
            "$pkg_name is listed in the environment but could not be loaded ($err). " *
            "Check HSL license / HSL_jll."))
    end
    return nothing
end

"""
    resolve_madnlp_type(type_name; extra_packages, hint)

Resolve a MadNLP type (linear solver, Hessian approximation, KKT system, ...) by
name, once the owning package is loaded.  `extra_packages` are searched before
`MadNLP` itself — that is how the HSL solver types are found inside `MadNLPHSL`.
MadNLP options such as `linear_solver` and `hessian_approximation` are typed
`::Type` on the MadNLP side, so passing the name as a `String` raises a
`MethodError` inside `MadNLP.set_options!`; always stamp the resolved type.
"""
function resolve_madnlp_type(
    type_name::String;
    extra_packages::Vector{String} = String[],
    hint::String = "",
)
    sym = Symbol(type_name)
    # MadNLPHSL typically `export`s the solver types into the caller's module on
    # `using`, and may also attach them under its own module or under MadNLP.
    candidates = copy(extra_packages)
    push!(candidates, "MadNLP")

    for pkg_name in candidates
        pkgid = Base.identify_package(pkg_name)
        pkgid === nothing && continue
        try
            mod = Base.require(pkgid)
            if isdefined(mod, sym)
                return getfield(mod, sym)
            end
        catch
            continue
        end
    end
    if isdefined(Main, sym)
        return getfield(Main, sym)
    end
    throw(ArgumentError(
        "MadNLP type \"$type_name\" is not defined after loading " *
        "$(join(candidates, ", ")). " * hint))
end

"""Resolve a MadNLP linear-solver type by name after the extension package is loaded."""
function resolve_madnlp_linear_solver_type(type_name::String)
    extra = String[]
    if type_name in ("Ma57Solver", "Ma97Solver", "Ma27Solver", "Ma86Solver", "Ma77Solver")
        push!(extra, "MadNLPHSL")
    end
    return resolve_madnlp_type(
        type_name;
        extra_packages = extra,
        hint = "Ensure MadNLPHSL matches your MadNLP version and that " *
               "`using MadNLPHSL` succeeds.",
    )
end

"""Configure MadNLP HSL linear solver when requested."""
function configure_madnlp_linear_solver!(
    model::JuMP.Model,
    solver_name::String,
    madnlp::MadNLPSolverConfig,
)::JuMP.Model
    type_name = madnlp_linear_solver_type_name(solver_name)
    if type_name !== nothing
        ensure_madnlp_linear_solver_package!(solver_name)
        ls_type = resolve_madnlp_linear_solver_type(type_name)
        JuMP.set_optimizer_attribute(model, "linear_solver", ls_type)
    end
    if madnlp.ma97_num_threads !== nothing
        JuMP.set_optimizer_attribute(model, "ma97_num_threads", madnlp.ma97_num_threads)
    end
    return model
end
