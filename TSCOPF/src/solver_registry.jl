# ==============================================================================
#  Optimizer factory registry (core: Ipopt + HiGHS; optional via package extensions)
# ==============================================================================
# Extensions `TSCOPFGurobiExt` / `TSCOPFMadNLPExt` call `register_optimizer_builder!`
# when the user loads Gurobi.jl or MadNLP.jl after `using TSCOPF`.
# ==============================================================================

const _OPTIMIZER_BUILDERS = Dict{String, Function}()
const _SOLVER_LOG_FILE_ATTR = Dict{String, String}()
const _GUROBI_AVAILABLE_FN = Ref{Union{Function, Nothing}}(nothing)

"""Register a JuMP model factory for `solver_name`."""
function register_optimizer_builder!(
    solver_name::String,
    builder::Function,
    log_file_attr::Union{Nothing, String} = nothing,
)
    _OPTIMIZER_BUILDERS[solver_name] = builder
    if log_file_attr !== nothing
        _SOLVER_LOG_FILE_ATTR[solver_name] = log_file_attr
    end
    return nothing
end

"""Register the Gurobi license probe used by `gurobi_available()` (Gurobi extension)."""
function register_gurobi_availability_check!(f::Function)
    _GUROBI_AVAILABLE_FN[] = f
    return nothing
end

function _require_solver_builder(solver_name::String)
    if haskey(_OPTIMIZER_BUILDERS, solver_name)
        return _OPTIMIZER_BUILDERS[solver_name]
    end
    if solver_name == "Gurobi"
        throw(ArgumentError(
            "solver_name \"Gurobi\" requires Gurobi.jl. " *
            "Add it to your environment and run `using Gurobi` before solving."))
    elseif is_madnlp_backend_solver(solver_name)
        throw(ArgumentError(
            "solver_name \"$solver_name\" requires MadNLP.jl. " *
            "Add it to your environment and run `using MadNLP` before solving" *
            (solver_name in MADNLP_HSL_LINEAR_SOLVERS ?
                " (also add MadNLPHSL + HSL_jll for HSL linear solvers)" : "") *
            "."))
    else
        throw(ArgumentError("Unknown solver_name \"$solver_name\" in Setup_Optim_Model."))
    end
end

function _build_ipopt_model(;
    ipopt::IpoptSolverConfig,
    silent::Bool,
    solver_name::String,
)::JuMP.Model
    model = JuMP.Model(Ipopt.Optimizer)
    apply_ipopt_options!(model, ipopt; silent=silent)
    configure_ipopt_linear_solver!(model, solver_name, ipopt)
    return model
end

function _build_highs_model(;
    highs::HiGHSSolverConfig,
    silent::Bool,
    solver_name::String,
)::JuMP.Model
    solver_name == "HiGHS" ||
        throw(ArgumentError("internal HiGHS builder called with solver_name=\"$solver_name\""))
    model = JuMP.Model(HiGHS.Optimizer)
    apply_highs_options!(model, highs; silent=silent)
    return model
end

function _register_core_solver_builders!()
    register_optimizer_builder!("Ipopt", _build_ipopt_model, "output_file")
    for hsl_name in IPOPT_HSL_LINEAR_SOLVERS
        register_optimizer_builder!(hsl_name, _build_ipopt_model, "output_file")
    end
    register_optimizer_builder!(IPOPT_PARDISO_SOLVER, _build_ipopt_model, "output_file")
    register_optimizer_builder!("HiGHS", _build_highs_model, "log_file")
    return nothing
end
