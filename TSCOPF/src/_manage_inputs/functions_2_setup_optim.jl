# Function to setup the Model
"""
    Setup_Optim_Model(solver_name; ipopt, highs, gurobi, madnlp, silent)

Build a JuMP model with the registered optimizer for `solver_name` and stamp
solver-specific options from the matching config struct (`IpoptSolverConfig`,
`HiGHSSolverConfig`, etc.).
"""
function Setup_Optim_Model(
    solver_name::String;
    ipopt::IpoptSolverConfig = IpoptSolverConfig(),
    highs::HiGHSSolverConfig = HiGHSSolverConfig(),
    gurobi::GurobiSolverConfig = GurobiSolverConfig(),
    madnlp::MadNLPSolverConfig = MadNLPSolverConfig(),
    silent::Bool = true,
)::JuMP.Model
    builder = _require_solver_builder(solver_name)
    if is_ipopt_backend_solver(solver_name)
        validate_ipopt_solver_config!(ipopt, solver_name)
        return builder(; ipopt=ipopt, silent=silent, solver_name=solver_name)
    elseif solver_name == "HiGHS"
        validate_highs_solver_config!(highs)
        return builder(; highs=highs, silent=silent, solver_name=solver_name)
    elseif solver_name == "Gurobi"
        validate_gurobi_solver_config!(gurobi)
        return builder(; gurobi=gurobi, silent=silent, solver_name=solver_name)
    elseif is_madnlp_backend_solver(solver_name)
        validate_madnlp_solver_config!(madnlp)
        return builder(; madnlp=madnlp, silent=silent, solver_name=solver_name)
    else
        throw(ArgumentError("Unknown solver_name \"$solver_name\" in Setup_Optim_Model."))
    end
end

# Point solver log output at an absolute path (replaces cd() before JuMP.optimize!).
"""
    set_solver_log_path!(model, solver_name, log_path)

Redirect solver text output to `log_path` (`output_file` for Ipopt/MadNLP/PARDISO,
`LogFile` for Gurobi, `log_file` for HiGHS).
"""
function set_solver_log_path!(model::Model, solver_name::String, log_path::String)::JuMP.Model
    attr = get(_SOLVER_LOG_FILE_ATTR, solver_name, nothing)
    if attr === nothing
        if is_ipopt_backend_solver(solver_name) || is_madnlp_backend_solver(solver_name)
            JuMP.set_optimizer_attribute(model, "output_file", log_path)
        end
    else
        JuMP.set_optimizer_attribute(model, attr, log_path)
    end
    return model
end

"""Unwrap `model` to the innermost MOI optimizer (past caching / bridge layers)."""
function innermost_optimizer(model::Model)
    opt = unsafe_backend(model)
    while opt isa MOI.Utilities.CachingOptimizer
        opt = MOI.Utilities.get_optimizer(opt)
    end
    while opt isa MOI.Bridges.LazyBridgeOptimizer
        opt = opt.inner
    end
    return opt
end

"""
    release_solver_backend!(model::Model)

Release native solver resources held by `model`, in particular Ipopt log-file
handles opened via `output_file`. Call only after all `JuMP.value` / `JuMP.dual`
queries needed from this solve are done; later queries may throw
`OptimizeNotCalled()`.

Uses `empty!(model)` so Ipopt's `MOI.empty!` drops the native `IpoptProblem`
(`inner = nothing`), then `GC.gc()` so finalizers run promptly on Windows.
An explicit `finalize` on any surviving Ipopt inner problem is attempted first.
Non-Ipopt backends are left to Julia GC (HiGHS `empty!` teardown can segfault on Windows).
"""
function release_solver_backend!(model::Model)::Nothing
    opt = innermost_optimizer(model)
    if opt isa Ipopt.Optimizer
        inner = opt.inner
        if inner !== nothing
            finalize(inner)
            opt.inner = nothing
        end
        empty!(model)
        GC.gc()
    end
    return nothing
end


# New configuration to use Gurobi in nonlinear/nonconvex problems
function Modify_Model_if_Gurobi_4_NonLinear_Problem(model::Model)::JuMP.Model
        JuMP.set_optimizer_attribute(model, "NonConvex", 2)
        JuMP.set_optimizer_attribute(model, "QCPDual", 1)
        JuMP.set_optimizer_attribute(model, "NumericFocus", 3)
        # MIPGap / NL barrier / OptimalityTarget come from GurobiSolverConfig

    return model
end
