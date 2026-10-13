module TSCOPFGurobiExt

using JuMP
using Gurobi
using TSCOPF

function _build_gurobi_model(;
    gurobi::GurobiSolverConfig,
    silent::Bool,
    solver_name::String,
)::JuMP.Model
    solver_name == "Gurobi" ||
        throw(ArgumentError("internal Gurobi builder called with solver_name=\"$solver_name\""))
    model = JuMP.Model(Gurobi.Optimizer)
    apply_gurobi_options!(model, gurobi; silent=silent)
    return model
end

function _gurobi_license_available()::Bool
    try
        env = Gurobi.Env()
        return env !== nothing
    catch
        return false
    end
end

function __init__()
    TSCOPF.register_optimizer_builder!("Gurobi", _build_gurobi_model, "LogFile")
    TSCOPF.register_gurobi_availability_check!(_gurobi_license_available)
    return nothing
end

end
