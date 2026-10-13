module TSCOPFMadNLPExt

using JuMP
using MadNLP
using TSCOPF

function _build_madnlp_model(;
    madnlp::MadNLPSolverConfig,
    silent::Bool,
    solver_name::String,
)::JuMP.Model
    is_madnlp_backend_solver(solver_name) ||
        throw(ArgumentError("internal MadNLP builder called with solver_name=\"$solver_name\""))
    model = JuMP.Model(MadNLP.Optimizer)
    apply_madnlp_options!(model, madnlp; silent=silent)
    configure_madnlp_linear_solver!(model, solver_name, madnlp)
    return model
end

function __init__()
    for name in MADNLP_BACKEND_SOLVERS
        TSCOPF.register_optimizer_builder!(name, _build_madnlp_model, "output_file")
    end
    return nothing
end

end
