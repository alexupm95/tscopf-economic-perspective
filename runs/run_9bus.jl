#=
================================================================================
 9-BUS RUNS OF THE PAPER (SECTION IV-A)
================================================================================
  Usage, from the root of the repository:

      julia --project=TSCOPF runs/run_9bus.jl <run>

      <run>    linear      INPUT_FILES/9bus_linear_H24_6_3           results/9bus/Results - 2026-10-08 112948
               colocated   INPUT_FILES/9bus_colocated_bus1_H24_6_3   results/9bus/Results - 2026-10-08 113323

  The model is `_9bus_config.jl`, with the time grid and the fault of the paper set
  below. A new `Results - <timestamp>` folder is written under new_runs/.
================================================================================
=#

using TSCOPF
using JuMP

const RUNS = Dict("linear" => "9bus_linear_H24_6_3", "colocated" => "9bus_colocated_bus1_H24_6_3")
length(ARGS) == 1 && haskey(RUNS, ARGS[1]) ||
    error("usage: julia --project=TSCOPF runs/run_9bus.jl <linear|colocated>")

const ROOT         = normpath(joinpath(@__DIR__, ".."))
const PATH_RESULTS = joinpath(ROOT, "new_runs")
mkpath(PATH_RESULTS)

CASE = RUNS[ARGS[1]]                     # read by _9bus_config.jl
include(joinpath(@__DIR__, "_9bus_config.jl"))   # defines cfg

# The fault window keeps the node at t_start_fault + clearing_time, so the fault
# lasts clearing_time + t_step: 0.095 s on a 5 ms grid is a 100 ms fault.
const TSTEP = 0.005
const CLEAR = 0.095

# psi_g export (Dispatch/CSV_duals/kkt_Pg_terms.csv), attached through the post-solve hook
include(joinpath(ROOT, "kkt", "kkt_pg_terms.jl"))
function kkt_post_solve_hook(; cfg, sys, path_names, opf_dict, dyn_model_dict, status_solved, kwargs...)
    (cfg.save_duals && status_solved) || return nothing
    Save_KKT_Pg_Terms!(path_names, cfg, sys, opf_dict, dyn_model_dict)
    return nothing
end

# Inductive solid short (−j1e6 pu added to the faulted bus), the package default
fault = FaultConfig(fault_type = SC, contingency_id = 2)

run_cfg = reconfigure(cfg;
    silent_solver = false, overwrite_results = false,
    save_duals = true, save_matrices = true, save_warmstart_dispatch = true,
    run_script = @__FILE__,
    post_solve_hook = kkt_post_solve_hook,
    transient = TSCOPF.reconfigure_transient(cfg.transient;
        simulation = TsSimulationConfig(; (f => f == :t_step ? TSTEP : f == :clearing_time ? CLEAR :
                                           getfield(cfg.transient.simulation, f)
                                           for f in fieldnames(TsSimulationConfig))...),
        dyn_model = DynModelConfig(; (f => f == :fault ? fault : getfield(cfg.transient.dyn_model, f)
                                      for f in fieldnames(DynModelConfig))...)))

validate_run_config!(run_cfg)
validate_dyn_config!(run_cfg)

println("Run : ", ARGS[1], "    case : ", CASE, "    t_step : ", TSTEP, "    clearing_time : ", CLEAR)
sys    = load_system(run_cfg, ROOT)
result = run_case!(run_cfg, sys, ROOT, PATH_RESULTS)
println("status : ", result.status, "    objective : ", result.obj_MVA)
println("folder : ", result.path_names[:pf_results_date])
