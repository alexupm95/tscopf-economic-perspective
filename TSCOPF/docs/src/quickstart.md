# Quick start

All commands run from the root of the repository, after [Install](install.md).

## The runs of the paper

| Case folder (`INPUT_FILES/`) | Command |
|---|---|
| `9bus_linear_H24_6_3` | `julia --project=TSCOPF runs/run_9bus.jl linear` |
| `9bus_colocated_bus1_H24_6_3` | `julia --project=TSCOPF runs/run_9bus.jl colocated` |
| `39bus_gfm` | `julia --project=TSCOPF runs/run_39bus.jl` |

Each solve writes a `Results - <timestamp>/` folder under `new_runs/` (not tracked). The
folders of the runs printed in the paper are under `results/`.

`runs/_9bus_config.jl` is a complete, working configuration: it sets every field of
`DispatchConfig`, `TsSimulationConfig`, `TsBuilderConfig`, `DynModelConfig` and the Ipopt
options for the 9-bus runs. `runs/run_39bus.jl` does the same for a mixed synchronous and
grid-forming fleet, and reads its settings from the environment variables listed in its header.
Start from one of them for a new study.

## Programmatic use

An ACOPF on the 9-bus case of the paper:

```julia
using TSCOPF

root = pwd()                                   # root of the repository
cfg = RunConfig(
    case = "9bus_linear_H24_6_3",
    load_factor = 1.0,
    dispatch = DispatchConfig(type_model = "ACOPF", cost_type = "linear", ineq_sbranch_upper = false),
    silent_solver = true,
    save_optim_matrices = false,
)
sys = load_system(cfg, root)
result = run_case!(cfg, sys, root, mkpath(joinpath(root, "new_runs")))

result.status          # MOI termination status
result.obj_MVA         # objective
result.RGEN.p_g        # generator dispatch [MW]
```

`load_system` reads `INPUT_FILES/<case>/` under its second argument.

!!! note "Check `result.status` before trusting duals"
    `run_case!` only writes output files on `OPTIMAL`, `LOCALLY_SOLVED`, or
    `ITERATION_LIMIT`. For nonconvex ACOPF/TSC-ACOPF, `LOCALLY_SOLVED` is a valid
    KKT point but not a certified global optimum — near-flat cost curvature or a
    poor warm start can land Ipopt in a different local solution, which changes
    the dual values even when the primal dispatch looks reasonable.

## Transient stability constraints

`trans_stab = true` adds the transient block, configured through `transient = TransientConfig(...)`:
the time grid and the angle tolerance in `TsSimulationConfig`, the machine and network models
and the fault in `DynModelConfig`. The run scripts above show the full set. See
[Running a case](running_a_case.md) for more example configurations and
[Parameter reference](parameter_reference.md) for `DynModelConfig`, fault types (`SC`, `GL`) and
builder toggles.

## Shadow prices

The run scripts attach the export of `kkt/kkt_pg_terms.jl` through `RunConfig.post_solve_hook`;
`runs/run_39bus.jl` adds `kkt/kkt_sensitivity.jl` with `J_KKT=2`. The notebooks in `kkt/` and
the scripts in `postprocess/` turn those exports into the prices and settlements of the paper;
`kkt/README.md` describes them.
