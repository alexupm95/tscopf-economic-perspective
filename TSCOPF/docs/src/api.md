# API reference

Exported types and functions from the `TSCOPF` module. This page is generated from
docstrings in `src/`; the narrative field-by-field reference lives in the
[Parameter reference](parameter_reference.md), [Running a case](running_a_case.md), and [Configuration map](configuration_map.md)
— read those first, come back here for exact signatures.

!!! note "Partial by design"
    Not every exported symbol has a docstring yet (`checkdocs = :none` in `docs/make.jl`
    reflects that). Internal helpers and one-off UC-toy functions are deliberately left
    off this page; if you need one, it is documented inline in `src/`.

## Orchestration

The entry-point trio: build a [`RunConfig`](@ref), load the network once with
[`load_system`](@ref), then call [`run_case!`](@ref) (repeatedly, for a sweep).

```@docs
RunConfig
reconfigure
reconfigure_transient
SystemData
load_system
run_case!
```

## Validation

```@docs
validate_run_config!
validate_dyn_config!
validate_dispatch_config!
```

## Steady-state dispatch (avenue 1)

```@docs
DispatchConfig
BoundEncoding
SusceptanceModel
build_opf_input_param
```

## Transient stability (avenue 2)

```@docs
TransientConfig
TsSimulationConfig
TsBuilderConfig
DynModelConfig
GenOrder
NetworkForm
MechPowerMode
SteadyStateHints
```

## Disturbances

```@docs
FaultConfig
FaultType
```

## Paths and I/O

```@docs
project_root
default_results_dir
input_files_dir
build_results_paths
build_results_path_names
results_folder_keys
Import_Matpower_Case
```

## Solver helpers

```@docs
IpoptSolverConfig
apply_ipopt_options!
validate_ipopt_solver_config!
HiGHSSolverConfig
apply_highs_options!
validate_highs_solver_config!
GurobiSolverConfig
apply_gurobi_options!
validate_gurobi_solver_config!
MadNLPSolverConfig
apply_madnlp_options!
validate_madnlp_solver_config!
Setup_Optim_Model
set_solver_log_path!
release_solver_backend!
gurobi_available
is_ipopt_backend_solver
hsl_jll_available
IPOPT_PARDISO_SOLVER
pardiso_available
configure_ipopt_linear_solver!
is_madnlp_backend_solver
MADNLP_BACKEND_SOLVERS
MADNLP_HSL_LINEAR_SOLVERS
madnlp_hsl_ext_available
configure_madnlp_linear_solver!
```
