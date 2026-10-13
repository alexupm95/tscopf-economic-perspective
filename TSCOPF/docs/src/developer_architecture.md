# Developer architecture

This page maps the source layout for maintainers. User-facing run options are in the [Parameter reference](parameter_reference.md) and [Configuration map](configuration_map.md).

## Transient-stability builders

The transient-stability layer keeps orchestration separate from low-level JuMP builders:

| File family | Responsibility |
|-------------|----------------|
| `functions_2_build_TS_model_common.jl` | Routes `DynModelConfig` to the selected dynamic model and export path. |
| `functions_2_build_TS_model_w_Kron.jl` | Orchestrates nonlinear Kron-reduced TSC-ACOPF blocks. |
| `functions_2_build_TS_model_w_Kron_Linear.jl` | Orchestrates Taylor-linearized Kron TSC-DCOPF blocks. |
| `functions_2_build_TS_model_w_FullBus.jl` | Orchestrates full-bus TSC-ACOPF blocks and disturbance handling. |
| `functions_4_TS_kron_variables.jl` | Kron-reduced dynamic variables (`E`, `δ`, `P_m`, time-indexed generator/COI variables). |
| `functions_4_TS_kron_eqconst.jl` | Kron equality constraints: initial conditions, COI equations, swing equations, nonlinear/Taylor `P_e`. |
| `functions_4_TS_kron_ineqconst.jl` | Kron stability and COI-relative inequality constraints. |
| `functions_4_TS_fullbus_variables.jl` | Full-bus time-indexed bus and generator variables. |
| `functions_4_TS_fullbus_eqconst.jl` | Full-bus KCL, `P_e`/`Q_e`, and swing equality constraints. |
| `functions_4_TS_fullbus_ineqconst.jl` | Full-bus voltage lower bounds and stability inequality dispatch. |
| `functions_4_TS_fullbus_helpers.jl` | Full-bus warm starts and sparse nodal-injection expressions. |
| `functions_2_build_TS_model_w_Kron_shared.jl` | Shared helper math used across Kron and FullBus paths. |

The `functions_2_*` files should remain orchestration-focused: choose fault windows, prepare admittance matrices, call variable builders, call equality/inequality builders, and store results in `dyn_model_dict`.

The `functions_4_TS_*` files should own individual JuMP variable or constraint families. This mirrors the steady-state `_common/functions_4_variables.jl`, `_common/functions_4_eqconst.jl`, and `_common/functions_4_ineqconst.jl` layout.

**Variable box limits.** Simple bounds can be encoded either as explicit `@constraint` inequalities (`BoundEncoding.CONSTRAINT`, default) or as JuMP variable bounds at creation (`BoundEncoding.VARIABLE`). Shared helpers live in `_common/bound_encoding.jl` and `_common/bound_builders.jl`; dual export uses `opf_dict[:meta][:bound_manifest]` / `dyn_model_dict[:meta][:bound_manifest]` when the VARIABLE path is selected so registry CSV names stay unchanged. `Export_OPF_Model` (`functions_2_save_dispatch_model.jl`) reads the same manifest for `model_details.txt`, writing equivalent ≤-form bound lines so VARIABLE runs document the full primal.

## Optional extensions

Core package load should stay open-solver and lightweight: Ipopt + HiGHS only. Optional packages register runtime callbacks through Julia package extensions:

| Extension | Trigger | Registration |
|-----------|---------|--------------|
| `TSCOPFGurobiExt` | `using Gurobi` | Gurobi optimizer builder and license probe. |
| `TSCOPFMadNLPExt` | `using MadNLP` | MadNLP optimizer builder (also `"MadNLP-ma57"` / `"MadNLP-ma97"`; loads `MadNLPHSL` on demand). |
| `TSCOPFPlotsExt` | `using Plots` | Transient-stability plot callbacks. |

Extension registration happens in each extension module's `__init__()` so precompilation does not lose registry side effects.

## Solver setup and log files

| Function | File | Role |
|----------|------|------|
| `Setup_Optim_Model` | `functions_2_setup_optim.jl` | Factory via `solver_registry.jl` (Ipopt, HiGHS; Gurobi/MadNLP via extensions). |
| `IpoptSolverConfig` | `IpoptSolverConfig.jl` | Ipopt tolerances, iteration cap, Hessian mode; PARDISO `pardiso_lib_path` when `solver_name = "Ipopt-pardiso"`. |
| `HiGHSSolverConfig` | `HiGHSSolverConfig.jl` | HiGHS LP/QP tolerances, iteration limits, solver mode, `raw_options`. |
| `GurobiSolverConfig` | `GurobiSolverConfig.jl` | Gurobi LP/MIP/NL barrier knobs + `raw_options`. |
| `MadNLPSolverConfig` | `MadNLPSolverConfig.jl` | MadNLP tolerances, Hessian mode, optional `ma97_num_threads`, `raw_options`. Linear solver via `solver_name` + `_common/madnlp_linear_solvers.jl`. |
| `set_solver_log_path!` | `functions_2_setup_optim.jl` | Redirect solver text output to an absolute path (`output_file` / `LogFile` / `log_file`). |
| `release_solver_backend!` | same | Tear down after exports: optional Ipopt `finalize`, `empty!(model)`, `GC.gc()`. Called at end of `run_case!` and explicit dual LP builders. |

`engine.jl` writes `solver_log.txt` for the main solve and `solver_log_warmstart.txt` when a DCOPF or FULL_BUS ACOPF pre-solve runs. See [Running a case — solver logs](running_a_case.md).

!!! tip "Adding a new optional backend"
    Follow the same shape: a `TSCOPFxExt` module gated on `[extensions]` /
    `[weakdeps]` in `Project.toml`, a registry function it populates from its own
    `__init__()`, and a `_ext_loaded()`-style query so core code can check
    availability (see `gurobi_available`, `plots_extension_loaded`) instead of
    `@isdefined`-guessing.

## Include order

`src/TSCOPF.jl` includes configuration and registry files first, then low-level transient-stability builders, then orchestration files. When adding new dynamic model variants, add low-level builder files before any `functions_2_build_*` file that calls them.

!!! warning "Never mutate the loaded network DataFrames"
    `SystemData.DBUS` / `DGEN` / `DGEN_DYN` / `DCIR` are read once by
    [`load_system`](api.md) and reused across a sweep. GL (generator/load trip)
    disturbances must `deepcopy` before scaling demand or flipping `g_status` —
    see `apply_gl_load_scaling!` in `FaultConfig.jl`. Mutating the originals in
    place silently corrupts every later scenario in the same sweep.
