# Numerical / physical constants live in `_common/constants.jl`.

# ==============================================================================
#  Project paths (A3 — no cwd assumptions)
# ==============================================================================

"""
    project_root() -> String

Repository / package root (`pkgdir(TSCOPF)` in a dev checkout with `Project.toml` at top level).
Use for default `path_main` when the entry script does not rely on `pwd()`.
"""
function project_root()::String
    root = pkgdir(@__MODULE__)
    root === nothing &&
        throw(ArgumentError("project_root(): TSCOPF is not loaded as an installed package; activate the project environment first."))
    return root
end

"""Default results tree: `<project_root>/RESULTS`."""
function default_results_dir()::String
    return joinpath(project_root(), "RESULTS")
end

"""
Bundled demo / case-study folder: `<project_root>/INPUT_FILES/<case>`.

Convenience for user scripts. It is **not** how the package locates inputs — a run reads
`<path_main>/INPUT_FILES/<case>`, where `path_main` is [`load_system`](@ref)'s second
argument, so case data rooted anywhere else is reached by passing that root. The test
suite does exactly this to read its frozen fixtures from `test/INPUT_FILES/`.
"""
function input_files_dir(case::String)::String
    return joinpath(project_root(), "INPUT_FILES", case)
end

# Function to clean the terminal
function Clean_Terminal()

    # If system == Windows
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # If system is based on Unix
    else
        Base.run(`clear`)
    end
    
end

# ==============================================================================
#  Results folder layout (simulation-dependent)
# ==============================================================================

"""
    _timestamped_results_dir(path_results) -> String

Timestamped run directory (`:` stripped for Windows paths), **claimed atomically**.

The name has second resolution. Two runs that finish inside the same second — routine
for ED/DC-OPF sweeps and for the test suite — would otherwise resolve to the same
folder and silently interleave their outputs: `Inputs/`, `Dispatch/CSV/` and the dual
exports become a union of both runs, with same-named files overwritten by whichever
finished last.

The directory is therefore claimed with `mkdir`, which **fails** if the path already
exists, and a ` (2)`, ` (3)`, … suffix is appended until the claim succeeds. Because
`mkdir` is atomic on both Windows and POSIX, this is also correct for concurrent
sweeps running in separate processes, not just sequential runs in one session.

The common case is unchanged: a run that gets its own second keeps the plain
`Results - yyyy-mm-dd HHMMSS` name, so existing frozen baselines and the paths quoted
in the docs stay valid. Suffixes are assigned in creation order, so chronological
sorting still holds within a second.
"""
function _timestamped_results_dir(path_results::String)
    dt = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    base = joinpath(path_results, replace("Results - " * dt, ":" => ""))
    mkpath(path_results)

    candidate = base
    n = 1
    while true
        try
            mkdir(candidate)
            return candidate
        catch err
            # Only "already exists" is retryable — a permission or disk error must
            # surface here rather than spin through 999 suffixes.
            (err isa Base.IOError && ispath(candidate)) || rethrow()
            n += 1
            n > 999 && error("Could not claim a unique results directory near $base " *
                             "(999 collisions). Check for a runaway loop writing to $path_results.")
            candidate = string(base, " (", n, ")")
        end
    end
end

"""All canonical path keys for a run; directories are not created here."""
function build_results_path_names(path_main::String,
                                  path_results::String,
                                  path_results_date::String)
    path_folder_dispatch                = joinpath(path_results_date, "Dispatch")
    path_folder_dispatch_CSV            = joinpath(path_folder_dispatch, "CSV")
    path_folder_dispatch_CSV_duals      = joinpath(path_folder_dispatch, "CSV_duals")
    path_folder_dispatch_dual           = joinpath(path_results_date, "Dispatch_Dual")
    path_folder_dispatch_dual_CSV       = joinpath(path_folder_dispatch_dual, "CSV")
    path_folder_dispatch_dual_CSV_duals = joinpath(path_folder_dispatch_dual, "CSV_duals")
    path_folder_dispatch_warmstart            = joinpath(path_results_date, "Dispatch_WarmStart")
    path_folder_dispatch_warmstart_CSV        = joinpath(path_folder_dispatch_warmstart, "CSV")
    path_folder_dispatch_warmstart_CSV_duals    = joinpath(path_folder_dispatch_warmstart, "CSV_duals")
    path_folder_bus_matrices            = joinpath(path_results_date, "Bus_Matrices")
    path_folder_TS                      = joinpath(path_results_date, "Transient_Stability")
    path_folder_TS_figures              = joinpath(path_folder_TS, "Figures")
    path_folder_TS_figures_duals        = joinpath(path_folder_TS, "Figures_Duals")
    path_folder_TS_CSV                  = joinpath(path_folder_TS, "CSV")
    path_folder_TS_CSV_duals            = joinpath(path_folder_TS, "CSV_duals")
    path_folder_inputs                  = joinpath(path_results_date, "Inputs")

    return OrderedDict(
        :pf_main                    => path_main,
        :pf_results                 => path_results,
        :pf_results_date           => path_results_date,
        :pf_inputs                 => path_folder_inputs,
        :pf_dispatch               => path_folder_dispatch,
        :pf_dispatch_CSV           => path_folder_dispatch_CSV,
        :pf_dispatch_CSV_duals     => path_folder_dispatch_CSV_duals,
        :pf_dispatch_dual           => path_folder_dispatch_dual,
        :pf_dispatch_dual_CSV       => path_folder_dispatch_dual_CSV,
        :pf_dispatch_dual_CSV_duals => path_folder_dispatch_dual_CSV_duals,
        :pf_dispatch_warmstart       => path_folder_dispatch_warmstart,
        :pf_dispatch_warmstart_CSV   => path_folder_dispatch_warmstart_CSV,
        :pf_dispatch_warmstart_CSV_duals => path_folder_dispatch_warmstart_CSV_duals,
        :pf_bus_matrices            => path_folder_bus_matrices,
        :pf_TS                      => path_folder_TS,
        :pf_TS_figures              => path_folder_TS_figures,
        :pf_TS_figures_duals        => path_folder_TS_figures_duals,
        :pf_TS_CSV                  => path_folder_TS_CSV,
        :pf_TS_CSV_duals            => path_folder_TS_CSV_duals,
    )
end

"""
    _has_usable_duals(model) -> Bool

`true` only when the solver returned a dual point it vouches for.

Deliberately **not** `JuMP.has_duals`, which is defined as `dual_status != NO_SOLUTION`
and therefore accepts `UNKNOWN_RESULT_STATUS` — exactly what Ipopt reports after an
`ITERATION_LIMIT`, where multipliers exist as iterates but certify nothing. Exporting
those into `CSV_duals/` would put numbers that look like prices next to numbers that are.
"""
function _has_usable_duals(model)::Bool
    return JuMP.dual_status(model) in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
end

"""
    _write_solve_status_header!(io, model, label)

Stamp the solver's termination / primal / dual status at the top of a dual report.

Dual values are only prices when the solve actually converged, and `run_case!`
deliberately accepts `ITERATION_LIMIT` and `ALMOST_LOCALLY_SOLVED` as "export the
primal anyway". Recording the status in the file itself means a CSV or TXT read
back months later cannot be mistaken for a clean optimum.
"""
function _write_solve_status_header!(io::IO, model, label::AbstractString)
    println(io, "======================================")
    println(io, "  $label")
    println(io, "======================================")
    println(io, "termination_status : ", JuMP.termination_status(model))
    println(io, "primal_status      : ", JuMP.primal_status(model))
    println(io, "dual_status        : ", JuMP.dual_status(model))
    println(io, "usable_duals       : ", _has_usable_duals(model))
    println(io)
    return nothing
end

"""
    results_folder_keys(cfg) -> Vector{Symbol}

Leaf directories to `mkpath` before a simulation run, based on what will actually
be written (steady-state dispatch, explicit dual, matrices, transient layer).
`cfg` is a `RunConfig` (defined in `engine.jl`; untyped here for include order).
"""
function results_folder_keys(cfg)
    keys_sym = [:pf_inputs, :pf_dispatch_CSV]
    cfg.save_duals && push!(keys_sym, :pf_dispatch_CSV_duals)
    if cfg.dispatch.solve_explicit_dual && cfg.dispatch.type_model in ("DCOPF", "ED")
        append!(keys_sym, [:pf_dispatch_dual_CSV, :pf_dispatch_dual_CSV_duals])
    end
    type_model = cfg.dispatch.type_model
    # Only ACOPF/DCOPF write Ybus/Bbus; ED and UC have no admittance dump.
    if type_model in ("ACOPF", "DCOPF") && (cfg.save_matrices || cfg.trans_stab)
        push!(keys_sym, :pf_bus_matrices)
    end
    if cfg.trans_stab
        push!(keys_sym, :pf_TS_CSV)
        cfg.save_duals && push!(keys_sym, :pf_TS_CSV_duals)
        # Figures_Duals/ holds the δ-COI dual SVGs, which are figures like any other:
        # they need `save_ts_plots` as well as `save_duals`. Creating the folder on
        # `save_duals` alone left an empty directory on every default run.
        cfg.save_ts_plots && push!(keys_sym, :pf_TS_figures)
        cfg.save_ts_plots && cfg.save_duals && push!(keys_sym, :pf_TS_figures_duals)
    end
    if warmstart_dispatch_applicable(cfg)
        push!(keys_sym, :pf_dispatch_warmstart_CSV)
        cfg.save_duals && push!(keys_sym, :pf_dispatch_warmstart_CSV_duals)
    end
    return keys_sym
end

"""`true` when FULL_BUS TSC-ACOPF coupling options apply to this run."""
function fullbus_tsc_acopf_applicable(cfg)::Bool
    cfg.trans_stab || return false
    cfg.dispatch.type_model == "ACOPF" || return false
    cfg.transient === nothing && return false
    return cfg.transient.dyn_model.network_form == FULL_BUS
end

"""
`true` when `save_warmstart_dispatch` applies to this run configuration.

Two TSC paths pre-solve a steady-state OPF before assembling the dynamics, and both
have something worth archiving: FULL_BUS ACOPF (the warm-start operating point that
seeds every coupling variable) and Kron TSC-DCOPF (the DC dispatch that fixes the
Taylor anchor `δ_ref`). Kron TSC-ACOPF solves once, jointly, so there is no pre-solve
to save.
"""
function warmstart_dispatch_applicable(cfg)::Bool
    cfg.save_warmstart_dispatch || return false
    fullbus_tsc_acopf_applicable(cfg) && return true
    return dcopf_tsc_presolve_applicable(cfg)
end

"""`true` when this run is a TSC-DCOPF whose DC pre-solve produces `δ_ref`."""
function dcopf_tsc_presolve_applicable(cfg)::Bool
    cfg.trans_stab || return false
    cfg.dispatch.type_model == "DCOPF" || return false
    cfg.transient === nothing && return false
    return cfg.transient.dyn_model.network_form == KRON_REDUCED
end

"""Remap dispatch path keys to `Dispatch_WarmStart/` (mirrors `dispatch_dual_path_names`)."""
function dispatch_warmstart_path_names(path_names::OrderedDict{Symbol, String})
    return OrderedDict{Symbol, String}(
        :pf_main => path_names[:pf_main],
        :pf_results_date => path_names[:pf_results_date],
        :pf_dispatch => path_names[:pf_dispatch_warmstart],
        :pf_dispatch_CSV => path_names[:pf_dispatch_warmstart_CSV],
        :pf_dispatch_CSV_duals => path_names[:pf_dispatch_warmstart_CSV_duals],
    )
end

function _mkpath_results_folders!(path_names::OrderedDict{Symbol, String},
                                  folder_keys::AbstractVector{Symbol})
    for k in folder_keys
        mkpath(path_names[k])
    end
    return path_names
end

"""
    build_results_paths(path_main, path_folder_results, cfg; overwrite_results=false)

Create the timestamped (or flat overwrite) results tree for `run_case!`.
Only subfolders required by `cfg` are created on disk; all path keys are
always returned so callers can reference optional locations safely.
"""
function build_results_paths(path_main::String,
                             path_folder_results::String,
                             cfg;
                             overwrite_results::Bool=cfg.overwrite_results)
    path_results_date = overwrite_results ?
        path_folder_results :
        _timestamped_results_dir(path_folder_results)
    path_names = build_results_path_names(path_main, path_folder_results, path_results_date)
    _mkpath_results_folders!(path_names, results_folder_keys(cfg))
    return path_names
end

"""
    build_import_results_paths(path_main, path_folder_results)

Timestamped run folder with `Inputs/` only — used by MATPOWER import and any
workflow that archives source CSVs before `run_case!` copies them again.
"""
function build_import_results_paths(path_main::String, path_folder_results::String)
    path_results_date = _timestamped_results_dir(path_folder_results)
    path_names = build_results_path_names(path_main, path_folder_results, path_results_date)
    _mkpath_results_folders!(path_names, [:pf_inputs])
    return path_names
end

"""Path keys for the explicit DC-OPF dual formulation (mirrors primal Dispatch layout)."""
function dispatch_dual_path_names(path_names::OrderedDict{Symbol, String})
    return OrderedDict{Symbol, String}(
        :pf_main                   => path_names[:pf_main],
        :pf_results_date           => path_names[:pf_results_date],
        :pf_dispatch               => path_names[:pf_dispatch_dual],
        :pf_dispatch_CSV           => path_names[:pf_dispatch_dual_CSV],
        :pf_dispatch_CSV_duals     => path_names[:pf_dispatch_dual_CSV_duals],
    )
end
