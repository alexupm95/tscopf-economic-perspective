#=
================================================================================
 functions_2_save_run_manifest.jl  —  machine-readable description of one run
================================================================================
 `run_manifest.toml` sits at the root of every results folder and answers, without
 any Julia in the loop, the question "what exactly produced these CSVs?".

 Why it exists. Three sources of run provenance already ship, and none of them is
 machine-readable:

   * `RunConfig.run_script` is a byte-for-byte copy of the script the user pointed
     at, so it contains only the fields that were *typed*. `main.jl` sets ~25 of
     the 150+ fields in the resolved config tree, and the defaults are not uniform
     (some `TsBuilderConfig` toggles default `true`, most default `false`). Reading
     it back means parsing Julia and mirroring every default.
   * `input_parameters.txt` is prose for humans, and it omits `load_factor`.
   * `Inputs/bus_data.csv` archives the *unscaled* demand — `load_system` applies
     `load_factor` after that copy is taken — so the demand actually solved cannot
     be recovered from a results folder at all without this file.

 The resolved facts that only exist inside the model (which generator the δ
 corridor was measured against, the SG/GFM partition, the time windows, which dual
 CSVs were written) are recorded in `[resolved]` and `[exports]`.

 Everything is written through `_toml_value`, so a config field added upstream
 appears here automatically: the `[run]`, `[dispatch]` and `[transient.*]` tables
 are reflected off the structs with `fieldnames`, never enumerated by hand.
================================================================================
=#

# --- key and value normalisation ------------------------------------------------

# Greek field names (`δ_tol`, `Δω_tol`, `bound_style_δ`) are idiomatic in this
# codebase but hostile to attribute access on the reading side, so they are
# transliterated on the way out. `Δω` is listed before `Δ` and `ω`: `replace` takes
# the first pattern that matches at a position, so the pair maps to a single
# readable `Delta_omega` instead of `Deltaomega`.
const _MANIFEST_KEY_TRANSLITERATION = [
    "Δω" => "Delta_omega",
    "δ" => "delta", "Δ" => "Delta", "ω" => "omega", "θ" => "theta",
    "λ" => "lambda", "φ" => "phi", "ν" => "nu", "μ" => "mu",
    "Γ" => "Gamma", "γ" => "gamma", "σ" => "sigma", "τ" => "tau",
    "ε" => "epsilon", "ρ" => "rho", "α" => "alpha", "β" => "beta",
    "η" => "eta", "ψ" => "psi", "Ω" => "Omega", "π" => "pi",
]

"""ASCII TOML key for a (possibly Greek) field name."""
function _manifest_key(name)::String
    return replace(String(name), _MANIFEST_KEY_TRANSLITERATION...)
end

"""
    _toml_value(x)

Coerce `x` into something `TOML.print` accepts, or `nothing` when the field has no
place in the manifest (the caller then omits the key entirely).

`nothing` is returned for `Nothing` and for `Function`: an unset optional field and
a callback such as `post_solve_hook` both have no serialisable content, and TOML has
no null literal to record them with. Non-finite floats become strings because TOML
has no `inf`/`nan` — a reader that expects a number there would rather see `"Inf"`
than silently read `0.0`. A `Complex` (the SC fault admittance) is `[re, im]`.
"""
function _toml_value(x)
    x === nothing && return nothing
    x isa Function && return nothing
    x isa Bool && return x
    x isa Integer && return Int(x)
    if x isa AbstractFloat
        isfinite(x) || return string(x)
        return Float64(x)
    end
    x isa AbstractString && return String(x)
    x isa Enum && return string(x)
    x isa Symbol && return String(x)
    # TOML has no complex type. A fault admittance is the one config field that
    # is complex; write it as [re, im] so a reader does not have to parse a string.
    if x isa Complex
        (isfinite(real(x)) && isfinite(imag(x))) || return string(x)
        return Any[Float64(real(x)), Float64(imag(x))]
    end
    if x isa Tuple || x isa AbstractVector
        vals = Any[_toml_value(v) for v in x]
        any(isnothing, vals) && return string(x)   # heterogeneous / unserialisable
        return vals
    end
    if x isa AbstractDict
        # Solver `raw_options` and friends: keep them as a real subtable rather than
        # stringifying the whole Dict, so a reader can see the options that were set.
        table = OrderedDict{String, Any}()
        for (k, v) in x
            converted = _toml_value(v)
            converted === nothing || (table[string(k)] = converted)
        end
        return table
    end
    return string(x)
end

"""
`true` for a nested configuration struct that should become its own TOML subtable.

Restricted to types defined in this module so that a field holding, say, a
`DataFrame` or a solver handle is stringified by `_toml_value` instead of being
walked field by field.
"""
function _is_config_struct(x)::Bool
    T = typeof(x)
    x isa Enum && return false
    return isstructtype(T) && parentmodule(T) === @__MODULE__
end

"""Reflect a configuration struct into an `OrderedDict` ready for `TOML.print`."""
function _config_table(obj)::OrderedDict{String, Any}
    table = OrderedDict{String, Any}()
    for field in fieldnames(typeof(obj))
        value = getfield(obj, field)
        key = _manifest_key(field)
        if _is_config_struct(value)
            table[key] = _config_table(value)
        else
            converted = _toml_value(value)
            converted === nothing || (table[key] = converted)
        end
    end
    return table
end

# --- resolved facts (not derivable from the config alone) -----------------------

"""Row of `DGEN_DYN` holding generator `gen`, or `nothing` (GFM units have none)."""
function _dyn_row_index(DGEN_DYN::DataFrame, gen::Int)
    if :id in propertynames(DGEN_DYN)
        return findfirst(==(gen), Int.(DGEN_DYN.id))
    end
    return (1 <= gen <= nrow(DGEN_DYN)) ? gen : nothing
end

"""Collect `column` of `DGEN_DYN` over `gens`, skipping units with no dynamic row."""
function _dyn_column(DGEN_DYN::DataFrame, gens::Vector{Int}, column::Symbol)
    column in propertynames(DGEN_DYN) || return Float64[]
    out = Float64[]
    for gen in gens
        row = _dyn_row_index(DGEN_DYN, gen)
        row === nothing || push!(out, Float64(DGEN_DYN[row, column]))
    end
    return out
end

function _resolved_table(cfg, sys, dyn_model_dict, dyn_parameters_dict)
    resolved = OrderedDict{String, Any}()

    # --- steady-state fleet and cost curve ---------------------------------------
    DGEN = sys.DGEN
    active = findall(==(1), DGEN.g_status)
    gen_ids = Int.(DGEN.id[active])
    resolved["gen_ids"] = gen_ids
    resolved["gen_bus"] = Int.(DGEN.bus[active])
    resolved["cost_type"] = String(cfg.dispatch.cost_type)
    resolved["c0"] = Float64.(DGEN.g_cost_0[active])
    resolved["c1"] = Float64.(DGEN.g_cost_1[active])
    resolved["c2"] = Float64.(DGEN.g_cost_2[active])

    dyn_model_dict === nothing && return resolved

    # --- dynamic fleet: SG / GFM partition ---------------------------------------
    # `meta[:sg_gens]` only exists on the mixed-fleet DQ path; everywhere else every
    # active unit is a machine (same fallback the debug exporter uses).
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    active_gen = Vector{Int}(get(dyn_model_dict, :active_gen, gen_ids))
    gfm_ids = Vector{Int}(get(meta, :gfm_gens, Int[]))
    sg_ids = Vector{Int}(get(meta, :sg_gens, Int[]))
    isempty(sg_ids) && (sg_ids = Int[g for g in active_gen if g ∉ gfm_ids])
    resolved["active_gen"] = active_gen
    resolved["sg_ids"] = sg_ids
    resolved["gfm_ids"] = gfm_ids

    # Converter equivalent inertia, when `gfm_in_coi` put the converters in the COI.
    # `H_eff = TPf/(2·mp)` is derived at build time from the system-base `DGFM` row, so it
    # appears in no CSV and in no config table — this is the only record of the weight the
    # COI mean actually used. Ordered like `gfm_ids`.
    if get(meta, :gfm_in_coi, false)
        params = get(meta, :gfm_params, nothing)
        resolved["gfm_in_coi"] = true
        if params !== nothing && !isempty(gfm_ids)
            weights = Float64[]
            for g in gfm_ids
                haskey(params, g) || continue
                p = params[g]
                p.mp > 0.0 || continue
                push!(weights, Float64(p.TPf) / (2.0 * Float64(p.mp)))
            end
            length(weights) == length(gfm_ids) && (resolved["gfm_H_eff"] = weights)
        end
    end

    # Overload-clip smoothing ε, per axis per converter, for the state clip and the output clip. Derived at build time from
    # Δt and the machine-base gains (`_gfm_pq_eps`, `_gfm_pq_eps_out`), so, like `gfm_H_eff`, no
    # CSV or config table carries it. Ordered like `gfm_ids`; only present when
    # `gfm_pq_overload` built the branches.
    pq_eps = get(meta, :gfm_pq_eps, nothing)
    if pq_eps !== nothing && all(g -> haskey(pq_eps, g), gfm_ids) && !isempty(gfm_ids)
        resolved["gfm_pq_eps_P"] = Float64[pq_eps[g].P for g in gfm_ids]
        resolved["gfm_pq_eps_Q"] = Float64[pq_eps[g].Q for g in gfm_ids]
        resolved["gfm_pq_eps_out_P"] = Float64[pq_eps[g].P_out for g in gfm_ids]
        resolved["gfm_pq_eps_out_Q"] = Float64[pq_eps[g].Q_out for g in gfm_ids]
        resolved["gfm_pq_eps_step_frac"] = _GFM_PQLIM_EPS_STEP_FRAC
        resolved["gfm_pq_out_err_tol"] = _GFM_PQLIM_OUT_ERR_TOL
    end

    # The δ corridor reference is resolved inside the builder (`:highest_H` ranks by
    # inertia at build time), so the config alone does not name it.
    ref_gen = get(meta, :δ_ref_gen_resolved, nothing)
    ref_gen === nothing || (resolved["delta_ref_gen"] = Int(ref_gen))

    # `bound_style_δ = :pairwise` resolves no reference; what a reader needs instead is which
    # unordered pairs the corridor actually covered, which depends on the post-disturbance
    # survivors and so is not derivable from the config either.
    let ineq = get(dyn_model_dict, :ineq_const, nothing)
        if ineq isa AbstractDict
            for key in (:ineq_const_δ_pair_tf_upper, :ineq_const_δ_pair_tf_lower,
                        :ineq_const_δ_pair_tpf_upper, :ineq_const_δ_pair_tpf_lower)
                haskey(ineq, key) || continue
                pairs = collect(keys(ineq[key]))
                resolved["delta_pair_count"] = length(pairs)
                # TOML has no tuple; a flat "i-j" string keeps the pairing readable.
                resolved["delta_pairs"] = String["$(p[1])-$(p[2])" for p in pairs]
                break
            end
        end
    end

    # Machine constants, SG only: GFM units have no row in DGEN_DYN and no H.
    if sys.DGEN_DYN !== nothing
        for (column, key) in ((:H, "H"), (:D, "D"), (:Xd_tr, "Xd_tr"))
            values = _dyn_column(sys.DGEN_DYN, sg_ids, column)
            isempty(values) || (resolved[key] = values)
        end
    end

    dyn_parameters_dict === nothing && return resolved

    # --- tolerances and the simulation timeline ----------------------------------
    common = get(dyn_parameters_dict, :common, OrderedDict{Symbol, Any}())
    for (field, key) in ((:δ_tol, "delta_tol"), (:δ_pair_tol, "delta_pair_tol"),
                         (:Δω_tol, "Delta_omega_tol"),
                         (:ω_syn, "omega_syn"), (:f_syn, "f_syn"), (:Δω_0, "Delta_omega_0"))
        haskey(common, field) || continue
        value = _toml_value(common[field])
        value === nothing || (resolved[key] = value)
    end

    time = get(dyn_parameters_dict, :time, OrderedDict{Symbol, Any}())
    for (field, key) in ((:t_step, "t_step"), (:t_start_sim, "t_start_sim"),
                         (:t_end_sim, "t_end_sim"), (:t_start_fault, "t_start_fault"),
                         (:t_clear_fault, "t_clear_fault"), (:clearing_time, "clearing_time"))
        haskey(time, field) || continue
        value = _toml_value(time[field])
        value === nothing || (resolved[key] = value)
    end
    # Window *lengths*, not the windows themselves: the dual CSVs are row-aligned
    # with the `t` column of Transient_Stability/CSV/electrical_power.csv, and these
    # counts are what a reader needs to split that axis into fault / post-fault.
    for (field, key) in ((:t_window_fault, "n_steps_fault"),
                         (:t_window_postf, "n_steps_postfault"),
                         (:t_window_total, "n_steps_total"))
        haskey(time, field) && (resolved[key] = length(time[field]))
    end

    return resolved
end

# --- which dual files this run actually wrote ------------------------------------

function _exports_table(opf_dict, dyn_model_dict)
    exports = OrderedDict{String, Any}()

    if opf_dict !== nothing
        files = String[spec.csv_file for spec in STEADY_STATE_DUAL_SPECS
                       if spec.csv_file !== nothing && spec_present(opf_dict, spec)]
        exports["dispatch_dual_files"] = unique(files)
    end

    if dyn_model_dict !== nothing
        # `:dual_registry` is already filtered to the families this run built
        # (`build_dual_registry!`), so no second presence check is needed.
        registry = get(dyn_model_dict, :dual_registry, DualRegistryEntry[])
        isempty(registry) && (registry = dual_spec_catalog(dyn_model_dict))
        files = String[entry.csv_file for entry in registry if entry.csv_file !== nothing]
        exports["ts_dual_files"] = unique(files)
    end

    return exports
end

# --- entry point -----------------------------------------------------------------

"""
    Save_Run_Manifest!(path_names, cfg, sys; kwargs...) -> String

Write `run_manifest.toml` at the root of the results folder and return its path.

Tables: `[run]` (every `RunConfig` scalar, `load_factor` included), `[dispatch]`,
`[transient.*]`, `[resolved]` (facts only the built model knows — δ reference
machine, SG/GFM partition, machine constants, timeline, cost curve), `[exports]`
(the dual CSV basenames this run wrote) and `[status]`.

Called on every run, including failed ones: a manifest whose `[status]` records a
non-optimal termination is exactly what is needed to tell an empty results folder
apart from a converged one.

Keyword arguments are all optional so steady-state-only runs (no `dyn_model_dict`)
and pre-solve failures (no `obj_MVA`) write the tables they can fill.
"""
function Save_Run_Manifest!(
    path_names,
    cfg,
    sys;
    opf_dict = nothing,
    dyn_model_dict = nothing,
    dyn_parameters_dict = nothing,
    status = nothing,
    obj_MVA = nothing,
    t_build = nothing,
    t_solve = nothing,
    filename::String = "run_manifest.toml",
)::String
    manifest = OrderedDict{String, Any}()

    run_table = _config_table(cfg)
    # Both are promoted to top-level tables below; leaving them nested as well would
    # duplicate the whole tree under [run].
    delete!(run_table, "dispatch")
    delete!(run_table, "transient")
    manifest["run"] = run_table
    manifest["dispatch"] = _config_table(cfg.dispatch)
    cfg.transient === nothing || (manifest["transient"] = _config_table(cfg.transient))

    manifest["resolved"] = _resolved_table(cfg, sys, dyn_model_dict, dyn_parameters_dict)
    manifest["exports"] = _exports_table(opf_dict, dyn_model_dict)

    status_table = OrderedDict{String, Any}()
    status === nothing || (status_table["termination_status"] = string(status))
    for (value, key) in ((obj_MVA, "objective_MVA"), (t_build, "t_build"), (t_solve, "t_solve"))
        converted = _toml_value(value)
        converted === nothing || (status_table[key] = converted)
    end
    manifest["status"] = status_table

    file_path = joinpath(path_names[:pf_results_date], filename)
    open(file_path, "w") do io
        TOML.print(io, manifest; sorted = true)
    end
    return file_path
end
