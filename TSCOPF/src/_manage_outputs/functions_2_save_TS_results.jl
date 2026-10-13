
# ===================================================================================
#                   PRINT THE OPTIMIZATION MODEL IN TXT FILE
# ===================================================================================

"""`true` when the dynamic sub-model is the DQ_4TH FULL_BUS machine path."""
function _is_dq_4th_model(dyn_model_dict::OrderedDict{Symbol, Any})::Bool
    return get(get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}()), :gen_order, "") == "DQ_4TH"
end

"""Print machine/control constraints not covered by the classical 2nd-order export.

Covers the dq machine (DQ_4TH only), the AVR (DQ only) and the turbine governor. The
governor runs on the **classical** FULL_BUS path too, so this block must not be gated on
`_is_dq_4th_model` — every family is selected by `haskey`, and the dq/AVR keys are simply
absent on a classical run.
"""

function _export_mech_power_mode_label(dyn_model_dict::OrderedDict{Symbol, Any})::String
    if !haskey(dyn_model_dict, :mech_power_mode)
        return "mech_power_mode: (not set)"
    end
    mode = dyn_model_dict[:mech_power_mode]
    if mode == USE_PG
        return "mech_power_mode: USE_PG (swing mechanical power = OPF P_g)"
    else
        return "mech_power_mode: USE_PM (swing mechanical power = explicit P_m)"
    end
end

"""
Model-metadata block printed once at the top of `dynamic_model_details.txt`.

Shared by the Kron, Kron-linear and FULL_BUS exports, so the constraint listings that
follow are never interrupted by configuration lines. The two corridor blocks name the
stability constraints that were actually built — including the resolved reference machine,
which is a run-time choice under `:highest_H` and appears nowhere else in the results tree
— and the ZIP splits (FULL_BUS only) describe the load model inside the nodal balances.
"""
function _export_dyn_meta_header(dyn_model_dict::OrderedDict{Symbol, Any})::Vector{String}
    lines = String[_export_mech_power_mode_label(dyn_model_dict)]
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    if get(meta, :constrain_δ, true)
        style = get(meta, :bound_style_δ, :coi_box)
        ref = get(meta, :δ_ref_gen_resolved, nothing)
        # The reference is resolved over the surviving machines, so record which one the
        # run actually used rather than only what was asked for.
        suffix = ref === nothing ? "" : " (reference: generator $ref)"
        push!(lines, "constrain_δ: true, bound_style_δ: $(style)$(suffix)")
    else
        push!(lines, "constrain_δ: false")
    end
    haskey(meta, :ode_first_step) && push!(lines, "ode_first_step: $(meta[:ode_first_step])")
    if get(meta, :constrain_Δω, false)
        # Signed pair — may be asymmetric (Δω_tol_pu_lower / Δω_tol_pu_upper).
        lo, hi = meta[:Δω_tol]
        style = get(meta, :bound_style_Δω, :coi_box)
        push!(lines,
            "constrain_Δω: true, bound_style_Δω: $(style) (Δω_tol = [$lo, +$hi] p.u.)")
    else
        push!(lines, "constrain_Δω: false")
    end
    haskey(meta, :zip_load_p) && push!(lines, "zip_load_p (Z,I,P): $(meta[:zip_load_p])")
    haskey(meta, :zip_load_q) && push!(lines, "zip_load_q (Z,I,P): $(meta[:zip_load_q])")
    # A per-load table overrides those splits bus by bus, so saying nothing here would
    # leave the dump claiming a system-wide load model the run did not use.
    if get(meta, :load_table_rows, 0) > 0
        push!(lines, "load_data.csv: $(meta[:load_table_rows]) loads " *
                     "(per-bus ZIP = demand-weighted mean of its loads)")
    end
    return lines
end

"""Metadata for pre-fault coupling exports (TXT header + CSV provenance)."""
function _prefault_coupling_metadata(dyn_model_dict::OrderedDict{Symbol, Any})
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    # Only the FULL_BUS builders record a coupling source; the Kron paths seed
    # constants and never reach the pre-fault starts export.
    init_src = get(meta, :coupling_init_source, "none")
    return (
        gen_order = string(get(meta, :gen_order, "CLASSICAL_2ND")),
        network_form = string(get(meta, :network_form, "KRON_REDUCED")),
        mech_power_mode_label = _export_mech_power_mode_label(dyn_model_dict),
        include_avr = Bool(get(meta, :include_avr, false)),
        include_governor = Bool(get(meta, :include_governor, false)),
        coupling_init_source = string(init_src),
    )
end

"""One row per active generator: steady-state variables that couple OPF to dynamics."""
function _prefault_coupling_table(
    dyn_model_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
)::DataFrame
    vars = dyn_model_dict[:vars]
    refs = get(dyn_model_dict, :refs, OrderedDict{Symbol, Any}())
    is_dq = _is_dq_4th_model(dyn_model_dict)
    has_avr = haskey(vars, :V_ref)
    has_gov = haskey(vars, :P_ref)
    P_mech = refs[:P_mech]
    P_g = refs[:P_g]
    gen_ids = sort(collect(keys(vars[:δ])))

    rows = NamedTuple[]
    for g in gen_ids
        push!(rows, (
            gen = "G$g",
            delta_rad = JuMP.value(vars[:δ][g]),
            E_pu = is_dq ? missing :
                (haskey(vars, :E) && haskey(vars[:E], g) ? JuMP.value(vars[:E][g]) : missing),
            E_fd_pu = is_dq && haskey(vars, :E_fd) && haskey(vars[:E_fd], g) ?
                JuMP.value(vars[:E_fd][g]) : missing,
            Ed_pu = is_dq && haskey(vars, :Ed) && haskey(vars[:Ed], g) ?
                JuMP.value(vars[:Ed][g]) : missing,
            Eq_pu = is_dq && haskey(vars, :Eq) && haskey(vars[:Eq], g) ?
                JuMP.value(vars[:Eq][g]) : missing,
            Id_pu = is_dq && haskey(vars, :Id) && haskey(vars[:Id], g) ?
                JuMP.value(vars[:Id][g]) : missing,
            Iq_pu = is_dq && haskey(vars, :Iq) && haskey(vars[:Iq], g) ?
                JuMP.value(vars[:Iq][g]) : missing,
            P_mech_pu = haskey(P_mech, g) ? JuMP.value(P_mech[g]) : missing,
            P_mech_MW = haskey(P_mech, g) ? JuMP.value(P_mech[g]) * base_MVA : missing,
            P_g_pu = JuMP.value(P_g[g]),
            P_g_MW = JuMP.value(P_g[g]) * base_MVA,
            V_ref_pu = has_avr && haskey(vars[:V_ref], g) ? JuMP.value(vars[:V_ref][g]) : missing,
            P_ref_pu = has_gov && haskey(vars[:P_ref], g) ? JuMP.value(vars[:P_ref][g]) : missing,
        ))
    end
    table = DataFrame(rows)
    # `E` is the classical internal voltage behind x'_d; the dq machine carries Ed/Eq
    # instead, so the column would be empty on every DQ row. Drop it rather than export
    # a blank column. (The TXT writer only reads `E_pu` on the non-DQ branch.)
    is_dq && select!(table, Not(:E_pu))
    return table
end

function _print_prefault_txt_section(
    io::IO,
    title::String,
    gens::AbstractVector,
    values::AbstractVector,
)
    println(io, "==========================")
    println(io, "         $title")
    println(io, "==========================")
    for (g, v) in zip(gens, values)
        println(io, "$g: ", v isa Missing ? "—" : v)
    end
    println(io, "\n")
    return nothing
end

"""Human-readable pre-fault coupling block for `dynamic_model_short_results.txt`."""
function _write_prefault_short_results_txt!(
    io::IO,
    table::DataFrame,
    metadata,
)
    println(io, "==========================")
    println(io, "   Pre-fault coupling (t=0)")
    println(io, "==========================")
    println(io, "gen_order: ", metadata.gen_order)
    println(io, "network_form: ", metadata.network_form)
    println(io, "coupling_init_source: ", metadata.coupling_init_source)
    println(io, metadata.mech_power_mode_label)
    metadata.include_avr && println(io, "include_avr: true")
    metadata.include_governor && println(io, "include_governor: true")
    println(io, "\n")

    gens = table.gen
    if metadata.gen_order == "DQ_4TH"
        _print_prefault_txt_section(io, "E_fd [p.u.]", gens, table.E_fd_pu)
    else
        _print_prefault_txt_section(io, "E [p.u.]", gens, table.E_pu)
    end
    _print_prefault_txt_section(io, "δ [deg]", gens, rad2deg.(table.delta_rad))
    _print_prefault_txt_section(io, "P_mech [MW]", gens, table.P_mech_MW)
    _print_prefault_txt_section(io, "P_g [MW]", gens, table.P_g_MW)

    if metadata.gen_order == "DQ_4TH"
        for (label, col) in (
            ("Ed [p.u.]", :Ed_pu), ("Eq [p.u.]", :Eq_pu),
            ("Id [p.u.]", :Id_pu), ("Iq [p.u.]", :Iq_pu),
        )
            _print_prefault_txt_section(io, label, gens, table[!, col])
        end
    end

    metadata.include_avr &&
        _print_prefault_txt_section(io, "V_ref [p.u.]", gens, table.V_ref_pu)
    if metadata.include_governor
        _print_prefault_txt_section(io, "P_ref [p.u.]", gens, table.P_ref_pu)
        println(io, "Governor equilibrium (t=0): P_valve_0 = P_mech, P_mech_0 = P_mech")
        println(io, "\n")
    end
    return nothing
end

"""Write machine-precision pre-fault coupling table to `CSV/prefault_coupling.csv`."""
function Save_Prefault_Coupling_CSV!(
    table::DataFrame,
    path_names::OrderedDict{Symbol, String},
)
    pf_ts_csv = path_names[:pf_TS_CSV]
    mkpath(pf_ts_csv)
    CSV.write(joinpath(pf_ts_csv, "prefault_coupling.csv"), table; delim=';')
    println("Pre-fault coupling variables saved to: ", pf_ts_csv)
    return nothing
end

# ===================================================================================
# Prefault coupling START values (Dispatch_WarmStart/, before joint optimize!)
# ===================================================================================

"""`JuMP.start_value` as Float64, or `missing` if unset."""
function _start_or_missing(v::JuMP.VariableRef)::Union{Float64, Missing}
    s = JuMP.start_value(v)
    return s === nothing ? missing : Float64(s)
end

"""Start value of `vars[sym][g]` when present; else `missing`."""
function _var_start(
    vars::OrderedDict{Symbol, Any},
    sym::Symbol,
    g::Int,
)::Union{Float64, Missing}
    haskey(vars, sym) || return missing
    d = vars[sym]
    (d isa AbstractDict && haskey(d, g)) || return missing
    return _start_or_missing(d[g])
end

"""Variable ref `vars[sym][g]` when present; else `nothing`."""
function _var_ref_or_nothing(
    vars::OrderedDict{Symbol, Any},
    sym::Symbol,
    g::Int,
)
    haskey(vars, sym) || return nothing
    d = vars[sym]
    (d isa AbstractDict && haskey(d, g)) || return nothing
    return d[g]
end

"""
One row per pre-fault machine: JuMP `start=` values that link OPF → dynamics.

Defensive over classical / DQ ± AVR/TG / GFM — absent families are `missing`.
"""
function _prefault_coupling_starts_table(
    dyn_model_dict::OrderedDict{Symbol, Any},
)::DataFrame
    vars = get(dyn_model_dict, :vars, OrderedDict{Symbol, Any}())
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    gfm_set = Set{Int}(Int.(get(meta, :gfm_gens, Int[])))

    empty_cols = (
        gen=Int[], kind=String[],
        delta_rad=Union{Float64, Missing}[],
        E=Union{Float64, Missing}[], E_fd=Union{Float64, Missing}[],
        Ed=Union{Float64, Missing}[], Eq=Union{Float64, Missing}[],
        Id=Union{Float64, Missing}[], Iq=Union{Float64, Missing}[],
        P_m=Union{Float64, Missing}[],
        V_ref=Union{Float64, Missing}[], P_ref=Union{Float64, Missing}[],
        P_meas=Union{Float64, Missing}[], Q_meas=Union{Float64, Missing}[],
        V_meas=Union{Float64, Missing}[], E_int=Union{Float64, Missing}[],
        V_set=Union{Float64, Missing}[],
    )
    (!haskey(vars, :δ) || !(vars[:δ] isa AbstractDict) || isempty(vars[:δ])) &&
        return DataFrame(empty_cols)

    rows = NamedTuple[]
    for g in sort(collect(keys(vars[:δ])))
        is_gfm = g in gfm_set ||
            (haskey(vars, :P_meas) && vars[:P_meas] isa AbstractDict && haskey(vars[:P_meas], g))
        push!(rows, (
            gen = g,
            kind = is_gfm ? "GFM" : "SG",
            delta_rad = _var_start(vars, :δ, g),
            E = _var_start(vars, :E, g),
            E_fd = _var_start(vars, :E_fd, g),
            Ed = _var_start(vars, :Ed, g),
            Eq = _var_start(vars, :Eq, g),
            Id = _var_start(vars, :Id, g),
            Iq = _var_start(vars, :Iq, g),
            P_m = _var_start(vars, :P_m, g),
            V_ref = _var_start(vars, :V_ref, g),
            P_ref = _var_start(vars, :P_ref, g),
            P_meas = _var_start(vars, :P_meas, g),
            Q_meas = _var_start(vars, :Q_meas, g),
            V_meas = _var_start(vars, :V_meas, g),
            E_int = _var_start(vars, :E_int, g),
            V_set = _var_start(vars, :V_set, g),
        ))
    end
    return DataFrame(rows)
end

"""Dump-style TXT section: only gens that have `sym`, labeled with JuMP variable names."""
function _print_prefault_starts_section!(
    io::IO,
    title::String,
    vars::OrderedDict{Symbol, Any},
    sym::Symbol,
    gen_ids::AbstractVector{Int},
)
    pairs = Tuple{String, Union{Float64, Missing}}[]
    for g in gen_ids
        v = _var_ref_or_nothing(vars, sym, g)
        v === nothing && continue
        push!(pairs, (string(JuMP.name(v)), _start_or_missing(v)))
    end
    isempty(pairs) && return nothing
    println(io, "==========================")
    println(io, "         $title")
    println(io, "==========================")
    for (i, (nm, val)) in enumerate(pairs)
        println(io, "$i: $nm = ", val isa Missing ? "—" : val)
    end
    println(io, "\n")
    return nothing
end

"""Human-readable pre-fault START snapshot (before joint TSC optimize!)."""
function _write_prefault_starts_txt!(
    io::IO,
    dyn_model_dict::OrderedDict{Symbol, Any},
    table::DataFrame,
)
    meta = _prefault_coupling_metadata(dyn_model_dict)
    vars = get(dyn_model_dict, :vars, OrderedDict{Symbol, Any}())
    println(io, "========================================")
    println(io, "  Prefault coupling START values (t=0)")
    println(io, "  (JuMP start= before joint optimize!)")
    println(io, "========================================")
    println(io, "gen_order: ", meta.gen_order)
    println(io, "network_form: ", meta.network_form)
    println(io, "coupling_init_source: ", meta.coupling_init_source)
    println(io, meta.mech_power_mode_label)
    meta.include_avr && println(io, "include_avr: true")
    meta.include_governor && println(io, "include_governor: true")
    println(io, "\n")

    isempty(table) && return nothing
    gen_ids = Int.(table.gen)

    # Order matches typical dynamic_model_details pre-fault listing.
    for (title, sym) in (
        ("E_fd", :E_fd),
        ("E (classical)", :E),
        ("Ed", :Ed),
        ("Eq", :Eq),
        ("Id", :Id),
        ("Iq", :Iq),
        ("δ", :δ),
        ("P_m / P_set", :P_m),
        ("P_meas", :P_meas),
        ("Q_meas", :Q_meas),
        ("V_meas", :V_meas),
        ("E_int", :E_int),
        ("V_set", :V_set),
        ("V_ref", :V_ref),
        ("P_ref", :P_ref),
    )
        _print_prefault_starts_section!(io, title, vars, sym, gen_ids)
    end
    return nothing
end

"""
Write pre-fault JuMP `start=` values to `Dispatch_WarmStart/` (not Transient_Stability).

Called after TS assembly when ACOPF warm-start was used, before joint `optimize!`.
"""
function Save_Prefault_Coupling_Starts!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    path_names::OrderedDict{Symbol, String},
)
    ws_dir = path_names[:pf_dispatch_warmstart]
    ws_csv = path_names[:pf_dispatch_warmstart_CSV]
    mkpath(ws_dir)
    mkpath(ws_csv)

    table = _prefault_coupling_starts_table(dyn_model_dict)
    CSV.write(joinpath(ws_csv, "prefault_coupling_starts.csv"), table; delim=';')
    open(joinpath(ws_dir, "prefault_coupling_starts.txt"), "w") do io
        _write_prefault_starts_txt!(io, dyn_model_dict, table)
    end
    println("Prefault coupling START values saved to: ", ws_dir)
    return nothing
end

# δCOI_tf is optional like the rest: a machine-referenced δ corridor keeps the COI as an
# expression, so there is no variable for this listing to show.


# Inertia H for generator `gen` from SG-only `DGEN_DYN` (0 for GFM / missing).
function _inertia_H(DGEN_DYN::DataFrame, gen::Int)::Float64
    if :id in propertynames(DGEN_DYN)
        i = findfirst(==(gen), Int.(DGEN_DYN.id))
        i === nothing && return 0.0
        return Float64(DGEN_DYN.H[i])
    end
    (1 <= gen <= nrow(DGEN_DYN)) || return 0.0
    return Float64(DGEN_DYN.H[gen])
end

# Total SG inertia over `gens`, or `nothing` for an all-GFM fleet.
#
# `nothing` rather than an error: an inertia-free fleet is a legitimate study case
# (100 % converter-based), and the H-weighted COI is simply undefined there. Every
# caller then drops the COI-weighted correction term, whose numerator is itself zero
# when every H is zero, so the remaining quantities stay correct. Throwing here used
# to abort *all* post-processing after a successful solve.
function _H_total(DGEN_DYN::DataFrame, gens)::Union{Nothing, Float64}
    s = sum(_inertia_H(DGEN_DYN, Int(g)) for g in gens; init=0.0)
    s > 0.0 || return nothing
    return s
end

# Build ΔωCOI and Δω_COI trajectories; use model ΔωCOI vars when constrained.
function _extract_ΔωCOI_trajectories(
    dyn_model_dict::OrderedDict{Symbol, Any},
    Δωt::OrderedDict{Int, Vector{Float64}},
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int},
    t_window_len::Int,
)::Tuple{Vector{Float64}, OrderedDict{Int, Vector{Float64}}}
    if haskey(dyn_model_dict[:vars], :ΔωCOI_tf)
        ΔωCOI_tf_values = [JuMP.value(v) for (_, v) in dyn_model_dict[:vars][:ΔωCOI_tf]]
        if haskey(dyn_model_dict[:vars], :ΔωCOI_tpf)
            ΔωCOI_tpf_values = [JuMP.value(v) for (_, v) in dyn_model_dict[:vars][:ΔωCOI_tpf]]
            ΔωCOIt = vcat(ΔωCOI_tf_values, ΔωCOI_tpf_values)
        else
            ΔωCOIt = ΔωCOI_tf_values
        end
    else
        ΔωCOIt = zeros(Float64, t_window_len)
        for (k, Δω_vector) in Δωt
            Hk = _inertia_H(DGEN_DYN, Int(k))
            Hk == 0.0 && continue  # GFM / no inertia — SG-only COI
            ΔωCOIt .+= Hk .* Δω_vector
        end
        # `nothing` = all-GFM fleet: the SG-weighted COI has no members, so the
        # reference stays at zero and Δω_COI degenerates to the absolute speed
        # deviation. The accumulator above is already all-zero in that case.
        H_total = _H_total(DGEN_DYN, active_gen)
        H_total === nothing || (ΔωCOIt ./= H_total)
    end
    Δω_COIt = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Δωt)
        Δω_COIt[k] = Δωt[k] .- ΔωCOIt
    end
    return ΔωCOIt, Δω_COIt
end

# =====================================================================================
# Exhaustive constraint reporting for dynamic_model_details.txt
# =====================================================================================
#
# The hand-written sections below name the families they print. That is fine for the
# ones they know about and useless for the ones they do not: the writer named 16 of
# roughly 90 inequality families, and every omission was SILENT. Worse than silent, in
# one case — the "Angle in Relation to the COI" header prints unconditionally and each
# side is then guarded by `haskey`, so a machine-referenced run
# (`bound_style_δ ∈ (:highest_H, :ref_gen)`, which stores under `δ_ref_*`) produced an
# empty section that reads as "this run has no angle corridor" when it has one. Same
# for `bound_style_Δω = :abs` and for the always-on FULL_BUS voltage floor
# `ineq_const_V_tf_lower`, which was invisible despite binding on the faulted bus.
#
# The three helpers here iterate `dyn_model_dict` instead of naming families, so a
# builder that gains a constraint shows up without anyone editing this file.

"""
Rows in one family container.

Three shapes occur in the model dictionaries and all three are counted here: `id → ref`,
`id → t → ref`, and `id → t → Vector{ref}` (the GFM dynamics blocks, which register a
whole vector of constraints per generator and step).
"""
function _family_row_count(container)::Int
    container isa AbstractDict || return 0
    n = 0
    for (_, v) in container
        if v isa AbstractDict
            for (_, w) in v
                n += w isa AbstractVector ? length(w) : 1
            end
        elseif v isa AbstractVector
            n += length(v)
        else
            n += 1
        end
    end
    return n
end

"""
Print one family in full, whichever container shape it uses.

The only per-family printer in this file. Everything a builder registers reaches the TXT
through here, so a new constraint or variable family needs no writer edit — only a title
in `_FAMILY_TITLES` if the raw symbol is not readable enough.
"""
function _dump_family!(io::IO, title::AbstractString, container)
    container isa AbstractDict || return nothing
    println(io, "-"^70)
    println(io, title, "   [", _family_row_count(container), " rows]")
    println(io, "-"^70)
    for (k1, v1) in container
        if v1 isa AbstractDict
            println(io, " ******* $(k1) ****** ")
            for (k2, v2) in v1
                if v2 isa AbstractVector
                    for (k, v3) in enumerate(v2)
                        println(io, "$(k2)#$(k): ", v3)
                    end
                else
                    println(io, "$(k2): ", v2)
                end
            end
            println(io)
        else
            println(io, "$(k1): ", v1)
        end
    end
    println(io, "\n")
    return nothing
end

"""
Census of every variable and constraint family in the model: name and row count, one line
each, sorted by name.

Cheap and, unlike any hand-maintained list, impossible to leave stale — it is derived
from the model that was actually built. Read it first: if a family you expect is not
here, the model did not build it; if it is here, the dump further down must contain it.
"""
function _export_constraint_census!(io::IO, dyn_model_dict::OrderedDict{Symbol, Any})
    println(io, "="^70)
    println(io, "Census — every variable and constraint family in the built model")
    println(io, "="^70)
    for (kind, key) in (("variable", :vars), ("equality", :eq_const),
                        ("inequality", :ineq_const))
        store = get(dyn_model_dict, key, nothing)
        store isa AbstractDict || continue
        fams = sort(collect(keys(store)); by = String)
        println(io, "-- $(kind): $(length(fams)) families, ",
                sum(_family_row_count(store[f]) for f in fams; init = 0), " rows --")
        for f in fams
            println(io, rpad(String(f), 44), lpad(_family_row_count(store[f]), 8), " rows")
        end
        println(io)
    end
    println(io, "\n")
    return nothing
end

"""
Which stability corridor this run actually built, named explicitly.

Without this an empty corridor section is ambiguous: it could mean the corridor was
switched off, or that it was built under a style this file does not print. Stating the
resolved style removes the ambiguity at the point of reading.
"""
function _export_corridor_style_header!(io::IO, dyn_model_dict::OrderedDict{Symbol, Any})
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    # Print a field only when the model actually carries it. A "?" placeholder would be
    # indistinguishable from a genuinely unset option, which is the ambiguity this whole
    # header exists to remove. `δ_tol` in particular lives in `dyn_parameters_dict`, not
    # in meta, so it is simply absent here — the bound itself is visible in the printed
    # corridor rows.
    function line(label, key, suffix = "")
        haskey(meta, key) || return nothing
        println(io, rpad(label, 17), ": ", meta[key], suffix)
        return nothing
    end
    style = get(meta, :bound_style_δ, nothing)
    println(io, "="^70)
    println(io, "Stability Corridors — configuration actually built")
    println(io, "="^70)
    line("constrain_δ", :constrain_δ)
    line("bound_style_δ", :bound_style_δ,
         "   (:coi_box / :swing_propagated store under δ_COI_*, " *
         ":highest_H / :ref_gen under δ_ref_*, :pairwise under δ_pair_*)")
    # Only the machine-referenced styles resolve a reference, and only they fall back to the
    # configured id. Printing that fallback under any other style would name a machine the
    # corridor never used — `δ_ref_gen_id` may still be set from an earlier run's config.
    if style ∈ (:highest_H, :ref_gen)
        let r = get(meta, :δ_ref_gen_resolved, get(meta, :δ_ref_gen_id, nothing))
            r === nothing || println(io, rpad("δ reference gen", 17), ": ", r,
                                     "   (this machine carries no corridor row of its own,",
                                     " so δ_ref_* has one row per OTHER machine)")
        end
    end
    # The pairwise band lives in meta (unlike δ_tol, which the builders keep in
    # `dyn_parameters_dict`), so it can be reported here — and it must be, because the
    # printed rows below give the corridor only implicitly.
    if style === :pairwise
        line("δ_pair_tol", :δ_pair_tol, " rad")
        ineq = get(dyn_model_dict, :ineq_const, nothing)
        if ineq isa AbstractDict
            for key in (:ineq_const_δ_pair_tf_upper, :ineq_const_δ_pair_tf_lower)
                haskey(ineq, key) || continue
                println(io, rpad("δ pairs bounded", 17), ": ", length(ineq[key]),
                        "   (every unordered pair of synchronised units, converters included)")
                break
            end
        end
    end
    line("constrain_Δω", :constrain_Δω)
    line("bound_style_Δω", :bound_style_Δω,
         "   (:coi_box stores under Δω_COI_*, :abs under Δω_abs_*)")
    line("Δω_tol", :Δω_tol, " pu")
    println(io, "\n")
    return nothing
end

# Families the hand-written sections already print in full. Anything NOT listed here is
# dumped by `_export_uncovered_ineq_const!`. A name that drifts out of this list causes
# a family to be printed TWICE — visible and harmless — whereas the previous design
# failed the other way, dropping families with no trace. Duplication is the safer error.
const _FAMILY_TITLES = Dict{Symbol, String}(
    # --- variables ----------------------------------------------------------------
    :E                 => "Pre-Fault — Internal EMF E",
    :E_fd              => "Pre-Fault — Field Voltage E_fd",
    :δ                 => "Pre-Fault — Rotor Angle δ",
    :P_m               => "Pre-Fault — Mechanical Power P_m",
    :Ed                => "Pre-Fault — Subtransient EMF Ed",
    :Eq                => "Pre-Fault — Subtransient EMF Eq",
    :Id                => "Pre-Fault — dq Current Id",
    :Iq                => "Pre-Fault — dq Current Iq",
    :V_ref             => "Pre-Fault — AVR Voltage Set-point V_ref",
    :P_ref             => "Pre-Fault — Governor Set-point P_ref",
    :Pe_tf             => "Electrical Power Pe_tf",
    :δ_tf              => "Rotor Angle δ_tf",
    :Δω_tf             => "Speed Deviation Δω_tf",
    :δCOI_tf           => "COI Angle δCOI_tf",
    :ΔωCOI_tf          => "COI Speed Deviation ΔωCOI_tf",
    :Pe_tpf            => "Electrical Power Pe_tpf",
    :δ_tpf             => "Rotor Angle δ_tpf",
    :Δω_tpf            => "Speed Deviation Δω_tpf",
    :δCOI_tpf          => "COI Angle δCOI_tpf",
    :ΔωCOI_tpf         => "COI Speed Deviation ΔωCOI_tpf",
    :V_tf              => "Bus Voltage Magnitude V_tf",
    :θ_tf              => "Bus Voltage Angle θ_tf",
    :Qe_tf             => "Generator Reactive Power Qe_tf",
    :Ed_tf             => "Subtransient EMF Ed_tf",
    :Eq_tf             => "Subtransient EMF Eq_tf",
    :Id_tf             => "dq Current Id_tf",
    :Iq_tf             => "dq Current Iq_tf",
    :Te_tf             => "Electrical Torque Te_tf",
    :E_LL_tf           => "AVR Lead-Lag E_LL_tf",
    :E_fd_unlim_tf     => "AVR E_fd_unlim_tf",
    :E_fd_tf           => "AVR E_fd_tf",
    :V_tpf             => "Bus Voltage Magnitude V_tpf",
    :θ_tpf             => "Bus Voltage Angle θ_tpf",
    :Qe_tpf            => "Generator Reactive Power Qe_tpf",
    :Ed_tpf            => "Subtransient EMF Ed_tpf",
    :Eq_tpf            => "Subtransient EMF Eq_tpf",
    :Id_tpf            => "dq Current Id_tpf",
    :Iq_tpf            => "dq Current Iq_tpf",
    :Te_tpf            => "Electrical Torque Te_tpf",
    :E_LL_tpf          => "AVR Lead-Lag E_LL_tpf",
    :E_fd_unlim_tpf    => "AVR E_fd_unlim_tpf",
    :E_fd_tpf          => "AVR E_fd_tpf",
    :Pv_raw_tf         => "Governor Valve Raw P_valve_raw_tf",
    :Pv_tf             => "Governor Valve P_valve_tf",
    :Pm_tf             => "Governor Mech Power P_mech_tf",
    :Pv_raw_tpf        => "Governor Valve Raw P_valve_raw_tpf",
    :Pv_tpf            => "Governor Valve P_valve_tpf",
    :Pm_tpf            => "Governor Mech Power P_mech_tpf",
    :P_meas            => "GFM P_meas",
    :Q_meas            => "GFM Q_meas",
    :V_meas            => "GFM V_meas",
    :E_int             => "GFM E_int",
    :V_set             => "GFM V_set",
    :P_meas_tf         => "GFM P_meas_tf",
    :Q_meas_tf         => "GFM Q_meas_tf",
    :V_meas_tf         => "GFM V_meas_tf",
    :E_int_raw_tf      => "GFM E_int_raw_tf",
    :E_int_tf          => "GFM E_int_tf",
    :E_droop_raw_tf    => "GFM E_droop_raw_tf",
    :E_droop_tf        => "GFM E_droop_tf",
    :P_meas_tpf        => "GFM P_meas_tpf",
    :Q_meas_tpf        => "GFM Q_meas_tpf",
    :V_meas_tpf        => "GFM V_meas_tpf",
    :E_int_raw_tpf     => "GFM E_int_raw_tpf",
    :E_int_tpf         => "GFM E_int_tpf",
    :E_droop_raw_tpf   => "GFM E_droop_raw_tpf",
    :E_droop_tpf       => "GFM E_droop_tpf",

    # --- equality constraints -----------------------------------------------------
    :eq_const_P_init          => "Initial Active Electrical Power of Generators",
    :eq_const_Q_init          => "Initial Reactive Electrical Power of Generators",
    :eq_const_Pm_init         => "Initial Mechanical Power of Generators",
    :eq_const_δCOI_tf         => "Angle of the COI (Fault Period)",
    :eq_const_δCOI_tpf        => "Angle of the COI (Post-Fault Period)",
    :eq_const_ΔωCOI_tf        => "Speed Deviation of COI (Fault Period)",
    :eq_const_ΔωCOI_tpf       => "Speed Deviation of COI (Post-Fault Period)",
    :eq_const_Pe_tf           => "Electrical Power (Fault Period)",
    :eq_const_Pe_tpf          => "Electrical Power (Post-Fault Period)",
    :eq_const_δ_tf            => "Angle — Swing Equation (Fault Period)",
    :eq_const_δ_tpf           => "Angle — Swing Equation (Post-Fault Period)",
    :eq_const_Δω_tf           => "Speed Deviation — Swing Equation (Fault Period)",
    :eq_const_Δω_tpf          => "Speed Deviation — Swing Equation (Post-Fault Period)",
    :eq_const_Ed_init         => "DQ Init — Subtransient EMF Ed",
    :eq_const_Eq_init         => "DQ Init — Subtransient EMF Eq",
    :eq_const_Vd_init         => "DQ Init — Stator Voltage Vd",
    :eq_const_Vq_init         => "DQ Init — Stator Voltage Vq",
    :eq_const_Efd_init        => "DQ Init — Exciter E_fd",
    :eq_const_Pref_init       => "Governor Init — Set-point P_ref",
    :eq_const_Te_tf           => "DQ Te = Ed·Id + Eq·Iq (Fault Period)",
    :eq_const_Te_tpf          => "DQ Te = Ed·Id + Eq·Iq (Post-Fault Period)",
    :eq_const_Vd_tf           => "DQ Stator Vd (Fault Period)",
    :eq_const_Vd_tpf          => "DQ Stator Vd (Post-Fault Period)",
    :eq_const_Vq_tf           => "DQ Stator Vq (Fault Period)",
    :eq_const_Vq_tpf          => "DQ Stator Vq (Post-Fault Period)",
    :eq_const_Ed_tf           => "DQ EMF Ed dynamics (Fault Period)",
    :eq_const_Ed_tpf          => "DQ EMF Ed dynamics (Post-Fault Period)",
    :eq_const_Eq_tf           => "DQ EMF Eq dynamics (Fault Period)",
    :eq_const_Eq_tpf          => "DQ EMF Eq dynamics (Post-Fault Period)",
    :eq_const_avr_leadlag_tf  => "AVR Lead-Lag (E_LL) (Fault Period)",
    :eq_const_avr_leadlag_tpf => "AVR Lead-Lag (E_LL) (Post-Fault Period)",
    :eq_const_E_fd_unlim_tf   => "AVR Exciter ODE (E_fd_unlim) (Fault Period)",
    :eq_const_E_fd_unlim_tpf  => "AVR Exciter ODE (E_fd_unlim) (Post-Fault Period)",
    :eq_const_E_fd_tf         => "AVR Field Softsat (E_fd) (Fault Period)",
    :eq_const_E_fd_tpf        => "AVR Field Softsat (E_fd) (Post-Fault Period)",
    :eq_const_gov_valve_tf    => "Governor Valve ODE (P_valve_raw) (Fault Period)",
    :eq_const_gov_valve_tpf   => "Governor Valve ODE (P_valve_raw) (Post-Fault Period)",
    :eq_const_gov_valve_limit_tf  => "Governor Valve Softsat (P_valve) (Fault Period)",
    :eq_const_gov_valve_limit_tpf => "Governor Valve Softsat (P_valve) (Post-Fault Period)",
    :eq_const_gov_mech_tf     => "Governor Mech Power ODE (P_mech) (Fault Period)",
    :eq_const_gov_mech_tpf    => "Governor Mech Power ODE (P_mech) (Post-Fault Period)",
    :eq_const_Qe_tf           => "Generator Qe (Fault Period)",
    :eq_const_Qe_tpf          => "Generator Qe (Post-Fault Period)",
    :eq_const_Pbalance_tf     => "Active Power Balance (Fault Period)",
    :eq_const_Pbalance_tpf    => "Active Power Balance (Post-Fault Period)",
    :eq_const_Qbalance_tf     => "Reactive Power Balance (Fault Period)",
    :eq_const_Qbalance_tpf    => "Reactive Power Balance (Post-Fault Period)",
    :eq_const_gfm_Pmeas_init  => "GFM Init — P_meas",
    :eq_const_gfm_Qmeas_init  => "GFM Init — Q_meas",
    :eq_const_gfm_Vmeas_init  => "GFM Init — V_meas",
    :eq_const_gfm_Vset_init   => "GFM Init — V_set equilibrium",
    :eq_const_gfm_tf          => "GFM dynamics — filter/droop/PI/limiter (Fault Period)",
    :eq_const_gfm_tpf         => "GFM dynamics — filter/droop/PI/limiter (Post-Fault Period)",

    # --- inequality constraints ---------------------------------------------------
    :ineq_const_δ_COI_tf_lower   => "δ Corridor vs COI — Fault Period, Lower Bound",
    :ineq_const_δ_COI_tf_upper   => "δ Corridor vs COI — Fault Period, Upper Bound",
    :ineq_const_δ_COI_tpf_lower  => "δ Corridor vs COI — Post-Fault Period, Lower Bound",
    :ineq_const_δ_COI_tpf_upper  => "δ Corridor vs COI — Post-Fault Period, Upper Bound",
    :ineq_const_δ_ref_tf_lower   => "δ Corridor vs Reference Machine — Fault Period, Lower Bound",
    :ineq_const_δ_ref_tf_upper   => "δ Corridor vs Reference Machine — Fault Period, Upper Bound",
    :ineq_const_δ_ref_tpf_lower  => "δ Corridor vs Reference Machine — Post-Fault Period, Lower Bound",
    :ineq_const_δ_ref_tpf_upper  => "δ Corridor vs Reference Machine — Post-Fault Period, Upper Bound",
    :ineq_const_δ_pair_tf_lower  => "δ Corridor, Pairwise Separation — Fault Period, Lower Bound",
    :ineq_const_δ_pair_tf_upper  => "δ Corridor, Pairwise Separation — Fault Period, Upper Bound",
    :ineq_const_δ_pair_tpf_lower => "δ Corridor, Pairwise Separation — Post-Fault Period, Lower Bound",
    :ineq_const_δ_pair_tpf_upper => "δ Corridor, Pairwise Separation — Post-Fault Period, Upper Bound",
    :ineq_const_Δω_COI_tf_lower  => "Δω Corridor vs COI — Fault Period, Lower Bound",
    :ineq_const_Δω_COI_tf_upper  => "Δω Corridor vs COI — Fault Period, Upper Bound",
    :ineq_const_Δω_COI_tpf_lower => "Δω Corridor vs COI — Post-Fault Period, Lower Bound",
    :ineq_const_Δω_COI_tpf_upper => "Δω Corridor vs COI — Post-Fault Period, Upper Bound",
    :ineq_const_Δω_abs_tf_lower  => "Δω Corridor, Absolute — Fault Period, Lower Bound",
    :ineq_const_Δω_abs_tf_upper  => "Δω Corridor, Absolute — Fault Period, Upper Bound",
    :ineq_const_Δω_abs_tpf_lower => "Δω Corridor, Absolute — Post-Fault Period, Lower Bound",
    :ineq_const_Δω_abs_tpf_upper => "Δω Corridor, Absolute — Post-Fault Period, Upper Bound",
    :ineq_const_E_lower          => "Initial Operating Point Box — E, Lower",
    :ineq_const_E_upper          => "Initial Operating Point Box — E, Upper",
    :ineq_const_E_fd_lower       => "Initial Operating Point Box — E_fd, Lower",
    :ineq_const_E_fd_upper       => "Initial Operating Point Box — E_fd, Upper",
    :ineq_const_δ_lower          => "Initial Operating Point Box — δ, Lower",
    :ineq_const_δ_upper          => "Initial Operating Point Box — δ, Upper",
    :ineq_const_P_m_lower        => "Initial Operating Point Box — P_m, Lower",
    :ineq_const_P_m_upper        => "Initial Operating Point Box — P_m, Upper",
    :ineq_const_V_tf_lower       => "Bus Voltage Floor (Fault Period)",
    :ineq_const_V_tpf_lower      => "Bus Voltage Floor (Post-Fault Period)",
)

"""Readable title for a family, falling back to the symbol so nothing can go unnamed."""
function _family_title(sym::Symbol)::String
    title = get(_FAMILY_TITLES, sym, nothing)
    return title === nothing ? String(sym) : string(title, "   [", sym, "]")
end

"""
Which period a family belongs to: 1 pre-fault, 2 fault, 3 post-fault.

Read off the name rather than from a list of families, so a new trajectory variable lands
in the right block on its own. `"_tpf"` does not contain `"_tf"`, so the two tests are
independent and their order does not matter.
"""
function _family_period(sym::Symbol)::Int
    name = String(sym)
    occursin("_tpf", name) && return 3
    occursin("_tf", name) && return 2
    return 1
end

"""Dump every family of one store, in build order, grouped and titled."""
function _export_family_section!(
    io::IO,
    dyn_model_dict::OrderedDict{Symbol, Any},
    store_key::Symbol,
    heading::AbstractString;
    by_period::Bool = false,
)
    store = get(dyn_model_dict, store_key, nothing)
    println(io, "="^70)
    println(io, heading)
    println(io, "="^70)
    if !(store isa AbstractDict) || isempty(store)
        println(io, "(none — the model built no family of this kind)\n\n")
        return nothing
    end
    # Insertion order is build order (these are OrderedDicts), which reads init → fault →
    # post-fault without anyone maintaining a sequence. Grouping by period only re-buckets
    # that order; it never sorts within a bucket.
    groups = by_period ?
        (("Pre-Fault Period", 1), ("Fault Period", 2), ("Post-Fault Period", 3)) :
        (("", 0),)
    for (label, period) in groups
        fams = period == 0 ? collect(keys(store)) :
               [f for f in keys(store) if _family_period(f) == period]
        isempty(fams) && continue
        if !isempty(label)
            println(io, "~"^70)
            println(io, label, "  (", length(fams), " families)")
            println(io, "~"^70, "\n")
        end
        for f in fams
            _dump_family!(io, _family_title(f), store[f])
        end
    end
    return nothing
end

"""
    Export_Dynamic_Model!(model, path_names, dyn_model_dict)

Write `dynamic_model_details.txt` (and `model_summary.txt`) for **every** dynamic path —
Kron classical, Kron linear, FULL_BUS classical, FULL_BUS DQ_4TH, with or without AVR,
governor and GFM converters.

Layout, top to bottom:

1. **Summary** — model configuration, the stability corridors actually built, and a census
   of every variable / equality / inequality family with its row count.
2. **Variables** — grouped pre-fault / fault / post-fault.
3. **Equality constraints** — every family, in build order.
4. **Inequality constraints** — every family, in build order.

Written **before** `optimize!`, so it shows the assembled model rather than a solution.

Nothing here names a family: the sections iterate `dyn_model_dict`'s own containers, so a
builder that gains a constraint or a variable shows up unaided. That is the point. The
three writers this replaces drove their sections off hand-written family lists with a
FULL_BUS appendix and a GFM appendix bolted on after the census — which both mis-ordered
the file (summary material at line 43 486 of 96 876) and let whole families go unprinted
when a name drifted out of a list.
"""
function Export_Dynamic_Model!(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    pf_ts = path_names[:pf_TS]
    mkpath(pf_ts)

    open(io -> show(io, model), joinpath(pf_ts, "model_summary.txt"), "w")

    open(joinpath(pf_ts, "dynamic_model_details.txt"), "w") do io
        # --- 1. summary -----------------------------------------------------------
        println(io, "="^70)
        println(io, "Dynamic Model Summary")
        println(io, "="^70)
        meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
        haskey(meta, :network_form) && println(io, "network_form: ", meta[:network_form])
        haskey(meta, :gen_order) && println(io, "gen_order: ", meta[:gen_order])
        for line in _export_dyn_meta_header(dyn_model_dict)
            println(io, line)
        end
        println(io, "\n")
        _export_corridor_style_header!(io, dyn_model_dict)
        _export_constraint_census!(io, dyn_model_dict)

        # --- 2. variables ---------------------------------------------------------
        _export_family_section!(io, dyn_model_dict, :vars,
                                "Variables — every family in the built model";
                                by_period = true)

        # Variables the dynamic model READS but does not own — under `mech_power_mode =
        # USE_PG` the swing has no `P_m` of its own and drives off the dispatch `P_g`, so
        # without this section such a run shows no mechanical power anywhere.
        refs = get(dyn_model_dict, :refs, nothing)
        if refs isa AbstractDict && !isempty(refs)
            _export_family_section!(io, dyn_model_dict, :refs,
                                    "References — dispatch variables the dynamic model reads")
        end

        # Expressions are neither variables nor constraints, but on the FULL_BUS paths Qe
        # enters the nodal balances as one, so it belongs with the quantities rather than
        # in a footnote. Absent where Qe is a variable.
        exprs = get(dyn_model_dict, :expressions, nothing)
        if exprs isa AbstractDict && !isempty(exprs)
            _export_family_section!(io, dyn_model_dict, :expressions,
                                    "Expressions — derived quantities used in the model";
                                    by_period = true)
        end

        # --- 3. equality constraints ----------------------------------------------
        _export_family_section!(io, dyn_model_dict, :eq_const,
                                "Equality Constraints — every family in the built model")

        # --- 4. inequality constraints --------------------------------------------
        _export_family_section!(io, dyn_model_dict, :ineq_const,
                                "Inequality Constraints — every family in the built model")
    end

    println("Dynamic Model successfully saved as TXT file in: ", pf_ts)
    return nothing
end

"""
    Export_Variable_Bounds!(model, outdir)

Write `variable_bounds.txt` listing every JuMP variable that has a VariableRef
lower and/or upper bound (including ±Inf). Unbounded vars are counted only in
the header. Does **not** list constraint-encoded boxes (`BoundEncoding::CONSTRAINT`).

Parity artefact with the reference `Transient_Stability/variable_bounds.txt`.
"""
function Export_Variable_Bounds!(model::Model, outdir::AbstractString)
    mkpath(outdir)
    vars = JuMP.all_variables(model)
    n_all = length(vars)

    rows = NamedTuple{(:name, :lb, :ub, :finite_lb, :finite_ub),
        Tuple{String, String, String, Bool, Bool}}[]
    n_unbounded = 0
    n_finite_box = 0
    n_one_sided = 0

    for v in vars
        has_lb = JuMP.has_lower_bound(v)
        has_ub = JuMP.has_upper_bound(v)
        if !has_lb && !has_ub
            n_unbounded += 1
            continue
        end
        lb = has_lb ? JuMP.lower_bound(v) : NaN
        ub = has_ub ? JuMP.upper_bound(v) : NaN
        flb = has_lb && isfinite(lb)
        fub = has_ub && isfinite(ub)
        if flb && fub
            n_finite_box += 1
        elseif flb || fub
            n_one_sided += 1
        end
        push!(rows, (
            name = String(JuMP.name(v)),
            lb = has_lb ? string(lb) : "—",
            ub = has_ub ? string(ub) : "—",
            finite_lb = flb,
            finite_ub = fub,
        ))
    end
    sort!(rows; by = r -> r.name)

    path = joinpath(outdir, "variable_bounds.txt")
    open(path, "w") do io
        println(io, "source: TSCOPF package")
        println(io, "artefact: JuMP VariableRef bounds only (not constraint-encoded boxes)")
        println(io, "num_variables              = $n_all")
        println(io, "num_with_any_var_bound     = $(length(rows))")
        println(io, "num_finite_box_lb_and_ub   = $n_finite_box")
        println(io, "num_one_sided_finite       = $n_one_sided")
        println(io, "num_unbounded_no_var_bound = $n_unbounded")
        println(io, "num_±Inf_only_or_partial   = $(length(rows) - n_finite_box - n_one_sided)")
        println(io)
        println(io, "name\tlb\tub\tfinite_lb\tfinite_ub")
        for r in rows
            println(io, "$(r.name)\t$(r.lb)\t$(r.ub)\t$(r.finite_lb)\t$(r.finite_ub)")
        end
    end
    println("Variable bounds saved to: ", path)
    return path
end


# ===================================================================================
#                  PRINT THE DUALS OF THE OPTIMIZATION PROBLEM
# ===================================================================================

# Function to obtain the duals of the Dynamic Simulation
function Save_Duals_Dynamic_Model_tsred(model::Model,
    path_names::OrderedDict{Symbol, String},
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any};
    save_ts_plots::Bool=false,
    )

    pf_ts = path_names[:pf_TS]
    pf_ts_duals_csv = path_names[:pf_TS_CSV_duals]
    pf_ts_figures_duals = path_names[:pf_TS_figures_duals]

    # `run_case!` treats ITERATION_LIMIT / ALMOST_LOCALLY_SOLVED as "solved enough" to
    # export primal trajectories, but a solve that stopped early carries no meaningful
    # dual certificate. Exporting it anyway would put numbers into the CSVs that look
    # like prices and are not. Ask JuMP directly rather than re-deriving from the status.
    if !_has_usable_duals(model)
        @warn "Dual solution unavailable (termination: $(JuMP.termination_status(model)), " *
              "dual status: $(JuMP.dual_status(model))). Skipping the TS dual export."
        return nothing
    end

    # =================================================================================================
    #                     Extract constraint duals (registry-driven, Phase 1.3)
    # =================================================================================================
    # At build time, `build_dual_registry!` stored in `dyn_model_dict[:dual_registry]`
    # the subset of `CLASSICAL_KRON_DUAL_SPECS` whose JuMP keys are present.
    # `extract_registered_duals` walks that list and returns a Dict keyed by
    # `export_name` (:dual_Pe, :dual_δ_COI_lower, …).
    #
    # Optional families (Pm_init for USE_PM, ΔωCOI for constrain_Δω, E bounds
    # for TSC-ACOPF, …) appear only when the builder created those constraints.
    #
    # Local variables below are kept so the TXT / CSV / XLSX writers further down
    # remain unchanged; each is `nothing` when the family was not in the model.
    _duals = extract_registered_duals(dyn_model_dict)
    _getdual(sym) = get(_duals, sym, nothing)

    # --- pre-fault initial-condition duals (λ on P/Q/Pm balance at t = 0) -------
    dual_Pe_init = _getdual(:dual_Pe_init)
    dual_Qe_init = _getdual(:dual_Qe_init)
    dual_Pm_init = _getdual(:dual_Pm_init)   # USE_PM only

    # --- COI reference equalities (λ on δ_COI and optional Δω_COI constraints) --
    dual_δCOI = _getdual(:dual_δCOI)
    dual_ΔωCOI = _getdual(:dual_ΔωCOI)      # constrain_Δω only

    # --- discretized swing / Pe equalities (per generator, fault + post-fault) ---
    dual_Pe = _getdual(:dual_Pe)
    dual_δ = _getdual(:dual_δ)
    dual_Δω = _getdual(:dual_Δω)

    # --- transient-stability inequality duals (δ and Δω w.r.t. COI) ------------
    dual_δ_COI_lower = _getdual(:dual_δ_COI_lower)
    dual_δ_COI_upper = _getdual(:dual_δ_COI_upper)
    dual_Δω_COI_lower = _getdual(:dual_Δω_COI_lower)
    dual_Δω_COI_upper = _getdual(:dual_Δω_COI_upper)

    # --- the other corridor conventions (only one δ family and one Δω family per run) ---
    # These reached the CSV/XLSX exports through the registry but had no local here, so the TXT
    # dump and the dual figures were blind to them: a `:highest_H` / `:ref_gen` run wrote no
    # corridor section and left `Figures_Duals/` empty.
    dual_δ_ref_lower = _getdual(:dual_δ_ref_lower)     # bound_style_δ ∈ (:highest_H, :ref_gen)
    dual_δ_ref_upper = _getdual(:dual_δ_ref_upper)
    dual_δ_pair_lower = _getdual(:dual_δ_pair_lower)  # bound_style_δ = :pairwise, keyed by (i,j)
    dual_δ_pair_upper = _getdual(:dual_δ_pair_upper)
    dual_Δω_abs_lower = _getdual(:dual_Δω_abs_lower)   # bound_style_Δω = :abs
    dual_Δω_abs_upper = _getdual(:dual_Δω_abs_upper)

    # Reference machine behind the δ_ref family — the id that has *no* corridor row, and so no
    # section or figure of its own below.
    δ_ref_gen = get(get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}()),
        :δ_ref_gen_resolved, nothing)

    # --- variable bound duals (explicit inequalities on E, δ, P_m) --------------
    dual_LB_E = _getdual(:dual_LB_E)        # TSC-ACOPF
    dual_UB_E = _getdual(:dual_UB_E)
    dual_LB_δ = _getdual(:dual_LB_δ)        # TSC-DCOPF linear path
    dual_UB_δ = _getdual(:dual_UB_δ)
    dual_LB_Pm = _getdual(:dual_LB_Pm)      # USE_PM only
    dual_UB_Pm = _getdual(:dual_UB_Pm)

    # --- FULL_BUS network duals (Qe trajectories, nodal KCL balances) ----------
    dual_Qe = _getdual(:dual_Qe)
    dual_Pbalance = _getdual(:dual_Pbalance)
    dual_Qbalance = _getdual(:dual_Qbalance)
    dual_LB_V = _getdual(:dual_LB_V)

    # --- turbine-governor duals (set-point + valve / mech ODEs) ----------------
    dual_Pref_init = _getdual(:dual_Pref_init)   # include_governor only
    dual_gov_valve = _getdual(:dual_gov_valve)
    dual_gov_mech = _getdual(:dual_gov_mech)

    # --- DQ_4TH machine duals (init, E_fd bounds, algebra, EMF dynamics) --------
    dual_Ed_init = _getdual(:dual_Ed_init)
    dual_Eq_init = _getdual(:dual_Eq_init)
    dual_Vd_init = _getdual(:dual_Vd_init)
    dual_Vq_init = _getdual(:dual_Vq_init)
    dual_LB_E_fd = _getdual(:dual_LB_E_fd)
    dual_UB_E_fd = _getdual(:dual_UB_E_fd)
    dual_Te = _getdual(:dual_Te)
    dual_Vd = _getdual(:dual_Vd)
    dual_Vq = _getdual(:dual_Vq)
    dual_Ed = _getdual(:dual_Ed)
    dual_Eq = _getdual(:dual_Eq)
    dual_Vref_init = _getdual(:dual_Vref_init)   # include_avr only
    dual_avr_leadlag = _getdual(:dual_avr_leadlag)   # lead-lag stage only
    dual_avr_E_fd = _getdual(:dual_avr_E_fd)
    dual_avr_E_fd_sat = _getdual(:dual_avr_E_fd_sat)


    # ========== WRITE TO TXT FILE ==========
    open(joinpath(pf_ts, "dynamic_model_duals.txt"), "w") do io
        _write_solve_status_header!(io, model, "Transient-stability dual solution")
        # Which machine the δ_ref sections are measured against. Without it the surviving
        # [G…] labels read as an arbitrary subset — the id that is missing IS the reference.
        if δ_ref_gen !== nothing
            println(io, "δ corridor reference: generator $(δ_ref_gen)")
            println(io)
        elseif dual_δ_pair_upper !== nothing
            # The pairwise corridor has no reference to name; what the reader needs instead
            # is how many pairs the [G<i>-G<j>] labels below are drawn from.
            println(io, "δ corridor: pairwise separation, " *
                        "$(length(dual_δ_pair_upper)) unordered pairs (no reference machine)")
            println(io)
        end

        function write_dual_active_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MW $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_reactive_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MVAr $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_apparent_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MVA $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_voltage_pu(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/p.u. $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_angle_rad(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/rad $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_speed_pu(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/p.u. $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_others(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t $val")
            end
            println(io)  # empty line between sections
        end

        dual_Pe_init !== nothing && write_dual_active_power(io, "dual_Pe_init", dual_Pe_init)
        dual_Qe_init !== nothing && write_dual_reactive_power(io, "dual_Qe_init", dual_Qe_init)
        dual_Pm_init !== nothing && write_dual_active_power(io, "dual_Pm_init", dual_Pm_init)
        dual_δCOI !== nothing && write_dual_angle_rad(io, "dual_eq_δCOI", dual_δCOI)
        dual_ΔωCOI !== nothing && write_dual_speed_pu(io, "dual_eq_ΔωCOI", dual_ΔωCOI)
        dual_Pe !== nothing && for (i, info) in dual_Pe write_dual_active_power(io, "dual_eq_Pe[G$i]", info) end
        dual_Qe !== nothing && for (i, info) in dual_Qe write_dual_reactive_power(io, "dual_eq_Qe[G$i]", info) end
        dual_Pbalance !== nothing && for (i, info) in dual_Pbalance write_dual_active_power(io, "dual_eq_Pbalance[bus$i]", info) end
        dual_Qbalance !== nothing && for (i, info) in dual_Qbalance write_dual_reactive_power(io, "dual_eq_Qbalance[bus$i]", info) end
        dual_LB_V !== nothing && for (i, info) in dual_LB_V write_dual_voltage_pu(io, "dual_ineq_V_lower[bus$i]", info) end
        dual_δ !== nothing && for (i, info) in dual_δ write_dual_angle_rad(io, "dual_eq_δ[G$i]", info) end
        dual_Δω !== nothing && for (i, info) in dual_Δω write_dual_speed_pu(io, "dual_eq_Δω[G$i]", info) end
        dual_δ_COI_lower !== nothing && for (i, info) in dual_δ_COI_lower write_dual_angle_rad(io, "dual_ineq_δ_COI_lower[G$i]", info) end
        dual_δ_COI_upper !== nothing && for (i, info) in dual_δ_COI_upper write_dual_angle_rad(io, "dual_ineq_δ_COI_upper[G$i]", info) end
        dual_Δω_COI_lower !== nothing && for (i, info) in dual_Δω_COI_lower write_dual_speed_pu(io, "dual_ineq_Δω_COI_lower[G$i]", info) end
        dual_Δω_COI_upper !== nothing && for (i, info) in dual_Δω_COI_upper write_dual_speed_pu(io, "dual_ineq_Δω_COI_upper[G$i]", info) end
        dual_δ_ref_lower !== nothing && for (i, info) in dual_δ_ref_lower write_dual_angle_rad(io, "dual_ineq_δ_ref_lower[G$i]", info) end
        dual_δ_ref_upper !== nothing && for (i, info) in dual_δ_ref_upper write_dual_angle_rad(io, "dual_ineq_δ_ref_upper[G$i]", info) end
        dual_δ_pair_lower !== nothing && for (k, info) in dual_δ_pair_lower write_dual_angle_rad(io, "dual_ineq_δ_pair_lower[G$(k[1])-G$(k[2])]", info) end
        dual_δ_pair_upper !== nothing && for (k, info) in dual_δ_pair_upper write_dual_angle_rad(io, "dual_ineq_δ_pair_upper[G$(k[1])-G$(k[2])]", info) end
        dual_Δω_abs_lower !== nothing && for (i, info) in dual_Δω_abs_lower write_dual_speed_pu(io, "dual_ineq_Δω_abs_lower[G$i]", info) end
        dual_Δω_abs_upper !== nothing && for (i, info) in dual_Δω_abs_upper write_dual_speed_pu(io, "dual_ineq_Δω_abs_upper[G$i]", info) end
        dual_LB_E !== nothing && write_dual_voltage_pu(io, "dual_LB_E", dual_LB_E)
        dual_UB_E !== nothing && write_dual_voltage_pu(io, "dual_UB_E", dual_UB_E)
        dual_LB_δ !== nothing && write_dual_angle_rad(io, "dual_LB_δ", dual_LB_δ)
        dual_UB_δ !== nothing && write_dual_angle_rad(io, "dual_UB_δ", dual_UB_δ)
        dual_LB_Pm !== nothing && write_dual_active_power(io, "dual_LB_Pm", dual_LB_Pm)
        dual_UB_Pm !== nothing && write_dual_active_power(io, "dual_UB_Pm", dual_UB_Pm)
        dual_Ed_init !== nothing && write_dual_voltage_pu(io, "dual_Ed_init", dual_Ed_init)
        dual_Eq_init !== nothing && write_dual_voltage_pu(io, "dual_Eq_init", dual_Eq_init)
        dual_Vd_init !== nothing && write_dual_voltage_pu(io, "dual_Vd_init", dual_Vd_init)
        dual_Vq_init !== nothing && write_dual_voltage_pu(io, "dual_Vq_init", dual_Vq_init)
        dual_LB_E_fd !== nothing && write_dual_voltage_pu(io, "dual_LB_E_fd", dual_LB_E_fd)
        dual_UB_E_fd !== nothing && write_dual_voltage_pu(io, "dual_UB_E_fd", dual_UB_E_fd)
        dual_Te !== nothing && for (i, info) in dual_Te write_dual_active_power(io, "dual_eq_Te[G$i]", info) end
        dual_Vd !== nothing && for (i, info) in dual_Vd write_dual_voltage_pu(io, "dual_eq_Vd[G$i]", info) end
        dual_Vq !== nothing && for (i, info) in dual_Vq write_dual_voltage_pu(io, "dual_eq_Vq[G$i]", info) end
        dual_Ed !== nothing && for (i, info) in dual_Ed write_dual_voltage_pu(io, "dual_eq_Ed[G$i]", info) end
        dual_Eq !== nothing && for (i, info) in dual_Eq write_dual_voltage_pu(io, "dual_eq_Eq[G$i]", info) end
        dual_Vref_init !== nothing && write_dual_voltage_pu(io, "dual_Vref_init", dual_Vref_init)
        dual_avr_leadlag !== nothing && for (i, info) in dual_avr_leadlag write_dual_voltage_pu(io, "dual_eq_E_LL[G$i]", info) end
        dual_avr_E_fd !== nothing && for (i, info) in dual_avr_E_fd write_dual_voltage_pu(io, "dual_eq_E_fd[G$i]", info) end
        dual_avr_E_fd_sat !== nothing && for (i, info) in dual_avr_E_fd_sat write_dual_voltage_pu(io, "dual_eq_E_fd_sat[G$i]", info) end
    end
    println("Duals of the dynamic model successfully saved as TXT file in: ", path_names[:pf_TS])

    # ========== SAVE IN A CSV FILE ==========
    # Registry-driven: every catalog entry that declares a `csv_file` and whose family
    # exists in this run is written here. Adding a `DualRegistryEntry` is therefore the
    # only step needed to export a new constraint family — this block never needs edits.
    # (It replaced a hand-maintained list that silently skipped registered families such
    # as the optional tf/tpf bound duals and the governor valve limiter.)
    begin
        written = String[]
        for entry in dual_spec_catalog(dyn_model_dict)
            entry.csv_file === nothing && continue
            val = _getdual(entry.export_name)
            (val === nothing || isempty(val)) && continue
            csv_path = joinpath(pf_ts_duals_csv, entry.csv_file)
            if entry.layout in (GEN_INDEXED, TIME_INDEXED, TIME_INDEXED_MERGE)
                # Flat vector → single column named after the file stem, matching the
                # historical headers (`dual_delta_COI` in dual_delta_COI.csv, …), plus a
                # leading index column so the rows are self-describing. Readers that select
                # by name are unaffected; readers that aligned two files by row position
                # can now check they are talking about the same units.
                col = Symbol(first(splitext(entry.csv_file)))
                df = DataFrame(col => val)
                ids = extract_dual_entry_ids(dyn_model_dict, entry)
                if entry.layout == GEN_INDEXED && ids !== nothing && length(ids) == length(val)
                    insertcols!(df, 1, :Gen_ID => ids)
                elseif entry.layout != GEN_INDEXED
                    # Time-indexed families merge the fault and post-fault windows, so the
                    # container keys restart mid-vector; a running row counter is the only
                    # index that stays monotonic across the join.
                    insertcols!(df, 1, :Index => collect(1:length(val)))
                end
                CSV.write(csv_path, df; delim = ';')
            else
                # Per-generator / per-bus trajectories → one column per key. The bus-indexed
                # families carry bus ids, and labelling them `Gen_` (the historical default)
                # made `dual_Pbalance.csv` read as if it were per generator.
                prefix = entry.layout == PER_BUS_TIME_MERGE ? "Bus" : "Gen"
                save_ordered_dict_to_csv(val, csv_path; id_prefix = prefix)
            end
            push!(written, entry.csv_file)
        end

        println("Duals of the dynamic model successfully saved as CSV file in: ",
            pf_ts_duals_csv, " (", length(written), " files)")
    end

    # ========== SAVE IN A XLSX FILE ==========
    # Map registry export names → keyword names expected by Save_Duals_2_Excel_tsred
    Save_Duals_2_Excel_tsred(path_names; duals_to_xlsx_kwargs(_duals, dyn_model_dict)...)

    # Corridor dual figures are gated on `save_ts_plots` like every other SVG; without the
    # gate a default run with Plots loaded produced figures nobody asked for, and a run
    # without Plots left `Figures_Duals/` empty.
    #
    # Every corridor family goes in, whichever `bound_style_δ` / `bound_style_Δω` built it:
    # the families are mutually exclusive per run, so at most one δ entry and one Δω entry
    # carry data and the rest are skipped inside the plotter. Feeding it only the δ-COI pair
    # meant a machine-referenced or `:abs` run produced no dual figure at all.
    if save_ts_plots
        mkpath(pf_ts_figures_duals)
        # Three mutually exclusive δ families now: COI-framed, reference-machine-framed and
        # pairwise. Pick whichever one this run built, and name the frame in the title —
        # these are different quantities and the figures must not be confusable.
        δ_lower_data = something(dual_δ_COI_lower, dual_δ_ref_lower, dual_δ_pair_lower, Some(nothing))
        δ_upper_data = something(dual_δ_COI_upper, dual_δ_ref_upper, dual_δ_pair_upper, Some(nothing))
        δ_frame = if dual_δ_pair_upper !== nothing || dual_δ_pair_lower !== nothing
            "pairwise"
        elseif δ_ref_gen === nothing
            "ref. COI"
        else
            "ref. G$(δ_ref_gen)"
        end
        corridor_dual_families = [
            (title = "δ corridor ($δ_frame) - lower", side = "lower", tag = "",        ylabel = "€/rad",
             data = δ_lower_data),
            (title = "δ corridor ($δ_frame) - upper", side = "upper", tag = "",        ylabel = "€/rad",
             data = δ_upper_data),
            (title = "Δω corridor - lower",           side = "lower", tag = "Domega ", ylabel = "€/p.u.",
             data = dual_Δω_COI_lower === nothing ? dual_Δω_abs_lower : dual_Δω_COI_lower),
            (title = "Δω corridor - upper",           side = "upper", tag = "Domega ", ylabel = "€/p.u.",
             data = dual_Δω_COI_upper === nothing ? dual_Δω_abs_upper : dual_Δω_COI_upper),
        ]
        invoke_save_dual_ts_constraint_svgs!(
            dyn_parameters_dict[:time][:t_window_total],
            corridor_dual_families,
            pf_ts_figures_duals,
        )
    end

end

"""
Column name for one key of a per-key dual trajectory.

Single-id families ("Gen", "Bus") keep their historical `<prefix>_<id>` header. The pairwise
angle corridor is keyed by a `(i, j)` tuple, and string-interpolating that raw would produce
`Gen_(1, 2)` — a header that reads as a generator and does not survive a round trip through
CSV readers. Dispatching on the key type keeps both cases in one place, since the CSV and
XLSX writers below must agree on the spelling.
"""
_dual_key_label(key::Tuple{Int, Int}, ::AbstractString) = "Pair_G$(key[1])_G$(key[2])"
_dual_key_label(key, prefix::AbstractString) = "$(prefix)_$(key)"

# A helper function to convert OrderedDict{key, Vector{Float64}} to a CSV
# `id_prefix` names what the keys are: "Gen" (default) for per-generator families,
# "Bus" for the nodal ones (PER_BUS_TIME_MERGE). Pair-keyed families ignore it.
function save_ordered_dict_to_csv(data_dict, filename; id_prefix::String = "Gen")
    # Convert the dictionary to a DataFrame
    # Column names will be "Gen_1", "Gen_2", … or "Pair_G1_G2" for the pairwise corridor.
    df = DataFrame()
    for (gen_id, values) in data_dict
        df[!, Symbol(_dual_key_label(gen_id, id_prefix))] = values
    end
    CSV.write(filename, df; delim = ';')
end

# Write registry-extracted duals to `Transient_Stability/OPF_Duals_Results.xlsx`.
# Accepts any `*_xlsx` kwargs from `duals_to_xlsx_kwargs` (vectors or OrderedDicts)
# so new registry entries (GFM, V_ref bounds, …) do not require signature edits.
function Save_Duals_2_Excel_tsred(path_names; kwargs...)
    file_path = joinpath(path_names[:pf_TS], "OPF_Duals_Results.xlsx")

    function to_df(data_dict)
        if isnothing(data_dict) || isempty(data_dict)
            return nothing
        end
        df = DataFrame()
        for (gen_id, values) in data_dict
            df[!, Symbol(_dual_key_label(gen_id, "Gen"))] = values
        end
        return df
    end

    function sheet_name(kw::Symbol)::String
        s = String(kw)
        endswith(s, "_xlsx") && (s = chop(s; tail=5))
        startswith(s, "dual_") && (s = chop(s; head=5, tail=0))
        return s * "_Duals"
    end

    sheets = Pair{String, DataFrame}[]
    for (kw, data) in kwargs
        data === nothing && continue
        if data isa AbstractDict
            df = to_df(data)
            df === nothing && continue
            push!(sheets, sheet_name(kw) => df)
        elseif data isa AbstractVector
            isempty(data) && continue
            push!(sheets, sheet_name(kw) => DataFrame(Value = data))
        else
            @warn "Skipping Excel dual sheet for $kw: unsupported type $(typeof(data))"
        end
    end

    if !isempty(sheets)
        XLSX.writetable(file_path, sheets..., overwrite = true)
        println("Duals successfully saved to: $file_path")
    else
        @warn "No dual data was provided. Excel file not created."
    end
    return nothing
end

# ===================================================================================
#                   MANAGE RESULTS FROM THE DYNAMIC SIMULATION
# ===================================================================================
# Function to save results across the whole time window (fault and post-fault)
# Accelerating power of each generator relative to the COI, expressed in pu of
# base_MVA. Keyed by generator index — which MUST match DGEN_DYN row indices —
# NOT by enumeration position: under a GL generator trip `active_gen` is
# non-contiguous, so a position-based H lookup would weight the wrong machine
# (review finding C2). Unit-tested in test/runtests_fast_unit.jl.
function _pacc_relative_to_coi(
    Pacc::OrderedDict{Int, Vector{Float64}},
    PaccCOI_total::Vector{Float64},
    DGEN_DYN::DataFrame,
    H_total::Union{Nothing, Float64},
    base_MVA::Float64,
)
    Pacc_COI = OrderedDict{Int, Vector{Float64}}()
    for (i, info) in Pacc
        # All-GFM fleet: every Hi is zero, so the COI term drops out and the
        # COI-relative accelerating power collapses to the absolute one.
        if H_total === nothing
            Pacc_COI[i] = info ./ base_MVA
            continue
        end
        Hi = _inertia_H(DGEN_DYN, Int(i))
        Pacc_COI[i] = (info .- ((Hi .* PaccCOI_total) ./ H_total)) ./ base_MVA
    end
    return Pacc_COI
end

"""
    _δ_relative_to_ref(δt, δ_ref_gen) -> OrderedDict | nothing

Rotor angles relative to the corridor's reference machine, in **radians**, keyed like `δt`.

Returns `nothing` on the COI-referenced styles (`δ_ref_gen === nothing`), where no such frame
exists. The reference machine's own entry is identically zero — it is `δ_ref − δ_ref` — which is
exactly what `bound_style_δ ∈ (:highest_H, :ref_gen)` bounds the other machines against, and why
that machine carries no corridor row of its own.
"""
function _δ_relative_to_ref(
    δt::OrderedDict{Int, Vector{Float64}},
    δ_ref_gen::Union{Nothing, Int},
)
    (δ_ref_gen === nothing || !haskey(δt, δ_ref_gen)) && return nothing
    δ_ref_trace = δt[δ_ref_gen]
    return OrderedDict{Int, Vector{Float64}}(k => v .- δ_ref_trace for (k, v) in δt)
end

"""
    _compute_pacc_trajectories(Pet, Pmt, δ_COIt, DGEN_DYN, active_gen, base_MVA)

Accelerating power (MW), COI-relative Pacc (pu then scaled to MW), and potential
energy from ∫Pacc_COI dδ. `Pmt` must be keyed like `Pet` with matching time lengths.
"""
function _compute_pacc_trajectories(
    Pet::OrderedDict{Int, Vector{Float64}},
    Pmt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    DGEN_DYN::DataFrame,
    active_gen,
    base_MVA::Float64,
)
    H_total = _H_total(DGEN_DYN, active_gen)
    Pacc = OrderedDict{Int, Vector{Float64}}()
    PaccCOI = OrderedDict{Int, Vector{Float64}}()
    first_gen = true
    for (i, Pe_vals) in Pet
        Pm_vals = Pmt[i]
        length(Pm_vals) == length(Pe_vals) ||
            throw(ArgumentError("Pm and Pe length mismatch for generator $i"))
        pacc_i = Pm_vals .- Pe_vals
        Pacc[i] = pacc_i
        if first_gen
            PaccCOI[1] = copy(pacc_i)
            first_gen = false
        else
            PaccCOI[1] .+= pacc_i
        end
    end
    Pacc_COI = _pacc_relative_to_coi(Pacc, PaccCOI[1], DGEN_DYN, H_total, base_MVA)

    Vpe = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Pacc_COI)
        haskey(δ_COIt, k) || continue
        # Classical PE is SG-only (GFM has H=0 and is omitted from Vke).
        _inertia_H(DGEN_DYN, Int(k)) == 0.0 && continue
        integral = Float64[]
        for (idx, _) in enumerate(Pacc_COI[k])
            if idx == 1
                push!(integral, 0.0)
            else
                push!(integral, -trapz(δ_COIt[k][1:idx], Pacc_COI[k][1:idx]))
            end
        end
        Vpe[k] = integral
    end
    for (i, info) in Pacc_COI
        Pacc_COI[i] = info .* base_MVA
    end
    return Pacc, PaccCOI, Pacc_COI, Vpe
end

# Function to save results across the whole time window for the dynamic model (fault and post-fault).
# Model-agnostic: handles both the linear (TSC-DCOPF) and nonlinear (TSC-ACOPF) reduced models,
# since both populate the same dyn_model_dict structure. (Renamed from
# Save_Results_Dynamic_Model_tsred_linear, which misleadingly implied linear-only — audit SF-2.)
function Save_Results_Dynamic_Model(model::Model,
    path_names::OrderedDict{Symbol, String},
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
    P_m::OrderedDict{Int, JuMP.VariableRef},
    DGEN_DYN::DataFrame;
    save_ts_plots::Bool=false,
    save_ts_debug_csv::Bool=false,
    )

    active_gen = dyn_model_dict[:active_gen]     # Indices of active generators

    pf_ts = path_names[:pf_TS]
    mkpath(pf_ts)

    coupling_meta = _prefault_coupling_metadata(dyn_model_dict)
    coupling_table = _prefault_coupling_table(dyn_model_dict, base_MVA)

    # ====================================
    # Save short results in TXT
    # ====================================
    open(joinpath(pf_ts, "dynamic_model_short_results.txt"), "w") do io
        _write_prefault_short_results_txt!(io, coupling_table, coupling_meta)
    end
    Save_Prefault_Coupling_CSV!(coupling_table, path_names)

    # ====================================
    #                P_e
    # ====================================
    Pe_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in dyn_model_dict[:vars][:Pe_tf])  # Values of Pe during the fault

    if haskey(dyn_model_dict[:vars], :Pe_tpf)
        Pe_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in dyn_model_dict[:vars][:Pe_tpf]) # Values of Pe in the post-fault

        Pet = OrderedDict{Int, Vector{Float64}}()
        for k in keys(Pe_tf_values)
            Pet[k] = vcat(Pe_tf_values[k], Pe_tpf_values[k])
        end
    else
        Pet = Pe_tf_values
    end

    Qet = _extract_gen_time_trajectory(
        dyn_model_dict, :Qe_tf, :Qe_tpf; scale=base_MVA)

    # ====================================
    #                 δ
    # ====================================
    δ_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in dyn_model_dict[:vars][:δ_tf])   # Values of δ during the fault
    # δ_COI is a variable when the corridor references it and an expression otherwise;
    # `JuMP.value` reads both, so `angle_rel_COI.csv` is produced either way.
    δCOI_tf_values  = [JuMP.value(v) for (i, v) in _coi_series(dyn_model_dict, :δCOI_tf)]  # Values of δ_COI during the fault

    if haskey(dyn_model_dict[:vars], :δ_tpf)
        δ_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in dyn_model_dict[:vars][:δ_tpf]) # Values of Pe in the post-fault
        δCOI_tpf_values = [JuMP.value(v) for (i, v) in _coi_series(dyn_model_dict, :δCOI_tpf)] # Values of δ_COI in the post-fault

        δt = OrderedDict{Int, Vector{Float64}}()
        δ_COIt = OrderedDict{Int, Vector{Float64}}()  # δ individual in relation to COI across the whole simulation
        δCOIt = vcat(δCOI_tf_values, δCOI_tpf_values) # δ do COI across the whole simulation

        for k in keys(δ_tf_values)
            δt[k] = vcat(δ_tf_values[k], δ_tpf_values[k])
            δ_COIt[k] = δt[k] .- δCOIt
        end

    else
        δt = δ_tf_values
        δ_COIt = OrderedDict{Int, Vector{Float64}}()  # δ individual in relation to COI across the whole simulation
        δCOIt = δCOI_tf_values # δ do COI across the whole simulation

        for k in keys(δ_tf_values)
            δ_COIt[k] = δt[k] .- δCOIt
        end
    end

    # ====================================
    #                 Δω
    # ====================================
    Δω_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in dyn_model_dict[:vars][:Δω_tf])  # Values of Δω during the fault

    if haskey(dyn_model_dict[:vars], :Δω_tpf)
        Δω_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in dyn_model_dict[:vars][:Δω_tpf]) # Values of Pe in the post-fault

        Δωt = OrderedDict{Int, Vector{Float64}}()

        for k in keys(Δω_tf_values)
            Δωt[k] = vcat(Δω_tf_values[k], Δω_tpf_values[k])
        end

    else
        Δωt = Δω_tf_values
    end

    ΔωCOIt, Δω_COIt = _extract_ΔωCOI_trajectories(
        dyn_model_dict, Δωt, DGEN_DYN, active_gen, length(dyn_parameters_dict[:time][:t_window_total]))


    # ====================================
    #              RoCoF
    # ====================================
    time_RoCoF = []
    RoCoF = OrderedDict{Int, Vector{Float64}}()  # RoCoF individual across the whole simulation
    for (i, info) in Δωt
        time_RoCoF, RoCoF[i] = Calculate_Derivative_CFD(dyn_parameters_dict[:time][:t_window_total], dyn_parameters_dict[:common][:f_syn] .* (1.0 .+ info))
    end

    time_RoCoF, RoCoFCOI = Calculate_Derivative_CFD(dyn_parameters_dict[:time][:t_window_total], dyn_parameters_dict[:common][:f_syn] .* (1.0 .+ ΔωCOIt))

    # ====================================
    #          Kinetic Energy
    # ====================================
    Vke = OrderedDict{Int, Vector{Float64}}()  # Vke individual across the whole simulation
    for (i, info) in Δω_COIt
        Hi = _inertia_H(DGEN_DYN, Int(i))
        Hi == 0.0 && continue  # GFM has no classical kinetic energy term
        Vke[i] = Hi .* dyn_parameters_dict[:common][:ω_syn] .* (info .^2)
    end

    # ====================================
    #        Accelerating Power
    # ====================================
    t_window_total = dyn_parameters_dict[:time][:t_window_total]
    Pmt = _mechanical_power_trajectories(
        dyn_model_dict, Pet, P_m, base_MVA, length(t_window_total))
    Pacc, PaccCOI, Pacc_COI, Vpe = _compute_pacc_trajectories(
        Pet, Pmt, δ_COIt, DGEN_DYN, active_gen, base_MVA)

    # ====================================
    #           SAVE FIGURES
    # ====================================
    # Resolved δ corridor reference (`:highest_H` / `:ref_gen` only). Both the figures and the
    # CSV writer below take it, and both derive the reference-relative angles from the same
    # helper so `angle_rel_ref.csv` and "Delta (ref G<k>) vs time.svg" cannot drift apart.
    δ_ref_gen = get(dyn_model_dict[:meta], :δ_ref_gen_resolved, nothing)
    if save_ts_plots
        Save_Dynamic_Results_Plots_tsred(
            t_window_total, Pmt, Pet, Qet, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt,
            time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe,
            base_MVA, dyn_parameters_dict[:common][:f_syn], dyn_parameters_dict[:common][:δ_tol],
            path_names;
            δ_ref_gen = δ_ref_gen,
            δ_reft = _δ_relative_to_ref(δt, δ_ref_gen))
    end
    t_clear_fault = get(dyn_parameters_dict[:time], :t_clear_fault, nothing)

    # ====================================
    #           SAVE CSV FILES
    # ====================================
    Save_Dynamic_Results_CSV_tsred(
        t_window_total, Pmt, Pet, Qet, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt,
        time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe,
        base_MVA, dyn_parameters_dict[:common][:f_syn], path_names;
        has_governor = haskey(dyn_model_dict[:vars], :Pm_tf),
        governed_gens = haskey(dyn_model_dict[:vars], :Pm_tf) ?
            collect(Int, keys(dyn_model_dict[:vars][:Pm_tf])) : Int[],
        δ_ref_gen = δ_ref_gen)

    # Pre-fault internal EMF magnitude. `:E` exists on the classical Kron path as well as
    # on classical FULL_BUS, but the writer used to sit inside `Save_FullBus_Network_Results!`
    # — so a Kron run solved for E and then exported nothing. (DQ_4TH has no `:E`; it uses
    # Ed/Eq, written by `Save_Dq_Generator_Results!`.)
    if haskey(dyn_model_dict[:vars], :E)
        E_container = dyn_model_dict[:vars][:E]
        CSV.write(joinpath(path_names[:pf_TS_CSV], "generator_internal_voltage.csv"),
            DataFrame(generator = collect(keys(E_container)),
                      E_pu = [JuMP.value(v) for (_, v) in E_container]); delim = ';')
    end

    # Governor trajectories (no-op unless include_governor produced Pm_tf).
    Save_Governor_Trajectories_CSV(dyn_model_dict, dyn_parameters_dict, base_MVA, path_names;
        t_clear_fault=t_clear_fault, save_ts_plots=save_ts_plots)

    if get(get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}()), :network_form, "") == "FULL_BUS"
        Save_FullBus_Network_Results!(
            dyn_model_dict, dyn_parameters_dict, base_MVA, path_names;
            t_clear_fault=t_clear_fault, save_ts_plots=save_ts_plots)
        if _is_dq_4th_model(dyn_model_dict)
            Save_Dq_Generator_Results!(
                dyn_model_dict, dyn_parameters_dict, base_MVA, path_names;
                t_clear_fault=t_clear_fault, save_ts_plots=save_ts_plots)
            Save_Gfm_Converter_Results!(
                dyn_model_dict, dyn_parameters_dict, path_names;
                t_clear_fault=t_clear_fault, save_ts_plots=save_ts_plots)
        end
    end

    if save_ts_debug_csv
        Save_TS_Debug_CSV!(dyn_model_dict, dyn_parameters_dict, DGEN_DYN, path_names)
    end

end

# Function to save the transient stability results in CSV files
function Save_Dynamic_Results_CSV_tsred(t_window_total::Vector{Float64},
    Pmt::OrderedDict{Int, Vector{Float64}},
    Pet::OrderedDict{Int, Vector{Float64}},
    Qet::Union{Nothing, OrderedDict{Int, Vector{Float64}}},
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    path_names::OrderedDict{Symbol, String};
    has_governor::Bool = false,
    governed_gens::Vector{Int} = Int[],
    δ_ref_gen::Union{Nothing, Int} = nothing,
    )

    pf_ts_csv = path_names[:pf_TS_CSV]
    
    # ========================================================
    # Names of the variables
    # ========================================================

    gen_names = ["G$(i)" for (i, inner) in Pet]

    delta_names = ["G$(i)" for (i, inner) in δt]
    push!(delta_names, "COI")

    omega_names = ["G$(i)" for (i, inner) in Δωt]
    push!(omega_names, "COI")

    # ========================================================
    # Creating the matrices
    # ========================================================
    # Electrical Power
    Pe_matrix = zeros(Float64, length(t_window_total), length(Pet))
    aux_count = 0
    for (gen_id, values) in Pet
        aux_count += 1
        Pe_matrix[:,aux_count] = values
    end
    df_Pe = DataFrame(hcat(t_window_total, Pe_matrix), vcat("t", gen_names))

    Pm_matrix = zeros(Float64, length(t_window_total), length(Pmt))
    aux_count = 0
    for (gen_id, values) in Pmt
        aux_count += 1
        Pm_matrix[:, aux_count] = values
    end
    pm_gen_names = ["G$(i)" for (i, _) in Pmt]
    df_Pm = DataFrame(hcat(t_window_total, Pm_matrix), vcat("t", pm_gen_names))

    df_Qe = nothing
    if Qet !== nothing
        Qe_matrix = zeros(Float64, length(t_window_total), length(Qet))
        aux_count = 0
        for (gen_id, values) in Qet
            aux_count += 1
            Qe_matrix[:, aux_count] = values
        end
        q_gen_names = ["G$(i)" for (i, _) in Qet]
        df_Qe = DataFrame(hcat(t_window_total, Qe_matrix), vcat("t", q_gen_names))
    end

    # Accelerating Power
    Pacc_matrix = zeros(Float64, length(t_window_total), length(Pacc)+1)
    aux_count = 0
    for (gen_id, values) in Pacc
        aux_count += 1
        Pacc_matrix[:,aux_count] = values
    end
    Pacc_matrix[:,end] = PaccCOI[1]
    df_Pacc = DataFrame(hcat(t_window_total, Pacc_matrix), vcat("t", delta_names))

    # Accelerating Power in relation to COI
    Pacc_COI_matrix = zeros(Float64, length(t_window_total), length(Pacc_COI))
    aux_count = 0
    for (gen_id, values) in Pacc_COI
        aux_count += 1
        Pacc_COI_matrix[:,aux_count] = values
    end
    df_Pacc_COI = DataFrame(hcat(t_window_total, Pacc_COI_matrix), vcat("t", gen_names))

    # Absolute values of δ
    δ_matrix = zeros(Float64, length(t_window_total), length(δt)+1)
    aux_count = 0
    for (gen_id, values) in δt
        aux_count += 1
        δ_matrix[:,aux_count] = rad2deg.(values)
    end
    δ_matrix[:,end] = rad2deg.(δCOIt)
    df_δ = DataFrame(hcat(t_window_total, δ_matrix), vcat("t", delta_names))

    # δ values in relation to the COI
    δ_COI_matrix = zeros(Float64, length(t_window_total), length(δ_COIt))
    aux_count = 0
    for (gen_id, values) in δ_COIt
        aux_count += 1
        δ_COI_matrix[:,aux_count] = rad2deg.(values)
    end
    df_δ_COI = DataFrame(hcat(t_window_total, δ_COI_matrix), vcat("t", gen_names))

    # δ values relative to the reference machine, when the run used one. This — not the
    # COI-relative trajectory — is the quantity the :highest_H / :ref_gen corridors bound,
    # so it is what the corridor plots need. `nothing` on the COI-referenced styles.
    df_δ_ref = nothing
    δ_reft = _δ_relative_to_ref(δt, δ_ref_gen)   # same helper the figures use
    if δ_reft !== nothing
        δ_ref_matrix = zeros(Float64, length(t_window_total), length(δ_reft))
        aux_count = 0
        for (gen_id, values) in δ_reft
            aux_count += 1
            δ_ref_matrix[:, aux_count] = rad2deg.(values)
        end
        df_δ_ref = DataFrame(hcat(t_window_total, δ_ref_matrix),
            vcat("t", ["G$(i)" for (i, _) in δ_reft]))
    end


    # Absolute values of Δω
    Δω_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0
    for (gen_id, values) in Δωt
        aux_count += 1
        Δω_matrix[:,aux_count] = values
    end
    Δω_matrix[:,end] = ΔωCOIt
    df_Δω = DataFrame(hcat(t_window_total, Δω_matrix), vcat("t", omega_names))

    # Δω values in relation to the COI
    Δω_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0
    for (gen_id, values) in Δω_COIt
        aux_count += 1
        Δω_COI_matrix[:,aux_count] = values
    end
    df_Δω_COI = DataFrame(hcat(t_window_total, Δω_COI_matrix), vcat("t", gen_names))

    # Frequency values
    f_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0
    for (gen_id, values) in Δωt
        aux_count += 1
        f_matrix[:,aux_count] = f_syn .* (1.0 .+ values)
    end
    f_matrix[:,end] = f_syn .* (1.0 .+ ΔωCOIt)
    df_f = DataFrame(hcat(t_window_total, f_matrix), vcat("t", omega_names))

    # Frequency values in relation to the COI
    f_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0
    for (gen_id, values) in Δω_COIt
        aux_count += 1
        f_COI_matrix[:,aux_count] = f_syn .* (1.0 .+ values)
    end
    df_f_COI = DataFrame(hcat(t_window_total, f_COI_matrix), vcat("t", gen_names))
    
    # RoCoF values
    RoCoF_matrix = zeros(Float64, length(time_RoCoF), length(RoCoF)+1)
    aux_count = 0
    for (gen_id, values) in RoCoF
        aux_count += 1
        RoCoF_matrix[:,aux_count] = values
    end
    RoCoF_matrix[:,end] = RoCoFCOI
    df_RoCoF = DataFrame(hcat(time_RoCoF, RoCoF_matrix), vcat("t", omega_names))

    # Vke values (SG inertia only — GFM keys omitted)
    Vke_matrix = zeros(Float64, length(t_window_total), length(Vke))
    aux_count = 0
    for (gen_id, values) in Vke
        aux_count += 1
        Vke_matrix[:,aux_count] = values
    end
    vke_names = ["G$(i)" for (i, _) in Vke]
    df_Vke = DataFrame(hcat(t_window_total, Vke_matrix), vcat("t", vke_names))

    # Vpe values (may omit GFM when δ_COI is SG-only)
    Vpe_matrix = zeros(Float64, length(t_window_total), length(Vpe))
    aux_count = 0
    for (gen_id, values) in Vpe
        aux_count += 1
        Vpe_matrix[:,aux_count] = values
    end
    vpe_names = ["G$(i)" for (i, _) in Vpe]
    df_Vpe = DataFrame(hcat(t_window_total, Vpe_matrix), vcat("t", vpe_names))

    # Total transient energy Ve = Vke + Vpe, over the generators that have both terms
    # (SG only: GFM units carry no classical kinetic energy). This is the quantity the
    # equal-area / energy-function reading of the swing curves rests on, and it used to
    # exist only inside "Total Energy vs time.svg" — computed in the plot and discarded
    # with it, so a run without `save_ts_plots` produced no record of it at all.
    ve_ids = [i for i in keys(Vpe) if haskey(Vke, i)]
    Ve_matrix = zeros(Float64, length(t_window_total), length(ve_ids))
    for (j, gen_id) in enumerate(ve_ids)
        Ve_matrix[:, j] = Vke[gen_id] .+ Vpe[gen_id]
    end
    df_Ve = DataFrame(hcat(t_window_total, Ve_matrix),
        vcat("t", ["G$(i)" for i in ve_ids]))

    # ========================================================
    # Saving the CSVs
    # ========================================================
    CSV.write(joinpath(pf_ts_csv, "electrical_power.csv"), df_Pe; delim=';')
    # With a governor, P_m(t) is the governor mechanical-power state, so for an SG-only
    # fleet this file would duplicate `governor_P_mech.csv` byte for byte and is skipped.
    # It is *not* a duplicate once GFM units are present: `governor_P_mech.csv` covers the
    # governed SGs only, while `Pmt` also carries the converters' constant P_set — which
    # otherwise reached no file at all. So suppress only when every column is governed.
    gov_covers_all = has_governor && length(governed_gens) == length(Pmt) &&
        all(g -> g in governed_gens, keys(Pmt))
    gov_covers_all || CSV.write(joinpath(pf_ts_csv, "mechanical_power.csv"), df_Pm; delim=';')
    df_Qe !== nothing && CSV.write(joinpath(pf_ts_csv, "electrical_reactive_power.csv"), df_Qe; delim=';')
    CSV.write(joinpath(pf_ts_csv, "accelerating_power.csv"), df_Pacc; delim=';')
    CSV.write(joinpath(pf_ts_csv, "accelerating_power_COI.csv"), df_Pacc_COI; delim=';')
    CSV.write(joinpath(pf_ts_csv, "angle_abs.csv"), df_δ; delim=';')
    CSV.write(joinpath(pf_ts_csv, "angle_rel_COI.csv"), df_δ_COI; delim=';')
    df_δ_ref === nothing ||
        CSV.write(joinpath(pf_ts_csv, "angle_rel_ref.csv"), df_δ_ref; delim=';')
    CSV.write(joinpath(pf_ts_csv, "speed_dev.csv"), df_Δω; delim=';')
    CSV.write(joinpath(pf_ts_csv, "speed_dev_rel_COI.csv"), df_Δω_COI; delim=';')
    CSV.write(joinpath(pf_ts_csv, "frequency.csv"), df_f; delim=';')
    CSV.write(joinpath(pf_ts_csv, "frequency_rel_COI.csv"), df_f_COI; delim=';')
    CSV.write(joinpath(pf_ts_csv, "RoCoF.csv"), df_RoCoF; delim=';')
    CSV.write(joinpath(pf_ts_csv, "kinetic_energy.csv"), df_Vke; delim=';')
    CSV.write(joinpath(pf_ts_csv, "potential_energy.csv"), df_Vpe; delim=';')
    CSV.write(joinpath(pf_ts_csv, "total_energy.csv"), df_Ve; delim=';')

    println("Results of the dynamic model successfully saved as CSV files in: ", pf_ts_csv)

end

# ===================================================================================
#                   FULL_BUS network trajectories (bus P/Q/V/θ)
# ===================================================================================

"""Merge fault-on and post-fault per-generator time series from `vars` or `expressions`."""
function _extract_gen_time_trajectory(
    dyn_model_dict::OrderedDict{Symbol, Any},
    key_tf::Symbol,
    key_tpf::Symbol;
    scale::Float64=1.0,
)
    for source in (:vars, :expressions)
        store = get(dyn_model_dict, source, nothing)
        store === nothing && continue
        haskey(store, key_tf) || continue
        # A family can exist as a container and still hold no rows — the SG-only dq/AVR/
        # governor families on a converter-only fleet are exactly that. Report it the same
        # way as a missing key so every caller's `=== nothing` guard covers both.
        isempty(store[key_tf]) && return nothing
        tf_vals = OrderedDict(
            gen => scale .* [JuMP.value(v) for (_, v) in inner]
            for (gen, inner) in store[key_tf])
        if haskey(store, key_tpf)
            tpf_vals = OrderedDict(
                gen => scale .* [JuMP.value(v) for (_, v) in inner]
                for (gen, inner) in store[key_tpf])
            return OrderedDict(k => vcat(tf_vals[k], tpf_vals[k]) for k in keys(tf_vals))
        end
        return tf_vals
    end
    return nothing
end

"""
    _mechanical_power_trajectories(dyn_model_dict, Pet, P_m, base_MVA, n_time)

Per-generator mechanical power [MW] over the merged fault+post-fault window.
Governor runs use `:Pm_tf`/`:Pm_tpf`; otherwise constant pre-fault `P_m` per generator.
"""
function _mechanical_power_trajectories(
    dyn_model_dict::OrderedDict{Symbol, Any},
    Pet::OrderedDict{Int, Vector{Float64}},
    P_m::OrderedDict{Int, JuMP.VariableRef},
    base_MVA::Float64,
    n_time::Int,
)::OrderedDict{Int, Vector{Float64}}
    vars = get(dyn_model_dict, :vars, nothing)
    if vars !== nothing && haskey(vars, :Pm_tf)
        traj = _extract_gen_time_trajectory(dyn_model_dict, :Pm_tf, :Pm_tpf; scale=base_MVA)
        if traj !== nothing
            # Governor trajectories cover SG only; GFM uses constant P_set (= P_m).
            Pmt = OrderedDict{Int, Vector{Float64}}()
            for gen_id in keys(Pet)
                if haskey(traj, gen_id)
                    Pmt[gen_id] = traj[gen_id]
                else
                    Pm_mw = JuMP.value(P_m[gen_id]) * base_MVA
                    Pmt[gen_id] = fill(Pm_mw, n_time)
                end
            end
            return Pmt
        end
    end
    Pmt = OrderedDict{Int, Vector{Float64}}()
    for gen_id in keys(Pet)
        Pm_mw = JuMP.value(P_m[gen_id]) * base_MVA
        Pmt[gen_id] = fill(Pm_mw, n_time)
    end
    return Pmt
end

"""
    Save_Governor_Trajectories_CSV(dyn_model_dict, dyn_parameters_dict, base_MVA, path_names)

Write the turbine-governor trajectories to CSV (and SVG when `save_ts_plots`) if a
governor is present in the run; no-op otherwise. Reuses `_extract_gen_time_trajectory`
for the fault + post-fault merge.

Produces `governor_P_mech.csv`, `governor_P_valve.csv`, and `governor_P_valve_raw.csv`
(columns `t` + `G<id>`, powers in MW) plus `governor_P_ref.csv` (the constant
per-generator set-point, in p.u.). `P_valve_raw` is the integrator state ahead of the
valve limiter, so `P_valve_raw − P_valve` is the active clamp under `GOV_SMOOTH`.
"""
function Save_Governor_Trajectories_CSV(
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
    path_names::OrderedDict{Symbol, String};
    t_clear_fault=nothing,
    save_ts_plots::Bool=false,
)
    haskey(dyn_model_dict[:vars], :Pm_tf) || return nothing  # no governor in this run

    pf_ts_csv = path_names[:pf_TS_CSV]
    mkpath(pf_ts_csv)
    pf_figures = path_names[:pf_TS_figures]
    save_ts_plots && mkpath(pf_figures)
    t = dyn_parameters_dict[:time][:t_window_total]

    for (label, key_tf, key_tpf, title, ylab, svg_name) in (
        ("P_mech", :Pm_tf, :Pm_tpf,
            "Governor Mechanical Power", "P_m (MW)", "governor_P_mech.svg"),
        ("P_valve", :Pv_tf, :Pv_tpf,
            "Governor Valve Output", "P_valve (MW)", "governor_P_valve.svg"),
        # Raw integrator state ahead of the valve limiter. Identical to P_valve under
        # GOV_NO_LIMIT / GOV_HARD_BOUND; under GOV_SMOOTH the difference is the clamp.
        ("P_valve_raw", :Pv_raw_tf, :Pv_raw_tpf,
            "Governor Valve (unsaturated)", "P_valve_raw (MW)", "governor_P_valve_raw.svg"),
    )
        traj = _extract_gen_time_trajectory(dyn_model_dict, key_tf, key_tpf; scale=base_MVA)
        traj === nothing && continue
        gen_names = ["G$(g)" for g in keys(traj)]
        mat = zeros(Float64, length(t), length(traj))
        for (j, (_, vals)) in enumerate(traj)
            mat[:, j] = vals
        end
        CSV.write(joinpath(pf_ts_csv, "governor_$(label).csv"),
            DataFrame(hcat(t, mat), vcat("t", gen_names)); delim = ';')
        if save_ts_plots
            _save_fullbus_svg(t, traj, title, ylab, pf_figures, svg_name;
                t_clear_fault=t_clear_fault, label_prefix="G")
        end
    end

    if haskey(dyn_model_dict[:vars], :P_ref)
        P_ref = dyn_model_dict[:vars][:P_ref]
        CSV.write(joinpath(pf_ts_csv, "governor_P_ref.csv"),
            DataFrame(gen = ["G$(g)" for g in keys(P_ref)],
                      P_ref_pu = [JuMP.value(v) for (_, v) in P_ref]); delim = ';')
    end

    println("Governor trajectories saved to: ", pf_ts_csv)
    return nothing
end

"""
    Save_Dq_Generator_Results!(dyn_model_dict, dyn_parameters_dict, base_MVA, path_names)

Write dq-machine trajectories (`E_fd`, `Ed/Eq/Id/Iq/Te`) to CSV and optional SVG plots.
Called from `Save_Results_Dynamic_Model` when `meta[:gen_order] == "DQ_4TH"`.
"""
function Save_Dq_Generator_Results!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
    path_names::OrderedDict{Symbol, String};
    t_clear_fault=nothing,
    save_ts_plots::Bool=false,
)
    pf_ts_csv = path_names[:pf_TS_CSV]
    mkpath(pf_ts_csv)
    t_window = dyn_parameters_dict[:time][:t_window_total]
    pf_figures = path_names[:pf_TS_figures]
    save_ts_plots && mkpath(pf_figures)

    if haskey(dyn_model_dict[:vars], :V_ref)
        CSV.write(joinpath(pf_ts_csv, "avr_V_ref.csv"),
            DataFrame(gen = ["G$(g)" for g in keys(dyn_model_dict[:vars][:V_ref])],
                      V_ref_pu = [JuMP.value(v) for (_, v) in dyn_model_dict[:vars][:V_ref]]);
            delim = ';')
    end

    if haskey(dyn_model_dict[:vars], :E_fd)
        CSV.write(joinpath(pf_ts_csv, "dq_E_fd.csv"),
            DataFrame(gen = ["G$(g)" for g in keys(dyn_model_dict[:vars][:E_fd])],
                      E_fd_pu = [JuMP.value(v) for (_, v) in dyn_model_dict[:vars][:E_fd]]);
            delim = ';')
    end

    prefault_rows = Vector{NamedTuple}()
    for (label, sym) in (
        ("Ed_pu", :Ed), ("Eq_pu", :Eq), ("Id_pu", :Id), ("Iq_pu", :Iq),
    )
        haskey(dyn_model_dict[:vars], sym) || continue
        for (g, v) in dyn_model_dict[:vars][sym]
            push!(prefault_rows, (gen = "G$g", variable = label, value = JuMP.value(v)))
        end
    end
    # Terminal dq voltages at t = 0 are expressions, not variables, so they live in
    # `:expressions` — same rows, same file, so the pre-fault algebraic state is complete.
    exprs_store = get(dyn_model_dict, :expressions, OrderedDict{Symbol, Any}())
    for (label, sym) in (("Vd_pu", :Vd_init), ("Vq_pu", :Vq_init))
        haskey(exprs_store, sym) || continue
        for (g, e) in exprs_store[sym]
            push!(prefault_rows, (gen = "G$g", variable = label, value = JuMP.value(e)))
        end
    end
    if !isempty(prefault_rows)
        CSV.write(joinpath(pf_ts_csv, "dq_prefault_algebraic.csv"), DataFrame(prefault_rows); delim = ';')
    end

    dq_series = (
        ("dq_Ed_pu", :Ed_tf, :Ed_tpf, 1.0, "Subtransient EMF Ed", "Ed (p.u.)", "dq_Ed.svg"),
        ("dq_Eq_pu", :Eq_tf, :Eq_tpf, 1.0, "Subtransient EMF Eq", "Eq (p.u.)", "dq_Eq.svg"),
        ("dq_Id_pu", :Id_tf, :Id_tpf, 1.0, "dq Current Id", "Id (p.u.)", "dq_Id.svg"),
        ("dq_Iq_pu", :Iq_tf, :Iq_tpf, 1.0, "dq Current Iq", "Iq (p.u.)", "dq_Iq.svg"),
        ("dq_Te_MW", :Te_tf, :Te_tpf, base_MVA, "Electrical Torque Te", "T_e (MW)", "dq_Te.svg"),
        ("dq_E_fd_pu", :E_fd_tf, :E_fd_tpf, 1.0, "Field Voltage E_fd", "E_fd (p.u.)", "dq_E_fd_traj.svg"),
        # Terminal dq voltages: JuMP expressions, picked up by `_extract_gen_time_trajectory`
        # from `:expressions` exactly like the nodal injections.
        ("dq_Vd_pu", :Vd_tf, :Vd_tpf, 1.0, "Terminal Voltage Vd", "Vd (p.u.)", "dq_Vd.svg"),
        ("dq_Vq_pu", :Vq_tf, :Vq_tpf, 1.0, "Terminal Voltage Vq", "Vq (p.u.)", "dq_Vq.svg"),
        # Pre-saturation exciter output (AVR only); equals E_fd wherever the softsat
        # clamp is inactive, so the gap between the two shows when the AVR is limiting.
        ("dq_E_fd_unlim_pu", :E_fd_unlim_tf, :E_fd_unlim_tpf, 1.0,
            "Unsaturated Field Voltage E_fd_unlim", "E_fd_unlim (p.u.)", "dq_E_fd_unlim.svg"),
        # SEXS lead-lag output feeding the gain-lag; covers only the generators that have
        # the stage, and the file is absent entirely on the Ta_exc = Tb_exc = 0 bypass.
        ("dq_E_LL_pu", :E_LL_tf, :E_LL_tpf, 1.0,
            "AVR Lead-Lag Output E_LL", "E_LL (p.u.)", "dq_E_LL.svg"),
    )
    for (fname, key_tf, key_tpf, scale, title, ylab, svg_name) in dq_series
        traj = _extract_gen_time_trajectory(dyn_model_dict, key_tf, key_tpf; scale=scale)
        traj === nothing && continue
        _write_fullbus_csv(joinpath(pf_ts_csv, "$(fname).csv"), t_window, traj, "G")
        if save_ts_plots
            _save_fullbus_svg(t_window, traj, title, ylab, pf_figures, svg_name;
                t_clear_fault=t_clear_fault, label_prefix="G")
        end
    end

    println("DQ generator trajectories saved as CSV in: ", pf_ts_csv,
        save_ts_plots ? " (figures in $pf_figures)" : "")
    return nothing
end

"""
    Save_Gfm_Converter_Results!(dyn_model_dict, dyn_parameters_dict, path_names; …)

Write grid-forming converter trajectories to CSV and optional SVG plots.

Covers the measured states (`P_meas`, `Q_meas`, `V_meas`), both Q–V PI states raw and
clipped, and the current loading `|I| / Imax` — the channel that says whether the limiter
is actually saturating (`1.0` = hard against `Imax`). Two overlay figures put each raw PI
state against its clipped counterpart, so the visible gap is the anti-windup clip at work.

No-op unless the run carried GFM units.
"""
function Save_Gfm_Converter_Results!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    path_names::OrderedDict{Symbol, String};
    t_clear_fault=nothing,
    save_ts_plots::Bool=false,
)
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    gfm_gens = Vector{Int}(get(meta, :gfm_gens, Int[]))
    isempty(gfm_gens) && return nothing
    params = get(meta, :gfm_params, nothing)
    params === nothing && return nothing

    pf_ts_csv = path_names[:pf_TS_CSV]
    mkpath(pf_ts_csv)
    t_window = dyn_parameters_dict[:time][:t_window_total]
    pf_figures = path_names[:pf_TS_figures]
    save_ts_plots && mkpath(pf_figures)

    # --- pre-fault converter set-points and t=0 filter states -------------------
    # `V_set` is the Q–V droop reference the PI integrator works against, and the four
    # t=0 scalars are the initial conditions every transient filter starts from. None of
    # them reached a results file before: they appeared only in the pre-solve `start=`
    # dump under Dispatch_WarmStart/, which records what the solver was *seeded* with,
    # not what it converged to.
    setpoint_cols = (
        ("V_set_pu", :V_set), ("E_int_pu", :E_int), ("P_meas_pu", :P_meas),
        ("Q_meas_pu", :Q_meas), ("V_meas_pu", :V_meas),
    )
    present_cols = [(name, sym) for (name, sym) in setpoint_cols
                    if haskey(dyn_model_dict[:vars], sym)]
    if !isempty(present_cols)
        df_set = DataFrame(gen = ["GFM$(g)" for g in gfm_gens])
        for (name, sym) in present_cols
            container = dyn_model_dict[:vars][sym]
            df_set[!, name] = [haskey(container, g) ? JuMP.value(container[g]) : missing
                               for g in gfm_gens]
        end
        CSV.write(joinpath(pf_ts_csv, "gfm_setpoints.csv"), df_set; delim = ';')
    end

    gfm_series = (
        ("gfm_P_meas", :P_meas_tf, :P_meas_tpf, "Filtered Active Power", "P_meas (p.u.)",
            "gfm_P_meas.svg"),
        ("gfm_Q_meas", :Q_meas_tf, :Q_meas_tpf, "Filtered Reactive Power", "Q_meas (p.u.)",
            "gfm_Q_meas.svg"),
        ("gfm_V_meas", :V_meas_tf, :V_meas_tpf, "Filtered Terminal Voltage", "V_meas (p.u.)",
            "gfm_V_meas.svg"),
        ("gfm_E_int", :E_int_tf, :E_int_tpf, "Q–V PI Integrator (clipped)", "E_int (p.u.)",
            "gfm_E_int.svg"),
        ("gfm_E_int_raw", :E_int_raw_tf, :E_int_raw_tpf, "Q–V PI Integrator (raw)",
            "E_int_raw (p.u.)", "gfm_E_int_raw.svg"),
        ("gfm_E_droop", :E_droop_tf, :E_droop_tpf, "Internal Voltage (clipped)",
            "E_droop (p.u.)", "gfm_E_droop.svg"),
        ("gfm_E_droop_raw", :E_droop_raw_tf, :E_droop_raw_tpf, "Internal Voltage (raw)",
            "E_droop_raw (p.u.)", "gfm_E_droop_raw.svg"),
    )
    traj_cache = Dict{String, OrderedDict{Int, Vector{Float64}}}()
    for (fname, key_tf, key_tpf, title, ylab, svg_name) in gfm_series
        traj = _extract_gen_time_trajectory(dyn_model_dict, key_tf, key_tpf)
        traj === nothing && continue
        traj_cache[fname] = traj
        _write_fullbus_csv(joinpath(pf_ts_csv, "$(fname).csv"), t_window, traj, "GFM")
        if save_ts_plots
            _save_fullbus_svg(t_window, traj, title, ylab, pf_figures, svg_name;
                t_clear_fault=t_clear_fault, label_prefix="GFM ")
        end
    end

    # |I| / Imax is not a stored variable: rebuild it from the dq currents, which the DQ
    # builder creates for every unit, restricted to the converter ids.
    Id = _extract_gen_time_trajectory(dyn_model_dict, :Id_tf, :Id_tpf)
    Iq = _extract_gen_time_trajectory(dyn_model_dict, :Iq_tf, :Iq_tpf)
    if Id !== nothing && Iq !== nothing
        loading = OrderedDict{Int, Vector{Float64}}()
        for gen in gfm_gens
            (haskey(Id, gen) && haskey(Iq, gen) && haskey(params, gen)) || continue
            Imax = params[gen].Imax
            Imax > 0 || continue
            loading[gen] = sqrt.(Id[gen] .^ 2 .+ Iq[gen] .^ 2) ./ Imax
        end
        if !isempty(loading)
            _write_fullbus_csv(joinpath(pf_ts_csv, "gfm_current_loading.csv"),
                t_window, loading, "GFM")
            if save_ts_plots
                _save_fullbus_svg(t_window, loading, "Converter Current Loading",
                    "|I| / Imax (-)", pf_figures, "gfm_current_loading.svg";
                    t_clear_fault=t_clear_fault, label_prefix="GFM ")
            end
        end
    end

    # Raw vs clipped overlays: the gap between the pair is the anti-windup clip acting.
    if save_ts_plots
        for (raw_key, clip_key, title, ylab, svg_name) in (
            ("gfm_E_int_raw", "gfm_E_int", "Q–V PI Integrator: raw vs clipped",
                "E_int (p.u.)", "gfm_E_int_clip.svg"),
            ("gfm_E_droop_raw", "gfm_E_droop", "Internal Voltage: raw vs clipped",
                "E_droop (p.u.)", "gfm_E_droop_clip.svg"),
        )
            (haskey(traj_cache, raw_key) && haskey(traj_cache, clip_key)) || continue
            _save_fullbus_svg(t_window, traj_cache[clip_key], title, ylab,
                pf_figures, svg_name; t_clear_fault=t_clear_fault,
                label_prefix="GFM ", overlay=traj_cache[raw_key],
                overlay_label_suffix=" raw")
        end
    end

    println("GFM converter trajectories saved as CSV in: ", pf_ts_csv,
        save_ts_plots ? " (figures in $pf_figures)" : "")
    return nothing
end

# ===================================================================================
# Per-(window, generator, step) diagnostic dumps  —  Transient_Stability/CSV/Debug/
# ===================================================================================
# Column layout mirrors the reference implementation's Debug/ folder so a parity check is
# a direct CSV diff. Each row carries the value, its distance to the relevant limit, and
# the dual, so "was this binding at t, and what did it cost?" is one line rather than a
# join across four files. Nothing here is stored by the builders: the limiter algebra and
# the filter coefficients are recomputed from solved primal values plus `meta[:gfm_params]`.

"""`(period_label, suffix, time_vector)` for each window present in the run."""
function _ts_debug_windows(dyn_parameters_dict::OrderedDict{Symbol, Any})
    t = dyn_parameters_dict[:time]
    out = Tuple{String, String, Vector{Float64}}[]
    haskey(t, :t_window_fault) && push!(out, ("fault", "tf", t[:t_window_fault]))
    if haskey(t, :t_window_postf) && !isempty(t[:t_window_postf])
        push!(out, ("postfault", "tpf", t[:t_window_postf]))
    end
    return out
end

"""Solved value of `container[gen][t]`, or `NaN` when the entry is absent."""
function _dbg_val(container, gen::Int, t::Int)
    container === nothing && return NaN
    haskey(container, gen) || return NaN
    inner = container[gen]
    haskey(inner, t) || return NaN
    return Float64(JuMP.value(inner[t]))
end

"""Dual of `dyn_model_dict[store][key][gen][t]`, or `NaN` when unavailable."""
function _dbg_dual(dyn_model_dict::OrderedDict{Symbol, Any}, store::Symbol,
                   key::Symbol, gen::Int, t::Int)
    s = get(dyn_model_dict, store, nothing)
    s === nothing && return NaN
    haskey(s, key) || return NaN
    fam = s[key]
    haskey(fam, gen) || return NaN
    haskey(fam[gen], t) || return NaN
    try
        return Float64(JuMP.dual(fam[gen][t]))
    catch
        return NaN
    end
end

_dbg_vars(dyn_model_dict, name::String, suffix::String) =
    get(dyn_model_dict[:vars], Symbol(name, "_", suffix), nothing)

"""
Filter states, coefficients and duals — one row per (window, gen, step, channel).

`α` / `β` are the coefficients of the update actually written into the model, which
depends on `meta[:gfm_integrator]` (and, under `:follow_ode_first_step`, on
`meta[:ode_first_step]` at `t = 1`):

| scheme | row | α | β |
|---|---|---|---|
| backward Euler | `x_t = α·x_prev + β·u_t` | `Tf/(Tf+Δt)` | `Δt/(Tf+Δt)` |
| trapezoidal | `x_t = α·x_prev + β·(u_t + u_prev)` | `(1−c)/(1+c)` | `c/(1+c)`, `c = Δt/(2Tf)` |

The trapezoidal row needs `u_prev`, which is deliberately not its own column (the header
is asserted verbatim by the parity test): it is `u_t` on the previous row of the same
`(period, gen, channel)` group.

`Tf` is the constant of *that row's channel* — `TPf`, `TQf` or `TVf` — not one shared
value per converter. The column keeps its short name because `channel` already says which
of the three it reports.
"""
function _write_gfm_filter_debug!(io_path, dyn_model_dict, windows, gfm_gens, params, Δt)
    rows = NamedTuple[]
    meta = dyn_model_dict[:meta]
    gfm_int = get(meta, :gfm_integrator, :backward_euler)
    ode_fs = get(meta, :ode_first_step, :trapezoidal)
    be_at(t) = gfm_int === :backward_euler ||
        (gfm_int === :follow_ode_first_step && t == 1 && ode_fs === :backward_euler)
    for (period, sfx, tvec) in windows
        # The fifth slot names the `gfm_params` field holding that channel's own time
        # constant (REGFM_A1 eq. 4–6); the three are independent, so α/β are resolved
        # per channel rather than once per converter.
        chans = (
            ("P", _dbg_vars(dyn_model_dict, "P_meas", sfx), _dbg_vars(dyn_model_dict, "Pe", sfx),
                Symbol("eq_const_gfm_filter_P_", sfx), :TPf),
            ("Q", _dbg_vars(dyn_model_dict, "Q_meas", sfx), _dbg_vars(dyn_model_dict, "Qe", sfx),
                Symbol("eq_const_gfm_filter_Q_", sfx), :TQf),
            ("V", _dbg_vars(dyn_model_dict, "V_meas", sfx), nothing,
                Symbol("eq_const_gfm_filter_V_", sfx), :TVf),
        )
        V_bus = _dbg_vars(dyn_model_dict, "V", sfx)
        prev_anchor = _gfm_debug_prev_anchor(dyn_model_dict, sfx)
        for gen in gfm_gens
            p = params[gen]
            for (ch, state, input, dual_key, Tf_field) in chans
                # `Tf` in the exported row is this channel's constant; the `channel`
                # column says which of TPf/TQf/TVf it is.
                Tf = Float64(getfield(p, Tf_field))
                c = Tf <= 1e-9 ? 0.0 : Δt / (2.0 * Tf)
                for t in eachindex(tvec)
                    α, β = if Tf <= 1e-9
                        0.0, 1.0                                  # pass-through bypass
                    elseif be_at(t)
                        Tf / (Tf + Δt), Δt / (Tf + Δt)
                    else
                        (1 - c) / (1 + c), c / (1 + c)
                    end
                    x_t = _dbg_val(state, gen, t)
                    x_prev = t == 1 ? prev_anchor(ch, gen) : _dbg_val(state, gen, t - 1)
                    u_t = ch == "V" ?
                        (V_bus === nothing ? NaN : _dbg_val(V_bus, p.bus, t)) :
                        _dbg_val(input, gen, t)
                    push!(rows, (
                        period = period, gen = gen, t_idx = t, t_s = tvec[t], channel = ch,
                        Tf = Tf, Δt = Δt, Tf_over_Δt = Tf / Δt, α = α, β = β,
                        x_t = x_t, x_prev = x_prev, u_t = u_t,
                        x_minus_u = x_t - u_t, x_minus_prev = x_t - x_prev,
                        dual_filter = _dbg_dual(dyn_model_dict, :eq_const, dual_key, gen, t),
                    ))
                end
            end
        end
    end
    isempty(rows) || CSV.write(io_path, DataFrame(rows); delim = ';')
    return length(rows)
end

"""`t = 1` previous value for a filter channel: pre-fault scalar, or last fault-on step."""
function _gfm_debug_prev_anchor(dyn_model_dict::OrderedDict{Symbol, Any}, sfx::String)
    vars = dyn_model_dict[:vars]
    return function (ch::String, gen::Int)
        key = ch == "P" ? :P_meas : ch == "Q" ? :Q_meas : :V_meas
        if sfx == "tf"
            c = get(vars, key, nothing)
            (c === nothing || !haskey(c, gen)) && return NaN
            return Float64(JuMP.value(c[gen]))
        end
        c = get(vars, Symbol(key, "_tf"), nothing)
        (c === nothing || !haskey(c, gen)) && return NaN
        return Float64(JuMP.value(last(c[gen])[2]))
    end
end

"""Current-limiter geometry, margins and duals — one row per (window, gen, step)."""
function _write_gfm_limiter_debug!(io_path, dyn_model_dict, windows, gfm_gens, params)
    rows = NamedTuple[]
    for (period, sfx, tvec) in windows
        V_bus = _dbg_vars(dyn_model_dict, "V", sfx)
        θ_bus = _dbg_vars(dyn_model_dict, "θ", sfx)
        δ_w = _dbg_vars(dyn_model_dict, "δ", sfx)
        Id_w = _dbg_vars(dyn_model_dict, "Id", sfx)
        Iq_w = _dbg_vars(dyn_model_dict, "Iq", sfx)
        E_dr = _dbg_vars(dyn_model_dict, "E_droop", sfx)
        for gen in gfm_gens
            p = params[gen]
            Xl, Imax = p.Xl, p.Imax
            Imax_Xl = Imax * Xl
            bypassed = Imax >= _GFM_IMAX_NO_LIMIT
            for t in eachindex(tvec)
                Vmag = _dbg_val(V_bus, p.bus, t)
                θ = _dbg_val(θ_bus, p.bus, t)
                δ = _dbg_val(δ_w, gen, t)
                Ed = _dbg_val(E_dr, gen, t)
                Id = _dbg_val(Id_w, gen, t)
                Iq = _dbg_val(Iq_w, gen, t)
                Vd = Vmag * sin(δ - θ)
                Vq = Vmag * cos(δ - θ)
                ΔVd = Ed - Vq
                ΔVq = Vd
                ΔVraw = sqrt(ΔVd^2 + ΔVq^2)
                Iraw = ΔVraw / Xl
                Iout = sqrt(Id^2 + Iq^2)
                den_Xl = bypassed ? ΔVraw : max(ΔVraw, Imax_Xl)
                scale = bypassed ? 1.0 : (den_Xl == 0 ? 1.0 : Imax_Xl / den_Xl)
                push!(rows, (
                    period = period, gen = gen, t_idx = t, t_s = tvec[t], bus = p.bus,
                    Xl = Xl, Imax = Imax, Vd = Vd, Vq = Vq, ΔVd = ΔVd, ΔVq = ΔVq,
                    ΔVraw = ΔVraw, ΔVlim = Imax_Xl, ΔV_margin = Imax_Xl - ΔVraw,
                    Iraw = Iraw, Iout = Iout, I_margin = Imax - Iout,
                    Iout_over_Imax = Imax == 0 ? NaN : Iout / Imax,
                    transition_distance_I = Iraw - Imax,
                    Iraw_Xl = ΔVraw, Iout_Xl = Xl * Iout, Imax_Xl = Imax_Xl,
                    den_Xl = den_Xl, scale = scale,
                    dual_c_Id = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_limiter_Id_", sfx), gen, t),
                    dual_c_Iq = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_limiter_Iq_", sfx), gen, t),
                    dual_c_Pe = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Pe_", sfx), gen, t),
                    dual_c_Qe = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Qe_", sfx), gen, t),
                ))
            end
        end
    end
    isempty(rows) || CSV.write(io_path, DataFrame(rows); delim = ';')
    return length(rows)
end

"""Q–V PI states, clip distances and duals — one row per (window, gen, step)."""
function _write_gfm_voltage_debug!(io_path, dyn_model_dict, windows, gfm_gens, params)
    rows = NamedTuple[]
    for (period, sfx, tvec) in windows
        V_meas = _dbg_vars(dyn_model_dict, "V_meas", sfx)
        Q_meas = _dbg_vars(dyn_model_dict, "Q_meas", sfx)
        E_int_raw = _dbg_vars(dyn_model_dict, "E_int_raw", sfx)
        E_int = _dbg_vars(dyn_model_dict, "E_int", sfx)
        E_dr_raw = _dbg_vars(dyn_model_dict, "E_droop_raw", sfx)
        E_dr = _dbg_vars(dyn_model_dict, "E_droop", sfx)
        V_set = get(dyn_model_dict[:vars], :V_set, nothing)
        for gen in gfm_gens
            p = params[gen]
            vset = (V_set === nothing || !haskey(V_set, gen)) ? NaN :
                Float64(JuMP.value(V_set[gen]))
            for t in eachindex(tvec)
                err = vset - p.mq * _dbg_val(Q_meas, gen, t) - _dbg_val(V_meas, gen, t)
                eir = _dbg_val(E_int_raw, gen, t); ei = _dbg_val(E_int, gen, t)
                edr = _dbg_val(E_dr_raw, gen, t); ed = _dbg_val(E_dr, gen, t)
                push!(rows, (
                    period = period, gen = gen, t_idx = t, t_s = tvec[t],
                    Emin = p.Emin, Emax = p.Emax, voltage_error = err,
                    E_int_raw = eir, E_int = ei,
                    E_int_to_Emax = p.Emax - ei, E_int_to_Emin = ei - p.Emin,
                    E_droop_raw = edr, E_droop = ed,
                    E_droop_to_Emax = p.Emax - ed, E_droop_to_Emin = ed - p.Emin,
                    dual_Eint_raw = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Eint_raw_", sfx), gen, t),
                    dual_Eint_clip = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Eint_clip_", sfx), gen, t),
                    dual_Edroop_raw = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Edroop_raw_", sfx), gen, t),
                    dual_Edroop_clip = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_gfm_Edroop_clip_", sfx), gen, t),
                ))
            end
        end
    end
    isempty(rows) || CSV.write(io_path, DataFrame(rows); delim = ';')
    return length(rows)
end

"""
COI angle/speed series for `key`, wherever the builder put it.

It is a JuMP variable when the corridor references it and a plain `AffExpr` otherwise
(see `_attach_δCOI!`). `JuMP.value` handles both, so the diagnostics do not care which
form a run produced — only that they find it.
"""
function _coi_series(dyn_model_dict::OrderedDict{Symbol, Any}, key::Symbol)
    series = get(dyn_model_dict[:vars], key, nothing)
    series === nothing || return series
    exprs = get(dyn_model_dict, :expressions, nothing)
    exprs === nothing && return nothing
    return get(exprs, key, nothing)
end

"""
SG swing states with their δ–COI (and Δω–COI) corridor margins.

The reference columns come first so a diff lines up; the package extras follow. `δ_util`
is the corridor analogue of `Iout_over_Imax` in the limiter file: `1.0` means the machine
sits on the stability corridor, and the matching dual is what holding it there costs.
`dual_accel` is `NaN` by construction — the package inlines the accelerating power into
the swing rows instead of carrying a separate `a` constraint, so the *value* is available
but there is no row to price.
"""
function _write_swing_debug!(io_path, dyn_model_dict, dyn_parameters_dict, DGEN_DYN,
                             windows, sg_gens)
    rows = NamedTuple[]
    meta = dyn_model_dict[:meta]
    Δt = dyn_parameters_dict[:time][:t_step]
    # δ_tol is the (lower, upper) corridor in **radians**, matching the δ variables here;
    # angle_rel_COI.csv reports degrees, so compare ratios, not raw values.
    δ_tol = get(dyn_parameters_dict[:common], :δ_tol, (NaN, NaN))
    Δω_tol = get(meta, :Δω_tol, (NaN, NaN))
    has_Δω_corridor = get(meta, :constrain_Δω, false)
    # Reference machine under bound_style_δ = :highest_H / :ref_gen; `nothing` on the
    # COI-referenced styles, where the δ_ref columns stay NaN.
    δ_ref_gen = get(meta, :δ_ref_gen_resolved, nothing)
    for (period, sfx, tvec) in windows
        δ_w = _dbg_vars(dyn_model_dict, "δ", sfx)
        Δω_w = _dbg_vars(dyn_model_dict, "Δω", sfx)
        Pe_w = _dbg_vars(dyn_model_dict, "Pe", sfx)
        Pm_w = _dbg_vars(dyn_model_dict, "Pm", sfx)
        # δCOI is a variable only when the δ corridor references it; a machine-referenced
        # run keeps it as an unconstrained expression so this column stays populated.
        δCOI = _coi_series(dyn_model_dict, Symbol("δCOI_", sfx))
        ΔωCOI = get(dyn_model_dict[:vars], Symbol("ΔωCOI_", sfx), nothing)
        P_m0 = get(dyn_model_dict[:vars], :P_m, nothing)
        for gen in sg_gens
            # By id: this is the SG-only `sys.DGEN_DYN`, not the builder's aligned copy.
            row = _dyn_row_index(DGEN_DYN, gen)
            H = Float64(DGEN_DYN.H[row]); D = Float64(DGEN_DYN.D[row])
            for t in eachindex(tvec)
                δ = _dbg_val(δ_w, gen, t)
                Δω = _dbg_val(Δω_w, gen, t)
                δ_prev = t == 1 ? NaN : _dbg_val(δ_w, gen, t - 1)
                Δω_prev = t == 1 ? NaN : _dbg_val(Δω_w, gen, t - 1)
                Pm = Pm_w !== nothing ? _dbg_val(Pm_w, gen, t) :
                    (P_m0 !== nothing && haskey(P_m0, gen) ?
                        Float64(JuMP.value(P_m0[gen])) : NaN)
                Pe = _dbg_val(Pe_w, gen, t)
                a = H == 0 ? NaN : (Pm - Pe - D * Δω) / (2 * H)
                coi = δCOI === nothing || !haskey(δCOI, t) ? NaN :
                    Float64(JuMP.value(δCOI[t]))
                δ_rel = δ - coi
                ω_coi = ΔωCOI === nothing || !haskey(ΔωCOI, t) ? NaN :
                    Float64(JuMP.value(ΔωCOI[t]))
                Δω_rel = Δω - ω_coi
                # Angle relative to the reference machine — the quantity the :highest_H /
                # :ref_gen corridors actually bound, and the reference's own row is 0.
                δ_ref_val = δ_ref_gen === nothing ? NaN :
                    _dbg_val(δ_w, δ_ref_gen, t)
                δ_rel_ref = δ - δ_ref_val
                push!(rows, (
                    period = period, gen = gen, t_idx = t, t_s = tvec[t],
                    H = H, D = D, Δt = Δt, a = a,
                    Δω_curr = Δω, Δω_prev = Δω_prev, Δω_step = Δω - Δω_prev,
                    Δω_step_over_Δt = (Δω - Δω_prev) / Δt,
                    δ_curr = δ, δ_prev = δ_prev, δ_step = δ - δ_prev,
                    δ_step_over_Δt = (δ - δ_prev) / Δt,
                    dual_accel = NaN,
                    dual_delta = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_δ_", sfx), gen, t),
                    dual_omega = _dbg_dual(dyn_model_dict, :eq_const,
                        Symbol("eq_const_Δω_", sfx), gen, t),
                    # --- package extras: distance to the stability corridor -------------
                    δ_rel_COI = δ_rel, δ_tol_lo = δ_tol[1], δ_tol_hi = δ_tol[2],
                    δ_margin_lo = δ_rel - δ_tol[1], δ_margin_hi = δ_tol[2] - δ_rel,
                    δ_util = _corridor_utilisation(δ_rel, δ_tol),
                    dual_δ_COI_lower = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_δ_COI_", sfx, "_lower"), gen, t),
                    dual_δ_COI_upper = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_δ_COI_", sfx, "_upper"), gen, t),
                    Δω_rel_COI = Δω_rel,
                    Δω_margin_lo = has_Δω_corridor ? Δω_rel - Δω_tol[1] : NaN,
                    Δω_margin_hi = has_Δω_corridor ? Δω_tol[2] - Δω_rel : NaN,
                    Δω_util = has_Δω_corridor ? _corridor_utilisation(Δω_rel, Δω_tol) : NaN,
                    dual_Δω_COI_lower = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_Δω_COI_", sfx, "_lower"), gen, t),
                    dual_Δω_COI_upper = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_Δω_COI_", sfx, "_upper"), gen, t),
                    # --- machine-referenced δ corridor and absolute Δω corridor ---------
                    # NaN unless the matching style built them; the columns are always
                    # present so runs under different styles stay diff-able side by side.
                    δ_ref_gen = δ_ref_gen === nothing ? 0 : δ_ref_gen,
                    δ_rel_ref = δ_rel_ref,
                    δ_ref_margin_lo = δ_rel_ref - δ_tol[1],
                    δ_ref_margin_hi = δ_tol[2] - δ_rel_ref,
                    δ_ref_util = _corridor_utilisation(δ_rel_ref, δ_tol),
                    dual_δ_ref_lower = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_δ_ref_", sfx, "_lower"), gen, t),
                    dual_δ_ref_upper = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_δ_ref_", sfx, "_upper"), gen, t),
                    Δω_abs_margin_lo = has_Δω_corridor ? Δω - Δω_tol[1] : NaN,
                    Δω_abs_margin_hi = has_Δω_corridor ? Δω_tol[2] - Δω : NaN,
                    dual_Δω_abs_lower = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_Δω_abs_", sfx, "_lower"), gen, t),
                    dual_Δω_abs_upper = _dbg_dual(dyn_model_dict, :ineq_const,
                        Symbol("ineq_const_Δω_abs_", sfx, "_upper"), gen, t),
                ))
            end
        end
    end
    isempty(rows) || CSV.write(io_path, DataFrame(rows); delim = ';')
    return length(rows)
end

"""
Fraction of the corridor a deviation uses: `1.0` = sitting on the limit.

Each side is measured against its own tolerance, so an asymmetric corridor
(`δ_tol_deg_lower` ≠ `δ_tol_deg_upper`) reports honestly on both sides.
"""
function _corridor_utilisation(x::Float64, tol::Tuple{Float64, Float64})::Float64
    (isnan(x) || isnan(tol[1]) || isnan(tol[2])) && return NaN
    lim = x < 0 ? abs(tol[1]) : abs(tol[2])
    lim == 0 && return NaN
    return abs(x) / lim
end

"""
    Save_TS_Debug_CSV!(dyn_model_dict, dyn_parameters_dict, DGEN_DYN, path_names)

Write the per-(window, generator, step) diagnostic dumps under
`Transient_Stability/CSV/Debug/`. Gated by `RunConfig.save_ts_debug_csv` (default `false`)
because the filter file alone is three rows per converter per step.
"""
function Save_TS_Debug_CSV!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    DGEN_DYN::DataFrame,
    path_names::OrderedDict{Symbol, String},
)
    windows = _ts_debug_windows(dyn_parameters_dict)
    isempty(windows) && return nothing
    pf_debug = joinpath(path_names[:pf_TS_CSV], "Debug")
    mkpath(pf_debug)

    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    gfm_gens = Vector{Int}(get(meta, :gfm_gens, Int[]))
    params = get(meta, :gfm_params, nothing)
    Δt = dyn_parameters_dict[:time][:t_step]
    written = String[]

    if !isempty(gfm_gens) && params !== nothing
        n = _write_gfm_filter_debug!(joinpath(pf_debug, "gfm_filter_debug.csv"),
            dyn_model_dict, windows, gfm_gens, params, Δt)
        n > 0 && push!(written, "gfm_filter_debug.csv")
        n = _write_gfm_limiter_debug!(joinpath(pf_debug, "gfm_limiter_debug.csv"),
            dyn_model_dict, windows, gfm_gens, params)
        n > 0 && push!(written, "gfm_limiter_debug.csv")
        n = _write_gfm_voltage_debug!(joinpath(pf_debug, "gfm_voltage_debug.csv"),
            dyn_model_dict, windows, gfm_gens, params)
        n > 0 && push!(written, "gfm_voltage_debug.csv")
    end

    # `meta[:sg_gens]` is stamped only on the GFM-aware paths, so on an SG-only run it may
    # be absent and every active unit is a machine. That fallback must NOT fire when the
    # run has converters: on a converter-only fleet `sg_gens` is legitimately empty, and
    # falling back to `active_gen` would hand `_write_swing_debug!` a list of converters to
    # look up `H`/`D` for — rows that do not exist in `DGEN_DYN`. Writing no file is the
    # correct outcome; `swing_debug.csv` is a machine diagnostic by construction.
    sg_gens = Vector{Int}(get(meta, :sg_gens, Int[]))
    if isempty(sg_gens) && isempty(Vector{Int}(get(meta, :gfm_ids, Int[])))
        sg_gens = Vector{Int}(get(dyn_model_dict, :active_gen, Int[]))
    end
    if !isempty(sg_gens)
        n = _write_swing_debug!(joinpath(pf_debug, "swing_debug.csv"),
            dyn_model_dict, dyn_parameters_dict, DGEN_DYN, windows, sg_gens)
        n > 0 && push!(written, "swing_debug.csv")
    end

    println("TS debug CSVs saved in: ", pf_debug, " (", join(written, ", "), ")")
    return nothing
end

"""Concatenate fault-on and post-fault per-bus time series."""
function _merge_bus_time_dict(
    tf::OrderedDict{Int, Vector{Float64}},
    tpf::Union{Nothing, OrderedDict{Int, Vector{Float64}}},
)
    out = OrderedDict{Int, Vector{Float64}}()
    for k in keys(tf)
        out[k] = tpf === nothing ? tf[k] : vcat(tf[k], tpf[k])
    end
    return out
end

"""Extract per-bus trajectories from nested JuMP variable containers."""
function _fullbus_values_from_vars(
    tf_vars::OrderedDict,
    tpf_vars::Union{Nothing, OrderedDict};
    scale::Float64=1.0,
    to_deg::Bool=false,
)
    tf_vals = OrderedDict(
        bus => scale .* [JuMP.value(v) for (_, v) in inner]
        for (bus, inner) in tf_vars)
    tpf_vals = tpf_vars === nothing ? nothing : OrderedDict(
        bus => scale .* [JuMP.value(v) for (_, v) in inner]
        for (bus, inner) in tpf_vars)
    merged = _merge_bus_time_dict(tf_vals, tpf_vals)
    to_deg && (merged = OrderedDict(k => rad2deg.(v) for (k, v) in merged))
    return merged
end

"""Extract per-bus nodal injection trajectories from stored JuMP expressions."""
function _fullbus_values_from_exprs(
    tf_exprs::OrderedDict,
    tpf_exprs::Union{Nothing, OrderedDict};
    scale::Float64=1.0,
)
    tf_vals = OrderedDict(
        bus => scale .* [JuMP.value(expr) for (_, expr) in inner]
        for (bus, inner) in tf_exprs)
    tpf_vals = tpf_exprs === nothing ? nothing : OrderedDict(
        bus => scale .* [JuMP.value(expr) for (_, expr) in inner]
        for (bus, inner) in tpf_exprs)
    return _merge_bus_time_dict(tf_vals, tpf_vals)
end

"""
    Save_FullBus_Network_Results!(dyn_model_dict, dyn_parameters_dict, base_MVA, path_names)

Save bus-level active/reactive power, voltage magnitude, and angle as CSV and SVG.
Called from `Save_Results_Dynamic_Model` when `meta[:network_form] == "FULL_BUS"`.
"""
function Save_FullBus_Network_Results!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    dyn_parameters_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
    path_names::OrderedDict{Symbol, String};
    t_clear_fault=nothing,
    save_ts_plots::Bool=false,
)
    t_window = dyn_parameters_dict[:time][:t_window_total]
    pf_ts_csv = path_names[:pf_TS_CSV]
    mkpath(pf_ts_csv)
    pf_figures = path_names[:pf_TS_figures]
    save_ts_plots && mkpath(pf_figures)

    V_t = _fullbus_values_from_vars(
        dyn_model_dict[:vars][:V_tf],
        get(dyn_model_dict[:vars], :V_tpf, nothing))
    θ_t = _fullbus_values_from_vars(
        dyn_model_dict[:vars][:θ_tf],
        get(dyn_model_dict[:vars], :θ_tpf, nothing);
        to_deg=true)
    exprs = dyn_model_dict[:expressions]
    P_bus = _fullbus_values_from_exprs(
        exprs[:P_inj_tf], get(exprs, :P_inj_tpf, nothing); scale=base_MVA)
    Q_bus = _fullbus_values_from_exprs(
        exprs[:Q_inj_tf], get(exprs, :Q_inj_tpf, nothing); scale=base_MVA)

    Pet = _extract_gen_time_trajectory(dyn_model_dict, :Pe_tf, :Pe_tpf; scale=base_MVA)
    Qet = _extract_gen_time_trajectory(dyn_model_dict, :Qe_tf, :Qe_tpf; scale=base_MVA)

    _write_fullbus_csv(joinpath(pf_ts_csv, "bus_active_power.csv"), t_window, P_bus, "Bus_")
    _write_fullbus_csv(joinpath(pf_ts_csv, "bus_reactive_power.csv"), t_window, Q_bus, "Bus_")
    _write_fullbus_csv(joinpath(pf_ts_csv, "bus_voltage_magnitude.csv"), t_window, V_t, "Bus_")
    _write_fullbus_csv(joinpath(pf_ts_csv, "bus_voltage_angle.csv"), t_window, θ_t, "Bus_")

    if Pet !== nothing
        _write_fullbus_csv(joinpath(pf_ts_csv, "generator_active_power.csv"), t_window, Pet, "G")
    end
    if Qet !== nothing
        _write_fullbus_csv(joinpath(pf_ts_csv, "generator_reactive_power.csv"), t_window, Qet, "G")
    end

    if save_ts_plots
        _save_fullbus_svg(t_window, P_bus, "Bus Active Power", "P (MW)", pf_figures,
            "Bus active power vs time.svg"; t_clear_fault=t_clear_fault)
        _save_fullbus_svg(t_window, Q_bus, "Bus Reactive Power", "Q (MVAr)", pf_figures,
            "Bus reactive power vs time.svg"; t_clear_fault=t_clear_fault)
        _save_fullbus_svg(t_window, V_t, "Bus Voltage Magnitude", "V (p.u.)", pf_figures,
            "Bus voltage magnitude vs time.svg"; t_clear_fault=t_clear_fault)
        _save_fullbus_svg(t_window, θ_t, "Bus Voltage Angle", "θ (deg)", pf_figures,
            "Bus voltage angle vs time.svg"; t_clear_fault=t_clear_fault)
        Pet !== nothing && _save_fullbus_svg(t_window, Pet, "Generator Active Power", "P_e (MW)", pf_figures,
            "Generator active power vs time.svg"; t_clear_fault=t_clear_fault, label_prefix="G")
        Qet !== nothing && _save_fullbus_svg(t_window, Qet, "Generator Reactive Power", "Q_e (MVAr)", pf_figures,
            "Generator reactive power vs time.svg"; t_clear_fault=t_clear_fault, label_prefix="G")
    end

    println("FULL_BUS network trajectories saved as CSV in: ", pf_ts_csv,
        save_ts_plots ? " (figures in $pf_figures)" : "")
    return nothing
end

"""
    _write_fullbus_csv(filename, t_window, data, prefix)

Write `t` + one column per id of `data`, columns ordered by sorted id and named
`"\$(prefix)\$(id)"` (e.g. `Bus_5`, `G2`, `GFM3`).

The header is built here, from the same sorted key vector used to lay out the
matrix, so a caller cannot pass names derived from a *different* dict and have
the columns silently mislabelled — which is what the earlier
`col_names::Vector{String}` signature allowed.
"""
function _write_fullbus_csv(
    filename::String,
    t_window::Vector{Float64},
    data::OrderedDict{Int, Vector{Float64}},
    prefix::String,
)
    ids = sort(collect(keys(data)))
    # Nothing to write rather than an `hcat` DimensionMismatch against the time column.
    # Reachable whenever a family spans no units — an SG-only family on a converter-only
    # fleet, or every unit of a family tripped by the disturbance.
    isempty(ids) && return nothing
    mat = hcat([data[i] for i in ids]...)
    df = DataFrame(hcat(t_window, mat), vcat("t", ["$(prefix)$(i)" for i in ids]))
    CSV.write(filename, df; delim=';')
    return nothing
end

