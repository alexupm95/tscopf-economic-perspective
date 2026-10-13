#=
================================================================================
 bound_builders.jl — shared variable-bound construction and dual extraction
================================================================================
=#

"""One registered variable-bound family for registry-driven dual export."""
struct BoundManifestEntry
    export_key::Symbol
    side::Symbol              # :lower | :upper
    vars::Any                 # OrderedDict{Int,VariableRef} or nested gen×time
    layout::Symbol            # :gen_indexed | :time_indexed | :per_gen_time | :per_bus_time
end

function ensure_model_meta!(dict::OrderedDict{Symbol, Any})
    if !haskey(dict, :meta)
        dict[:meta] = OrderedDict{Symbol, Any}()
    end
    meta = dict[:meta]
    if !haskey(meta, :bound_manifest)
        meta[:bound_manifest] = BoundManifestEntry[]
    end
    return meta
end

"""Initialize a steady-state `opf_dict` with vars/eq/ineq/meta containers."""
function init_opf_dict!(extras::Pair{Symbol, <:Any}...)
    opf_dict = OrderedDict{Symbol, Any}()
    opf_dict[:vars] = OrderedDict{Symbol, Any}()
    opf_dict[:eq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    opf_dict[:ineq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    opf_dict[:meta] = OrderedDict{Symbol, Any}()
    ensure_model_meta!(opf_dict)
    for (k, v) in extras
        opf_dict[k] = v
    end
    return opf_dict
end

function _limit_at(lim, idx::Int)
    lim isa AbstractVector && return lim[idx]
    return lim
end

function _var_start(start, idx::Int)
    start === nothing && return nothing
    start isa AbstractVector && return start[idx]
    return start
end

"""JuMP lower/upper kwargs for one side; non-finite limits → unbounded on that side."""
function _jump_bound_kw(lim)
    lim === nothing && return nothing
    isfinite(lim) || return nothing
    return Float64(lim)
end

"""
    create_bounded_variable!(model; base_name, start, min_lim, max_lim, encoding)

Create one scalar decision variable. Returns `VariableRef` only.
"""
function create_bounded_variable!(
    model::JuMP.Model;
    base_name::String,
    start=nothing,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
)
    if encoding == VARIABLE && (min_lim !== nothing || max_lim !== nothing)
        lb = _jump_bound_kw(min_lim)
        ub = _jump_bound_kw(max_lim)
        if start !== nothing
            if lb !== nothing && ub !== nothing
                return JuMP.@variable(model, lower_bound = lb, upper_bound = ub,
                    base_name = base_name, start = start)
            elseif lb !== nothing
                return JuMP.@variable(model, lower_bound = lb,
                    base_name = base_name, start = start)
            elseif ub !== nothing
                return JuMP.@variable(model, upper_bound = ub,
                    base_name = base_name, start = start)
            end
        end
        if lb !== nothing && ub !== nothing
            return JuMP.@variable(model, lower_bound = lb, upper_bound = ub, base_name = base_name)
        elseif lb !== nothing
            return JuMP.@variable(model, lower_bound = lb, base_name = base_name)
        elseif ub !== nothing
            return JuMP.@variable(model, upper_bound = ub, base_name = base_name)
        end
    end
    if start !== nothing
        return JuMP.@variable(model, base_name = base_name, start = start)
    end
    return JuMP.@variable(model, base_name = base_name)
end

function register_bound_manifest!(
    meta::OrderedDict{Symbol, Any},
    export_key::Symbol,
    side::Symbol,
    vars::Any,
    layout::Symbol,
)
    if !haskey(meta, :bound_manifest)
        meta[:bound_manifest] = BoundManifestEntry[]
    end
    push!(meta[:bound_manifest], BoundManifestEntry(export_key, side, vars, layout))
    return nothing
end

"""
    build_indexed_scalar_vars!(
        model, ids, var_name; bounded, min_lim, max_lim, start, encoding, meta,
        export_key_lower, export_key_upper,
    ) -> NamedTuple{(:vars, :lower, :upper)}

Fixed-shape return for all `var_*!` steady-state / pre-fault scalar builders.
"""
function build_indexed_scalar_vars!(
    model::JuMP.Model,
    ids::AbstractVector{Int},
    var_name::String;
    bounded::Bool=true,
    min_lim=nothing,
    max_lim=nothing,
    start=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta::Union{Nothing, OrderedDict{Symbol, Any}}=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    vars = OrderedDict{Int, JuMP.VariableRef}()
    lower = nothing
    upper = nothing

    use_var_bounds = bounded && encoding == VARIABLE &&
        (min_lim !== nothing || max_lim !== nothing)

    for (idx, id) in enumerate(ids)
        v = create_bounded_variable!(
            model;
            base_name = var_name * "[$id]",
            start = _var_start(start, idx),
            min_lim = bounded ? _limit_at(min_lim, idx) : nothing,
            max_lim = bounded ? _limit_at(max_lim, idx) : nothing,
            encoding = use_var_bounds ? VARIABLE : CONSTRAINT,
        )
        vars[id] = v
    end

    if bounded && !use_var_bounds
        if min_lim !== nothing
            lower = OrderedDict{Int, JuMP.ConstraintRef}()
            for (idx, id) in enumerate(ids)
                lb = _limit_at(min_lim, idx)
                isfinite(lb) || continue
                lower[id] = JuMP.@constraint(model, lb - vars[id] ≤ 0.0)
            end
            isempty(lower) && (lower = nothing)
        end
        if max_lim !== nothing
            upper = OrderedDict{Int, JuMP.ConstraintRef}()
            for (idx, id) in enumerate(ids)
                ub = _limit_at(max_lim, idx)
                isfinite(ub) || continue
                upper[id] = JuMP.@constraint(model, vars[id] - ub ≤ 0.0)
            end
            isempty(upper) && (upper = nothing)
        end
    elseif use_var_bounds && meta !== nothing
        if min_lim !== nothing && export_key_lower !== nothing &&
           any(isfinite, min_lim isa AbstractVector ? min_lim : (min_lim,))
            register_bound_manifest!(meta, export_key_lower, :lower, vars, :gen_indexed)
        end
        if max_lim !== nothing && export_key_upper !== nothing &&
           any(isfinite, max_lim isa AbstractVector ? max_lim : (max_lim,))
            register_bound_manifest!(meta, export_key_upper, :upper, vars, :gen_indexed)
        end
    end

    return (vars = vars, lower = lower, upper = upper)
end

"""
    store_scalar_var_bounds!(
        dict, var_key, result, lower_key, upper_key,
    )

Write variables and bound containers (ineq_const or bound_manifest) into a model dict.
"""
function store_scalar_var_bounds!(
    dict::OrderedDict{Symbol, Any},
    var_key::Symbol,
    result::NamedTuple,
    lower_key::Symbol,
    upper_key::Symbol,
)
    dict[:vars][var_key] = result.vars
    if result.lower !== nothing
        dict[:ineq_const][lower_key] = result.lower
    end
    if result.upper !== nothing
        dict[:ineq_const][upper_key] = result.upper
    end
    return result.vars
end

"""
    extract_var_side_dual(var, side::Symbol) -> Float64

Dual multiplier for a variable bound (`:lower` or `:upper`), normalized to match
explicit ≤-form constraint duals used elsewhere in TSCOPF export.
"""
function extract_var_side_dual(var::JuMP.VariableRef, side::Symbol)::Float64
    if side == :lower
        JuMP.has_lower_bound(var) || return 0.0
        # JuMP reports LowerBoundRef duals with opposite sign to explicit
        # ≤-form rows `lb - x ≤ 0` used in the CONSTRAINT encoding.
        return -JuMP.dual(JuMP.LowerBoundRef(var))
    elseif side == :upper
        JuMP.has_upper_bound(var) || return 0.0
        return JuMP.dual(JuMP.UpperBoundRef(var))
    end
    throw(ArgumentError("side must be :lower or :upper, got $side"))
end

"""Extract duals from a gen-indexed variable container for one bound side."""
function extract_gen_indexed_var_bound_duals(
    vars::OrderedDict{Int, JuMP.VariableRef},
    side::Symbol,
)::Tuple{Vector{Int}, Vector{Float64}}
    ids = collect(keys(vars))
    vals = [extract_var_side_dual(vars[id], side) for id in ids]
    return ids, vals
end

"""Extract duals from a time-indexed variable container for one bound side."""
function extract_time_indexed_var_bound_duals(
    vars::OrderedDict{Int, JuMP.VariableRef},
    side::Symbol,
)::Tuple{Vector{Int}, Vector{Float64}}
    ids = collect(keys(vars))
    vals = [extract_var_side_dual(vars[t], side) for t in ids]
    return ids, vals
end

"""Extract duals from a per-gen×time variable container for one bound side."""
function extract_per_gen_time_var_bound_duals(
    vars,
    side::Symbol,
)::OrderedDict{Int, Vector{Float64}}
    out = OrderedDict{Int, Vector{Float64}}()
    for (gen, inner) in vars
        out[gen] = [extract_var_side_dual(inner[t], side) for t in collect(keys(inner))]
    end
    return out
end

"""Extract duals from a per-bus×time variable container for one bound side."""
extract_per_bus_time_var_bound_duals(vars, side::Symbol) =
    extract_per_gen_time_var_bound_duals(vars, side)

function bound_manifest_entry_present(
    dict::OrderedDict{Symbol, Any},
    export_key::Symbol,
)::Bool
    meta = get(dict, :meta, nothing)
    meta === nothing && return false
    manifest = get(meta, :bound_manifest, BoundManifestEntry[])
    return any(e -> e.export_key == export_key, manifest)
end

function find_bound_manifest_entry(
    dict::OrderedDict{Symbol, Any},
    export_key::Symbol,
)::Union{Nothing, BoundManifestEntry}
    meta = get(dict, :meta, nothing)
    meta === nothing && return nothing
    for entry in get(meta, :bound_manifest, BoundManifestEntry[])
        entry.export_key == export_key && return entry
    end
    return nothing
end

"""
    extract_bound_manifest_duals(dict, export_key) -> (ids, vals) or nested dict

Returns the same shape as the corresponding constraint-based extractor.
"""
function extract_bound_manifest_duals(dict::OrderedDict{Symbol, Any}, export_key::Symbol)
    entry = find_bound_manifest_entry(dict, export_key)
    entry === nothing && throw(ArgumentError("bound manifest entry $export_key not found"))
    if entry.layout == :gen_indexed
        return extract_gen_indexed_var_bound_duals(entry.vars, entry.side)
    elseif entry.layout == :time_indexed
        return extract_time_indexed_var_bound_duals(entry.vars, entry.side)
    elseif entry.layout == :per_gen_time
        return extract_per_gen_time_var_bound_duals(entry.vars, entry.side)
    elseif entry.layout == :per_bus_time
        return extract_per_bus_time_var_bound_duals(entry.vars, entry.side)
    end
    throw(ArgumentError("Unhandled bound manifest layout: $(entry.layout)"))
end
