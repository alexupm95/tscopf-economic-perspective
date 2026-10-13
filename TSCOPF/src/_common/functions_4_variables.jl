# ===========================================================================
# ------------------------------- Variables ---------------------------------
#
# All bounded scalar builders return NamedTuple (vars, lower, upper):
#   - CONSTRAINT encoding: lower/upper are ConstraintRef OrderedDicts or nothing
#   - VARIABLE encoding: bounds on @variable; lower/upper nothing; manifest in meta
# ===========================================================================

# Function to create the variable voltage magnitude
function var_voltage_magnitude!(
    model::Model,
    bus_id::AbstractVector{Int},
    var_name::String;
    bounded::Bool=true,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(bus_id)), var_name;
        bounded = bounded,
        min_lim = min_lim,
        max_lim = max_lim,
        start = 1.0,
        encoding = encoding,
        meta = meta,
        export_key_lower = export_key_lower,
        export_key_upper = export_key_upper,
    )
end

# Function to create the variable voltage angle (it can be bounded or unbouded)
function var_voltage_angle!(
    model::Model,
    bus_id::AbstractVector{Int},
    var_name::String;
    bounded::Bool=true,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(bus_id)), var_name;
        bounded = bounded,
        min_lim = min_lim,
        max_lim = max_lim,
        start = 0.0,
        encoding = encoding,
        meta = meta,
        export_key_lower = export_key_lower,
        export_key_upper = export_key_upper,
    )
end

# Function to create the variable active power generated (it can be bounded or unbouded)
function var_gen_power_active!(
    model::Model,
    gen_id::AbstractVector{Int},
    var_name::String;
    bounded::Bool=true,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(gen_id)), var_name;
        bounded = bounded,
        min_lim = min_lim,
        max_lim = max_lim,
        encoding = encoding,
        meta = meta,
        export_key_lower = export_key_lower,
        export_key_upper = export_key_upper,
    )
end

# Function to create the variable reactive power generated (it can be bounded or unbouded)
function var_gen_power_reactive!(
    model::Model,
    gen_id::AbstractVector{Int},
    var_name::String;
    bounded::Bool=true,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(gen_id)), var_name;
        bounded = bounded,
        min_lim = min_lim,
        max_lim = max_lim,
        encoding = encoding,
        meta = meta,
        export_key_lower = export_key_lower,
        export_key_upper = export_key_upper,
    )
end

"""Branch-flow scalar builder (custom base_name includes bus pair)."""
function _build_branch_flow_vars!(
    model::Model,
    branch_id::AbstractVector{Int},
    from_bus::AbstractVector{Int},
    to_bus::AbstractVector{Int},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
    flip_buses::Bool=false,
)
    vars = OrderedDict{Int, JuMP.VariableRef}()
    lower = nothing
    upper = nothing
    use_var_bounds = bounded && encoding == VARIABLE &&
        (min_lim !== nothing || max_lim !== nothing)

    for (idx, branch) in enumerate(branch_id)
        fb = flip_buses ? to_bus[idx] : from_bus[idx]
        tb = flip_buses ? from_bus[idx] : to_bus[idx]
        base_name = var_name * "[$branch, ($fb, $tb)]"
        v = create_bounded_variable!(
            model;
            base_name = base_name,
            min_lim = bounded ? _limit_at(min_lim, idx) : nothing,
            max_lim = bounded ? _limit_at(max_lim, idx) : nothing,
            encoding = use_var_bounds ? VARIABLE : CONSTRAINT,
        )
        vars[branch] = v
    end

    if bounded && !use_var_bounds
        if min_lim !== nothing
            lower = OrderedDict{Int, JuMP.ConstraintRef}()
            for (idx, branch) in enumerate(branch_id)
                lb = _limit_at(min_lim, idx)
                isfinite(lb) || continue
                lower[branch] = JuMP.@constraint(model, lb - vars[branch] ≤ 0.0)
            end
            isempty(lower) && (lower = nothing)
        end
        if max_lim !== nothing
            upper = OrderedDict{Int, JuMP.ConstraintRef}()
            for (idx, branch) in enumerate(branch_id)
                ub = _limit_at(max_lim, idx)
                isfinite(ub) || continue
                upper[branch] = JuMP.@constraint(model, vars[branch] - ub ≤ 0.0)
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

# Function to create the variable active power flow from bus i to bus k
function var_powerflow_ik_active!(
    model::Model,
    branch_id::AbstractVector{Int},
    from_bus::AbstractVector{Int},
    to_bus::AbstractVector{Int},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return _build_branch_flow_vars!(
        model, branch_id, from_bus, to_bus, var_name;
        bounded = bounded, min_lim = min_lim, max_lim = max_lim,
        encoding = encoding, meta = meta,
        export_key_lower = export_key_lower, export_key_upper = export_key_upper,
    )
end

function var_powerflow_ik_reactive!(
    model::Model,
    branch_id::AbstractVector{Int},
    from_bus::AbstractVector{Int},
    to_bus::AbstractVector{Int},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return _build_branch_flow_vars!(
        model, branch_id, from_bus, to_bus, var_name;
        bounded = bounded, min_lim = min_lim, max_lim = max_lim,
        encoding = encoding, meta = meta,
        export_key_lower = export_key_lower, export_key_upper = export_key_upper,
    )
end

function var_powerflow_ki_active!(
    model::Model,
    branch_id::AbstractVector{Int},
    from_bus::AbstractVector{Int},
    to_bus::AbstractVector{Int},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return _build_branch_flow_vars!(
        model, branch_id, from_bus, to_bus, var_name;
        bounded = bounded, min_lim = min_lim, max_lim = max_lim,
        encoding = encoding, meta = meta,
        export_key_lower = export_key_lower, export_key_upper = export_key_upper,
        flip_buses = true,
    )
end

function var_powerflow_ki_reactive!(
    model::Model,
    branch_id::AbstractVector{Int},
    from_bus::AbstractVector{Int},
    to_bus::AbstractVector{Int},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return _build_branch_flow_vars!(
        model, branch_id, from_bus, to_bus, var_name;
        bounded = bounded, min_lim = min_lim, max_lim = max_lim,
        encoding = encoding, meta = meta,
        export_key_lower = export_key_lower, export_key_upper = export_key_upper,
        flip_buses = true,
    )
end
