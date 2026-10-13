#=
================================================================================
 functions_4_TS_kron_variables.jl — Kron TS variable builders
================================================================================
=#

# ===================================================================================
#                           Decision Variables (shared)
# ===================================================================================

function var_kron_gen_rotor_angle!(
    model::JuMP.Model,
    gen_ids::Vector{Int64},
    var_name::String;
    bounded::Bool=true,
    min_lim::Real=-π,
    max_lim::Real=π,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(gen_ids)), var_name;
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

function var_kron_gen_mech_power!(
    model::JuMP.Model,
    gen_ids::Vector{Int64},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(gen_ids)), var_name;
        bounded = bounded,
        min_lim = min_lim,
        max_lim = max_lim,
        encoding = encoding,
        meta = meta,
        export_key_lower = export_key_lower,
        export_key_upper = export_key_upper,
    )
end

# Function to create generic variables for generators for the fault and post-fault period (variables without bounds)
function var_kron_gen_time_generic!(
    model::JuMP.Model,
    gen_ids::Vector{Int64},
    var_name::String,
    time_window::Vector{Float64},
)
    var_dict = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()

    for gen in gen_ids
        var_dict[gen] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            var_dict[gen][t] = JuMP.@variable(model, base_name = var_name * "[$gen, $t]")
        end
    end

    return var_dict
end

# Function to create generic variables for COI (variables without bounds)
function var_kron_COI_time_generic!(
    model::JuMP.Model,
    var_name::String,
    time_window::Vector{Float64},
)
    var_dict = OrderedDict{Int, JuMP.VariableRef}()

    for t in eachindex(time_window)
        var_dict[t] = JuMP.@variable(model, base_name = var_name * "[$t]")
    end

    return var_dict
end

function var_tsred_gen_voltage_magnitude!(
    model::JuMP.Model,
    gen_ids::Vector{Int64},
    var_name::String;
    bounded::Bool=false,
    min_lim=nothing,
    max_lim=nothing,
    encoding::BoundEncoding=CONSTRAINT,
    meta=nothing,
    export_key_lower::Union{Nothing, Symbol}=nothing,
    export_key_upper::Union{Nothing, Symbol}=nothing,
)
    return build_indexed_scalar_vars!(
        model, collect(Int.(gen_ids)), var_name;
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
