#=
================================================================================
 functions_4_TS_fullbus_ineqconst.jl — FullBus TS inequality builders
================================================================================
=#

# The δ-COI `bound_style_δ` dispatchers (`_add_δ_bounds_fault!` / `_postf!`) live in
# functions_4_TS_kron_ineqconst.jl next to the two constraint families they choose
# between; every network form calls the same pair.

"""Attach forced lower bounds on bus voltage magnitudes (FullBus windows)."""
function attach_fullbus_V_lower_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    V_bus,
    V_min::Float64,
    manifest_key::Symbol,
)
    meta = dyn_model_dict[:meta]
    if bound_encoding_from_meta(meta) == VARIABLE
        for bus in keys(V_bus)
            for t in keys(V_bus[bus])
                JuMP.set_lower_bound(V_bus[bus][t], V_min)
            end
        end
        register_bound_manifest!(meta, manifest_key, :lower, V_bus, :per_bus_time)
        return nothing
    end
    dyn_model_dict[:ineq_const][manifest_key] =
        ineq_const_fullbus_V_lower!(model, V_bus, V_min)
    return nothing
end

"""Explicit ≤-form lower bounds on bus voltage magnitude: V ≥ min_lim."""
function ineq_const_fullbus_V_lower!(
    model::JuMP.Model,
    V_bus::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    min_lim::Float64=0.0,
)
    lower = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for bus in keys(V_bus)
        lower[bus] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in keys(V_bus[bus])
            lower[bus][t] = JuMP.@constraint(model, min_lim - V_bus[bus][t] ≤ 0.0)
        end
    end
    return lower
end
