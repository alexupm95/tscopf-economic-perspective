#=
================================================================================
 functions_4_TS_fullbus_variables.jl — FullBus TS variable builders
================================================================================
=#

# ===================================================================================
# Variable creators — one function per family
# ===================================================================================
# Each creator returns a nested OrderedDict keyed [entity][t] of JuMP variables for the
# given time window. `suffix` ("tf"/"tpf") disambiguates fault-on vs post-fault names, and
# `start=...` seeds each variable with the corresponding steady-state value for warm starts.

"""Per-bus, per-step voltage magnitude variables, seeded at the steady-state `val_V[bus]`."""
function var_fullbus_bus_voltage_magnitude_time!(
    model::JuMP.Model,
    bus_ids::Vector{Int},
    time_window::Vector{Float64},
    val_V::Dict;
    suffix::String,
)
    V_bus = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for bus in bus_ids
        V_bus[bus] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            V_bus[bus][t] = JuMP.@variable(model,
                base_name="V_$suffix[$bus,$t]", start=val_V[bus])
        end
    end
    return V_bus
end
"""Per-bus, per-step voltage angle variables, seeded at the steady-state `val_θ[bus]`."""
function var_fullbus_bus_angle_time!(
    model::JuMP.Model,
    bus_ids::Vector{Int},
    time_window::Vector{Float64},
    val_θ::Dict;
    suffix::String,
)
    θ_bus = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for bus in bus_ids
        θ_bus[bus] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            θ_bus[bus][t] = JuMP.@variable(model,
                base_name="θ_$suffix[$bus,$t]", start=val_θ[bus])
        end
    end
    return θ_bus
end

"""Per-generator, per-step active electrical power, seeded at the dispatched `val_Pg[gen]`."""
function var_fullbus_gen_Pe_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    val_Pg::Dict;
    suffix::String,
)
    Pe = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for gen in active_gen
        Pe[gen] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            Pe[gen][t] = JuMP.@variable(model,
                base_name="Pe_$suffix[$gen,$t]", start=val_Pg[gen])
        end
    end
    return Pe
end

"""Per-generator, per-step reactive electrical power, seeded at the dispatched `val_Qg[gen]`."""
function var_fullbus_gen_Qe_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    val_Qg::Dict;
    suffix::String,
)
    Qe = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for gen in active_gen
        Qe[gen] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            Qe[gen][t] = JuMP.@variable(model,
                base_name="Qe_$suffix[$gen,$t]", start=val_Qg[gen])
        end
    end
    return Qe
end

"""Per-generator, per-step rotor angle, seeded from the reference angle `δ_ref[gen]`
(the pre-fault δ variable; its `start_value` is read via `_opf_scalar_hint`)."""
function var_fullbus_gen_rotor_angle_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    δ_ref::OrderedDict{Int, JuMP.VariableRef};
    suffix::String,
)
    δ = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for gen in active_gen
        δ[gen] = OrderedDict{Int, JuMP.VariableRef}()
        d0 = _opf_scalar_hint(δ_ref[gen], 0.0)  # reuse the pre-fault δ hint as the seed
        for t in eachindex(time_window)
            δ[gen][t] = JuMP.@variable(model,
                base_name="δ_$suffix[$gen,$t]", start=d0)
        end
    end
    return δ
end

"""Per-generator, per-step speed deviation Δω, seeded at 0 (synchronous speed)."""
function var_fullbus_gen_speed_dev_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    time_window::Vector{Float64};
    suffix::String,
)
    Δω = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    for gen in active_gen
        Δω[gen] = OrderedDict{Int, JuMP.VariableRef}()
        for t in eachindex(time_window)
            Δω[gen][t] = JuMP.@variable(model,
                base_name="Δω_$suffix[$gen,$t]", start=0.0)
        end
    end
    return Δω
end
