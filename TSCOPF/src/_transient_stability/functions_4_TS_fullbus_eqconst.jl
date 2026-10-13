#=
================================================================================
 functions_4_TS_fullbus_eqconst.jl — FullBus TS equality builders
================================================================================
=#

# ===================================================================================
# Network equality constraints — nodal KCL + classical Pe/Qe
# ===================================================================================

"""
Classical-generator active electrical power behind the transient reactance Xd_tr:

    Pe_g = (E_g · V_busg / Xd_tr) · sin(δ_g − θ_busg)

i.e. the power transferred from the internal EMF node (E∠δ) to the terminal bus (V∠θ)
across Xd_tr, for every active generator and time step.
"""
function eq_const_fullbus_gen_Pe!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    E::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    δ::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64},
)
    eq_Pe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for g in active_gen
        eq_Pe[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        inv_Xd = 1.0 / DGEN_DYN.Xd_tr[g]  # 1/Xd_tr = transfer susceptance behind the machine
        bus_g = DGEN.bus[g]
        for t in eachindex(time_window)
            eq_Pe[g][t] = JuMP.@constraint(model,
                Pe[g][t] == inv_Xd * E[g] * V_t[bus_g][t] * sin(δ[g][t] - θ_t[bus_g][t]))
        end
    end
    return eq_Pe
end

"""
Classical-generator reactive electrical power behind Xd_tr:

    Qe_g = (E_g · V_busg / Xd_tr) · cos(δ_g − θ_busg) − V_busg² / Xd_tr

(the standard E-behind-reactance reactive injection), for every active gen and time step.
"""
function eq_const_fullbus_gen_Qe!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    E::OrderedDict{Int, JuMP.VariableRef},
    Qe::OrderedDict,
    δ::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64},
)
    eq_Qe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for g in active_gen
        eq_Qe[g] = OrderedDict{Int, JuMP.ConstraintRef}()
        inv_Xd = 1.0 / DGEN_DYN.Xd_tr[g]
        bus_g = DGEN.bus[g]
        for t in eachindex(time_window)
            eq_Qe[g][t] = JuMP.@constraint(model,
                Qe[g][t] == inv_Xd * E[g] * V_t[bus_g][t] * cos(δ[g][t] - θ_t[bus_g][t])
                - inv_Xd * V_t[bus_g][t]^2)
        end
    end
    return eq_Qe
end

"""
Active-power nodal balance (KCL) at every bus and time step:

    network injection (terms_Pb) = generator injection (ΣPe) − ZIP load demand.

`ZIP` here is the **active-demand** split, one row per bus (`nBUS × 3`, built by
`zip_by_bus` from `DynModelConfig.zip_load_p` and the optional `load_data.csv`); the
reactive balance below receives its own matrix (`zip_load_q`), so P and Q need not share
a load model. Coefficient order (Z, I, P): ZIP[i,1]=impedance, ZIP[i,2]=current,
ZIP[i,3]=power (each row sums to 1). The balance is scaled through by V[i]² (the
steady-state / nominal voltage) to avoid dividing by V_t, so the load term
p_d·[Z·(V_t/V_i)² + I·(V_t/V_i) + P]·V_i² expands to
ZIP[i,1]·V_t² + ZIP[i,2]·V_t·V_i + ZIP[i,3]·V_i². Buses with zero load skip the V[i]²
scaling (the special-cased `== inj_P` / `== 0` forms) to keep the constraint simpler and
better conditioned.
"""
function eq_const_fullbus_Pbalance!(
    model::JuMP.Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    terms_Pb::OrderedDict,
    V_t::OrderedDict,
    time_window::Vector{Float64},
    ZIP::AbstractMatrix{Float64},
)
    eq_Pb = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for i in 1:nBUS
        gen_ids = bus_gen_circ_dict[i][:gen_ids]
        eq_Pb[i] = OrderedDict{Int, JuMP.ConstraintRef}()
        p_d = DBUS.p_d[i] / base_MVA  # per-unit active load at this bus
        for t in eachindex(time_window)
            if isempty(gen_ids)
                # Load-only / passive bus: injection equals (negative) ZIP demand.
                if isapprox(p_d, 0.0; atol=1e-12)
                    eq_Pb[i][t] = JuMP.@constraint(model, terms_Pb[i][t] == 0.0)
                else
                    eq_Pb[i][t] = JuMP.@constraint(model,
                        terms_Pb[i][t] * V[i]^2 ==
                        -p_d * (ZIP[i, 1] * V_t[i][t]^2 + ZIP[i, 2] * V_t[i][t] * V[i] + ZIP[i, 3] * V[i]^2))
                end
            else
                # Generator bus: network injection = generator injection − ZIP demand.
                inj_P = sum(Pe[g][t] for g in gen_ids)
                if isapprox(p_d, 0.0; atol=1e-12)
                    eq_Pb[i][t] = JuMP.@constraint(model, terms_Pb[i][t] == inj_P)
                else
                    eq_Pb[i][t] = JuMP.@constraint(model,
                        terms_Pb[i][t] * V[i]^2 == inj_P * V[i]^2 -
                        p_d * (ZIP[i, 1] * V_t[i][t]^2 + ZIP[i, 2] * V_t[i][t] * V[i] + ZIP[i, 3] * V[i]^2))
                end
            end
        end
    end
    return eq_Pb
end

"""
Reactive-power nodal balance — the Q-counterpart of `eq_const_fullbus_Pbalance!`.

Identical structure and V[i]²-scaled ZIP treatment, using the reactive demand q_d and the
reactive injection expressions terms_Qb / Qe. `ZIP` here is the per-bus
**reactive-demand** split (`DynModelConfig.zip_load_q` + `load_data.csv`), independent of
the active one — a load can be a different share of its bus's MVAr than of its MW.
"""
function eq_const_fullbus_Qbalance!(
    model::JuMP.Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    Qe::OrderedDict,
    terms_Qb::OrderedDict,
    V_t::OrderedDict,
    time_window::Vector{Float64},
    ZIP::AbstractMatrix{Float64},
)
    eq_Qb = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for i in 1:nBUS
        gen_ids = bus_gen_circ_dict[i][:gen_ids]
        eq_Qb[i] = OrderedDict{Int, JuMP.ConstraintRef}()
        q_d = DBUS.q_d[i] / base_MVA  # per-unit reactive load at this bus
        for t in eachindex(time_window)
            if isempty(gen_ids)
                if isapprox(q_d, 0.0; atol=1e-12)
                    eq_Qb[i][t] = JuMP.@constraint(model, terms_Qb[i][t] == 0.0)
                else
                    eq_Qb[i][t] = JuMP.@constraint(model,
                        terms_Qb[i][t] * V[i]^2 ==
                        -q_d * (ZIP[i, 1] * V_t[i][t]^2 + ZIP[i, 2] * V_t[i][t] * V[i] + ZIP[i, 3] * V[i]^2))
                end
            else
                inj_Q = sum(Qe[g][t] for g in gen_ids)
                if isapprox(q_d, 0.0; atol=1e-12)
                    eq_Qb[i][t] = JuMP.@constraint(model, terms_Qb[i][t] == inj_Q)
                else
                    eq_Qb[i][t] = JuMP.@constraint(model,
                        terms_Qb[i][t] * V[i]^2 == inj_Q * V[i]^2 -
                        q_d * (ZIP[i, 1] * V_t[i][t]^2 + ZIP[i, 2] * V_t[i][t] * V[i] + ZIP[i, 3] * V[i]^2))
                end
            end
        end
    end
    return eq_Qb
end

# ===================================================================================
# Swing equations (trapezoidal; optional BE at t=1 via ode_first_step)
# ===================================================================================
# Classical 2nd-order swing model per generator:
#   dδ/dt = ω_syn · Δω
#   2H · dΔω/dt = P_m − P_e − D · Δω
# Default `:trapezoidal` uses x_t − x_{t-1} = (Δt/2)(f_t + f_{t-1}) at every step.
# With `ode_first_step=:backward_euler`, the first step of each fault/post-fault window
# matches the reference implementation (BE; remainder stays trapezoidal). BE at t=1 uses
# only current P_mech and Pe (no average with P_g / Pe_ant / P_mech₀).

"""
Rotor-angle update for the fault-on window.

Trapezoidal: δ_t − δ_{t-1} = ω_syn·(Δt/2)·(Δω_t + Δω_{t-1}).
Backward Euler at t=1 (`ode_first_step=:backward_euler`): δ_t − δ_0 = ω_syn·Δt·Δω_t.
"""
function eq_const_fullbus_δ_swing_fault!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    δ::OrderedDict,
    Δω::OrderedDict,
    δ_0::OrderedDict,
    Δω_0::Float64,
    time_window::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_δ = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        eq_δ[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ_0[gen] - ω_syn * Δt * Δω[gen][t] == 0.0)
            elseif t == 1
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ_0[gen] - ω_syn * (Δt / 2) * (Δω[gen][t] + Δω_0) == 0.0)
            else
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ[gen][t - 1] - ω_syn * (Δt / 2) * (Δω[gen][t] + Δω[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_δ
end

"""
Speed-deviation update for the fault-on window (constant mechanical power).

Trapezoidal:
    (1 + DΔt/4H)·Δω_t − (1 − DΔt/4H)·Δω_{t-1} − (Δt/4H)·(2·P_m − P_e,t − P_e,{t-1}) = 0.
Backward Euler at t=1:
    (1 + DΔt/2H)·Δω_t − Δω_0 − (Δt/2H)·(P_m − P_e,t) = 0.
"""
function eq_const_fullbus_Δω_swing_fault!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_mech::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    Δω::OrderedDict,
    Δω_0::Float64,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        H, D = DGEN_DYN.H[gen], DGEN_DYN.D[gen]  # inertia constant and damping
        eq_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (2H)) * Δω[gen][t] - Δω_0
                    - (Δt / (2H)) * (P_mech[gen] - Pe[gen][t]) == 0.0)
            elseif t == 1
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω_0
                    - (Δt / (4H)) * (2 * P_mech[gen] - Pe[gen][t] - P_g[gen]) == 0.0)
            else
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω[gen][t - 1]
                    - (Δt / (4H)) * (2 * P_mech[gen] - Pe[gen][t] - Pe[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_Δω
end

"""
Post-fault rotor-angle swing update; t=1 links to the last fault-on step (δ_ant, Δω_ant).

Same trapezoidal / optional-BE-at-t=1 split as the fault-on angle update.
"""
function eq_const_fullbus_δ_swing_postf!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    δ::OrderedDict,
    Δω::OrderedDict,
    δ_ant::OrderedDict{Int64, JuMP.VariableRef},
    Δω_ant::OrderedDict{Int64, JuMP.VariableRef},
    time_window::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_δ = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        eq_δ[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ_ant[gen] - ω_syn * Δt * Δω[gen][t] == 0.0)
            elseif t == 1
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ_ant[gen] - ω_syn * (Δt / 2) * (Δω[gen][t] + Δω_ant[gen]) == 0.0)
            else
                eq_δ[gen][t] = JuMP.@constraint(model,
                    δ[gen][t] - δ[gen][t - 1] - ω_syn * (Δt / 2) * (Δω[gen][t] + Δω[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_δ
end

"""
Post-fault speed-deviation swing update (constant mechanical power).

Trapezoidal t=1 uses Δω_ant and Pe_ant; BE t=1 uses only current Pe (reference-style).
"""
function eq_const_fullbus_Δω_swing_postf!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_mech::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    Δω::OrderedDict,
    Δω_ant::OrderedDict{Int64, JuMP.VariableRef},
    Pe_ant::OrderedDict{Int64, JuMP.VariableRef},
    Pe_tf::OrderedDict,
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        H, D = DGEN_DYN.H[gen], DGEN_DYN.D[gen]  # inertia constant and damping
        eq_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (2H)) * Δω[gen][t] - Δω_ant[gen]
                    - (Δt / (2H)) * (P_mech[gen] - Pe[gen][t]) == 0.0)
            elseif t == 1
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω_ant[gen]
                    - (Δt / (4H)) * (2 * P_mech[gen] - Pe[gen][t] - Pe_ant[gen]) == 0.0)
            else
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω[gen][t - 1]
                    - (Δt / (4H)) * (2 * P_mech[gen] - Pe[gen][t] - Pe[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_Δω
end

# ===================================================================================
# Swing equations with a TIME-VARYING mechanical power (turbine governor)
# ===================================================================================
# When `include_governor` is on, the governor produces a per-step mechanical power
# `P_mech[gen][t]` instead of the constant `P_m`. The trapezoidal average that the
# constant-P_m methods collapse into `2*P_mech[gen]` becomes the genuine pair
# `(P_mech[gen][t] + P_mech[gen][t-1])`, with the t=1 step using `P_mech₀` (the
# pre-fault equilibrium P_m for the fault window, or the last fault-on P_mech for the
# post-fault window). BE at t=1 uses only current P_mech[t] and Pe[t].

"""Fault-on Δω swing with a governor-driven, time-varying mechanical power.

`P_mech` is per (gen, t); `P_mech₀` is the pre-fault equilibrium mechanical power
supplying the trapezoidal t=1 previous step (alongside `P_g` for Pe and `Δω_0`)."""
function eq_const_fullbus_Δω_swing_fault!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_mech::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    P_mech₀::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    Δω::OrderedDict,
    Δω_0::Float64,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        H, D = DGEN_DYN.H[gen], DGEN_DYN.D[gen]
        eq_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (2H)) * Δω[gen][t] - Δω_0
                    - (Δt / (2H)) * (P_mech[gen][t] - Pe[gen][t]) == 0.0)
            elseif t == 1
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω_0
                    - (Δt / (4H)) * ((P_mech[gen][t] + P_mech₀[gen]) - Pe[gen][t] - P_g[gen]) == 0.0)
            else
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω[gen][t - 1]
                    - (Δt / (4H)) * ((P_mech[gen][t] + P_mech[gen][t - 1]) - Pe[gen][t] - Pe[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_Δω
end

"""Post-fault Δω swing with a governor-driven, time-varying mechanical power.

Trapezoidal t=1 anchors: Δω_ant, Pe_ant, P_mech₀. BE t=1: current P_mech and Pe only."""
function eq_const_fullbus_Δω_swing_postf!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_mech::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    P_mech₀::OrderedDict{Int, JuMP.VariableRef},
    Pe::OrderedDict,
    Δω::OrderedDict,
    Δω_ant::OrderedDict{Int64, JuMP.VariableRef},
    Pe_ant::OrderedDict{Int64, JuMP.VariableRef},
    Pe_tf::OrderedDict,
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        H, D = DGEN_DYN.H[gen], DGEN_DYN.D[gen]
        eq_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            if t == 1 && ode_first_step === :backward_euler
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (2H)) * Δω[gen][t] - Δω_ant[gen]
                    - (Δt / (2H)) * (P_mech[gen][t] - Pe[gen][t]) == 0.0)
            elseif t == 1
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω_ant[gen]
                    - (Δt / (4H)) * ((P_mech[gen][t] + P_mech₀[gen]) - Pe[gen][t] - Pe_ant[gen]) == 0.0)
            else
                eq_Δω[gen][t] = JuMP.@constraint(model,
                    (1 + D * Δt / (4H)) * Δω[gen][t] - (1 - D * Δt / (4H)) * Δω[gen][t - 1]
                    - (Δt / (4H)) * ((P_mech[gen][t] + P_mech[gen][t - 1]) - Pe[gen][t] - Pe[gen][t - 1]) == 0.0)
            end
        end
    end
    return eq_Δω
end
