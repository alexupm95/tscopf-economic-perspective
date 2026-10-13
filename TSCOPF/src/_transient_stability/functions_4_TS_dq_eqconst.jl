#=
================================================================================
 functions_4_TS_dq_eqconst.jl — DQ 4th-order machine equality builders
================================================================================
 1. Pre-fault steady state (t=0): subtransient emf defs + stator + terminal P/Q
 2. Transient algebra: Te, Vd/Vq, Pe/Qe (optional (1+Δω) via `speed_in_algebra`)
 3. Nodal KCL: generator injection from Id/Iq (not classical Pe on Xd′)
 4. EMF ODEs: trapezoidal Ed/Eq; constant or AVR-driven E_fd
================================================================================
=#

# ===================================================================================
# Pre-fault steady-state (t = 0)
# ===================================================================================

"""
Pre-fault dq steady state at t = 0, coupled to the dispatch (V, θ, P_g, Q_g).

Pins Ed, Eq, Id, Iq, E_fd, and δ so the machine sits on the solved ACOPF point.
No (1+Δω) here — pre-fault is synchronous (Δω = 0).

Returns the `NamedTuple` `(; eq_Ed, eq_Eq, eq_Vd, eq_Vq, eq_P, eq_Q, Vd_init, Vq_init)`; the
last two are the t = 0 terminal dq projections, kept for export in the same spirit as the
transient `Vd`/`Vq` expressions (see `eq_const_dq_machine_algebra!`). Named fields rather
than a bare 8-tuple: two adjacent `ConstraintRef` dictionaries are interchangeable to the
compiler, so a positional slip here would build a silently wrong model instead of erroring.
"""
function eq_const_dq_init_steady_state!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef},
    E_fd::OrderedDict{Int, JuMP.VariableRef},
    δ::OrderedDict{Int, JuMP.VariableRef},
    Ed::OrderedDict{Int, JuMP.VariableRef},
    Eq::OrderedDict{Int, JuMP.VariableRef},
    Id::OrderedDict{Int, JuMP.VariableRef},
    Iq::OrderedDict{Int, JuMP.VariableRef},
)
    eq_Ed = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Eq = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Vd = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Vq = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_P = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Q = OrderedDict{Int, JuMP.ConstraintRef}()
    Vd_init = OrderedDict{Int, Any}()
    Vq_init = OrderedDict{Int, Any}()
    for gen in active_gen
        bus = DGEN.bus[gen]
        Xd = DGEN_DYN.Xd[gen]
        Xq = DGEN_DYN.Xq[gen]
        Xd_tr = DGEN_DYN.Xd_tr[gen]
        Xq_tr = DGEN_DYN.Xq_tr[gen]
        Ra = DGEN_DYN.Ra[gen]
        eq_Ed[gen] = JuMP.@constraint(model, Ed[gen] == (Xq - Xq_tr) * Iq[gen])
        eq_Eq[gen] = JuMP.@constraint(model, Eq[gen] + (Xd - Xd_tr) * Id[gen] == E_fd[gen])
        eq_Vd[gen] = JuMP.@constraint(model,
            V[bus] * sin(δ[gen] - θ[bus]) - Ed[gen] + Ra * Id[gen] - Xq_tr * Iq[gen] == 0)
        eq_Vq[gen] = JuMP.@constraint(model,
            V[bus] * cos(δ[gen] - θ[bus]) - Eq[gen] + Ra * Iq[gen] + Xd_tr * Id[gen] == 0)
        eq_P[gen] = JuMP.@constraint(model,
            P_g[gen] == V[bus] * sin(δ[gen] - θ[bus]) * Id[gen] +
                V[bus] * cos(δ[gen] - θ[bus]) * Iq[gen])
        eq_Q[gen] = JuMP.@constraint(model,
            Q_g[gen] == V[bus] * cos(δ[gen] - θ[bus]) * Id[gen] -
                V[bus] * sin(δ[gen] - θ[bus]) * Iq[gen])
        Vd_init[gen] = JuMP.@expression(model, V[bus] * sin(δ[gen] - θ[bus]))
        Vq_init[gen] = JuMP.@expression(model, V[bus] * cos(δ[gen] - θ[bus]))
    end
    return (; eq_Ed, eq_Eq, eq_Vd, eq_Vq, eq_P, eq_Q, Vd_init, Vq_init)
end

# ===================================================================================
# Transient dq algebra (stator + Te/Pe/Qe)
# ===================================================================================

"""
Transient dq stator algebra and electrical power at each time step.

`Te = Ed·Id + Eq·Iq`; `Pe = rot·Te` with `rot = dq_rotor_scale(Δω, speed_in_algebra)`.
Stator Vd/Vq links use the same `rot` factor on Ed/Eq when speed is included.

Returns the `NamedTuple` `(; eq_Pe, eq_Qe, eq_Vd, eq_Vq, eq_Te, Vd_out, Vq_out)`. The last two are the
terminal-voltage projections `V·sin(δ−θ)` / `V·cos(δ−θ)` per (gen, t): they were already
built here and thrown away, so the machine's own dq voltages — the axis the whole 4th-order
model is written in — could not be exported. Collecting them costs nothing (the same
`@expression` objects, kept rather than discarded) and adds no constraint to the solve.
"""
function eq_const_dq_machine_algebra!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    Pe::OrderedDict,
    Qe::OrderedDict,
    δ::OrderedDict,
    Δω::OrderedDict,
    Ed::OrderedDict,
    Eq::OrderedDict,
    Id::OrderedDict,
    Iq::OrderedDict,
    Te::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64};
    speed_in_algebra::Bool=true,
)
    eq_Pe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_Qe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_Vd = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_Vq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_Te = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    Vd_out = OrderedDict{Int, OrderedDict{Int, Any}}()
    Vq_out = OrderedDict{Int, OrderedDict{Int, Any}}()
    for gen in active_gen
        bus_g = DGEN.bus[gen]
        Xd_tr = DGEN_DYN.Xd_tr[gen]
        Xq_tr = DGEN_DYN.Xq_tr[gen]
        Ra = DGEN_DYN.Ra[gen]
        eq_Pe[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_Qe[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_Vd[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_Vq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_Te[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        Vd_out[gen] = OrderedDict{Int, Any}()
        Vq_out[gen] = OrderedDict{Int, Any}()
        for t in eachindex(time_window)
            Vd_expr = JuMP.@expression(model, V_t[bus_g][t] * sin(δ[gen][t] - θ_t[bus_g][t]))
            Vq_expr = JuMP.@expression(model, V_t[bus_g][t] * cos(δ[gen][t] - θ_t[bus_g][t]))
            Vd_out[gen][t] = Vd_expr
            Vq_out[gen][t] = Vq_expr
            rot = dq_rotor_scale(Δω[gen][t], speed_in_algebra)
            eq_Te[gen][t] = JuMP.@constraint(model, Te[gen][t] == Ed[gen][t] * Id[gen][t] + Eq[gen][t] * Iq[gen][t])
            eq_Vd[gen][t] = JuMP.@constraint(model,
                Vd_expr - rot * Ed[gen][t] + Ra * Id[gen][t] - Xq_tr * Iq[gen][t] == 0)
            eq_Vq[gen][t] = JuMP.@constraint(model,
                Vq_expr - rot * Eq[gen][t] + Ra * Iq[gen][t] + Xd_tr * Id[gen][t] == 0)
            eq_Pe[gen][t] = JuMP.@constraint(model, Pe[gen][t] == rot * Te[gen][t])
            eq_Qe[gen][t] = JuMP.@constraint(model,
                Qe[gen][t] == Vq_expr * Id[gen][t] - Vd_expr * Iq[gen][t])
        end
    end
    return (; eq_Pe, eq_Qe, eq_Vd, eq_Vq, eq_Te, Vd_out, Vq_out)
end

"""
Nodal KCL with dq current injection (replaces classical ΣPe balance).

`ZIP` is the **active-demand** (Z, I, P) split, one row per bus (`nBUS × 3`, from
`zip_by_bus`: `DynModelConfig.zip_load_p` plus the optional `load_data.csv`);
`eq_const_dq_Qbalance!` takes its own matrix, so P and Q need not share a load model.
"""
function eq_const_dq_Pbalance!(
    model::JuMP.Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    δ::OrderedDict,
    Id::OrderedDict,
    Iq::OrderedDict,
    terms_Pb::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64},
    ZIP::AbstractMatrix{Float64},
)
    eq_Pb = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for i in 1:nBUS
        gen_ids = bus_gen_circ_dict[i][:gen_ids]
        eq_Pb[i] = OrderedDict{Int, JuMP.ConstraintRef}()
        p_d = DBUS.p_d[i] / base_MVA
        for t in eachindex(time_window)
            if isempty(gen_ids)
                if isapprox(p_d, 0.0; atol=1e-12)
                    eq_Pb[i][t] = JuMP.@constraint(model, terms_Pb[i][t] == 0.0)
                else
                    eq_Pb[i][t] = JuMP.@constraint(model,
                        terms_Pb[i][t] * V[i]^2 ==
                        -p_d * (ZIP[i, 1] * V_t[i][t]^2 + ZIP[i, 2] * V_t[i][t] * V[i] + ZIP[i, 3] * V[i]^2))
                end
            else
                inj_P = JuMP.@expression(model, sum(
                    V_t[i][t] * Id[g][t] * sin(δ[g][t] - θ_t[i][t]) +
                    V_t[i][t] * Iq[g][t] * cos(δ[g][t] - θ_t[i][t]) for g in gen_ids))
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

"""Reactive counterpart; `ZIP` is the per-bus **reactive-demand** split (`zip_load_q`)."""
function eq_const_dq_Qbalance!(
    model::JuMP.Model,
    DBUS::DataFrame,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    δ::OrderedDict,
    Id::OrderedDict,
    Iq::OrderedDict,
    terms_Qb::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64},
    ZIP::AbstractMatrix{Float64},
)
    eq_Qb = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for i in 1:nBUS
        gen_ids = bus_gen_circ_dict[i][:gen_ids]
        eq_Qb[i] = OrderedDict{Int, JuMP.ConstraintRef}()
        q_d = DBUS.q_d[i] / base_MVA
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
                inj_Q = JuMP.@expression(model, sum(
                    V_t[i][t] * Id[g][t] * cos(δ[g][t] - θ_t[i][t]) -
                    V_t[i][t] * Iq[g][t] * sin(δ[g][t] - θ_t[i][t]) for g in gen_ids))
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
# EMF dynamics (first step: trap or BE via ode_first_step; trap thereafter)
# ===================================================================================

"""
Ed / Eq dynamics. First step of each window uses `ode_first_step`
(`:trapezoidal` default, or `:backward_euler`). Remaining steps use trapezoidal.
With `E_fd_time=nothing`, `E_fd` is constant (no AVR). With AVR,
pass `E_fd_tf`/`E_fd_tpf` via `E_fd_time` and anchor `E_fd_prev0` at the
window boundary.
"""
function eq_const_dq_emf_dynamics!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    E_fd::OrderedDict{Int, JuMP.VariableRef},
    Ed_prev::OrderedDict,
    Eq_prev::OrderedDict,
    Id_prev::OrderedDict,
    Iq_prev::OrderedDict,
    Ed::OrderedDict,
    Eq::OrderedDict,
    Id::OrderedDict,
    Iq::OrderedDict,
    time_window::Vector{Float64},
    Δt::Float64;
    anchor_first_step::Bool=true,
    E_fd_time::Union{Nothing, OrderedDict}=nothing,
    E_fd_prev0::Union{Nothing, OrderedDict}=nothing,
    ode_first_step::Symbol=:trapezoidal,
)
    eq_Ed = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    eq_Eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        Xd = DGEN_DYN.Xd[gen]
        Xd_tr = DGEN_DYN.Xd_tr[gen]
        Xq = DGEN_DYN.Xq[gen]
        Xq_tr = DGEN_DYN.Xq_tr[gen]
        Td = DGEN_DYN.Td[gen]
        Tq = DGEN_DYN.Tq[gen]
        eq_Ed[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        eq_Eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            use_be = t == 1 && anchor_first_step && ode_first_step === :backward_euler
            if use_be
                # --- BACKWARD EULER FOR STEP 1 (reference GFM path) ---
                # Ignores discontinuous previous Id/Iq (and previous E_fd).
                Ed_p = Ed_prev[gen]
                Eq_p = Eq_prev[gen]
                eq_Ed[gen][t] = JuMP.@constraint(model,
                    Ed[gen][t] * (1 + Δt / Tq) - Ed_p -
                    (Δt / Tq) * (Xq - Xq_tr) * Iq[gen][t] == 0.0)
                E_fd_term = E_fd_time === nothing ? E_fd[gen] : E_fd_time[gen][t]
                eq_Eq[gen][t] = JuMP.@constraint(model,
                    Eq[gen][t] * (1 + Δt / Td) - Eq_p -
                    (Δt / Td) * (E_fd_term - (Xd - Xd_tr) * Id[gen][t]) == 0.0)
            else
                # --- TRAPEZOIDAL (t=1 with anchors, or remainder) ---
                if t == 1 && anchor_first_step
                    Ed_p = Ed_prev[gen]
                    Eq_p = Eq_prev[gen]
                    Id_p = Id_prev[gen]
                    Iq_p = Iq_prev[gen]
                elseif t == 1
                    Ed_p = Ed[gen][t]
                    Eq_p = Eq[gen][t]
                    Id_p = Id[gen][t]
                    Iq_p = Iq[gen][t]
                else
                    Ed_p = Ed[gen][t - 1]
                    Eq_p = Eq[gen][t - 1]
                    Id_p = Id[gen][t - 1]
                    Iq_p = Iq[gen][t - 1]
                end
                eq_Ed[gen][t] = JuMP.@constraint(model,
                    Ed[gen][t] * (1 + Δt / (2 * Tq)) - Ed_p * (1 - Δt / (2 * Tq)) -
                    (Δt / (2 * Tq)) * (Xq - Xq_tr) * (Iq[gen][t] + Iq_p) == 0.0)
                if E_fd_time === nothing
                    E_fd_term = 2 * E_fd[gen]
                else
                    E_fd_curr = E_fd_time[gen][t]
                    E_fd_p = (t == 1) ? E_fd_prev0[gen] : E_fd_time[gen][t - 1]
                    E_fd_term = E_fd_curr + E_fd_p
                end
                eq_Eq[gen][t] = JuMP.@constraint(model,
                    Eq[gen][t] * (1 + Δt / (2 * Td)) - Eq_p * (1 - Δt / (2 * Td)) -
                    (Δt / (2 * Td)) * (E_fd_term - (Xd - Xd_tr) * (Id[gen][t] + Id_p)) == 0.0)
            end
        end
    end
    return eq_Ed, eq_Eq
end
