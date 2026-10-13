#=
================================================================================
 functions_4_TS_governor.jl — optional TGOV1 turbine-governor control layer
================================================================================
 Network-agnostic primary frequency control for the classical 2nd-order machine.
 The governor sees only the per-generator speed deviation Δω and produces a
 *time-varying* mechanical power P_mech that replaces the constant P_m inside the
 swing equation.  Currently wired into the FULL_BUS path only (it relies on the
 mandatory ACOPF warm start).  The AVR layer (later PR) will reuse the same
 variable / discretization conventions.

 Governor model (IEEE TGOV1, matching the DGEN_DYN columns R, T1, T2, T3):

     Valve:  T1·dPv/dt = (P_ref − Δω)/R − Pv
     Mech:   T3·dPm/dt + Pm = T2·dPv/dt + Pv          (turbine lead-lag, on Pv alone)

 The mech row is the *direct* discretization of the turbine transfer function
 (1 + s·T2)/(1 + s·T3) acting on the valve output.  It used to be written in the
 substituted state-space form obtained by eliminating dPv/dt with the valve ODE,

     T3·dPm/dt = (1 − T2/T1)·Pv + (T2/T1)·(P_ref − Δω)/R − Pm,          [removed]

 which is algebraically exact only while the valve is *unsaturated*: once a limiter
 clamps Pv, dPv/dt is no longer (u − Pv)/T1 and the (T2/T1)·u term smuggles the raw,
 unclamped droop signal past the limiter into the turbine.  The direct form above
 touches nothing but the limited Pv, so saturation propagates correctly.  For
 GOV_NO_LIMIT / GOV_HARD_BOUND (where Pv ≡ Pv_raw and the valve equality still holds)
 the two forms are identical row-by-row, trapezoidal and backward Euler alike.

 Initialization at the Δω = 0 equilibrium (P_m = P_g) gives, consistently,
     P_ref = R·P_m,   Pv₀ = P_m,   Pm₀ = P_m.

 Discretization is trapezoidal (module convention), the first step of each window
 anchored to the previous-window equilibrium/state.  If Ipopt struggles, the t=1
 step can be switched to backward Euler (the reference form) — not done here.

 Valve saturation is selectable via `GovernorLimiter`:
   GOV_NO_LIMIT   — Pv is the raw ODE state, unbounded.
   GOV_SMOOTH     — Pv_unlim is the ODE state; Pv = smooth min/max clamp (nonconvex).
   GOV_HARD_BOUND — Pv is the ODE state + explicit ≤-form bounds (can be infeasible).
 Only the valve treatment differs between modes; the mech ODE always consumes the
 (possibly-limited) valve output Pv, and nothing else.

 The valve integrator itself is **anti-windup**: its previous-step term is the limited
 output Pv, not the raw state Pv_raw (the AVR exciter has always used this convention —
 see `eq_const_avr_exciter!`).  With raw feedback the state is a free integrator, so
 under GOV_SMOOTH it ramps past the clamp for as long as the droop signal demands and
 must unwind that excess before the output can leave saturation.  Feeding the clamped
 value back caps the raw state one step's travel beyond the limit instead.  For
 GOV_NO_LIMIT / GOV_HARD_BOUND (both encodings) Pv *is* Pv_raw — the same container —
 so those paths are unchanged row-for-row.
================================================================================
=#

# `_GOV_SMOOTH_RHO` (GOV_SMOOTH sqrt min/max smoothing) lives in
# `_transient_stability/numerics.jl`.

"""Index `x[gen]` when `x` is a per-generator dict, else return the scalar `x`.

Lets the ODE builders take the t=1 "previous" value either as a per-generator
container (post-fault window: last fault-on states) or as a shared scalar
(fault window: Δω_prev = 0 at the equilibrium)."""
_gov_at(x::AbstractDict, gen::Int) = x[gen]
_gov_at(x, ::Int) = x

# ===================================================================================
# Variables
# ===================================================================================

"""Constant governor set-point `P_ref[g]` (one scalar per active generator).

Warm-started at `R·P_m` so the pre-fault equilibrium `Pv₀ = Pm₀ = P_m` holds."""
function var_gov_setpoint!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_m::OrderedDict{Int, JuMP.VariableRef},
)
    P_ref = OrderedDict{Int, JuMP.VariableRef}()
    for gen in active_gen
        P_ref[gen] = JuMP.@variable(model, base_name = "P_ref[$gen]")
        pm_seed = _gov_start(P_m[gen])
        JuMP.set_start_value(P_ref[gen], DGEN_DYN.R[gen] * pm_seed)
    end
    return P_ref
end

"""`JuMP.start_value(v)` when set, else `0.0` — used to seed governor states."""
_gov_start(v::JuMP.VariableRef)::Float64 = (s = JuMP.start_value(v); isnothing(s) ? 0.0 : Float64(s))

"""Time-indexed governor state (valve or mech), warm-started flat at the P_m seed."""
function var_gov_state_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    var_name::String,
    time_window::Vector{Float64},
    P_m::OrderedDict{Int, JuMP.VariableRef},
)
    var = var_kron_gen_time_generic!(model, active_gen, var_name, time_window)
    for gen in active_gen
        seed = _gov_start(P_m[gen])
        for t in eachindex(time_window)
            JuMP.set_start_value(var[gen][t], seed)
        end
    end
    return var
end

# ===================================================================================
# Set-point initialization: P_ref = R·P_m
# ===================================================================================

"""Pin the governor set-point to the equilibrium: `P_ref[g] = R·P_m[g]`."""
function eq_const_gov_setpoint_init!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    P_ref::OrderedDict{Int, JuMP.VariableRef},
    P_m::OrderedDict{Int, JuMP.VariableRef},
)
    eq = OrderedDict{Int, JuMP.ConstraintRef}()
    for gen in active_gen
        eq[gen] = JuMP.@constraint(model, P_ref[gen] - DGEN_DYN.R[gen] * P_m[gen] == 0.0)
    end
    return eq
end

# ===================================================================================
# Valve ODE (trapezoidal): T1·dPv/dt = (P_ref − Δω)/R − Pv
# ===================================================================================
# `Pv_raw` is the integrated valve state (= Pv itself for NO_LIMIT/HARD_BOUND, or the
# pre-saturation Pv_unlim for SMOOTH); `Pv_sat` is the limited output that the previous-step
# term feeds back (anti-windup). `Pv_prev0`/`Δω_prev0` supply the t=1 anchor: fault window →
# (P_m, 0.0); post-fault window → (last fault-on *limited* Pv, last fault-on Δω).

"""Valve-state update for one window: first step via `ode_first_step`, then trap.

The integrated state is `Pv_raw`, but the previous-step anchor is the **saturated** output
`Pv_sat` (the exciter convention, see `eq_const_avr_exciter!`) — that is what makes the valve
integrator anti-windup. `Pv_sat === Pv_raw` for every mode except `GOV_SMOOTH`, where the two
differ by the clamp."""
function eq_const_gov_valve!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pv_raw::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    Pv_sat::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    P_ref::OrderedDict{Int, JuMP.VariableRef},
    Δω::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    Pv_prev0::OrderedDict{Int, JuMP.VariableRef},
    Δω_prev0,
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        R, T1 = DGEN_DYN.R[gen], DGEN_DYN.T1[gen]
        c = Δt / (2 * T1)
        eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            Pv_prev = t == 1 ? Pv_prev0[gen] : Pv_sat[gen][t - 1]
            if t == 1 && ode_first_step === :backward_euler
                # --- BACKWARD EULER FOR STEP 1 (reference GFM path) ---
                eq[gen][t] = JuMP.@constraint(model,
                    Pv_raw[gen][t] * (1 + Δt / T1) - Pv_prev
                    - (Δt / (R * T1)) * (P_ref[gen] - Δω[gen][t]) == 0.0)
            elseif t == 1
                # --- TRAPEZOIDAL AT t=1 (SG pin default) ---
                Δω_prev = _gov_at(Δω_prev0, gen)
                eq[gen][t] = JuMP.@constraint(model,
                    Pv_raw[gen][t] * (1 + c) - Pv_prev * (1 - c)
                    - (Δt / (2 * R * T1)) * (2 * P_ref[gen] - Δω[gen][t] - Δω_prev) == 0.0)
            else
                # --- TRAPEZOIDAL FOR REMAINDER ---
                Δω_prev = Δω[gen][t - 1]
                eq[gen][t] = JuMP.@constraint(model,
                    Pv_raw[gen][t] * (1 + c) - Pv_prev * (1 - c)
                    - (Δt / (2 * R * T1)) * (2 * P_ref[gen] - Δω[gen][t] - Δω_prev) == 0.0)
            end
        end
    end
    return eq
end

# ===================================================================================
# Valve saturation — three selectable treatments (see GovernorLimiter)
# ===================================================================================

"""No limiter: the valve output is the raw ODE state itself (`Pv = Pv_raw`)."""
function apply_gov_valve_limit_none!(
    ::JuMP.Model,
    ::Vector{Int64},
    Pv_raw::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    ::Vector{Float64},
    ::Dict{Int, Float64},
    ::Dict{Int, Float64},
)
    return Pv_raw, nothing, nothing  # output = raw; no extra constraints/bounds
end

"""Smooth anti-windup: `Pv = smoothmax(smoothmin(Pv_raw, Pmax), Pmin)` via sqrt.

Introduces a separate limited-output variable `Pv` tied to the raw state `Pv_raw`
(= Pv_unlim) by one equality per (gen, t). Physically saturating but nonconvex."""
function apply_gov_valve_limit_smooth!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    Pv_raw::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    p_min::Dict{Int, Float64},
    p_max::Dict{Int, Float64};
    out_name::String="P_valve",
)
    ρ = _GOV_SMOOTH_RHO
    Pv = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        Pmax = p_max[gen]
        Pmin = p_min[gen]
        Pv[gen] = OrderedDict{Int, JuMP.VariableRef}()
        eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            # Window-specific name (P_valve_tf / P_valve_tpf) so dumps do not collide.
            Pv[gen][t] = JuMP.@variable(model, base_name = "$out_name[$gen, $t]",
                start = _gov_start(Pv_raw[gen][t]))
            u = Pv_raw[gen][t]
            # Smooth min(u, Pmax) then smooth max(·, Pmin).
            P_aux = JuMP.@expression(model, (u + Pmax - sqrt((u - Pmax)^2 + ρ)) / 2)
            eq[gen][t] = JuMP.@constraint(model,
                Pv[gen][t] == (P_aux + Pmin + sqrt((P_aux - Pmin)^2 + ρ)) / 2)
        end
    end
    return Pv, eq, nothing
end

"""Hard bounds: valve output equals the raw state, plus `p_min ≤ Pv ≤ p_max`
(dual-friendly, but over-constrains the ODE and can be infeasible).

Honours `bound_encoding`, like every other bound family: `CONSTRAINT` emits explicit
≤-form rows, `VARIABLE` stamps JuMP bounds on `Pv_raw` itself and leaves the dual export
to the bound manifest (registered by `_store_gov_limit!`, which owns the key names).

Returns `(Pv_output, extra_eq, extra_ineq)` with `extra_ineq` either the
`(lower, upper)` constraint pair or the `:variable_bounds` marker — the third slot of
the contract shared with the `none` and `smooth` limiters."""
function apply_gov_valve_limit_hard!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    Pv_raw::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    p_min::Dict{Int, Float64},
    p_max::Dict{Int, Float64};
    encoding::BoundEncoding=CONSTRAINT,
)
    if encoding == VARIABLE
        for gen in active_gen, t in eachindex(time_window)
            JuMP.set_lower_bound(Pv_raw[gen][t], p_min[gen])
            JuMP.set_upper_bound(Pv_raw[gen][t], p_max[gen])
        end
        return Pv_raw, nothing, :variable_bounds
    end

    lower = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    upper = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        Pmax = p_max[gen]
        Pmin = p_min[gen]
        lower[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        upper[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            lower[gen][t] = JuMP.@constraint(model, Pmin - Pv_raw[gen][t] ≤ 0.0)
            upper[gen][t] = JuMP.@constraint(model, Pv_raw[gen][t] - Pmax ≤ 0.0)
        end
    end
    return Pv_raw, nothing, (lower, upper)  # output = raw state, bounded in place
end

"""Dispatch the valve limiter on `mode`; returns `(Pv_output, extra_eq, extra_ineq)`.

`extra_eq`/`extra_ineq` are `nothing` unless the mode adds constraints (SMOOTH → the
clamp equality in `extra_eq`; HARD_BOUND → the `(lower, upper)` pair, or the
`:variable_bounds` marker under the VARIABLE encoding, in `extra_ineq`)."""
function apply_gov_valve_limit!(
    model::JuMP.Model,
    mode::GovernorLimiter,
    active_gen::Vector{Int64},
    Pv_raw::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    p_min::Dict{Int, Float64},
    p_max::Dict{Int, Float64};
    out_name::String="P_valve",
    encoding::BoundEncoding=CONSTRAINT,
)
    if mode == GOV_NO_LIMIT
        return apply_gov_valve_limit_none!(model, active_gen, Pv_raw, time_window, p_min, p_max)
    elseif mode == GOV_SMOOTH
        return apply_gov_valve_limit_smooth!(
            model, active_gen, Pv_raw, time_window, p_min, p_max; out_name=out_name)
    else # GOV_HARD_BOUND
        return apply_gov_valve_limit_hard!(
            model, active_gen, Pv_raw, time_window, p_min, p_max; encoding=encoding)
    end
end

# ===================================================================================
# Mech ODE (trapezoidal): T3·dPm/dt + Pm = T2·dPv/dt + Pv
# ===================================================================================
# `Pv` is the (possibly-limited) valve output — the only signal the turbine sees, which
# is what keeps the valve limiter effective (see the file header). Neither P_ref, Δω, R
# nor T1 enter here. t=1 anchors: fault → (Pm=P_m, Pv=P_m); post-fault → (last fault-on
# Pm, last fault-on *limited* Pv).

"""Mechanical-power update for one window: first step via `ode_first_step`, then trap.

Discretizes the turbine lead-lag `T3·dPm/dt + Pm = T2·dPv/dt + Pv` directly, both sides
integrated over the step, then divided through by `2·T3` (`T3` for backward Euler) so the
row keeps the `(1 + c)` normalization used by the valve ODE — this fixes the dual scale of
the `dual_gov_mech` family.  `Pv_prev0` must be the *limited* valve output of the previous
window, not the raw integrator state."""
function eq_const_gov_mech!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pm::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    Pv::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    Pm_prev0::OrderedDict{Int, JuMP.VariableRef},
    Pv_prev0::OrderedDict{Int, JuMP.VariableRef},
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        T2, T3 = DGEN_DYN.T2[gen], DGEN_DYN.T3[gen]
        c    = Δt / (2 * T3)      # trapezoidal half-step, module convention
        lead = T2 / T3            # turbine lead/lag ratio
        eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            Pm_prev = t == 1 ? Pm_prev0[gen] : Pm[gen][t - 1]
            Pv_prev = t == 1 ? Pv_prev0[gen] : Pv[gen][t - 1]
            if t == 1 && ode_first_step === :backward_euler
                # --- BACKWARD EULER FOR STEP 1 (reference GFM path) ---
                eq[gen][t] = JuMP.@constraint(model,
                    Pm[gen][t] * (1 + Δt / T3) - Pm_prev
                    - Pv[gen][t] * (lead + Δt / T3) + Pv_prev * lead == 0.0)
            else
                # --- TRAPEZOIDAL FOR ALL t ---
                eq[gen][t] = JuMP.@constraint(model,
                    Pm[gen][t] * (1 + c) - Pm_prev * (1 - c)
                    - Pv[gen][t] * (lead + c) + Pv_prev * (lead - c) == 0.0)
            end
        end
    end
    return eq
end

# ===================================================================================
# Orchestrators — called from the FULL_BUS builder before the Δω swing constraint
# ===================================================================================

"""Per-generator governor valve `(p_min, p_max)` from resolved `var_limit_specs[:gov_valve]`."""
function _gov_valve_limit_dicts(
    active_gen::Vector{Int64},
    specs::OrderedDict{Symbol, Tuple{Any, Any}},
)::Tuple{Dict{Int, Float64}, Dict{Int, Float64}}
    gmin, gmax = specs[:gov_valve]
    p_min = Dict{Int, Float64}()
    p_max = Dict{Int, Float64}()
    for (idx, gen) in enumerate(active_gen)
        p_min[gen] = gmin isa AbstractVector ? gmin[idx] : Float64(gmin)
        p_max[gen] = gmax isa AbstractVector ? gmax[idx] : Float64(gmax)
    end
    return p_min, p_max
end

"""
Attach the governor to the **fault-on** window and return the time-varying `Pm_tf`.

Creates the set-point `P_ref` (pinned to `R·P_m`) on first use, the valve state
`Pv_raw_tf` (+ limited output `Pv_tf`), and the mechanical power `Pm_tf`, wiring the
trapezoidal valve/mech ODEs with the pre-fault equilibrium (`P_m`, Δω=0) as the t=1
anchor. All families are stored in `dyn_model_dict` for dual export / results.
"""
function Attach_Governor_fault!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    base_MVA::Float64,
    Δω_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    Δt::Float64,
    limiter::GovernorLimiter,
)
    P_m = dyn_model_dict[:vars][:P_m]
    specs = dyn_model_dict[:meta][:var_limit_specs]
    p_min, p_max = _gov_valve_limit_dicts(active_gen, specs)

    # Anchor the governor equilibrium to the dispatch: without an explicit P_m = P_g pin the
    # governor states (P_ref, valve, mech) form a shift-invariant system and P_m drifts to a
    # spurious level decoupled from the OPF. Classical / DQ builders already stamp
    # `eq_const_Pm_init` whenever `mech_power_mode = USE_PM`; keep this fallback for safety
    # if a path attached the governor without that pin.
    if !haskey(dyn_model_dict[:eq_const], :eq_const_Pm_init)
        P_g = dyn_model_dict[:refs][:P_g]
        dyn_model_dict[:eq_const][:eq_const_Pm_init] =
            eq_const_kron_initial_mechanical_power!(model, P_m, P_g, active_gen)
    end

    # Constant set-point P_ref = R·P_m (created once, reused by the post-fault window).
    P_ref = var_gov_setpoint!(model, active_gen, DGEN_DYN, P_m)
    dyn_model_dict[:vars][:P_ref] = P_ref
    dyn_model_dict[:eq_const][:eq_const_Pref_init] =
        eq_const_gov_setpoint_init!(model, active_gen, DGEN_DYN, P_ref, P_m)
    attach_gov_prefault_bounds!(model, dyn_model_dict)

    # Raw valve state + limited output. The limiter is built *before* the valve ODE because
    # that ODE feeds the limited output back as its previous-step term (anti-windup), so the
    # `Pv_tf` container has to exist first — same ordering as `Attach_Avr_fault!`.
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)
    Pv_raw_tf = var_gov_state_time!(model, active_gen, "P_valve_raw_tf", time_window, P_m)
    dyn_model_dict[:vars][:Pv_raw_tf] = Pv_raw_tf

    Pv_tf, limit_eq, limit_ineq = apply_gov_valve_limit!(
        model, limiter, active_gen, Pv_raw_tf, time_window, p_min, p_max;
        out_name="P_valve_tf", encoding=bound_encoding_from_meta(dyn_model_dict[:meta]))
    dyn_model_dict[:vars][:Pv_tf] = Pv_tf
    _store_gov_limit!(dyn_model_dict, "tf", limit_eq, limit_ineq, Pv_raw_tf)

    # t=1 anchors on the pre-fault equilibrium, where Pv_raw = Pv = Pm = P_m (Δω = 0).
    dyn_model_dict[:eq_const][:eq_const_gov_valve_tf] = eq_const_gov_valve!(
        model, active_gen, DGEN_DYN, Pv_raw_tf, Pv_tf, P_ref, Δω_tf, P_m, 0.0, time_window, Δt;
        ode_first_step=ode_fs)

    # Mechanical power. The turbine consumes the *limited* valve output Pv_tf only; the
    # t=1 anchor is the pre-fault equilibrium, where Pm = Pv = P_m.
    Pm_tf = var_gov_state_time!(model, active_gen, "P_mech_tf", time_window, P_m)
    dyn_model_dict[:vars][:Pm_tf] = Pm_tf
    dyn_model_dict[:eq_const][:eq_const_gov_mech_tf] = eq_const_gov_mech!(
        model, active_gen, DGEN_DYN, Pm_tf, Pv_tf, P_m, P_m, time_window, Δt;
        ode_first_step=ode_fs)

    return Pm_tf
end

"""
Attach the governor to the **post-fault** window and return `Pm_tpf`.

Reuses `P_ref` from the fault-on call; the t=1 anchor is the last fault-on governor
state (the *limited* valve output and the mech state) and last fault-on Δω, mirroring
the swing/EMF continuity.
"""
function Attach_Governor_postf!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    base_MVA::Float64,
    Δω_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    Δt::Float64,
    limiter::GovernorLimiter,
)
    P_m   = dyn_model_dict[:vars][:P_m]
    P_ref = dyn_model_dict[:vars][:P_ref]
    specs = dyn_model_dict[:meta][:var_limit_specs]
    p_min, p_max = _gov_valve_limit_dicts(active_gen, specs)

    # Last fault-on states = post-fault t=1 anchors. Both the valve and the turbine take the
    # *limited* `Pv_out_last`: the valve because its previous-step term is the saturated output
    # everywhere else in the window, the turbine because it never sees the raw state at all.
    # Anchoring the valve on the raw state instead would leave a one-step discontinuity in the
    # recurrence exactly at the window seam under GOV_SMOOTH.
    Pv_out_last = OrderedDict(g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pv_tf])
    Pm_last     = OrderedDict(g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pm_tf])
    Δω_last     = OrderedDict(g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Δω_tf])

    Pv_raw_tpf = var_gov_state_time!(model, active_gen, "P_valve_raw_tpf", time_window, P_m)
    dyn_model_dict[:vars][:Pv_raw_tpf] = Pv_raw_tpf
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)

    Pv_tpf, limit_eq, limit_ineq = apply_gov_valve_limit!(
        model, limiter, active_gen, Pv_raw_tpf, time_window, p_min, p_max;
        out_name="P_valve_tpf", encoding=bound_encoding_from_meta(dyn_model_dict[:meta]))
    dyn_model_dict[:vars][:Pv_tpf] = Pv_tpf
    _store_gov_limit!(dyn_model_dict, "tpf", limit_eq, limit_ineq, Pv_raw_tpf)

    dyn_model_dict[:eq_const][:eq_const_gov_valve_tpf] = eq_const_gov_valve!(
        model, active_gen, DGEN_DYN, Pv_raw_tpf, Pv_tpf, P_ref, Δω_tpf,
        Pv_out_last, Δω_last, time_window, Δt; ode_first_step=ode_fs)

    # The turbine anchor is `Pv_out_last` — the *limited* valve output of the last fault-on
    # step, never the raw integrator state. Under GOV_SMOOTH the two differ exactly by the
    # clamp, and feeding the raw state here would leak the unsaturated signal across the seam.
    Pm_tpf = var_gov_state_time!(model, active_gen, "P_mech_tpf", time_window, P_m)
    dyn_model_dict[:vars][:Pm_tpf] = Pm_tpf
    dyn_model_dict[:eq_const][:eq_const_gov_mech_tpf] = eq_const_gov_mech!(
        model, active_gen, DGEN_DYN, Pm_tpf, Pv_tpf,
        Pm_last, Pv_out_last, time_window, Δt; ode_first_step=ode_fs)

    return Pm_tpf
end

"""Store the optional valve-limiter constraints (SMOOTH equality or HARD bounds) by window.

`limit_ineq` is `nothing` (no limiter), the `(lower, upper)` constraint pair produced by
the CONSTRAINT encoding, or the `:variable_bounds` marker — in which case the bounds live
on `Pv_raw` itself and the export goes through the bound manifest, exactly as for the δ,
Pe and dq box families. The registry rows `ineq_const_gov_valve_<win>_<side>` serve both
encodings, so the exported column names do not depend on the choice."""
function _store_gov_limit!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    win::String,
    limit_eq,
    limit_ineq,
    Pv_raw,
)
    if limit_eq !== nothing
        dyn_model_dict[:eq_const][Symbol("eq_const_gov_valve_limit_$win")] = limit_eq
    end
    key_lower = Symbol("ineq_const_gov_valve_$(win)_lower")
    key_upper = Symbol("ineq_const_gov_valve_$(win)_upper")
    if limit_ineq === :variable_bounds
        meta = dyn_model_dict[:meta]
        register_bound_manifest!(meta, key_lower, :lower, Pv_raw, :per_gen_time)
        register_bound_manifest!(meta, key_upper, :upper, Pv_raw, :per_gen_time)
    elseif limit_ineq !== nothing
        lower, upper = limit_ineq
        dyn_model_dict[:ineq_const][key_lower] = lower
        dyn_model_dict[:ineq_const][key_upper] = upper
    end
    return nothing
end
