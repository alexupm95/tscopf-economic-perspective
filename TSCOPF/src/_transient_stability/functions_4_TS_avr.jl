#=
================================================================================
 functions_4_TS_avr.jl — optional SEXS AVR (exciter) for DQ_4TH
================================================================================
 Voltage regulator for the 4th-order dq machine.  The exciter sees the generator
 terminal voltage `V` and drives the field voltage `E_fd` that feeds the Eq EMF
 ODE.  Wired into the DQ FULL_BUS path only (`include_avr=true`).

 AVR model — the SEXS chain (ANDES `sexs.py`: `LeadLag` → `LagAntiWindup`):

     u        = V_ref − V_terminal                              (summing junction)
     E_LL     = (1 + T_a·s)/(1 + T_b·s) · u                     (lead-lag, optional)
     E_fd_unlim = K_exc/(1 + T_exc·s) · E_LL                    (gain-lag)
     E_fd     = smooth clamp(E_fd_unlim; E_min, E_max)          (field limits)

 The gain-lag is a first-order lag, NOT an integrator; the −E_fd_unlim self-decay
 term is what makes the pre-fault link below the steady state of the same ODE.
 The lead-lag has unity DC gain, so it does not disturb that link:

     Pre-fault:  E_fd = K_exc · (V_ref − V_terminal)
     Transient:  T_b · dE_LL/dt + E_LL = T_a · du/dt + u
                 T_exc · d(E_fd_unlim)/dt = K_exc · E_LL − E_fd_unlim

 The lead-lag stage is skipped per generator when `Ta_exc` and `Tb_exc` are BOTH
 zero, in which case `E_LL ≡ u` and the chain collapses to the gain-lag alone —
 no variables, no rows.  `Tb_exc = 0` with `Ta_exc > 0` is a bare differentiator
 (improper) and is rejected by `validate_dyn_data!`.

 Discretization is trapezoidal (module convention).  The t=1 step of each window
 is anchored to the previous-window **saturated** E_fd (pre-fault scalar for the
 fault window; last fault-on E_fd_tf for post-fault), and to the previous-window
 lead-lag output / input for the lead-lag row.
================================================================================
=#

# `_AVR_SMOOTH_RHO` (sqrt min/max smoothing) lives in `_transient_stability/numerics.jl`.

"""Index `x[gen]` when `x` is a per-generator dict, else return the scalar `x`."""
_avr_at(x::AbstractDict, gen::Int) = x[gen]
_avr_at(x, ::Int) = x

_avr_start(v::JuMP.VariableRef)::Float64 = (s = JuMP.start_value(v); isnothing(s) ? 0.0 : Float64(s))

"""
Generators whose SEXS lead-lag stage is actually built — those with `Ta_exc` or
`Tb_exc` non-zero.

`Ta_exc = Tb_exc = 0` is an exact pass-through (`E_LL ≡ u`), so building the block
would only add one variable and one linear row per (gen, t) with `E_LL = u` as their
unique solution.  Skipping it keeps legacy AVR cases identical row-for-row to the
pre-lead-lag model — see the ANDES `zero_out=True` flag on the same block.
"""
function _avr_leadlag_gens(active_gen::Vector{Int64}, DGEN_DYN::DataFrame)::Vector{Int64}
    return Int64[gen for gen in active_gen
                 if !(iszero(DGEN_DYN.Ta_exc[gen]) && iszero(DGEN_DYN.Tb_exc[gen]))]
end

"""
Warn when the trapezoidal lead-lag row alternates in sign at the chosen step size.

The per-step factor on the previous output is `(2·T_b − Δt)/(2·T_b + Δt)`, whose magnitude
is below 1 for any `T_b > 0` — the recurrence never grows — but which turns negative once
`Δt > 2·T_b`, showing up as step-to-step ringing on `E_LL`.  Same trade-off, and the same
treatment (warn, do not fail), as the GFM measurement filters.

`T_b = 0` is exempt: there the factor is exactly `−1` and the row is the pass-through
`y_t + y_{t−1} = u_t + u_{t−1}`, whose alternating mode carries zero amplitude.
"""
function warn_avr_leadlag_step!(ll_gens::Vector{Int64}, DGEN_DYN::DataFrame, Δt::Float64)
    ringing = Int64[g for g in ll_gens
                    if DGEN_DYN.Tb_exc[g] > 0.0 && Δt > 2.0 * DGEN_DYN.Tb_exc[g]]
    isempty(ringing) && return nothing
    Tb_min = minimum(DGEN_DYN.Tb_exc[g] for g in ringing)
    @warn "AVR lead-lag is trapezoidal at t_step=$(Δt) s against a fastest lag " *
          "Tb_exc=$(Tb_min) s (generators $(ringing)). Above Δt = 2·Tb_exc the update " *
          "oscillates step-to-step; reduce t_step or raise Tb_exc."
    return nothing
end

# ===================================================================================
# Variables
# ===================================================================================

"""Constant voltage set-point `V_ref[g]` (one scalar per active SG).

Warm-started at `val_V[bus] + E_fd_start/K_exc` (reference AVR), using the
ACOPF warm-start bus voltage dictionary — not `JuMP.start_value(V)`, which is
still ~1.0 when AVR init runs (before OPF re-stamp). Falls back to
`start_value(V)` only if `val_V` lacks that bus.
"""
function var_avr_V_ref!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    E_fd::OrderedDict{Int, JuMP.VariableRef},
    val_V::Dict,
)
    V_ref = OrderedDict{Int, JuMP.VariableRef}()
    for gen in active_gen
        bus = DGEN.bus[gen]
        K_exc = DGEN_DYN.K_exc[gen]
        start_Efd = _avr_start(E_fd[gen])
        v_seed = haskey(val_V, bus) ? Float64(val_V[bus]) : _avr_start(V[bus])
        V_ref[gen] = JuMP.@variable(model, base_name = "V_ref[$gen]")
        JuMP.set_start_value(V_ref[gen], v_seed + start_Efd / K_exc)
    end
    return V_ref
end

"""
Time-indexed lead-lag output `E_LL`, warm-started at the pre-fault input `E_fd/K_exc`.

At the equilibrium the block passes its input through, so `E_LL = u = V_ref − V`, which
the steady-state gain link makes equal to `E_fd/K_exc`.
"""
function var_avr_leadlag_time!(
    model::JuMP.Model,
    ll_gens::Vector{Int64},
    var_name::String,
    time_window::Vector{Float64},
    DGEN_DYN::DataFrame,
    E_fd::OrderedDict{Int, JuMP.VariableRef},
)
    E_LL = var_kron_gen_time_generic!(model, ll_gens, var_name, time_window)
    for gen in ll_gens
        seed = _avr_start(E_fd[gen]) / DGEN_DYN.K_exc[gen]
        for t in eachindex(time_window)
            JuMP.set_start_value(E_LL[gen][t], seed)
        end
    end
    return E_LL
end

"""Time-indexed pre-saturation field voltage `E_fd_unlim`, warm-started at pre-fault `E_fd`."""
function var_avr_E_fd_unlim_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    var_name::String,
    time_window::Vector{Float64},
    E_fd::OrderedDict{Int, JuMP.VariableRef},
)
    E_fd_unlim = var_kron_gen_time_generic!(model, active_gen, var_name, time_window)
    for gen in active_gen
        seed = _avr_start(E_fd[gen])
        for t in eachindex(time_window)
            JuMP.set_start_value(E_fd_unlim[gen][t], seed)
        end
    end
    return E_fd_unlim
end

# ===================================================================================
# Pre-fault link: E_fd = K_exc · (V_ref − V_terminal)
# ===================================================================================

"""Pin pre-fault field voltage to the exciter steady-state gain."""
function eq_const_avr_Efd_init!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    V_ref::OrderedDict{Int, JuMP.VariableRef},
    E_fd::OrderedDict{Int, JuMP.VariableRef},
)
    eq = OrderedDict{Int, JuMP.ConstraintRef}()
    for gen in active_gen
        bus = DGEN.bus[gen]
        K_exc = DGEN_DYN.K_exc[gen]
        eq[gen] = JuMP.@constraint(model, E_fd[gen] - K_exc * (V_ref[gen] - V[bus]) == 0.0)
    end
    return eq
end

# ===================================================================================
# Lead-lag stage: T_b·dE_LL/dt + E_LL = T_a·du/dt + u,  u = V_ref − V_terminal
# ===================================================================================

"""
Lead-lag update for one window (fault or post-fault).

Discretizes `T_b·dy/dt + y = T_a·du/dt + u` directly, derivative of the input included,
then normalizes so the coefficient on `y_t` is exactly 1:

    trapezoidal   y_t − y_{t−1}·(2T_b−Δt)/n − u_t·(2T_a+Δt)/n + u_{t−1}·(2T_a−Δt)/n = 0,
                  n = 2·T_b + Δt
    backward Euler y_t − y_{t−1}·T_b/n − u_t·(T_a+Δt)/n + u_{t−1}·T_a/n = 0,
                  n = T_b + Δt

Normalizing by `n` rather than by `2·T_b` — the `(1+c)` form the lag ODEs use — is what
keeps the row well scaled as `T_b` gets small: `n` is strictly positive for any `Δt > 0`,
where `2·T_b` would head for zero and take the coefficients with it.

That is robustness, not a supported configuration: every generator reaching this builder
has `T_b > 0`. `_avr_leadlag_gens` filters out `T_a = T_b = 0` (the exact pass-through,
bypassed outright), and `validate_avr_data!` rejects `T_b = 0` with `T_a > 0` (a bare
differentiator, whose `(2T_a ± Δt)/Δt` coefficients diverge as the step shrinks). Do not
read the scaling argument above as licence to drop that validation — it covers the
arithmetic, not the modelling.

Normalizing to a unit `y_t` coefficient also keeps `dual_avr_leadlag` comparable across
machines with different `T_b`; the raw `(2T_b ± Δt)` form scales each row by its own `T_b`.

`u_prev0` / `E_LL_prev0` supply the t=1 anchor: fault window → the pre-fault equilibrium
(where both equal `V_ref − V`); post-fault → the last fault-on input and output.
"""
function eq_const_avr_leadlag!(
    model::JuMP.Model,
    ll_gens::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    E_LL::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    V_ref::OrderedDict{Int, JuMP.VariableRef},
    V_terminal::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    E_LL_prev0::AbstractDict,
    u_prev0::AbstractDict,
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
)
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in ll_gens
        bus = DGEN.bus[gen]
        Ta = DGEN_DYN.Ta_exc[gen]
        Tb = DGEN_DYN.Tb_exc[gen]
        n_tr = 2 * Tb + Δt
        n_be = Tb + Δt
        y_tr, u_tr, up_tr = (2 * Tb - Δt) / n_tr, (2 * Ta + Δt) / n_tr, (2 * Ta - Δt) / n_tr
        y_be, u_be, up_be = Tb / n_be, (Ta + Δt) / n_be, Ta / n_be
        eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            y_prev = t == 1 ? E_LL_prev0[gen] : E_LL[gen][t - 1]
            u_curr = JuMP.@expression(model, V_ref[gen] - V_terminal[bus][t])
            u_prev = t == 1 ? u_prev0[gen] :
                JuMP.@expression(model, V_ref[gen] - V_terminal[bus][t - 1])
            if t == 1 && ode_first_step === :backward_euler
                # --- BACKWARD EULER FOR STEP 1 (reference GFM path) ---
                eq[gen][t] = JuMP.@constraint(model,
                    E_LL[gen][t] - y_be * y_prev - u_be * u_curr + up_be * u_prev == 0.0)
            else
                # --- TRAPEZOIDAL FOR ALL OTHER STEPS ---
                eq[gen][t] = JuMP.@constraint(model,
                    E_LL[gen][t] - y_tr * y_prev - u_tr * u_curr + up_tr * u_prev == 0.0)
            end
        end
    end
    return eq
end

# ===================================================================================
# Exciter ODE (first step: trap or BE via ode_first_step; trap thereafter)
# ===================================================================================

"""
Exciter update for one window (fault or post-fault).

First step uses `ode_first_step` (`:trapezoidal` default, or `:backward_euler`);
remaining steps use trapezoidal. The ODE integrates `E_fd_unlim`.  The
previous-step anchor uses the **saturated** output `E_fd_sat` (the reference
convention), not the raw unlimited state.

`V_terminal` is the per-step generator-bus voltage from `V_tf` / `V_tpf`.
`E_fd_prev0` / `V_prev0` supply the t=1 anchor: fault window → pre-fault
`(E_fd, V)`; post-fault → last fault-on `(E_fd_tf, V_tf)`.

When a generator has a lead-lag stage, `E_LL` carries its output and the lag is driven
by that instead of by the raw `V_ref − V`: the trapezoidal input average becomes
`E_LL[t] + E_LL[t−1]` in place of `(V_ref − V_t) + (V_ref − V_{t−1})`, which is what the
`(2·V_ref − V_curr − V_prev)` term spells out.  `E_LL_prev0` is its t=1 anchor.
Generators absent from `E_LL` (the `Ta_exc = Tb_exc = 0` bypass) keep the original rows.
"""
function eq_const_avr_exciter!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    E_fd_unlim::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    E_fd_sat::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    V_ref::OrderedDict{Int, JuMP.VariableRef},
    V_terminal::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    E_fd_prev0::OrderedDict{Int, JuMP.VariableRef},
    V_prev0,
    time_window::Vector{Float64},
    Δt::Float64;
    ode_first_step::Symbol=:trapezoidal,
    E_LL=nothing,
    E_LL_prev0=nothing,
)
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        bus = DGEN.bus[gen]
        K_exc = DGEN_DYN.K_exc[gen]
        T_exc = DGEN_DYN.T_exc[gen]
        c = Δt / (2 * T_exc)
        use_ll = E_LL !== nothing && haskey(E_LL, gen)
        eq[gen] = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
        for t in eachindex(time_window)
            E_fd_prev = t == 1 ? E_fd_prev0[gen] : E_fd_sat[gen][t - 1]
            V_curr = V_terminal[bus][t]
            if use_ll
                # --- LEAD-LAG OUTPUT DRIVES THE LAG (V_ref / V enter via E_LL only) ---
                u_curr = E_LL[gen][t]
                u_prev = t == 1 ? E_LL_prev0[gen] : E_LL[gen][t - 1]
                if t == 1 && ode_first_step === :backward_euler
                    eq[gen][t] = JuMP.@constraint(model,
                        E_fd_unlim[gen][t] * (1 + Δt / T_exc) - E_fd_prev
                        - ((K_exc*Δt) / T_exc) * u_curr == 0.0)
                else
                    eq[gen][t] = JuMP.@constraint(model,
                        E_fd_unlim[gen][t] * (1 + c) - E_fd_prev * (1 - c)
                        - (K_exc * c) * (u_curr + u_prev) == 0.0)
                end
            elseif t == 1 && ode_first_step === :backward_euler
                # --- BACKWARD EULER FOR STEP 1 (reference GFM path) ---
                eq[gen][t] = JuMP.@constraint(model,
                    E_fd_unlim[gen][t] * (1 + Δt / T_exc) - E_fd_prev
                    - ((K_exc*Δt) / T_exc) * (V_ref[gen] - V_curr) == 0.0)
            elseif t == 1
                # --- TRAPEZOIDAL AT t=1 (SG pin default) ---
                V_prev = _avr_at(V_prev0, bus)
                eq[gen][t] = JuMP.@constraint(model,
                    E_fd_unlim[gen][t] * (1 + c) - E_fd_prev * (1 - c)
                    - (K_exc * c) * (2 * V_ref[gen] - V_curr - V_prev) == 0.0)
            else
                # --- TRAPEZOIDAL FOR REMAINDER ---
                V_prev = V_terminal[bus][t - 1]
                eq[gen][t] = JuMP.@constraint(model,
                    E_fd_unlim[gen][t] * (1 + c) - E_fd_prev * (1 - c)
                    - (K_exc * c) * (2 * V_ref[gen] - V_curr - V_prev) == 0.0)
            end
        end
    end
    return eq
end

# ===================================================================================
# Field-voltage saturation (smooth sqrt min/max)
# ===================================================================================

"""Per-generator field-voltage `(E_min, E_max)` from resolved `var_limit_specs[:E]`."""
function _avr_field_limit_dicts(
    active_gen::Vector{Int64},
    specs::OrderedDict{Symbol, Tuple{Any, Any}},
)::Tuple{Dict{Int, Float64}, Dict{Int, Float64}}
    emin, emax = specs[:E]
    e_min = Dict{Int, Float64}()
    e_max = Dict{Int, Float64}()
    for (idx, gen) in enumerate(active_gen)
        e_min[gen] = emin isa AbstractVector ? emin[idx] : Float64(emin)
        e_max[gen] = emax isa AbstractVector ? emax[idx] : Float64(emax)
    end
    return e_min, e_max
end

"""
Smooth anti-windup: `E_fd = smoothmax(smoothmin(E_fd_unlim, E_max), E_min)` via sqrt.

Introduces a separate limited-output variable `E_fd` tied to the raw state
`E_fd_unlim` by one equality per (gen, t).  Physically saturating but nonconvex.
"""
function apply_avr_field_limit_smooth!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    E_fd_unlim::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    e_min::Dict{Int, Float64},
    e_max::Dict{Int, Float64};
    out_name::String="E_fd",
)
    ρ = _AVR_SMOOTH_RHO
    E_fd = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    eq = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for gen in active_gen
        Emax = e_max[gen]
        Emin = e_min[gen]
        E_fd[gen] = OrderedDict{Int, JuMP.VariableRef}()
        eq[gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        for t in eachindex(time_window)
            # Window-specific name (E_fd_tf / E_fd_tpf) so dumps do not collide.
            E_fd[gen][t] = JuMP.@variable(model, base_name = "$out_name[$gen, $t]",
                start = _avr_start(E_fd_unlim[gen][t]))
            u = E_fd_unlim[gen][t]
            E_aux = JuMP.@expression(model, (u + Emax - sqrt((u - Emax)^2 + ρ)) / 2)
            eq[gen][t] = JuMP.@constraint(model,
                E_fd[gen][t] == (E_aux + Emin + sqrt((E_aux - Emin)^2 + ρ)) / 2)
        end
    end
    return E_fd, eq
end

# ===================================================================================
# Orchestrators — called from the DQ FULL_BUS builder
# ===================================================================================

"""
Attach AVR pre-fault variables and the steady-state exciter link.

Call from `Define_Initial_Condition_4_dq!` when `include_avr=true`.
`val_V` seeds `V_ref` the reference way: `V_ref ← val_V[bus] + E_fd/K_exc`.
"""
function Attach_Avr_init!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    val_V::Dict,
)
    E_fd = dyn_model_dict[:vars][:E_fd]
    V_ref = var_avr_V_ref!(model, active_gen, DGEN, DGEN_DYN, V, E_fd, val_V)
    dyn_model_dict[:vars][:V_ref] = V_ref
    dyn_model_dict[:eq_const][:eq_const_Efd_init] =
        eq_const_avr_Efd_init!(model, active_gen, DGEN, DGEN_DYN, V, V_ref, E_fd)
    return nothing
end

"""
Attach the exciter to the **fault-on** window; returns saturated time-varying `E_fd_tf`.

The t=1 anchor uses pre-fault `(E_fd, V_terminal)`.  When any generator carries a
lead-lag stage its output `E_LL_tf` is built first and anchored at the pre-fault
equilibrium, where the block passes through and `E_LL = V_ref − V`.
"""
function Attach_Avr_fault!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    V_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    Δt::Float64,
)
    E_fd = dyn_model_dict[:vars][:E_fd]
    V_ref = dyn_model_dict[:vars][:V_ref]
    specs = dyn_model_dict[:meta][:var_limit_specs]
    e_min, e_max = _avr_field_limit_dicts(active_gen, specs)
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)

    # Lead-lag stage (skipped entirely when every generator has Ta_exc = Tb_exc = 0, so
    # the registry's presence filter drops `dual_avr_leadlag` on those runs).
    ll_gens = _avr_leadlag_gens(active_gen, DGEN_DYN)
    E_LL_tf = nothing
    E_LL_prev0 = nothing
    if !isempty(ll_gens)
        warn_avr_leadlag_step!(ll_gens, DGEN_DYN, Δt)
        E_LL_tf = var_avr_leadlag_time!(model, ll_gens, "E_LL_tf", time_window, DGEN_DYN, E_fd)
        dyn_model_dict[:vars][:E_LL_tf] = E_LL_tf
        # Pre-fault equilibrium: the block passes through, so y_0 = u_0 = V_ref − V.
        E_LL_prev0 = OrderedDict(
            gen => JuMP.@expression(model, V_ref[gen] - V[DGEN.bus[gen]]) for gen in ll_gens)
        dyn_model_dict[:eq_const][:eq_const_avr_leadlag_tf] = eq_const_avr_leadlag!(
            model, ll_gens, DGEN, DGEN_DYN, E_LL_tf, V_ref, V_tf, E_LL_prev0, E_LL_prev0,
            time_window, Δt; ode_first_step=ode_fs)
    end

    E_fd_unlim_tf = var_avr_E_fd_unlim_time!(model, active_gen, "E_fd_unlim_tf", time_window, E_fd)
    dyn_model_dict[:vars][:E_fd_unlim_tf] = E_fd_unlim_tf

    E_fd_tf, sat_eq = apply_avr_field_limit_smooth!(
        model, active_gen, E_fd_unlim_tf, time_window, e_min, e_max; out_name="E_fd_tf")
    dyn_model_dict[:vars][:E_fd_tf] = E_fd_tf
    dyn_model_dict[:eq_const][:eq_const_E_fd_tf] = sat_eq

    V_prev0 = OrderedDict(bus => V[bus] for bus in keys(V_tf))
    dyn_model_dict[:eq_const][:eq_const_E_fd_unlim_tf] = eq_const_avr_exciter!(
        model, active_gen, DGEN, DGEN_DYN, E_fd_unlim_tf, E_fd_tf, V_ref, V_tf, E_fd, V_prev0,
        time_window, Δt; ode_first_step=ode_fs, E_LL=E_LL_tf, E_LL_prev0=E_LL_prev0)
    return E_fd_tf
end

"""
Attach the exciter to the **post-fault** window; returns saturated `E_fd_tpf`.

Reuses `V_ref` from init; t=1 anchors to the last fault-on saturated exciter state
and terminal voltage.  The lead-lag row anchors on the last fault-on `(E_LL_tf, V_tf)`
pair, mirroring the exciter/EMF continuity across the seam.
"""
function Attach_Avr_postfault!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    V_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    V_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    time_window::Vector{Float64},
    Δt::Float64,
)
    E_fd = dyn_model_dict[:vars][:E_fd]
    V_ref = dyn_model_dict[:vars][:V_ref]
    specs = dyn_model_dict[:meta][:var_limit_specs]
    e_min, e_max = _avr_field_limit_dicts(active_gen, specs)
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)

    E_fd_last = OrderedDict(g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:E_fd_tf])
    V_last = OrderedDict(bus => last(inner).second for (bus, inner) in V_tf)

    ll_gens = _avr_leadlag_gens(active_gen, DGEN_DYN)
    E_LL_tpf = nothing
    E_LL_last = nothing
    if !isempty(ll_gens)
        E_LL_tf = dyn_model_dict[:vars][:E_LL_tf]
        E_LL_tpf = var_avr_leadlag_time!(model, ll_gens, "E_LL_tpf", time_window, DGEN_DYN, E_fd)
        dyn_model_dict[:vars][:E_LL_tpf] = E_LL_tpf
        E_LL_last = OrderedDict(g => last(E_LL_tf[g]).second for g in ll_gens)
        u_last = OrderedDict(
            gen => JuMP.@expression(model, V_ref[gen] - V_last[DGEN.bus[gen]]) for gen in ll_gens)
        dyn_model_dict[:eq_const][:eq_const_avr_leadlag_tpf] = eq_const_avr_leadlag!(
            model, ll_gens, DGEN, DGEN_DYN, E_LL_tpf, V_ref, V_tpf, E_LL_last, u_last,
            time_window, Δt; ode_first_step=ode_fs)
    end

    E_fd_unlim_tpf = var_avr_E_fd_unlim_time!(model, active_gen, "E_fd_unlim_tpf", time_window, E_fd)
    dyn_model_dict[:vars][:E_fd_unlim_tpf] = E_fd_unlim_tpf

    E_fd_tpf, sat_eq = apply_avr_field_limit_smooth!(
        model, active_gen, E_fd_unlim_tpf, time_window, e_min, e_max; out_name="E_fd_tpf")
    dyn_model_dict[:vars][:E_fd_tpf] = E_fd_tpf
    dyn_model_dict[:eq_const][:eq_const_E_fd_tpf] = sat_eq

    dyn_model_dict[:eq_const][:eq_const_E_fd_unlim_tpf] = eq_const_avr_exciter!(
        model, active_gen, DGEN, DGEN_DYN, E_fd_unlim_tpf, E_fd_tpf, V_ref, V_tpf, E_fd_last, V_last,
        time_window, Δt; ode_first_step=ode_fs, E_LL=E_LL_tpf, E_LL_prev0=E_LL_last)
    return E_fd_tpf
end
