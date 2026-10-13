#=
================================================================================
 functions_4_TS_gfm.jl  —  grid-forming inverter builders
================================================================================
 GFM dynamics live in a separate `gfm_dynamic_data.csv` → in-memory `DGFM`.
 SG dynamics stay in `gen_dynamic_data*.csv` → `DGEN_DYN`.

 Phase G0: masks, meta registration, Attach stubs.
 Phase G1: pre-fault init + steady-state ACOPF GFM limits.
 Phase G2+: transient filters · droop · PI · current limiter.
================================================================================
=#

# ===================================================================================
# Module constants
# ===================================================================================
# Structural constants of the GFM builders. None of them is configuration: the
# user-facing GFM limits live in `TsBoundLimitsConfig` (fields `gfm_*`) and the device
# data in `gfm_dynamic_data.csv`.

# --- Equality-family registries ----------------------------------------------------
# Named equality families, stored per (gen, t) so `DynDualRegistry` can address them.
# The flat `eq_const_gfm_<window>` vector is kept alongside for the text appendix.
# The overload families are *not* in `_GFM_EQ_FAMILIES`: the heavy test iterates that
# const and expects every family on the default (knob-off) solve. They are stamped
# only when `gfm_pq_overload` is on.
const _GFM_EQ_FAMILIES = (
    :filter_P, :filter_Q, :filter_V, :droop, :delta,
    :Eint_raw, :Eint_clip, :Edroop_raw, :Edroop_clip,
    :limiter_Id, :limiter_Iq, :Pe, :Qe,
)

const _GFM_PQ_EQ_FAMILIES = (
    :Pmax_raw, :Pmax_clip, :Pmin_raw, :Pmin_clip,
    :Qmax_raw, :Qmax_clip, :Qmin_raw, :Qmin_clip,
)

# --- Numerical constants -----------------------------------------------------------
# The smoothing widths (`_GFM_EPS_*`, `_GFM_PQLIM_*`) and the current-limiter bypass
# threshold (`_GFM_IMAX_NO_LIMIT`) live in `_transient_stability/numerics.jl`, with every
# other numerical constant of the dynamic models. Restart Julia after editing one.

"""Set of generator ids present in `DGFM` (empty if `DGFM` is `nothing`)."""
function gfm_id_set(DGFM::Union{Nothing, DataFrame})::Set{Int}
    DGFM === nothing && return Set{Int}()
    return Set{Int}(Int.(DGFM.id))
end

"""Active generators that are synchronous (not listed in `DGFM`)."""
function sg_active_gens(
    active_gen::AbstractVector{Int},
    DGFM::Union{Nothing, DataFrame},
)::Vector{Int}
    gfm = gfm_id_set(DGFM)
    return Int[g for g in active_gen if g ∉ gfm]
end

"""Active generators that are grid-forming (listed in `DGFM`)."""
function gfm_active_gens(
    active_gen::AbstractVector{Int},
    DGFM::Union{Nothing, DataFrame},
)::Vector{Int}
    gfm = gfm_id_set(DGFM)
    return Int[g for g in active_gen if g ∈ gfm]
end

"""Row index in `DGFM` for generator `gen` (throws if missing)."""
function dgfm_row(DGFM::DataFrame, gen::Int)::Int
    i = findfirst(==(gen), Int.(DGFM.id))
    i === nothing && throw(ArgumentError("Generator id $gen not found in DGFM."))
    return i
end

"""
Stamp GFM membership into `dyn_model_dict[:meta]`.

Also records the per-unit converter parameters under `meta[:gfm_params]`. The save layer
has no access to `DGFM` — it reconstructs the limiter algebra and the filter coefficients
from solved primal values, and needs `Xl`, `Imax`, the three filter constants
`TPf`/`TQf`/`TVf`, `Emin`/`Emax` and the gains to do it. Plain `Float64`s, so they survive
`release_solver_backend!` unlike anything JuMP.
"""
function register_gfm_meta!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    DGFM::Union{Nothing, DataFrame},
)
    haskey(dyn_model_dict, :meta) || (dyn_model_dict[:meta] = OrderedDict{Symbol, Any}())
    ids = sort!(collect(gfm_id_set(DGFM)))
    dyn_model_dict[:meta][:allow_gfm] = !isempty(ids)
    dyn_model_dict[:meta][:gfm_ids] = ids
    dyn_model_dict[:meta][:n_gfm] = length(ids)
    if DGFM !== nothing && !isempty(ids)
        params = OrderedDict{Int, NamedTuple}()
        for gen in ids
            r = dgfm_row(DGFM, gen)
            params[gen] = (
                bus = Int(DGFM.bus[r]),
                Xl = Float64(DGFM.Xl[r]), Imax = Float64(DGFM.Imax[r]),
                TPf = Float64(DGFM.TPf[r]),
                TQf = Float64(DGFM.TQf[r]),
                TVf = Float64(DGFM.TVf[r]),
                Emin = Float64(DGFM.Emin[r]), Emax = Float64(DGFM.Emax[r]),
                mp = Float64(DGFM.mp[r]), mq = Float64(DGFM.mq[r]),
                Kpv = Float64(DGFM.Kpv[r]), Kiv = Float64(DGFM.Kiv[r]),
                Pmin = Float64(DGFM.Pmin[r]), Pmax = Float64(DGFM.Pmax[r]),
                Qmin = Float64(DGFM.Qmin[r]), Qmax = Float64(DGFM.Qmax[r]),
                Kppmax = Float64(DGFM.Kppmax[r]), Kipmax = Float64(DGFM.Kipmax[r]),
                Kpqmax = Float64(DGFM.Kpqmax[r]), Kiqmax = Float64(DGFM.Kiqmax[r]),
            )
        end
        dyn_model_dict[:meta][:gfm_params] = params
    end
    return dyn_model_dict
end

"""
    gfm_machine_warmstart(v, θ, P_g, Q_g, Xl) -> (δ, E, Id, Iq)

Map ACOPF terminal (V, θ, P, Q) to GFM internal angle/voltage behind `Xl`
(purely inductive coupling, PSS/E-style).
"""
function gfm_machine_warmstart(
    v::Real,
    θ::Real,
    P_g::Real,
    Q_g::Real,
    Xl::Real,
)::NTuple{4, Float64}
    V_c = Float64(v) * exp(1im * Float64(θ))
    I_c = (Float64(P_g) - 1im * Float64(Q_g)) / conj(V_c)
    V_int = V_c + I_c * (1im * Float64(Xl))
    δ = angle(V_int)
    E = abs(V_int)
    I_mag = abs(I_c)
    I_ang = angle(I_c)
    Id = I_mag * sin(δ - I_ang)
    Iq = I_mag * cos(δ - I_ang)
    return (Float64(δ), Float64(E), Float64(Id), Float64(Iq))
end

# ===================================================================================
# Steady-state ACOPF GFM inequalities (Phase G1)
# ===================================================================================

"""
Attach GFM current and internal-voltage limits to a solved/building ACOPF model.

    `P² + Q² ≤ (V · 1.0 · mach_base_MVA / base_MVA)²`
    `(V + Q·Xl/V)² + (P·Xl/V)² ≤ Emax²`

with `Xl` already on system base after GFM conversion.

**Current-rating convention (decision D4).** The steady-state limit is the converter
*nameplate* current, 1.0 pu on machine base, deliberately **not** the dynamic CSV `Imax`.
A dispatch is a sustained operating point and cannot bank on overload. The transient
builders use `DGFM.Imax` instead (typically 1.2 pu): a converter tolerates brief
overcurrent up to the thermal limit of its semiconductor junctions, and that headroom is
available during a fault but not in the pre-fault schedule. The two windows therefore use
different current ratings by design — see `docs/src/model/09_grid_forming.md`.
"""
function attach_gfm_acopf_limits!(
    model::Model,
    opf_dict::OrderedDict{Symbol, Any},
    DGEN::DataFrame,
    DGFM::DataFrame,
    base_MVA::Float64,
)
    isempty(DGFM) && return nothing
    V = opf_dict[:vars][:V]
    P_g = opf_dict[:vars][:P_g]
    Q_g = opf_dict[:vars][:Q_g]
    haskey(opf_dict, :ineq_const) ||
        (opf_dict[:ineq_const] = OrderedDict{Symbol, Any}())

    ineq_I = OrderedDict{Int, JuMP.ConstraintRef}()
    ineq_E = OrderedDict{Int, JuMP.ConstraintRef}()
    for r in 1:nrow(DGFM)
        gen = Int(DGFM.id[r])
        DGEN.g_status[gen] == 1 || continue
        haskey(P_g, gen) || continue
        bus = Int(DGEN.bus[gen])
        # Reference: hardcoded 1.0 on machine base → system-base current rating
        I_acopf = Float64(DGFM.mach_base_MVA[r]) / base_MVA
        Xl = Float64(DGFM.Xl[r])
        Emax = Float64(DGFM.Emax[r])
        ineq_I[gen] = @constraint(model,
            P_g[gen]^2 + Q_g[gen]^2 - (V[bus] * I_acopf)^2 ≤ 0.0)
        ineq_E[gen] = @constraint(model,
            (V[bus] + Q_g[gen] * (Xl / V[bus]))^2 +
            (P_g[gen] * (Xl / V[bus]))^2 - Emax^2 ≤ 0.0)
    end
    opf_dict[:ineq_const][:ineq_const_gfm_Imax] = ineq_I
    opf_dict[:ineq_const][:ineq_const_gfm_Emax] = ineq_E
    return ineq_I, ineq_E
end

# ===================================================================================
# Pre-fault GFM init (Phase G1)
# ===================================================================================

"""
    Attach_GFM_init!(model, dyn_model_dict, gfm_gens, DGEN, DGFM, V, θ, P_g, Q_g, hints)

Create pre-fault GFM variables and steady-state equalities (reference GRFORM init).
Merges `δ`, `Ed`, `Eq`, `Id`, `Iq`, and `P_m` (as `P_set`) into `dyn_model_dict[:vars]`
alongside any SG entries already created.
"""
function Attach_GFM_init!(
    model::Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    gfm_gens::AbstractVector{Int},
    DGEN::DataFrame,
    DGFM::DataFrame,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef},
    val_V::Dict,
    val_θ::Dict,
    val_Pg::Dict,
    val_Qg::Dict;
)
    register_gfm_meta!(dyn_model_dict, DGFM)
    isempty(gfm_gens) && return dyn_model_dict

    haskey(dyn_model_dict, :vars) || (dyn_model_dict[:vars] = OrderedDict{Symbol, Any}())
    haskey(dyn_model_dict, :eq_const) || (dyn_model_dict[:eq_const] = OrderedDict{Symbol, Any}())
    vars = dyn_model_dict[:vars]
    eqs = dyn_model_dict[:eq_const]
    encoding = bound_encoding_from_meta(dyn_model_dict[:meta])
    specs = _gfm_limit_specs(dyn_model_dict, DGFM, gfm_gens)
    δ_box_on = _gfm_box_enabled(dyn_model_dict[:meta], :gfm_δ)
    δ_gfm_only = OrderedDict{Int, JuMP.VariableRef}()
    # `P_m` (exported as `P_set`) is pinned to the dispatch by `P_m == P_g`. The dispatch
    # box itself belongs to `DGEN` (`pg_min`/`pg_max`). `DGFM.Pmin`/`Pmax` are device
    # ratings — comparator inputs to the overload PIs — and must not box the same variable.

    δ = get!(vars, :δ, OrderedDict{Int, JuMP.VariableRef}())
    Ed = get!(vars, :Ed, OrderedDict{Int, JuMP.VariableRef}())
    Eq = get!(vars, :Eq, OrderedDict{Int, JuMP.VariableRef}())
    Id = get!(vars, :Id, OrderedDict{Int, JuMP.VariableRef}())
    Iq = get!(vars, :Iq, OrderedDict{Int, JuMP.VariableRef}())
    P_m = get!(vars, :P_m, OrderedDict{Int, JuMP.VariableRef}())

    P_meas = OrderedDict{Int, JuMP.VariableRef}()
    Q_meas = OrderedDict{Int, JuMP.VariableRef}()
    V_meas = OrderedDict{Int, JuMP.VariableRef}()
    E_int = OrderedDict{Int, JuMP.VariableRef}()
    V_set = OrderedDict{Int, JuMP.VariableRef}()

    eq_Ed = get!(eqs, :eq_const_Ed_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Eq = get!(eqs, :eq_const_Eq_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Vd = get!(eqs, :eq_const_Vd_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Vq = get!(eqs, :eq_const_Vq_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_P = get!(eqs, :eq_const_P_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Q = get!(eqs, :eq_const_Q_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Pm = get!(eqs, :eq_const_Pm_init, OrderedDict{Int, JuMP.ConstraintRef}())
    eq_Pmeas = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Qmeas = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Vmeas = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_Vset = OrderedDict{Int, JuMP.ConstraintRef}()

    for (i, gen) in enumerate(gfm_gens)
        r = dgfm_row(DGFM, gen)
        bus = Int(DGEN.bus[gen])
        Xl = Float64(DGFM.Xl[r])
        mq = Float64(DGFM.mq[r])
        δ_lo, δ_hi = _gfm_lim(specs, :gfm_δ, i)

        v_val = Float64(val_V[bus])
        th_val = Float64(val_θ[bus])
        pg_val = Float64(val_Pg[gen])
        qg_val = Float64(val_Qg[gen])
        start_δ, start_E, start_Id, start_Iq = gfm_machine_warmstart(
            v_val, th_val, pg_val, qg_val, Xl)

        δ[gen] = create_bounded_variable!(model; base_name = "δ_gfm[$gen]",
            start = start_δ,
            min_lim = δ_box_on ? δ_lo : nothing, max_lim = δ_box_on ? δ_hi : nothing,
            encoding = encoding)
        δ_gfm_only[gen] = δ[gen]
        P_m[gen] = JuMP.@variable(model, base_name = "P_set[$gen]", start = pg_val)
        V_set[gen] = JuMP.@variable(model, base_name = "V_set[$gen]",
            start = v_val + mq * qg_val)
        P_meas[gen] = JuMP.@variable(model, base_name = "P_meas[$gen]", start = pg_val)
        Q_meas[gen] = JuMP.@variable(model, base_name = "Q_meas[$gen]", start = qg_val)
        V_meas[gen] = JuMP.@variable(model, base_name = "V_meas[$gen]", start = v_val)
        E_int[gen] = JuMP.@variable(model, base_name = "E_int[$gen]", start = start_E)
        Ed[gen] = JuMP.@variable(model, base_name = "Ed_gfm[$gen]", start = 0.0)
        Eq[gen] = JuMP.@variable(model, base_name = "Eq_gfm[$gen]", start = start_E)
        Id[gen] = JuMP.@variable(model, base_name = "Id_gfm[$gen]", start = start_Id)
        Iq[gen] = JuMP.@variable(model, base_name = "Iq_gfm[$gen]", start = start_Iq)

        # Terminal injection (same map as SG / reference GFM init)
        eq_P[gen] = JuMP.@constraint(model,
            P_g[gen] == V[bus] * sin(δ[gen] - θ[bus]) * Id[gen] +
                V[bus] * cos(δ[gen] - θ[bus]) * Iq[gen])
        eq_Q[gen] = JuMP.@constraint(model,
            Q_g[gen] == V[bus] * cos(δ[gen] - θ[bus]) * Id[gen] -
                V[bus] * sin(δ[gen] - θ[bus]) * Iq[gen])
        eq_Pm[gen] = JuMP.@constraint(model, P_m[gen] == P_g[gen])

        eq_Ed[gen] = JuMP.@constraint(model, Ed[gen] == 0.0)
        eq_Eq[gen] = JuMP.@constraint(model, Eq[gen] == E_int[gen])
        eq_Vd[gen] = JuMP.@constraint(model,
            V[bus] * sin(δ[gen] - θ[bus]) - Ed[gen] - Xl * Iq[gen] == 0)
        eq_Vq[gen] = JuMP.@constraint(model,
            V[bus] * cos(δ[gen] - θ[bus]) - Eq[gen] + Xl * Id[gen] == 0)

        eq_Pmeas[gen] = JuMP.@constraint(model, P_meas[gen] == P_g[gen])
        eq_Qmeas[gen] = JuMP.@constraint(model, Q_meas[gen] == Q_g[gen])
        eq_Vmeas[gen] = JuMP.@constraint(model, V_meas[gen] == V[bus])
        eq_Vset[gen] = JuMP.@constraint(model,
            V_set[gen] - mq * Q_meas[gen] - V_meas[gen] == 0)
    end

    vars[:P_meas] = P_meas
    vars[:Q_meas] = Q_meas
    vars[:V_meas] = V_meas
    vars[:E_int] = E_int
    vars[:V_set] = V_set
    eqs[:eq_const_gfm_Pmeas_init] = eq_Pmeas
    eqs[:eq_const_gfm_Qmeas_init] = eq_Qmeas
    eqs[:eq_const_gfm_Vmeas_init] = eq_Vmeas
    eqs[:eq_const_gfm_Vset_init] = eq_Vset
    attach_gfm_scalar_box!(model, dyn_model_dict, δ_gfm_only,
        specs[:gfm_δ][1], specs[:gfm_δ][2],
        :ineq_const_gfm_δ_lower, :ineq_const_gfm_δ_upper, :gfm_δ)
    dyn_model_dict[:meta][:gfm_phase] = :G1_init
    return dyn_model_dict
end

# ===================================================================================
# Transient GFM dynamics (Phase G2) — filters, droop, Q–V PI, current limiter
# ===================================================================================

function smooth_max_expr(model::Model, x, y; eps::Float64 = 1e-3)
    return JuMP.@expression(model, 0.5 * (x + y + sqrt((x - y)^2 + eps^2)))
end

function smooth_min_expr(model::Model, x, y; eps::Float64 = 1e-3)
    return JuMP.@expression(model, 0.5 * (x + y - sqrt((x - y)^2 + eps^2)))
end

function smooth_clip_expr(model::Model, x, lo, hi; eps::Float64 = 1e-3)
    clipped_hi = smooth_min_expr(model, x, hi; eps = eps)
    return smooth_max_expr(model, clipped_hi, lo; eps = eps)
end

function _gfm_safe_start(x, lo::Float64, hi::Float64, fallback::Float64)
    if x === nothing || !isfinite(Float64(x))
        return clamp(fallback, lo, hi)
    end
    return clamp(Float64(x), lo, hi)
end

function _gfm_start_of(x, fallback::Float64)
    try
        sx = JuMP.start_value(x)
        return sx === nothing ? fallback : Float64(sx)
    catch
        return fallback
    end
end

function _gfm_last_var(d::OrderedDict)
    return last(d)[2]
end

# ===================================================================================
# GFM variable boxes — resolved from TsBoundLimitsConfig + per-unit DGFM data
# ===================================================================================

"""
    resolve_gfm_bound_limits(limits, DGFM, gfm_gens) -> OrderedDict{Symbol, Tuple}

Per-unit `(min, max)` vectors for every GFM variable box, ordered like `gfm_gens`.

Margins and floors come from the `gfm_*` fields of `TsBoundLimitsConfig`; the per-unit
numbers they act on (`Emin`, `Emax`, `Imax`, `Pmin`, `Pmax`) come from `DGFM`. Shipped
defaults reproduce the values the reference port hard-coded, so a run configured with
`TsBoundLimitsConfig()` is unchanged by this indirection.

Keys `:gfm_δ`, `:gfm_P_meas`, `:gfm_Q_meas`, `:gfm_V_meas`, `:gfm_E_raw`, `:gfm_E_clip`,
`:gfm_I`, `:gfm_pq_raw`, `:gfm_pq_clip_hi`, `:gfm_pq_clip_lo` are boxes.
`:gfm_E_raw_start` / `:gfm_E_clip_start` are the tighter ranges used only to clamp
warm-start values — never emitted as bounds.

`:gfm_P_meas` and `:gfm_Q_meas` are solver guard-rails, not device ratings. The
ratings live on `DGFM` (`Pmin`/`Pmax`/`Qmin`/`Qmax`) and are comparator inputs to
the overload PIs. The Q measurement box is derived from the current rating: the
converter trades P against Q inside `|I|`, and that guard-rail is deliberately
wider than the `Qmin`/`Qmax` band.

These boxes are solver guard-rails, not physics: the physical converter limit is the
current-limiter equality. A converged solution sitting on one of them is a signal to
inspect the case, not a meaningful shadow price.
"""
function resolve_gfm_bound_limits(
    limits::TsBoundLimitsConfig,
    DGFM::DataFrame,
    gfm_gens::AbstractVector{Int},
)::OrderedDict{Symbol, Tuple{Any, Any}}
    n = length(gfm_gens)
    P_lo = Vector{Float64}(undef, n); P_hi = Vector{Float64}(undef, n)
    Q_lo = Vector{Float64}(undef, n); Q_hi = Vector{Float64}(undef, n)
    I_lo = Vector{Float64}(undef, n); I_hi = Vector{Float64}(undef, n)
    Eraw_lo = Vector{Float64}(undef, n); Eraw_hi = Vector{Float64}(undef, n)
    Eclip_lo = Vector{Float64}(undef, n); Eclip_hi = Vector{Float64}(undef, n)
    Eraw_s_lo = Vector{Float64}(undef, n); Eraw_s_hi = Vector{Float64}(undef, n)
    Eclip_s_lo = Vector{Float64}(undef, n); Eclip_s_hi = Vector{Float64}(undef, n)

    scale = limits.gfm_PQ_bound_scale
    offset = limits.gfm_PQ_bound_offset_pu

    for (i, gen) in enumerate(gfm_gens)
        r = dgfm_row(DGFM, gen)
        Emin = Float64(DGFM.Emin[r]); Emax = Float64(DGFM.Emax[r])
        Imax = Float64(DGFM.Imax[r])
        Pmin = Float64(DGFM.Pmin[r]); Pmax = Float64(DGFM.Pmax[r])

        # P box: widest of the floor and a margin over the set-point range.
        # A solver guard-rail, not the device rating — `DGFM.Pmin`/`Pmax` feed the
        # overload PI, and `DGEN` owns the dispatch box.
        P_bound = max(limits.gfm_P_meas_floor_pu - offset,
            scale * max(abs(Pmin) - offset, offset + abs(Pmax), 1.0))
        # Q box: derived from the current rating. The converter does have a reactive
        # nameplate (`DGFM.Qmin`/`Qmax`, the overload-PI comparators); this guard-rail
        # is wider than that band on purpose, so the PI — not the box — is what binds.
        Q_bound = max(limits.gfm_Q_meas_floor_pu - offset,
            scale * max(abs(Imax), 1.0))
        # Current box: a guard-rail around the limiter, clamped into [floor, ceiling].
        I_bound = min(limits.gfm_I_ceiling_pu,
            max(limits.gfm_I_floor_pu, scale * max(abs(Imax), 1.0)))

        P_lo[i] = -P_bound; P_hi[i] = P_bound
        Q_lo[i] = -Q_bound; Q_hi[i] = Q_bound
        I_lo[i] = -I_bound; I_hi[i] = I_bound

        # Raw (pre-clip) PI states get the widest band; the clipped states sit just
        # outside [Emin, Emax] so the smooth clip, not the box, is what saturates them.
        Eraw_s_lo[i] = Emin - limits.gfm_E_raw_extra_pu
        Eraw_s_hi[i] = Emax + limits.gfm_E_raw_extra_pu
        Eraw_lo[i] = Eraw_s_lo[i] - limits.gfm_E_raw_slack_pu
        Eraw_hi[i] = Eraw_s_hi[i] + limits.gfm_E_raw_slack_pu
        Eclip_lo[i] = Emin - limits.gfm_E_clip_slack_pu
        Eclip_hi[i] = Emax + limits.gfm_E_clip_slack_pu
        Eclip_s_lo[i] = Emin; Eclip_s_hi[i] = Emax
    end

    specs = OrderedDict{Symbol, Tuple{Any, Any}}()
    specs[:gfm_δ] = (fill(limits.gfm_δ_min_rad, n), fill(limits.gfm_δ_max_rad, n))
    specs[:gfm_P_meas] = (P_lo, P_hi)
    specs[:gfm_Q_meas] = (Q_lo, Q_hi)
    specs[:gfm_V_meas] = (fill(limits.gfm_V_meas_min_pu, n), fill(limits.gfm_V_meas_max_pu, n))
    specs[:gfm_E_raw] = (Eraw_lo, Eraw_hi)
    specs[:gfm_E_clip] = (Eclip_lo, Eclip_hi)
    specs[:gfm_I] = (I_lo, I_hi)
    specs[:gfm_E_raw_start] = (Eraw_s_lo, Eraw_s_hi)
    specs[:gfm_E_clip_start] = (Eclip_s_lo, Eclip_s_hi)
    # Same numbers on every converter: the overload states are deviations from 0,
    # not from a per-unit rating. Upper clip lives at or below 0, lower clip at or
    # above 0, each with `gfm_pq_clip_slack_pu` of room past the rail.
    span = limits.gfm_pq_raw_span_pu
    slack = limits.gfm_pq_clip_slack_pu
    specs[:gfm_pq_raw] = (fill(-span, n), fill(span, n))
    specs[:gfm_pq_clip_hi] = (fill(-span, n), fill(slack, n))
    specs[:gfm_pq_clip_lo] = (fill(-slack, n), fill(span, n))
    return specs
end

"""
Resolved GFM specs for this build.

Reads `meta[:var_limit_specs]` when the full pipeline registered them; otherwise
resolves from `TsBoundLimitsConfig()` defaults so the attachers stay callable from unit
tests that assemble a bare `dyn_model_dict`.
"""
function _gfm_limit_specs(
    dyn_model_dict::OrderedDict{Symbol, Any},
    DGFM::DataFrame,
    gfm_gens::AbstractVector{Int},
)::OrderedDict{Symbol, Tuple{Any, Any}}
    specs = get(dyn_model_dict[:meta], :var_limit_specs, nothing)
    (specs !== nothing && haskey(specs, :gfm_I)) && return specs
    return resolve_gfm_bound_limits(TsBoundLimitsConfig(), DGFM, gfm_gens)
end

"""Per-unit `(min, max)` from a resolved spec pair, by position in `gfm_gens`."""
_gfm_lim(specs, key::Symbol, i::Int) = _limit_for_gen(specs[key][1], specs[key][2], i)

"""
Is one GFM box switched on?

`TsBuilderConfig` defaults every `bound_gfm_*` toggle to **true** — these are
guard-rails the nonconvex limiter needs — so a missing entry (bare `dyn_model_dict`
in a unit test) also reads as on.
"""
function _gfm_box_enabled(meta::OrderedDict{Symbol, Any}, bound_key::Symbol)::Bool
    toggles = get(meta, :var_bounds, nothing)
    toggles === nothing && return true
    return get(toggles, bound_key, true)
end

"""
    attach_gfm_time_box!(model, dyn_model_dict, vars, min_lim, max_lim, key_lower, key_upper)

Materialise one GFM per-gen×time box under the model's `bound_encoding`.

`CONSTRAINT` (default) emits explicit ≤-rows into `dyn_model_dict[:ineq_const]`;
`VARIABLE` sets JuMP bounds and registers the container in the bound manifest. Either
way the dual is exported under the same key, because `extract_dual_entry` falls back to
the manifest when the constraint key is absent.
"""
function attach_gfm_time_box!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    vars,
    min_lim,
    max_lim,
    key_lower::Symbol,
    key_upper::Symbol,
    bound_key::Symbol,
)
    isempty(vars) && return nothing
    meta = dyn_model_dict[:meta]
    _gfm_box_enabled(meta, bound_key) || return nothing
    haskey(dyn_model_dict, :ineq_const) ||
        (dyn_model_dict[:ineq_const] = OrderedDict{Symbol, Any}())

    if bound_encoding_from_meta(meta) == VARIABLE
        apply_var_bounds_kron_gen_time!(vars, min_lim, max_lim)
        _use_lower_limit(min_lim) &&
            register_bound_manifest!(meta, key_lower, :lower, vars, :per_gen_time)
        _use_upper_limit(max_lim) &&
            register_bound_manifest!(meta, key_upper, :upper, vars, :per_gen_time)
        return nothing
    end

    lower, upper = ineq_const_kron_gen_time_bounds!(model, vars, min_lim, max_lim)
    lower !== nothing && (dyn_model_dict[:ineq_const][key_lower] = lower)
    upper !== nothing && (dyn_model_dict[:ineq_const][key_upper] = upper)
    return nothing
end

"""Pre-fault counterpart of `attach_gfm_time_box!` for per-generator scalars."""
function attach_gfm_scalar_box!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    vars::OrderedDict{Int, JuMP.VariableRef},
    min_lim,
    max_lim,
    key_lower::Symbol,
    key_upper::Symbol,
    bound_key::Symbol,
)
    isempty(vars) && return nothing
    meta = dyn_model_dict[:meta]
    _gfm_box_enabled(meta, bound_key) || return nothing
    haskey(dyn_model_dict, :ineq_const) ||
        (dyn_model_dict[:ineq_const] = OrderedDict{Symbol, Any}())

    if bound_encoding_from_meta(meta) == VARIABLE
        for (i, gen) in enumerate(keys(vars))
            min_g, max_g = _limit_for_gen(min_lim, max_lim, i)
            _use_lower_limit(min_g) && JuMP.set_lower_bound(vars[gen], Float64(min_g))
            _use_upper_limit(max_g) && JuMP.set_upper_bound(vars[gen], Float64(max_g))
        end
        _use_lower_limit(min_lim) &&
            register_bound_manifest!(meta, key_lower, :lower, vars, :gen_indexed)
        _use_upper_limit(max_lim) &&
            register_bound_manifest!(meta, key_upper, :upper, vars, :gen_indexed)
        return nothing
    end

    lower, upper = ineq_const_kron_gen_scalar_bounds!(model, vars, min_lim, max_lim)
    lower !== nothing && (dyn_model_dict[:ineq_const][key_lower] = lower)
    upper !== nothing && (dyn_model_dict[:ineq_const][key_upper] = upper)
    return nothing
end

"""GFM-only view of a shared per-gen×time container (`Id_tf`, `Iq_tf`, …)."""
function _gfm_view(vars, gfm_gens::AbstractVector{Int})
    out = OrderedDict{Int, Any}()
    for gen in gfm_gens
        haskey(vars, gen) && (out[gen] = vars[gen])
    end
    return out
end

"""
First-order measurement lag `Tf·ẋ = u − x`, discretised either way.

Backward Euler (`x_t = α·x_prev + β·u_t`, `α = Tf/(Tf+Δt)`, `β = Δt/(Tf+Δt)`) is
L-stable: it damps monotonically at any `Δt/Tf`. Trapezoidal is second-order accurate
but only A-stable — its per-step factor `(1−c)/(1+c)` with `c = Δt/(2Tf)` turns negative
once `Δt > 2·Tf`, which shows up as step-to-step ringing on the filtered signal. The
`Tf ≤ 1e-9` bypass comes first in both branches; it also guards the `2·Tf` division.

`u_prev` is only read on the trapezoidal branch.
"""
function _add_gfm_measurement_filter!(
    model::Model,
    x_t,
    x_prev,
    u_t,
    u_prev,
    Tf::Float64,
    Δt::Float64;
    backward_euler::Bool,
)
    if Tf <= 1e-9
        return JuMP.@constraint(model, x_t == u_t)
    end
    if backward_euler
        α_filter = Tf / (Tf + Δt)
        β_filter = Δt / (Tf + Δt)
        return JuMP.@constraint(model, x_t == α_filter * x_prev + β_filter * u_t)
    end
    # Trapezoidal, written in the same (1+c) / (1−c) / c shape as the AVR exciter and
    # governor rows so the three read alike:
    #   Tf·(x_t − x_prev) = (Δt/2)·[(u_t − x_t) + (u_prev − x_prev)]
    c = Δt / (2.0 * Tf)
    return JuMP.@constraint(model,
        x_t * (1 + c) - x_prev * (1 - c) - c * (u_t + u_prev) == 0.0)
end

function _add_gfm_current_limiter!(
    model::Model,
    Vmag,
    θbus,
    δgfm,
    E_droop,
    Xl::Float64,
    Imax::Float64,
    Id,
    Iq,
)
    Vd = JuMP.@expression(model, Vmag * sin(δgfm - θbus))
    Vq = JuMP.@expression(model, Vmag * cos(δgfm - θbus))
    ΔVd = JuMP.@expression(model, E_droop - Vq)
    ΔVq = JuMP.@expression(model, Vd)
    Imax_Xl = Imax * Xl
    bypassed = Imax >= _GFM_IMAX_NO_LIMIT

    if bypassed
        c_Id = JuMP.@constraint(model, Xl * Id == ΔVd)
        c_Iq = JuMP.@constraint(model, Xl * Iq == ΔVq)
        scale = JuMP.@expression(model, 1.0)
        Iraw_Xl = JuMP.@expression(model, sqrt(ΔVd^2 + ΔVq^2 + (_GFM_EPS_NORM_V * Xl)^2))
        Iout_Xl = JuMP.@expression(model, Xl * sqrt(Id^2 + Iq^2 + _GFM_EPS_NORM_V^2))
        den_Xl = Iraw_Xl
        Iraw = JuMP.@expression(model, Iraw_Xl / Xl)
        Iout = JuMP.@expression(model, Iout_Xl / Xl)
    else
        eps_lim_V = max(1e-8, _GFM_EPS_LIM_I * Xl)
        Iraw_Xl = JuMP.@expression(model, sqrt(ΔVd^2 + ΔVq^2 + (_GFM_EPS_NORM_V * Xl)^2))
        den_Xl = smooth_max_expr(model, Iraw_Xl, Imax_Xl; eps = eps_lim_V)
        scale = JuMP.@expression(model, Imax_Xl / den_Xl)
        c_Id = JuMP.@constraint(model, Xl * Id == scale * ΔVd)
        c_Iq = JuMP.@constraint(model, Xl * Iq == scale * ΔVq)
        Iraw = JuMP.@expression(model, Iraw_Xl / Xl)
        Iout = JuMP.@expression(model, sqrt(Id^2 + Iq^2 + _GFM_EPS_NORM_V^2))
        Iout_Xl = JuMP.@expression(model, Xl * Iout)
    end

    return (
        Vd = Vd, Vq = Vq, ΔVd = ΔVd, ΔVq = ΔVq,
        Iraw_Xl = Iraw_Xl, Iout_Xl = Iout_Xl, Imax_Xl = Imax_Xl, den_Xl = den_Xl,
        Iraw = Iraw, Iout = Iout, scale = scale, bypassed = bypassed,
        constraints = [c_Id, c_Iq], c_Id = c_Id, c_Iq = c_Iq,
    )
end

"""
Record and announce which GFM units had their current limiter bypassed.

The bypass swaps the clamp for unsaturated stator algebra, so the exported
`dual_gfm_limiter_*` rows mean something different for those units — worth one log line
rather than a silent change of model.
"""
function _record_gfm_limiter_bypass!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    bypassed_ids::Vector{Int},
    window::String,
)
    isempty(bypassed_ids) && return nothing
    meta = dyn_model_dict[:meta]
    known = get!(meta, :gfm_limiter_bypassed, Int[])
    union!(known, bypassed_ids)
    sort!(known)
    @info "GFM current limiter bypassed (Imax ≥ $_GFM_IMAX_NO_LIMIT pu): " *
          "gens $(bypassed_ids) in the $window window use unsaturated stator algebra."
    return nothing
end

"""Refresh `meta[:sg_gens]` / `meta[:gfm_gens]` from the current active set."""
function resolve_sg_gfm_gens!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::AbstractVector{Int},
    DGFM::Union{Nothing, DataFrame},
)::Tuple{Vector{Int64}, Vector{Int64}}
    sg = Vector{Int64}(sg_active_gens(active_gen, DGFM))
    gfm = Vector{Int64}(gfm_active_gens(active_gen, DGFM))
    # A converter-only fleet is opt-in. `validate_all_gfm_fleet!` has already rejected the
    # combinations that cannot work (`:swing_propagated`, `:highest_H`, a COI-referenced
    # corridor without `gfm_in_coi`, AVR/governor) before the warm-start solve, so by the
    # time we get here the flag alone is enough.
    if isempty(sg) && !get(dyn_model_dict[:meta], :allow_all_gfm, false)
        throw(ArgumentError(
            "DQ FULL_BUS requires at least one active SG (all-GFM fleets unsupported); " *
            "set DynModelConfig.allow_all_gfm=true to permit a converter-only fleet."))
    end
    dyn_model_dict[:meta][:sg_gens] = sg
    dyn_model_dict[:meta][:gfm_gens] = gfm
    return sg, gfm
end

"""
    _gfm_pq_eps(ki_sys, Δt, mach_base_MVA, sys_base_MVA) -> Float64

State-clip ε = c · Δt · Ki · E_ref (`_GFM_PQLIM_EPS_STEP_FRAC`), floored at
`_GFM_PQLIM_EPS_FLOOR`.

`ki_sys` is the system-base gain that `apply_gfm_base_conversion!` leaves on `DGFM`
(machine gain × S_sys/S_mach). E_ref is 1 pu on the **converter's** base, which is
S_mach/S_sys in system pu, so the product is the machine-base gain. The clip state
carries Δω (P axis) or voltage (Q axis), and neither depends on the MVA base, so the
base ratio does not cancel. Measuring the error against the converter's own rating
gives a small and a large unit with the same per-unit tuning the same ε.
"""
function _gfm_pq_eps(ki_sys::Float64, Δt::Float64,
                     mach_base_MVA::Float64, sys_base_MVA::Float64)
    ki_mach = ki_sys * mach_base_MVA / sys_base_MVA
    return max(_GFM_PQLIM_EPS_STEP_FRAC * Δt * abs(ki_mach), _GFM_PQLIM_EPS_FLOOR)
end

"""
    _gfm_pq_eps_out(kp_sys, mach_base_MVA, sys_base_MVA) -> Float64

Output-clip ε = Kp · tol (`_GFM_PQLIM_OUT_ERR_TOL`), floored at
`_GFM_PQLIM_EPS_FLOOR`. `kp_sys` is put back on the converter's base for the same
reason as in `_gfm_pq_eps`, so `tol` is a fraction of the unit's own rating.
"""
function _gfm_pq_eps_out(kp_sys::Float64, mach_base_MVA::Float64, sys_base_MVA::Float64)
    kp_mach = kp_sys * mach_base_MVA / sys_base_MVA
    return max(abs(kp_mach) * _GFM_PQLIM_OUT_ERR_TOL, _GFM_PQLIM_EPS_FLOOR)
end

"""
One REGFM_A1 Fig. 3 overload branch.

`upper` clamps the integrator state and the branch output from above at 0
(`smooth_min`); the lower branch clamps from below (`smooth_max`). The integrator
advances from the *clipped* previous state — the same anti-windup `E_int` uses.
Integrating the raw error and clipping only the output has no finite equilibrium:
inside the band the error is strictly positive and the state runs to +∞.

`eps` smooths the state clip and `eps_out` the output clip; see
`_GFM_PQLIM_EPS_STEP_FRAC` and `_GFM_PQLIM_OUT_ERR_TOL` for why they differ.
"""
function _add_gfm_overload_branch!(
    model::Model,
    x_raw,
    x_clip,
    x_prev,
    e_t,
    e_prev,
    ki::Float64,
    kp::Float64,
    Δt::Float64;
    backward_euler::Bool,
    upper::Bool,
    eps::Float64,
    eps_out::Float64,
)
    c_raw = if backward_euler
        JuMP.@constraint(model, x_raw - x_prev == Δt * ki * e_t)
    else
        JuMP.@constraint(model, x_raw - x_prev == (Δt / 2.0) * ki * (e_t + e_prev))
    end
    if upper
        c_clip = JuMP.@constraint(model,
            x_clip == smooth_min_expr(model, x_raw, 0.0; eps = eps))
    else
        c_clip = JuMP.@constraint(model,
            x_clip == smooth_max_expr(model, x_raw, 0.0; eps = eps))
    end
    u = _gfm_overload_output(model, kp, e_t, x_clip; upper = upper, eps = eps_out)
    return (c_raw = c_raw, c_clip = c_clip, u = u)
end

"""Branch output `u` from an already-clipped state. Used by the post-fault anchor."""
function _gfm_overload_output(model::Model, kp::Float64, e, x_clip;
                              upper::Bool, eps::Float64)
    if upper
        return smooth_min_expr(model, kp * e + x_clip, 0.0; eps = eps)
    end
    return smooth_max_expr(model, kp * e + x_clip, 0.0; eps = eps)
end

"""
Clipped previous state and previous error for one overload branch.

`t = 1` of the fault window starts from `x(0) = 0`. `t = 1` of the post-fault
window starts from the last fault-window *clipped* state. Later steps use the
clipped state at `t - 1`, never the raw one.
"""
function _gfm_pq_prev(t::Int, suffix::String, vars, gen::Int, clip_tf_key::Symbol,
                      clip_w, meas_w, limit::Float64, anchor_meas)
    if t == 1
        x_prev = suffix == "tf" ? 0.0 : _gfm_last_var(vars[clip_tf_key][gen])
        return x_prev, limit - anchor_meas
    end
    return clip_w[gen][t - 1], limit - meas_w[gen][t - 1]
end

function _gfm_pq_state!(model::Model, encoding::BoundEncoding, base::String,
                        suffix::String, gen::Int, t::Int, lo, hi, on::Bool)
    return create_bounded_variable!(model;
        base_name = "$(base)_$(suffix)[$gen,$t]",
        start = 0.0,
        min_lim = on ? lo : nothing,
        max_lim = on ? hi : nothing,
        encoding = encoding)
end

"""
    _attach_gfm_window!(model, dyn_model_dict, gfm_gens, DGEN, DGFM, Δt, ω_syn,
                        time_window, suffix, anchor, input_anchor)

Shared body of both GFM transient windows (measurement filters, P–f droop and δ
integration, Q–V PI with smooth clip, angle-preserving current limiter, Pe/Qe injection).

`suffix` is `"tf"` or `"tpf"` and selects the shared time-series containers.
`anchor(gen)` returns the `(P, Q, V, E, δ)` **state** references the window starts from —
the pre-fault scalars for `tf`, the last fault-window step for `tpf`. They serve both as
the `t = 1` previous values and as the warm-start seeds, which is the only structural
difference between the two windows. `input_anchor(gen)` returns the matching
`(Pe, Qe, V_bus)` **inputs** one step before the window, needed only by the trapezoidal
filter branch at `t = 1`.

Two independent dials govern the discretisation:

- `meta[:ode_first_step]` drives the δ integration exactly as it drives the SG swing,
  dq EMF, AVR and governor: `:backward_euler` makes `t = 1` of each window BE, everything
  after stays trapezoidal (`:trapezoidal`, the default, is trapezoidal throughout).
- `meta[:gfm_integrator]` drives the three measurement filters and the `E_int`
  integrator, and is GFM-only: `:backward_euler` (the default) is BE at every step —
  the reference implementation's scheme, and the only one whose per-iteration cost
  stays flat as the horizon grows; `:trapezoidal` is trapezoidal at every step; and
  `:follow_ode_first_step` applies the SG rule above. See `DynModelConfig` for the
  measured cost/accuracy trade-off behind that default.
"""
function _attach_gfm_window!(
    model::Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    gfm_gens::AbstractVector{Int},
    DGEN::DataFrame,
    DGFM::DataFrame,
    Δt::Float64,
    ω_syn::Float64,
    time_window::Vector{Float64},
    suffix::String,
    anchor::Function,
    input_anchor::Function,
)
    vars = dyn_model_dict[:vars]
    eqs = dyn_model_dict[:eq_const]
    meta = dyn_model_dict[:meta]
    encoding = bound_encoding_from_meta(meta)
    specs = _gfm_limit_specs(dyn_model_dict, DGFM, gfm_gens)
    ode_fs = get(meta, :ode_first_step, :trapezoidal)
    gfm_int = get(meta, :gfm_integrator, :backward_euler)
    pq_on = get(meta, :gfm_pq_overload, false)
    eq_fams = pq_on ? (_GFM_EQ_FAMILIES..., _GFM_PQ_EQ_FAMILIES...) : _GFM_EQ_FAMILIES

    # Single resolver for the filters and the Q–V PI, so the two families can never
    # drift apart. δ deliberately does not go through it — it stays on `ode_fs`.
    be_at(t::Int) = gfm_int === :backward_euler ||
        (gfm_int === :follow_ode_first_step && t == 1 && ode_fs === :backward_euler)

    # Trapezoidal is only A-stable: its per-step factor (1 − Δt/2Tf)/(1 + Δt/2Tf) turns
    # negative once Δt > 2·Tf, which shows up as a ringing filter output feeding the
    # droop → δ → limiter chain. Backward Euler is L-stable and never does this. Warned
    # once per run (the fault window is built first).
    # The criterion is per filter, so the binding constant is the smallest of all three
    # channels across all converters, not one shared value per unit.
    if suffix == "tf" && gfm_int !== :backward_euler && !isempty(gfm_gens)
        Tf_min = minimum(
            min(Float64(DGFM.TPf[dgfm_row(DGFM, g)]),
                Float64(DGFM.TQf[dgfm_row(DGFM, g)]),
                Float64(DGFM.TVf[dgfm_row(DGFM, g)]))
            for g in gfm_gens)
        if Tf_min > 1e-9 && Δt > 2.0 * Tf_min
            @warn "GFM measurement filters are trapezoidal at t_step=$(Δt) s against a " *
                  "fastest filter Tf=$(Tf_min) s. Above Δt = 2·Tf the trapezoidal update " *
                  "oscillates step-to-step; use gfm_integrator=:backward_euler or reduce t_step."
        end
    end

    # Per-family toggles (default true). A disabled box must also skip the JuMP
    # bound at creation, or `VARIABLE` encoding would keep it regardless.
    box_on(name::String) = _gfm_box_enabled(meta, Symbol("gfm_", name, "_", suffix))
    on_P = box_on("P_meas");        on_Q = box_on("Q_meas");   on_V = box_on("V_meas")
    on_Eir = box_on("E_int_raw");   on_Ei = box_on("E_int")
    on_Edr = box_on("E_droop_raw"); on_Ed = box_on("E_droop")
    # One switch per window for every overload branch. Missing key reads as on.
    on_pq_raw = pq_on && _gfm_box_enabled(meta, Symbol("gfm_pq_raw_", suffix))
    on_pq_clip = pq_on && _gfm_box_enabled(meta, Symbol("gfm_pq_clip_", suffix))

    shared(name::String) = vars[Symbol(name, "_", suffix)]
    V_w = shared("V"); θ_w = shared("θ")
    Pe_w = shared("Pe"); Qe_w = shared("Qe")
    δ_w = shared("δ"); Δω_w = shared("Δω")
    Id_w = shared("Id"); Iq_w = shared("Iq")
    V_set = vars[:V_set]
    P_set = vars[:P_m]

    P_meas_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Q_meas_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    V_meas_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    E_int_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    E_int_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    E_droop_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    E_droop_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Pmax_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Pmax_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Pmin_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Pmin_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Qmax_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Qmax_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Qmin_raw_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    x_Qmin_w = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()

    eq_flat = OrderedDict{Int, OrderedDict{Int, Vector{JuMP.ConstraintRef}}}()
    eq_named = OrderedDict{Symbol, OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}}(
        fam => OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
        for fam in eq_fams)
    bypassed = Int[]

    for (i, gen) in enumerate(gfm_gens)
        r = dgfm_row(DGFM, gen)
        bus = Int(DGEN.bus[gen])
        TPf = Float64(DGFM.TPf[r]); TQf = Float64(DGFM.TQf[r]); TVf = Float64(DGFM.TVf[r])
        mp = Float64(DGFM.mp[r]); mq = Float64(DGFM.mq[r])
        Kpv = Float64(DGFM.Kpv[r]); Kiv = Float64(DGFM.Kiv[r])
        Emax = Float64(DGFM.Emax[r]); Emin = Float64(DGFM.Emin[r])
        Imax = Float64(DGFM.Imax[r]); Xl = Float64(DGFM.Xl[r])
        # Gains and ratings are only read when the knob is on, so a hand-built
        # frame that predates the columns still builds the historical loops.
        Pmax_lim = 0.0; Pmin_lim = 0.0; Qmax_lim = 0.0; Qmin_lim = 0.0
        Kpp = 0.0; Kip = 0.0; Kpq = 0.0; Kiq = 0.0
        praw_lo = praw_hi = phi_lo = phi_hi = plo_lo = plo_hi = 0.0
        eps_P = eps_Q = eps_out_P = eps_out_Q = _GFM_PQLIM_EPS_FLOOR
        if pq_on
            Pmax_lim = Float64(DGFM.Pmax[r]); Pmin_lim = Float64(DGFM.Pmin[r])
            Qmax_lim = Float64(DGFM.Qmax[r]); Qmin_lim = Float64(DGFM.Qmin[r])
            Kpp = Float64(DGFM.Kppmax[r]); Kip = Float64(DGFM.Kipmax[r])
            Kpq = Float64(DGFM.Kpqmax[r]); Kiq = Float64(DGFM.Kiqmax[r])
            praw_lo, praw_hi = _gfm_lim(specs, :gfm_pq_raw, i)
            phi_lo, phi_hi = _gfm_lim(specs, :gfm_pq_clip_hi, i)
            plo_lo, plo_hi = _gfm_lim(specs, :gfm_pq_clip_lo, i)
            # Per-axis ε: state clip from the integrator step (`_GFM_PQLIM_EPS_STEP_FRAC`),
            # output clip from the proportional gain (`_GFM_PQLIM_OUT_ERR_TOL`).
            # Upper and lower branches share the axis gains, so they share ε too. A bare
            # `dyn_model_dict` without `:base_MVA` treats the converter as system-rated.
            sys_base = Float64(get(meta, :base_MVA, DGFM.mach_base_MVA[r]))
            mach_base = Float64(DGFM.mach_base_MVA[r])
            eps_P = _gfm_pq_eps(Kip, Δt, mach_base, sys_base)
            eps_Q = _gfm_pq_eps(Kiq, Δt, mach_base, sys_base)
            eps_out_P = _gfm_pq_eps_out(Kpp, mach_base, sys_base)
            eps_out_Q = _gfm_pq_eps_out(Kpq, mach_base, sys_base)
            eps_rec = get!(meta, :gfm_pq_eps, OrderedDict{Int, NamedTuple}())
            eps_rec[gen] = (P = eps_P, Q = eps_Q, P_out = eps_out_P, Q_out = eps_out_Q)
            for d in (x_Pmax_raw_w, x_Pmax_w, x_Pmin_raw_w, x_Pmin_w,
                      x_Qmax_raw_w, x_Qmax_w, x_Qmin_raw_w, x_Qmin_w)
                d[gen] = OrderedDict{Int, JuMP.VariableRef}()
            end
        end

        P_lo, P_hi = _gfm_lim(specs, :gfm_P_meas, i)
        Q_lo, Q_hi = _gfm_lim(specs, :gfm_Q_meas, i)
        V_lo, V_hi = _gfm_lim(specs, :gfm_V_meas, i)
        Eraw_lo, Eraw_hi = _gfm_lim(specs, :gfm_E_raw, i)
        Eclip_lo, Eclip_hi = _gfm_lim(specs, :gfm_E_clip, i)
        Eraw_s_lo, Eraw_s_hi = _gfm_lim(specs, :gfm_E_raw_start, i)
        Eclip_s_lo, Eclip_s_hi = _gfm_lim(specs, :gfm_E_clip_start, i)

        for d in (P_meas_w, Q_meas_w, V_meas_w,
                  E_int_raw_w, E_int_w, E_droop_raw_w, E_droop_w)
            d[gen] = OrderedDict{Int, JuMP.VariableRef}()
        end
        eq_flat[gen] = OrderedDict{Int, Vector{JuMP.ConstraintRef}}()
        for fam in eq_fams
            eq_named[fam][gen] = OrderedDict{Int, JuMP.ConstraintRef}()
        end

        P_anchor, Q_anchor, V_anchor, E_anchor, δ_anchor = anchor(gen)
        # Filter *inputs* one step before the window. In `tf` these coincide with the
        # state anchors above and the trapezoidal correction term vanishes identically:
        # `eq_const_gfm_{Pmeas,Qmeas,Vmeas}_init` pin `P_meas = P_g`, `Q_meas = Q_g`,
        # `V_meas = V[bus]`, and `eq_const_P_init` / `eq_const_Q_init` make `P_g`/`Q_g`
        # the pre-fault `Pe`/`Qe`. In `tpf` they are the genuinely distinct last
        # fault-window `Pe_tf` / `Qe_tf` / `V_tf` values.
        Pe_anchor, Qe_anchor, Vbus_anchor = input_anchor(gen)
        P0 = _gfm_start_of(P_anchor, 0.0)
        Q0 = _gfm_start_of(Q_anchor, 0.0)
        V0 = _gfm_start_of(V_anchor, 1.0)
        E0 = _gfm_start_of(E_anchor, 1.0)

        # Trapezoidal `t = 1` needs the step *before* the window. Both follow from the
        # anchors: the droop law and the voltage error, evaluated one step back. In the
        # fault window each is identically zero (the init equalities pin
        # `P_meas = P_g = P_set` and `V_set - mq·Q_meas - V_meas = 0`, and every
        # overload branch starts at `x(0) = 0`, so its output is zero inside the
        # band); in the post-fault window they reproduce the last fault-window `Δω_tf`
        # and error exactly, because those are defined by the same equations.
        # δ follows `ode_fs`; the Q–V PI and the overload integrators follow the GFM
        # resolver, so the anchors are built under separate conditions.
        Δω_anchor = ode_fs === :backward_euler ? nothing :
            JuMP.@expression(model, mp * (P_set[gen] - P_anchor))
        err_anchor = be_at(1) ? nothing :
            JuMP.@expression(model, V_set[gen] - mq * Q_anchor - V_anchor)
        # Post-fault trapezoidal `t = 1` must carry the last fault-window overload
        # correction. The fault window does not: `x(0) = 0` and the pre-fault point
        # sits inside the band, so every `u` is zero there.
        if pq_on && suffix == "tpf" && ode_fs !== :backward_euler
            u_Pmax_a = _gfm_overload_output(model, Kpp, Pmax_lim - P_anchor,
                _gfm_last_var(vars[:x_Pmax_tf][gen]); upper = true, eps = eps_out_P)
            u_Pmin_a = _gfm_overload_output(model, Kpp, Pmin_lim - P_anchor,
                _gfm_last_var(vars[:x_Pmin_tf][gen]); upper = false, eps = eps_out_P)
            Δω_anchor = JuMP.@expression(model,
                mp * (P_set[gen] - P_anchor) + u_Pmax_a + u_Pmin_a)
        end
        if pq_on && suffix == "tpf" && !be_at(1)
            u_Qmax_a = _gfm_overload_output(model, Kpq, Qmax_lim - Q_anchor,
                _gfm_last_var(vars[:x_Qmax_tf][gen]); upper = true, eps = eps_out_Q)
            u_Qmin_a = _gfm_overload_output(model, Kpq, Qmin_lim - Q_anchor,
                _gfm_last_var(vars[:x_Qmin_tf][gen]); upper = false, eps = eps_out_Q)
            err_anchor = JuMP.@expression(model,
                V_set[gen] - mq * Q_anchor + u_Qmax_a + u_Qmin_a - V_anchor)
        end

        # Voltage error one step back, carried across the `t` loop so the trapezoidal
        # Q–V PI does not rebuild the previous step's expression.
        err_prev = err_anchor

        for t in eachindex(time_window)
            eq_flat[gen][t] = JuMP.ConstraintRef[]
            # Every constraint lands in both the flat window vector (text appendix) and
            # its named family (dual registry).
            function record!(fam::Symbol, cref)
                push!(eq_flat[gen][t], cref)
                eq_named[fam][gen][t] = cref
                return cref
            end

            # Start values stay clamped to the resolved box even when a toggle is
            # off — that is a warm-start heuristic, not a constraint.
            P_meas_w[gen][t] = create_bounded_variable!(model;
                base_name = "P_meas_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(P0, P_lo, P_hi, 0.0),
                min_lim = on_P ? P_lo : nothing, max_lim = on_P ? P_hi : nothing,
                encoding = encoding)
            Q_meas_w[gen][t] = create_bounded_variable!(model;
                base_name = "Q_meas_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(Q0, Q_lo, Q_hi, 0.0),
                min_lim = on_Q ? Q_lo : nothing, max_lim = on_Q ? Q_hi : nothing,
                encoding = encoding)
            V_meas_w[gen][t] = create_bounded_variable!(model;
                base_name = "V_meas_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(V0, V_lo, V_hi, 1.0),
                min_lim = on_V ? V_lo : nothing, max_lim = on_V ? V_hi : nothing,
                encoding = encoding)
            E_int_raw_w[gen][t] = create_bounded_variable!(model;
                base_name = "E_int_raw_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(E0, Eraw_s_lo, Eraw_s_hi, 1.0),
                min_lim = on_Eir ? Eraw_lo : nothing, max_lim = on_Eir ? Eraw_hi : nothing,
                encoding = encoding)
            E_int_w[gen][t] = create_bounded_variable!(model;
                base_name = "E_int_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(E0, Eclip_s_lo, Eclip_s_hi, 1.0),
                min_lim = on_Ei ? Eclip_lo : nothing, max_lim = on_Ei ? Eclip_hi : nothing,
                encoding = encoding)
            E_droop_raw_w[gen][t] = create_bounded_variable!(model;
                base_name = "E_droop_raw_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(E0, Eraw_s_lo, Eraw_s_hi, 1.0),
                min_lim = on_Edr ? Eraw_lo : nothing, max_lim = on_Edr ? Eraw_hi : nothing,
                encoding = encoding)
            E_droop_w[gen][t] = create_bounded_variable!(model;
                base_name = "E_droop_$(suffix)[$gen,$t]",
                start = _gfm_safe_start(E0, Eclip_s_lo, Eclip_s_hi, 1.0),
                min_lim = on_Ed ? Eclip_lo : nothing, max_lim = on_Ed ? Eclip_hi : nothing,
                encoding = encoding)

            if t == 1
                P_prev = P_anchor; Q_prev = Q_anchor; V_prev = V_anchor
                E_prev = E_anchor; δ_prev = δ_anchor
                Pe_prev = Pe_anchor; Qe_prev = Qe_anchor; Vbus_prev = Vbus_anchor
            else
                P_prev = P_meas_w[gen][t - 1]; Q_prev = Q_meas_w[gen][t - 1]
                V_prev = V_meas_w[gen][t - 1]; E_prev = E_int_w[gen][t - 1]
                δ_prev = δ_w[gen][t - 1]
                Pe_prev = Pe_w[gen][t - 1]; Qe_prev = Qe_w[gen][t - 1]
                Vbus_prev = V_w[bus][t - 1]
            end

            filter_be = be_at(t)
            record!(:filter_P, _add_gfm_measurement_filter!(
                model, P_meas_w[gen][t], P_prev, Pe_w[gen][t], Pe_prev, TPf, Δt;
                backward_euler = filter_be))
            record!(:filter_Q, _add_gfm_measurement_filter!(
                model, Q_meas_w[gen][t], Q_prev, Qe_w[gen][t], Qe_prev, TQf, Δt;
                backward_euler = filter_be))
            record!(:filter_V, _add_gfm_measurement_filter!(
                model, V_meas_w[gen][t], V_prev, V_w[bus][t], Vbus_prev, TVf, Δt;
                backward_euler = filter_be))

            if pq_on
                # Same gains on the upper and lower branch of each axis. `u` stays an
                # expression inside the droop row and the voltage error — no variable,
                # no dual. The raw and clipped states are both variables: an expression
                # for the clip would nest across every step of the window.
                u_ol = Dict{Symbol, Any}()
                for b in (
                    (tag = :Pmax, fam_raw = :Pmax_raw, fam_clip = :Pmax_clip,
                     raw_w = x_Pmax_raw_w, clip_w = x_Pmax_w, clip_key = :x_Pmax_tf,
                     meas_w = P_meas_w, anchor_meas = P_anchor, limit = Pmax_lim,
                     kp = Kpp, ki = Kip, eps = eps_P, eps_out = eps_out_P, upper = true,
                     clip_lo = phi_lo, clip_hi = phi_hi),
                    (tag = :Pmin, fam_raw = :Pmin_raw, fam_clip = :Pmin_clip,
                     raw_w = x_Pmin_raw_w, clip_w = x_Pmin_w, clip_key = :x_Pmin_tf,
                     meas_w = P_meas_w, anchor_meas = P_anchor, limit = Pmin_lim,
                     kp = Kpp, ki = Kip, eps = eps_P, eps_out = eps_out_P, upper = false,
                     clip_lo = plo_lo, clip_hi = plo_hi),
                    (tag = :Qmax, fam_raw = :Qmax_raw, fam_clip = :Qmax_clip,
                     raw_w = x_Qmax_raw_w, clip_w = x_Qmax_w, clip_key = :x_Qmax_tf,
                     meas_w = Q_meas_w, anchor_meas = Q_anchor, limit = Qmax_lim,
                     kp = Kpq, ki = Kiq, eps = eps_Q, eps_out = eps_out_Q, upper = true,
                     clip_lo = phi_lo, clip_hi = phi_hi),
                    (tag = :Qmin, fam_raw = :Qmin_raw, fam_clip = :Qmin_clip,
                     raw_w = x_Qmin_raw_w, clip_w = x_Qmin_w, clip_key = :x_Qmin_tf,
                     meas_w = Q_meas_w, anchor_meas = Q_anchor, limit = Qmin_lim,
                     kp = Kpq, ki = Kiq, eps = eps_Q, eps_out = eps_out_Q, upper = false,
                     clip_lo = plo_lo, clip_hi = plo_hi),
                )
                    b.raw_w[gen][t] = _gfm_pq_state!(model, encoding,
                        "x_$(b.tag)_raw", suffix, gen, t, praw_lo, praw_hi, on_pq_raw)
                    b.clip_w[gen][t] = _gfm_pq_state!(model, encoding,
                        "x_$(b.tag)", suffix, gen, t, b.clip_lo, b.clip_hi, on_pq_clip)
                    x_prev, e_prev = _gfm_pq_prev(t, suffix, vars, gen, b.clip_key,
                        b.clip_w, b.meas_w, b.limit, b.anchor_meas)
                    e_t = b.limit - b.meas_w[gen][t]
                    br = _add_gfm_overload_branch!(model,
                        b.raw_w[gen][t], b.clip_w[gen][t], x_prev, e_t, e_prev,
                        b.ki, b.kp, Δt; backward_euler = be_at(t), upper = b.upper,
                        eps = b.eps, eps_out = b.eps_out)
                    record!(b.fam_raw, br.c_raw)
                    record!(b.fam_clip, br.c_clip)
                    u_ol[b.tag] = br.u
                end
                record!(:droop, JuMP.@constraint(model,
                    Δω_w[gen][t] == mp * (P_set[gen] - P_meas_w[gen][t]) +
                        u_ol[:Pmax] + u_ol[:Pmin]))
            else
                record!(:droop, JuMP.@constraint(model,
                    Δω_w[gen][t] == mp * (P_set[gen] - P_meas_w[gen][t])))
            end
            if t == 1 && ode_fs === :backward_euler
                record!(:delta, JuMP.@constraint(model,
                    δ_w[gen][t] - δ_prev == Δt * ω_syn * Δω_w[gen][t]))
            elseif t == 1
                record!(:delta, JuMP.@constraint(model,
                    δ_w[gen][t] - δ_prev ==
                    (Δt / 2.0) * ω_syn * (Δω_w[gen][t] + Δω_anchor)))
            else
                record!(:delta, JuMP.@constraint(model,
                    δ_w[gen][t] - δ_prev ==
                    (Δt / 2.0) * ω_syn * (Δω_w[gen][t] + Δω_w[gen][t - 1])))
            end

            voltage_error = pq_on ?
                JuMP.@expression(model,
                    V_set[gen] - mq * Q_meas_w[gen][t] + u_ol[:Qmax] + u_ol[:Qmin] -
                        V_meas_w[gen][t]) :
                JuMP.@expression(model,
                    V_set[gen] - mq * Q_meas_w[gen][t] - V_meas_w[gen][t])
            # `E_prev` stays the *clipped* `E_int`, the anti-windup convention, under
            # both schemes; only the integrand average changes.
            if be_at(t)
                record!(:Eint_raw, JuMP.@constraint(model,
                    E_int_raw_w[gen][t] - E_prev == Δt * Kiv * voltage_error))
            else
                record!(:Eint_raw, JuMP.@constraint(model,
                    E_int_raw_w[gen][t] - E_prev ==
                    (Δt / 2.0) * Kiv * (voltage_error + err_prev)))
            end
            err_prev = voltage_error
            record!(:Eint_clip, JuMP.@constraint(model,
                E_int_w[gen][t] == smooth_clip_expr(
                    model, E_int_raw_w[gen][t], Emin, Emax; eps = _GFM_EPS_E)))
            record!(:Edroop_raw, JuMP.@constraint(model,
                E_droop_raw_w[gen][t] == E_int_w[gen][t] + Kpv * voltage_error))
            record!(:Edroop_clip, JuMP.@constraint(model,
                E_droop_w[gen][t] == smooth_clip_expr(
                    model, E_droop_raw_w[gen][t], Emin, Emax; eps = _GFM_EPS_E)))

            lim = _add_gfm_current_limiter!(
                model, V_w[bus][t], θ_w[bus][t], δ_w[gen][t], E_droop_w[gen][t],
                Xl, Imax, Id_w[gen][t], Iq_w[gen][t])
            record!(:limiter_Id, lim.c_Id)
            record!(:limiter_Iq, lim.c_Iq)
            t == 1 && lim.bypassed && push!(bypassed, gen)
            record!(:Pe, JuMP.@constraint(model,
                Pe_w[gen][t] == lim.Vd * Id_w[gen][t] + lim.Vq * Iq_w[gen][t]))
            record!(:Qe, JuMP.@constraint(model,
                Qe_w[gen][t] == lim.Vq * Id_w[gen][t] - lim.Vd * Iq_w[gen][t]))
        end
    end

    for (name, container) in (
        ("P_meas", P_meas_w), ("Q_meas", Q_meas_w), ("V_meas", V_meas_w),
        ("E_int_raw", E_int_raw_w), ("E_int", E_int_w),
        ("E_droop_raw", E_droop_raw_w), ("E_droop", E_droop_w),
    )
        vars[Symbol(name, "_", suffix)] = container
    end
    if pq_on
        for (name, container) in (
            ("x_Pmax_raw", x_Pmax_raw_w), ("x_Pmax", x_Pmax_w),
            ("x_Pmin_raw", x_Pmin_raw_w), ("x_Pmin", x_Pmin_w),
            ("x_Qmax_raw", x_Qmax_raw_w), ("x_Qmax", x_Qmax_w),
            ("x_Qmin_raw", x_Qmin_raw_w), ("x_Qmin", x_Qmin_w),
        )
            vars[Symbol(name, "_", suffix)] = container
        end
    end

    # Boxes: explicit ≤-rows under CONSTRAINT, JuMP bounds + manifest under VARIABLE.
    # Id/Iq are shared with the SG path, so only the GFM ids are touched here.
    box_targets = (
        (P_meas_w, :gfm_P_meas, "P_meas"),
        (Q_meas_w, :gfm_Q_meas, "Q_meas"),
        (V_meas_w, :gfm_V_meas, "V_meas"),
        (E_int_raw_w, :gfm_E_raw, "E_int_raw"),
        (E_int_w, :gfm_E_clip, "E_int"),
        (E_droop_raw_w, :gfm_E_raw, "E_droop_raw"),
        (E_droop_w, :gfm_E_clip, "E_droop"),
        (_gfm_view(Id_w, gfm_gens), :gfm_I, "Id"),
        (_gfm_view(Iq_w, gfm_gens), :gfm_I, "Iq"),
    )
    for (container, key, name) in box_targets
        min_lim, max_lim = specs[key]
        attach_gfm_time_box!(model, dyn_model_dict, container, min_lim, max_lim,
            Symbol("ineq_const_gfm_", name, "_", suffix, "_lower"),
            Symbol("ineq_const_gfm_", name, "_", suffix, "_upper"),
            Symbol("gfm_", name, "_", suffix))
    end
    if pq_on
        # Upper clips are held at or below 0, lower clips at or above 0. The raw
        # columns share one symmetric box. `kind` selects the shared toggle.
        for (container, key, name, kind) in (
            (x_Pmax_raw_w, :gfm_pq_raw, "Pmax_raw", "raw"),
            (x_Pmax_w, :gfm_pq_clip_hi, "Pmax_clip", "clip"),
            (x_Pmin_raw_w, :gfm_pq_raw, "Pmin_raw", "raw"),
            (x_Pmin_w, :gfm_pq_clip_lo, "Pmin_clip", "clip"),
            (x_Qmax_raw_w, :gfm_pq_raw, "Qmax_raw", "raw"),
            (x_Qmax_w, :gfm_pq_clip_hi, "Qmax_clip", "clip"),
            (x_Qmin_raw_w, :gfm_pq_raw, "Qmin_raw", "raw"),
            (x_Qmin_w, :gfm_pq_clip_lo, "Qmin_clip", "clip"),
        )
            min_lim, max_lim = specs[key]
            attach_gfm_time_box!(model, dyn_model_dict, container, min_lim, max_lim,
                Symbol("ineq_const_gfm_", name, "_", suffix, "_lower"),
                Symbol("ineq_const_gfm_", name, "_", suffix, "_upper"),
                Symbol("gfm_pq_", kind, "_", suffix))
        end
    end

    eqs[Symbol("eq_const_gfm_", suffix)] = eq_flat
    for fam in eq_fams
        eqs[Symbol("eq_const_gfm_", fam, "_", suffix)] = eq_named[fam]
    end
    _record_gfm_limiter_bypass!(dyn_model_dict, bypassed, suffix)
    dyn_model_dict[:meta][:gfm_phase] = suffix == "tf" ? :G2_fault : :G2_postf
    return dyn_model_dict
end

"""
    Attach_GFM_fault!(model, dyn_model_dict, gfm_gens, DGEN, DGFM, Δt, ω_syn, time_window)

Fault-window GFM dynamics. Requires shared time-series vars (`Pe_tf`, `δ_tf`, …) already
created for `gfm_gens`, and the pre-fault scalars from `Attach_GFM_init!`.
"""
function Attach_GFM_fault!(
    model::Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    gfm_gens::AbstractVector{Int},
    DGEN::DataFrame,
    DGFM::DataFrame,
    Δt::Float64,
    ω_syn::Float64,
    time_window::Vector{Float64},
)
    register_gfm_meta!(dyn_model_dict, DGFM)
    isempty(gfm_gens) && return dyn_model_dict
    vars = dyn_model_dict[:vars]
    anchor = gen -> (vars[:P_meas][gen], vars[:Q_meas][gen], vars[:V_meas][gen],
                     vars[:E_int][gen], vars[:δ][gen])
    # Pre-fault filter inputs. The init equalities make them the state anchors exactly:
    # `P_meas = P_g = Pe_0`, `Q_meas = Q_g = Qe_0`, `V_meas = V[bus]`.
    input_anchor = gen -> (vars[:P_meas][gen], vars[:Q_meas][gen], vars[:V_meas][gen])
    return _attach_gfm_window!(model, dyn_model_dict, gfm_gens, DGEN, DGFM,
        Δt, ω_syn, time_window, "tf", anchor, input_anchor)
end

"""
    Attach_GFM_postf!(model, dyn_model_dict, gfm_gens, DGEN, DGFM, Δt, ω_syn, time_window)

Post-fault GFM dynamics; the `t = 1` previous values come from the last fault step.
"""
function Attach_GFM_postf!(
    model::Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    gfm_gens::AbstractVector{Int},
    DGEN::DataFrame,
    DGFM::DataFrame,
    Δt::Float64,
    ω_syn::Float64,
    time_window::Vector{Float64},
)
    register_gfm_meta!(dyn_model_dict, DGFM)
    isempty(gfm_gens) && return dyn_model_dict
    vars = dyn_model_dict[:vars]
    anchor = gen -> (
        _gfm_last_var(vars[:P_meas_tf][gen]), _gfm_last_var(vars[:Q_meas_tf][gen]),
        _gfm_last_var(vars[:V_meas_tf][gen]), _gfm_last_var(vars[:E_int_tf][gen]),
        _gfm_last_var(vars[:δ_tf][gen]))
    # Last fault-window filter *inputs*; unlike `tf` these differ from the state anchors,
    # so the trapezoidal `t = 1` correction term is genuinely non-zero here.
    input_anchor = gen -> (
        _gfm_last_var(vars[:Pe_tf][gen]), _gfm_last_var(vars[:Qe_tf][gen]),
        _gfm_last_var(vars[:V_tf][Int(DGEN.bus[gen])]))
    return _attach_gfm_window!(model, dyn_model_dict, gfm_gens, DGEN, DGFM,
        Δt, ω_syn, time_window, "tpf", anchor, input_anchor)
end

"""
    validate_gfm_dispatch_band!(DGEN, DGFM, base_MVA)

Each converter's dispatch box must sit inside the device band the overload PIs
initialise against. `x(0) = 0` is an equilibrium only inside that band, so a
pre-fault point outside it is rejected whether or not `gfm_pq_overload` is on.

`DGEN` limits are MW (`pg_min`/`pg_max`/`qg_min`/`qg_max`). `DGFM` limits are
already system pu after `apply_gfm_base_conversion!`. `DGEN` owns the dispatch
range; this check does not add a second box on `Q_g`.
"""
function validate_gfm_dispatch_band!(
    DGEN::DataFrame,
    DGFM::DataFrame,
    base_MVA::Float64,
)
    base_MVA > 0.0 || throw(ArgumentError(
        "validate_gfm_dispatch_band!: base_MVA must be positive (got $base_MVA)."))
    atol = 1.0e-8
    for r in 1:nrow(DGFM)
        gen = Int(DGFM.id[r])
        i = findfirst(==(gen), Int.(DGEN.id))
        i === nothing && throw(ArgumentError(
            "GFM generator $gen is listed in DGFM but missing from DGEN."))
        bands = (
            ("P", DGEN.pg_min[i] / base_MVA, DGEN.pg_max[i] / base_MVA,
             Float64(DGFM.Pmin[r]), Float64(DGFM.Pmax[r])),
            ("Q", DGEN.qg_min[i] / base_MVA, DGEN.qg_max[i] / base_MVA,
             Float64(DGFM.Qmin[r]), Float64(DGFM.Qmax[r])),
        )
        for (axis, dlo, dhi, glo, ghi) in bands
            if dlo < glo - atol || dhi > ghi + atol
                throw(ArgumentError(
                    "Generator $gen dispatch $axis box [$(dlo), $(dhi)] pu sits " *
                    "outside the DGFM device band [$(glo), $(ghi)] pu. DGEN owns " *
                    "the dispatch range and must lie inside the converter rating " *
                    "the overload PIs initialise against."))
            end
        end
    end
    return nothing
end
