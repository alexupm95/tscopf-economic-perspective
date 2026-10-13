#=
================================================================================
 ts_bound_limits.jl — resolve TsBoundLimitsConfig and attach optional time bounds
================================================================================
=#

"""
Scalar lower/upper pair for optional transient variable boxes.

Defaults `(-Inf, Inf)` mean **inactive / unbounded** — not a numeric placeholder.
At attach time:
- `VARIABLE` encoding: JuMP omits that side (native unbounded).
- `CONSTRAINT` encoding: no ≤-row is created for that side (Inf never enters an
  explicit inequality RHS; solvers reject that). Toggle on + Inf = no-op until
  finite limits are set.
"""
Base.@kwdef struct TsBoundLimitPair
    min::Float64 = -Inf
    max::Float64 =  Inf
end

"""
    TsBoundLimitsConfig

User-visible transient bound limits (`TsBuilderConfig.limits`).

Physical defaults: `E` [0, 2] pu, `δ` [-π, π] rad, `P_m` from `DGEN` (`P_m_source`),
FullBus `V_bus_min_pu` (forced lower on `V_tf` / `V_tpf`), governor valve
`[gov_valve_min_pu, pg_max/base_MVA]` when `include_governor` (`gov_valve_max_source`).
Optional tf/tpf pairs default to `(-Inf, Inf)` (inactive until you set finite
limits **and** enable the matching `bound_*` toggle). Same Inf policy for both
`CONSTRAINT` and `VARIABLE` encodings — see `TsBoundLimitPair`.
"""
Base.@kwdef struct TsBoundLimitsConfig
    E_min_pu::Float64 = 0.0
    E_max_pu::Float64 = 2.0
    δ_min_rad::Float64 = -π
    δ_max_rad::Float64 = π
    P_m_source::Symbol = :dgen_pg_limits
    V_bus_min_pu::Float64 = 0.0
    gov_valve_min_pu::Float64 = 0.0
    gov_valve_max_source::Symbol = :dgen_pg_limits
  # Classical / DQ shared optional tf/tpf boxes
    δ_tf::TsBoundLimitPair = TsBoundLimitPair()
    Δω_tf::TsBoundLimitPair = TsBoundLimitPair()
    Pe_tf::TsBoundLimitPair = TsBoundLimitPair()
    Qe_tf::TsBoundLimitPair = TsBoundLimitPair()
    δCOI_tf::TsBoundLimitPair = TsBoundLimitPair()
    δ_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Δω_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Pe_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Qe_tpf::TsBoundLimitPair = TsBoundLimitPair()
    δCOI_tpf::TsBoundLimitPair = TsBoundLimitPair()
  # DQ pre-fault algebraic states
    Ed_min_pu::Float64 = -2.0
    Ed_max_pu::Float64 =  2.0
    Eq_min_pu::Float64 = -2.0
    Eq_max_pu::Float64 =  2.0
    Id_min_pu::Float64 = -5.0
    Id_max_pu::Float64 =  5.0
    Iq_min_pu::Float64 = -5.0
    Iq_max_pu::Float64 =  5.0
  # AVR / governor pre-fault set-points
    V_ref_min_pu::Float64 = 0.8
    V_ref_max_pu::Float64 = 1.2
    P_ref_source::Symbol = :dgen_pg_limits
  # DQ tf/tpf machine states
    Ed_tf::TsBoundLimitPair = TsBoundLimitPair()
    Eq_tf::TsBoundLimitPair = TsBoundLimitPair()
    Id_tf::TsBoundLimitPair = TsBoundLimitPair()
    Iq_tf::TsBoundLimitPair = TsBoundLimitPair()
    Te_tf::TsBoundLimitPair = TsBoundLimitPair()
    Ed_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Eq_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Id_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Iq_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Te_tpf::TsBoundLimitPair = TsBoundLimitPair()
  # AVR / governor tf/tpf (research boxes; physics saturation unchanged)
    E_fd_unlim_tf::TsBoundLimitPair = TsBoundLimitPair()
    E_fd_unlim_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Pv_tf::TsBoundLimitPair = TsBoundLimitPair()
    Pv_tpf::TsBoundLimitPair = TsBoundLimitPair()
    Pm_tf::TsBoundLimitPair = TsBoundLimitPair()
    Pm_tpf::TsBoundLimitPair = TsBoundLimitPair()
  # GFM boxes (allow_gfm; DQ FULL_BUS). Margins and floors, not final numbers: the box
  # is derived per unit from DGFM (Emin/Emax, Imax, Pmin/Pmax) by
  # `resolve_gfm_bound_limits`. Unlike the optional tf/tpf pairs above these are always
  # built — the nonconvex limiter needs the guard-rails to converge — so the defaults
  # below reproduce the values the reference port hard-coded.
    gfm_δ_min_rad::Float64          = -π
    gfm_δ_max_rad::Float64          =  π
    gfm_V_meas_min_pu::Float64      = 0.0
    gfm_V_meas_max_pu::Float64      = 2.5
    gfm_E_raw_extra_pu::Float64     = 0.25   # band around [Emin, Emax] for the raw PI states
    gfm_E_raw_slack_pu::Float64     = 0.5    # further slack so the box never pre-empts the clip
    gfm_E_clip_slack_pu::Float64    = 0.05   # slack outside the smooth-clip range
    gfm_PQ_bound_scale::Float64     = 1.5    # margin factor in the P/Q/I box formulas
    gfm_PQ_bound_offset_pu::Float64 = 0.05   # offset in the same formulas
    gfm_P_meas_floor_pu::Float64    = 3.0
    gfm_Q_meas_floor_pu::Float64    = 3.0
    gfm_I_floor_pu::Float64         = 5.0
    gfm_I_ceiling_pu::Float64       = 20.0
    # Overload-PI guard-rails. The raw state sits off the clip rail inside the band,
    # and a free raw column next to a saturated clip is the geometry that drives the
    # clip Jacobian soft. Symmetric raw box; the clipped state is held near 0, with
    # a small slack on the feasible side of the rail (upper clip ≤ 0, lower clip ≥ 0).
    gfm_pq_raw_span_pu::Float64     = 20.0
    gfm_pq_clip_slack_pu::Float64   = 0.05
end

const _TS_LIMIT_PAIR_FIELDS = (
    :δ_tf, :Δω_tf, :Pe_tf, :Qe_tf, :δCOI_tf,
    :δ_tpf, :Δω_tpf, :Pe_tpf, :Qe_tpf, :δCOI_tpf,
    :Ed_tf, :Eq_tf, :Id_tf, :Iq_tf, :Te_tf,
    :Ed_tpf, :Eq_tpf, :Id_tpf, :Iq_tpf, :Te_tpf,
    :E_fd_unlim_tf, :E_fd_unlim_tpf, :Pv_tf, :Pv_tpf, :Pm_tf, :Pm_tpf,
)

function copy_ts_bound_limits(l::TsBoundLimitsConfig)::TsBoundLimitsConfig
    return TsBoundLimitsConfig(; (f => getfield(l, f) for f in fieldnames(TsBoundLimitsConfig))...)
end

function _validate_limit_pair!(name::Symbol, pair::TsBoundLimitPair)
    pair.min < pair.max ||
        throw(ArgumentError("TsBoundLimitsConfig: require $(name).min < $(name).max."))
    return nothing
end

function validate_ts_bound_limits!(limits::TsBoundLimitsConfig)
    limits.E_min_pu < limits.E_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require E_min_pu < E_max_pu."))
    limits.δ_min_rad < limits.δ_max_rad ||
        throw(ArgumentError("TsBoundLimitsConfig: require δ_min_rad < δ_max_rad."))
    limits.P_m_source == :dgen_pg_limits ||
        throw(ArgumentError("TsBoundLimitsConfig: P_m_source must be :dgen_pg_limits."))
    limits.P_ref_source == :dgen_pg_limits ||
        throw(ArgumentError("TsBoundLimitsConfig: P_ref_source must be :dgen_pg_limits."))
    limits.gov_valve_max_source == :dgen_pg_limits ||
        throw(ArgumentError("TsBoundLimitsConfig: gov_valve_max_source must be :dgen_pg_limits."))
    limits.Ed_min_pu < limits.Ed_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require Ed_min_pu < Ed_max_pu."))
    limits.Eq_min_pu < limits.Eq_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require Eq_min_pu < Eq_max_pu."))
    limits.Id_min_pu < limits.Id_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require Id_min_pu < Id_max_pu."))
    limits.Iq_min_pu < limits.Iq_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require Iq_min_pu < Iq_max_pu."))
    limits.V_ref_min_pu < limits.V_ref_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require V_ref_min_pu < V_ref_max_pu."))
    limits.gfm_δ_min_rad < limits.gfm_δ_max_rad ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_δ_min_rad < gfm_δ_max_rad."))
    limits.gfm_V_meas_min_pu < limits.gfm_V_meas_max_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_V_meas_min_pu < gfm_V_meas_max_pu."))
    limits.gfm_I_floor_pu <= limits.gfm_I_ceiling_pu ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_I_floor_pu ≤ gfm_I_ceiling_pu."))
    for name in (:gfm_E_raw_extra_pu, :gfm_E_raw_slack_pu, :gfm_E_clip_slack_pu,
                 :gfm_PQ_bound_offset_pu, :gfm_P_meas_floor_pu, :gfm_Q_meas_floor_pu,
                 :gfm_I_floor_pu)
        getfield(limits, name) >= 0.0 ||
            throw(ArgumentError("TsBoundLimitsConfig: require $(name) ≥ 0."))
    end
    limits.gfm_PQ_bound_scale > 0.0 ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_PQ_bound_scale > 0."))
    limits.gfm_pq_raw_span_pu > 0.0 ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_pq_raw_span_pu > 0."))
    limits.gfm_pq_clip_slack_pu >= 0.0 ||
        throw(ArgumentError("TsBoundLimitsConfig: require gfm_pq_clip_slack_pu ≥ 0."))
    for name in _TS_LIMIT_PAIR_FIELDS
        _validate_limit_pair!(name, getfield(limits, name))
    end
    return nothing
end

"""
True when a resolved limit side should be materialised.

Scalars: `isfinite`. Vectors: any finite element (per-gen loops still skip
non-finite entries via `_limit_for_gen`). `nothing` / all-Inf → inactive —
never write `±Inf` into CONSTRAINT rows or JuMP VARIABLE bounds.
"""
function _use_lower_limit(lim)::Bool
    lim === nothing && return false
    lim isa AbstractVector && return any(isfinite, lim)
    return isfinite(lim)
end

"""True when a resolved limit side should be materialised (see `_use_lower_limit`)."""
function _use_upper_limit(lim)::Bool
    lim === nothing && return false
    lim isa AbstractVector && return any(isfinite, lim)
    return isfinite(lim)
end

"""
    resolve_ts_bound_limits(limits, active_gen, DGEN, base_MVA) -> OrderedDict

Build per-family `(min, max)` specs for builders. Vectors are per-generator where needed;
scalars are broadcast inside attach helpers.
"""
function resolve_ts_bound_limits(
    limits::TsBoundLimitsConfig,
    active_gen::Vector{Int},
    DGEN::DataFrame,
    base_MVA::Float64,
)::OrderedDict{Symbol, Tuple{Any, Any}}
    n = length(active_gen)
    specs = OrderedDict{Symbol, Tuple{Any, Any}}()
    specs[:E] = (fill(limits.E_min_pu, n), fill(limits.E_max_pu, n))
    specs[:δ] = (limits.δ_min_rad, limits.δ_max_rad)
    if limits.P_m_source == :dgen_pg_limits
        specs[:P_m] = (
            [DGEN.pg_min[gen] / base_MVA for gen in active_gen],
            [DGEN.pg_max[gen] / base_MVA for gen in active_gen],
        )
    end
    if limits.P_ref_source == :dgen_pg_limits
        specs[:P_ref] = (
            [DGEN.pg_min[gen] / base_MVA for gen in active_gen],
            [DGEN.pg_max[gen] / base_MVA for gen in active_gen],
        )
    end
    if limits.gov_valve_max_source == :dgen_pg_limits
        specs[:gov_valve] = (
            fill(limits.gov_valve_min_pu, n),
            [DGEN.pg_max[gen] / base_MVA for gen in active_gen],
        )
    end
    specs[:V_bus_min] = (limits.V_bus_min_pu, nothing)
    specs[:Ed] = (limits.Ed_min_pu, limits.Ed_max_pu)
    specs[:Eq] = (limits.Eq_min_pu, limits.Eq_max_pu)
    specs[:Id] = (limits.Id_min_pu, limits.Id_max_pu)
    specs[:Iq] = (limits.Iq_min_pu, limits.Iq_max_pu)
    specs[:V_ref] = (limits.V_ref_min_pu, limits.V_ref_max_pu)
    for key in _TS_LIMIT_PAIR_FIELDS
        pair = getfield(limits, key)
        specs[key] = (pair.min, pair.max)
    end
    return specs
end

function _limit_for_gen(min_lim, max_lim, gen_idx::Int)
    min_g = min_lim isa AbstractVector ? min_lim[gen_idx] : min_lim
    max_g = max_lim isa AbstractVector ? max_lim[gen_idx] : max_lim
    return min_g, max_g
end

"""True when `vars` holds JuMP variables (not algebraic expressions)."""
function _time_series_is_variables(vars)::Bool
    isempty(vars) && return false
    first_inner = first(values(vars))
    if first_inner isa JuMP.VariableRef
        return true
    elseif first_inner isa OrderedDict && !isempty(first_inner)
        return first(values(first_inner)) isa JuMP.VariableRef
    end
    return false
end

"""Apply JuMP variable bounds on an existing per-gen×time series."""
function apply_var_bounds_kron_gen_time!(
    vars,
    min_lim,
    max_lim,
)
    for (idx, gen) in enumerate(keys(vars))
        min_g, max_g = _limit_for_gen(min_lim, max_lim, idx)
        for t in keys(vars[gen])
            _use_lower_limit(min_g) && JuMP.set_lower_bound(vars[gen][t], Float64(min_g))
            _use_upper_limit(max_g) && JuMP.set_upper_bound(vars[gen][t], Float64(max_g))
        end
    end
    return nothing
end

"""Apply JuMP variable bounds on an existing COI time series."""
function apply_var_bounds_kron_coi_time!(
    vars::OrderedDict{Int, JuMP.VariableRef},
    min_lim::Float64,
    max_lim::Float64,
)
    for t in keys(vars)
        _use_lower_limit(min_lim) && JuMP.set_lower_bound(vars[t], min_lim)
        _use_upper_limit(max_lim) && JuMP.set_upper_bound(vars[t], max_lim)
    end
    return nothing
end

"""Explicit ≤-form bounds on per-generator time-indexed variables or expressions."""
function ineq_const_kron_gen_time_bounds!(
    model::JuMP.Model,
    vars,
    min_lim,
    max_lim,
)
    lower = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    upper = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()
    for (idx, gen) in enumerate(keys(vars))
        min_g, max_g = _limit_for_gen(min_lim, max_lim, idx)
        for t in keys(vars[gen])
            if _use_lower_limit(min_g)
                haskey(lower, gen) || (lower[gen] = OrderedDict{Int, JuMP.ConstraintRef}())
                lower[gen][t] = JuMP.@constraint(model, min_g - vars[gen][t] ≤ 0.0)
            end
            if _use_upper_limit(max_g)
                haskey(upper, gen) || (upper[gen] = OrderedDict{Int, JuMP.ConstraintRef}())
                upper[gen][t] = JuMP.@constraint(model, vars[gen][t] - max_g ≤ 0.0)
            end
        end
    end
    isempty(lower) && (lower = nothing)
    isempty(upper) && (upper = nothing)
    return lower, upper
end

"""Explicit ≤-form bounds on per-generator scalars."""
function ineq_const_kron_gen_scalar_bounds!(
    model::JuMP.Model,
    vars::OrderedDict{Int, JuMP.VariableRef},
    min_lim,
    max_lim,
)
    lower = OrderedDict{Int, JuMP.ConstraintRef}()
    upper = OrderedDict{Int, JuMP.ConstraintRef}()
    for (idx, gen) in enumerate(keys(vars))
        min_g, max_g = _limit_for_gen(min_lim, max_lim, idx)
        _use_lower_limit(min_g) && (lower[gen] = JuMP.@constraint(model, min_g - vars[gen] ≤ 0.0))
        _use_upper_limit(max_g) && (upper[gen] = JuMP.@constraint(model, vars[gen] - max_g ≤ 0.0))
    end
    isempty(lower) && (lower = nothing)
    isempty(upper) && (upper = nothing)
    return lower, upper
end

"""Explicit ≤-form bounds on a single COI time series (`t => variable`)."""
function ineq_const_kron_coi_time_bounds!(
    model::JuMP.Model,
    vars::OrderedDict{Int, JuMP.VariableRef},
    min_lim::Float64,
    max_lim::Float64,
)
    lower = OrderedDict{Int, JuMP.ConstraintRef}()
    upper = OrderedDict{Int, JuMP.ConstraintRef}()
    for t in keys(vars)
        _use_lower_limit(min_lim) && (lower[t] = JuMP.@constraint(model, min_lim - vars[t] ≤ 0.0))
        _use_upper_limit(max_lim) && (upper[t] = JuMP.@constraint(model, vars[t] - max_lim ≤ 0.0))
    end
    isempty(lower) && (lower = nothing)
    isempty(upper) && (upper = nothing)
    return lower, upper
end

"""Attach optional explicit bounds when `meta[:var_bounds][bound_key]` is true."""
function attach_kron_gen_time_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    vars,
    bound_key::Symbol,
    limit_key::Symbol,
    ineq_key_lower::Symbol,
    ineq_key_upper::Symbol,
)
    get(dyn_model_dict[:meta][:var_bounds], bound_key, false) || return nothing
    specs = dyn_model_dict[:meta][:var_limit_specs]
    min_lim, max_lim = specs[limit_key]
    !_use_lower_limit(min_lim) && !_use_upper_limit(max_lim) && return nothing
    meta = dyn_model_dict[:meta]
    encoding = bound_encoding_from_meta(meta)

    if encoding == VARIABLE && _time_series_is_variables(vars)
        apply_var_bounds_kron_gen_time!(vars, min_lim, max_lim)
        _use_lower_limit(min_lim) &&
            register_bound_manifest!(meta, ineq_key_lower, :lower, vars, :per_gen_time)
        _use_upper_limit(max_lim) &&
            register_bound_manifest!(meta, ineq_key_upper, :upper, vars, :per_gen_time)
        return nothing
    end

    if encoding == VARIABLE && !_time_series_is_variables(vars)
        @info "Bound encoding VARIABLE requested for $bound_key but target is an expression; using CONSTRAINT bounds."
    end

    lower, upper = ineq_const_kron_gen_time_bounds!(model, vars, min_lim, max_lim)
    lower !== nothing && (dyn_model_dict[:ineq_const][ineq_key_lower] = lower)
    upper !== nothing && (dyn_model_dict[:ineq_const][ineq_key_upper] = upper)
    return nothing
end

"""Attach optional explicit bounds on pre-fault per-generator scalars."""
function attach_kron_gen_scalar_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    vars::OrderedDict{Int, JuMP.VariableRef},
    bound_key::Symbol,
    limit_key::Symbol,
    ineq_key_lower::Symbol,
    ineq_key_upper::Symbol,
)
    get(dyn_model_dict[:meta][:var_bounds], bound_key, false) || return nothing
    specs = dyn_model_dict[:meta][:var_limit_specs]
    min_lim, max_lim = specs[limit_key]
    !_use_lower_limit(min_lim) && !_use_upper_limit(max_lim) && return nothing
    meta = dyn_model_dict[:meta]
    encoding = bound_encoding_from_meta(meta)

    if encoding == VARIABLE
        for (idx, gen) in enumerate(keys(vars))
            min_g, max_g = _limit_for_gen(min_lim, max_lim, idx)
            _use_lower_limit(min_g) && JuMP.set_lower_bound(vars[gen], Float64(min_g))
            _use_upper_limit(max_g) && JuMP.set_upper_bound(vars[gen], Float64(max_g))
        end
        _use_lower_limit(min_lim) &&
            register_bound_manifest!(meta, ineq_key_lower, :lower, vars, :gen_indexed)
        _use_upper_limit(max_lim) &&
            register_bound_manifest!(meta, ineq_key_upper, :upper, vars, :gen_indexed)
        return nothing
    end

    lower, upper = ineq_const_kron_gen_scalar_bounds!(model, vars, min_lim, max_lim)
    lower !== nothing && (dyn_model_dict[:ineq_const][ineq_key_lower] = lower)
    upper !== nothing && (dyn_model_dict[:ineq_const][ineq_key_upper] = upper)
    return nothing
end

"""Attach optional explicit bounds on a COI scalar time series."""
function attach_kron_coi_time_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    vars::OrderedDict{Int, JuMP.VariableRef},
    bound_key::Symbol,
    limit_key::Symbol,
    ineq_key_lower::Symbol,
    ineq_key_upper::Symbol,
)
    get(dyn_model_dict[:meta][:var_bounds], bound_key, false) || return nothing
    specs = dyn_model_dict[:meta][:var_limit_specs]
    min_lim, max_lim = specs[limit_key]
    !_use_lower_limit(min_lim) && !_use_upper_limit(max_lim) && return nothing
    meta = dyn_model_dict[:meta]
    encoding = bound_encoding_from_meta(meta)

    if encoding == VARIABLE
        apply_var_bounds_kron_coi_time!(vars, min_lim, max_lim)
        _use_lower_limit(min_lim) &&
            register_bound_manifest!(meta, ineq_key_lower, :lower, vars, :time_indexed)
        _use_upper_limit(max_lim) &&
            register_bound_manifest!(meta, ineq_key_upper, :upper, vars, :time_indexed)
        return nothing
    end

    lower, upper = ineq_const_kron_coi_time_bounds!(model, vars, min_lim, max_lim)
    lower !== nothing && (dyn_model_dict[:ineq_const][ineq_key_lower] = lower)
    upper !== nothing && (dyn_model_dict[:ineq_const][ineq_key_upper] = upper)
    return nothing
end

"""Register resolved limit specs on `ts_input_param` and `dyn_model_dict[:meta]`."""
function register_ts_bound_limit_specs!(
    ts_input_param::OrderedDict{Symbol, Any},
    dyn_model_dict::OrderedDict{Symbol, Any},
    limits::TsBoundLimitsConfig,
    active_gen::Vector{Int},
    DGEN::DataFrame,
    base_MVA::Float64;
    DGFM::Union{Nothing, DataFrame}=nothing,
    gfm_gens::AbstractVector{Int}=Int[],
)
    specs = resolve_ts_bound_limits(limits, active_gen, DGEN, base_MVA)
    # GFM boxes are per-unit and derived from DGFM; they share the spec dict so the
    # transient windows read one source (see `resolve_gfm_bound_limits`).
    if DGFM !== nothing && !isempty(gfm_gens)
        merge!(specs, resolve_gfm_bound_limits(limits, DGFM, gfm_gens))
    end
    ts_input_param[:var_limit_specs] = specs
    dyn_model_dict[:meta][:var_limit_specs] = specs
    return specs
end

"""
SG-only view of a shared per-generator container.

The SG limit families are calibrated for machine states and their spec vectors are
indexed by the SG generator order, while `vars[:Id_tf]` and friends hold every active
unit. GFM units own their own boxes (`attach_gfm_time_box!`), and `Δω` on a GFM row is a
droop output rather than a swing state, so applying an SG box there would be wrong twice
over.
"""
function _sg_only_view(dyn_model_dict::OrderedDict{Symbol, Any}, vars)
    gfm = get(dyn_model_dict[:meta], :gfm_gens, Int[])
    isempty(gfm) && return vars
    gfm_set = Set{Int}(gfm)
    out = OrderedDict{Int, Any}()
    for (gen, inner) in vars
        gen in gfm_set || (out[gen] = inner)
    end
    return out
end

function _qe_time_series(dyn_model_dict::OrderedDict{Symbol, Any}, key::Symbol)
    vars = dyn_model_dict[:vars]
    haskey(vars, key) && return vars[key]
    exprs = get(dyn_model_dict, :expressions, nothing)
    return exprs !== nothing && haskey(exprs, key) ? exprs[key] : nothing
end

function _attach_if_present(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    var_key::Symbol,
    bound_key::Symbol,
    limit_key::Symbol,
    ineq_lo::Symbol,
    ineq_hi::Symbol,
)
    vars = dyn_model_dict[:vars]
    haskey(vars, var_key) || return nothing
    attach_kron_gen_time_var_bounds!(model, dyn_model_dict,
        _sg_only_view(dyn_model_dict, vars[var_key]),
        bound_key, limit_key, ineq_lo, ineq_hi)
    return nothing
end

"""Wire optional explicit bounds for fault-on time-indexed variables."""
function attach_fault_tf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :δ_tf, :δ_tf, :δ_tf,
        :ineq_const_δ_tf_lower, :ineq_const_δ_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Δω_tf, :Δω_tf, :Δω_tf,
        :ineq_const_Δω_tf_lower, :ineq_const_Δω_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Pe_tf, :Pe_tf, :Pe_tf,
        :ineq_const_Pe_tf_lower, :ineq_const_Pe_tf_upper)
    qe_tf = _qe_time_series(dyn_model_dict, :Qe_tf)
    qe_tf !== nothing && attach_kron_gen_time_var_bounds!(
        model, dyn_model_dict, _sg_only_view(dyn_model_dict, qe_tf), :Qe_tf, :Qe_tf,
        :ineq_const_Qe_tf_lower, :ineq_const_Qe_tf_upper)
    haskey(dyn_model_dict[:vars], :δCOI_tf) && attach_kron_coi_time_var_bounds!(
        model, dyn_model_dict, dyn_model_dict[:vars][:δCOI_tf], :δCOI_tf, :δCOI_tf,
        :ineq_const_δCOI_tf_lower, :ineq_const_δCOI_tf_upper)
    attach_dq_fault_tf_var_bounds!(model, dyn_model_dict)
    attach_avr_tf_var_bounds!(model, dyn_model_dict)
    attach_gov_tf_var_bounds!(model, dyn_model_dict)
    return nothing
end

"""Wire optional explicit bounds for post-fault time-indexed variables."""
function attach_postfault_tpf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :δ_tpf, :δ_tpf, :δ_tpf,
        :ineq_const_δ_tpf_lower, :ineq_const_δ_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Δω_tpf, :Δω_tpf, :Δω_tpf,
        :ineq_const_Δω_tpf_lower, :ineq_const_Δω_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Pe_tpf, :Pe_tpf, :Pe_tpf,
        :ineq_const_Pe_tpf_lower, :ineq_const_Pe_tpf_upper)
    qe_tpf = _qe_time_series(dyn_model_dict, :Qe_tpf)
    qe_tpf !== nothing && attach_kron_gen_time_var_bounds!(
        model, dyn_model_dict, _sg_only_view(dyn_model_dict, qe_tpf), :Qe_tpf, :Qe_tpf,
        :ineq_const_Qe_tpf_lower, :ineq_const_Qe_tpf_upper)
    haskey(dyn_model_dict[:vars], :δCOI_tpf) && attach_kron_coi_time_var_bounds!(
        model, dyn_model_dict, dyn_model_dict[:vars][:δCOI_tpf], :δCOI_tpf, :δCOI_tpf,
        :ineq_const_δCOI_tpf_lower, :ineq_const_δCOI_tpf_upper)
    attach_dq_postfault_tpf_var_bounds!(model, dyn_model_dict)
    attach_avr_tpf_var_bounds!(model, dyn_model_dict)
    attach_gov_tpf_var_bounds!(model, dyn_model_dict)
    return nothing
end

"""Optional pre-fault boxes on DQ algebraic states (`Ed`, `Eq`, `Id`, `Iq`)."""
function attach_dq_prefault_algebraic_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    vars = dyn_model_dict[:vars]
    haskey(vars, :Ed) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:Ed], :Ed, :Ed,
        :ineq_const_Ed_lower, :ineq_const_Ed_upper)
    haskey(vars, :Eq) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:Eq], :Eq, :Eq,
        :ineq_const_Eq_lower, :ineq_const_Eq_upper)
    haskey(vars, :Id) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:Id], :Id, :Id,
        :ineq_const_Id_lower, :ineq_const_Id_upper)
    haskey(vars, :Iq) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:Iq], :Iq, :Iq,
        :ineq_const_Iq_lower, :ineq_const_Iq_upper)
    return nothing
end

"""Optional AVR pre-fault box on `V_ref`."""
function attach_avr_prefault_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    vars = dyn_model_dict[:vars]
    haskey(vars, :V_ref) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:V_ref], :V_ref, :V_ref,
        :ineq_const_V_ref_lower, :ineq_const_V_ref_upper)
    return nothing
end

"""Optional governor pre-fault box on `P_ref`."""
function attach_gov_prefault_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    vars = dyn_model_dict[:vars]
    haskey(vars, :P_ref) && attach_kron_gen_scalar_var_bounds!(
        model, dyn_model_dict, vars[:P_ref], :P_ref, :P_ref,
        :ineq_const_P_ref_lower, :ineq_const_P_ref_upper)
    return nothing
end

"""DQ-only optional fault-window boxes."""
function attach_dq_fault_tf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :Ed_tf, :Ed_tf, :Ed_tf,
        :ineq_const_Ed_tf_lower, :ineq_const_Ed_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Eq_tf, :Eq_tf, :Eq_tf,
        :ineq_const_Eq_tf_lower, :ineq_const_Eq_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Id_tf, :Id_tf, :Id_tf,
        :ineq_const_Id_tf_lower, :ineq_const_Id_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Iq_tf, :Iq_tf, :Iq_tf,
        :ineq_const_Iq_tf_lower, :ineq_const_Iq_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Te_tf, :Te_tf, :Te_tf,
        :ineq_const_Te_tf_lower, :ineq_const_Te_tf_upper)
    return nothing
end

"""DQ-only optional post-fault boxes."""
function attach_dq_postfault_tpf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :Ed_tpf, :Ed_tpf, :Ed_tpf,
        :ineq_const_Ed_tpf_lower, :ineq_const_Ed_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Eq_tpf, :Eq_tpf, :Eq_tpf,
        :ineq_const_Eq_tpf_lower, :ineq_const_Eq_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Id_tpf, :Id_tpf, :Id_tpf,
        :ineq_const_Id_tpf_lower, :ineq_const_Id_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Iq_tpf, :Iq_tpf, :Iq_tpf,
        :ineq_const_Iq_tpf_lower, :ineq_const_Iq_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Te_tpf, :Te_tpf, :Te_tpf,
        :ineq_const_Te_tpf_lower, :ineq_const_Te_tpf_upper)
    return nothing
end

"""AVR optional research boxes on `E_fd_unlim_*` (physics saturation unchanged)."""
function attach_avr_tf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :E_fd_unlim_tf, :E_fd_unlim_tf, :E_fd_unlim_tf,
        :ineq_const_E_fd_unlim_tf_lower, :ineq_const_E_fd_unlim_tf_upper)
    return nothing
end

function attach_avr_tpf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :E_fd_unlim_tpf, :E_fd_unlim_tpf, :E_fd_unlim_tpf,
        :ineq_const_E_fd_unlim_tpf_lower, :ineq_const_E_fd_unlim_tpf_upper)
    return nothing
end

"""Governor optional research boxes on valve/mechanical trajectories (`GOV_NO_LIMIT` friendly)."""
function attach_gov_tf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :Pv_tf, :Pv_tf, :Pv_tf,
        :ineq_const_Pv_tf_lower, :ineq_const_Pv_tf_upper)
    _attach_if_present(model, dyn_model_dict, :Pm_tf, :Pm_tf, :Pm_tf,
        :ineq_const_Pm_tf_lower, :ineq_const_Pm_tf_upper)
    return nothing
end

function attach_gov_tpf_var_bounds!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    _attach_if_present(model, dyn_model_dict, :Pv_tpf, :Pv_tpf, :Pv_tpf,
        :ineq_const_Pv_tpf_lower, :ineq_const_Pv_tpf_upper)
    _attach_if_present(model, dyn_model_dict, :Pm_tpf, :Pm_tpf, :Pm_tpf,
        :ineq_const_Pm_tpf_lower, :ineq_const_Pm_tpf_upper)
    return nothing
end

"""Store pre-fault scalar bounds from a var builder NamedTuple into `dyn_model_dict`."""
function store_ts_scalar_var_bounds!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    var_key::Symbol,
    result::NamedTuple,
    lower_key::Symbol,
    upper_key::Symbol,
)
    dyn_model_dict[:vars][var_key] = result.vars
    if result.lower !== nothing
        dyn_model_dict[:ineq_const][lower_key] = result.lower
    end
    if result.upper !== nothing
        dyn_model_dict[:ineq_const][upper_key] = result.upper
    end
    return result.vars
end
