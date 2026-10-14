#=
================================================================================
 DynModelConfig.jl  —  runtime selection of dynamic-model complexity
================================================================================
 Phase 0 (migration scaffolding). Defaults reproduce the current integrated
 behaviour: 2nd-order classical generator, Kron-reduced network, P_g in swing.

 Phase 1.1 wires `mech_power_mode` into Kron builders; Phase 1.2 routes builders
 via `dynamic_gen_model` in `AbstractDynamicGenModel.jl`.  Validation runs in
 `run_case!` so invalid combinations fail fast.
================================================================================
=#

# --- enumerations -------------------------------------------------------------

"""
    GenOrder

Generator swing-equation order. `CLASSICAL_2ND` (rotor angle + speed
deviation) is implemented on Kron and FULL_BUS paths. `DQ_4TH` (dq-axis,
4th order, FULL_BUS only) implements the machine core with optional AVR and
optional TGOV1 governor on the DQ path.
"""
@enum GenOrder       CLASSICAL_2ND DQ_4TH

"""
    NetworkForm

Transient network representation. `KRON_REDUCED` eliminates non-generator
buses via Kron reduction; `FULL_BUS` keeps the full sparse nodal `Ybus` with
explicit bus voltage/angle states and requires a mandatory ACOPF warm start.
"""
@enum NetworkForm    KRON_REDUCED FULL_BUS

"""
    MechPowerMode

Mechanical power used in the swing equation. `USE_PG` (default) uses the
dispatch decision `P_g` directly; `USE_PM` introduces an explicit `P_m`
variable and requires `bound_style_δ = :coi_box`. `FULL_BUS` requires `USE_PM`.
"""
@enum MechPowerMode  USE_PG USE_PM

"""
    GovernorLimiter

Valve/gate saturation treatment for the optional turbine governor
(`include_governor`). `GOV_NO_LIMIT` leaves the valve state unbounded (linear
governor, cleanest duals); `GOV_SMOOTH` applies the reference sqrt-based smooth
min/max anti-windup clamp (physically saturating, nonconvex); `GOV_HARD_BOUND`
adds explicit ≤-form bounds on the valve state (clean bound duals but can render
a strong transient infeasible). Only consulted when `include_governor = true`.
"""
@enum GovernorLimiter  GOV_NO_LIMIT GOV_SMOOTH GOV_HARD_BOUND

# --- configuration struct -----------------------------------------------------

"""
    DynModelConfig

Physics, network form, and disturbance selection for a transient-stability
run (`RunConfig.transient.dyn_model`).

Core knobs: `gen_order` (machine model), `network_form` (Kron-reduced vs
full-bus nodal), `mech_power_mode` (`P_g` vs explicit `P_m` in the swing
equation), `zip_load_p` / `zip_load_q` (independent constant-Z/I/P demand
splits for `FULL_BUS`, each `(Z, I, P)` and each summing to 1), the two
stability corridors described below, and `fault::FaultConfig`.

The two ZIP vectors are independent so the active and reactive demand can follow
different TSO conventions — e.g. REE's constant-current active / constant-admittance
reactive is `zip_load_p = (0.0, 1.0, 0.0)` with `zip_load_q = (1.0, 0.0, 0.0)`.

# Stability corridors

The rotor-angle and speed corridors are **independent**: each has its own on/off
toggle and its own reference, and a transient run needs at least one of them.

Rotor angle — `constrain_δ` (default `true`) with `bound_style_δ`:

| `bound_style_δ` | Corridor |
|---|---|
| `:swing_propagated` | `δ_g − δ_COI` band with one swing step substituted for `δ_g[t]` |
| `:coi_box` | plain box on `δ_g − δ_COI` |
| `:highest_H` | plain box on `δ_g − δ_ref`, `ref` = largest-`H` surviving synchronous machine |
| `:ref_gen` | same, `ref` = `δ_ref_gen_id` (required for this style) |
| `:pairwise` | box on `δ_i − δ_j` for every unordered pair of synchronised units |

The half-widths come from `TsSimulationConfig.δ_tol_deg*`, except `:pairwise`, which
reads its own symmetric `TsSimulationConfig.δ_pair_tol_deg` (default 180°). The
machine-referenced styles follow the common convention of bounding the n−1 remaining
machines against one live reference unit; the reference itself carries no row.
`:pairwise` has no reference: it states the loss-of-synchronism criterion directly, so
its duals are invariant to the choice of angle frame, at the cost of n(n−1)/2 rows per
time step instead of n−1.

`:swing_propagated` is inertia-weighted and spans synchronous machines only;
`:highest_H`, `:ref_gen` and `:pairwise` span grid-forming converters too. `:coi_box`
spans machines only by default and converters as well under `gfm_in_coi = true`, which
also weights them into the COI reference itself.

Speed — `constrain_Δω` (default `false`) with `bound_style_Δω`:

| `bound_style_Δω` | Corridor |
|---|---|
| `:coi_box` | box on `Δω_g − Δω_COI` |
| `:abs` | box on the raw `Δω_g` of the swing equation; no COI is formed at all |

`Δω_tol_pu` sets a symmetric half-width, `Δω_tol_pu_lower` / `Δω_tol_pu_upper` set
independent below/above limits (positive magnitudes, p.u.). Both styles read the same
three fields.

Validation rules enforced in `validate_dyn_config!`: `USE_PM` and `DQ_4TH` reject
`:swing_propagated` (that row substitutes the dispatch power into the angle band, which
is inconsistent with an explicit `P_m`); `FULL_BUS` requires `USE_PM` plus a solved
ACOPF warm start; TSC-DCOPF + `FULL_BUS` is not implemented. The data-dependent
reference-generator checks live in `validate_δ_reference!`, also pre-warm-start.

See the parameter reference, section 4.
"""
Base.@kwdef struct DynModelConfig
    # Machine / network complexity (Phase 1+ will dispatch on these)
    gen_order::GenOrder           = CLASSICAL_2ND
    network_form::NetworkForm     = KRON_REDUCED
    mech_power_mode::MechPowerMode = USE_PG

    # DQ_4TH only: include (1+Δω) in stator/Pe algebraic equations (RMS-style); false neglects it.
    dq_speed_dev_in_algebra::Bool = true

    # Optional control layers (dq full-network only, Phase 5–6)
    include_avr::Bool             = false
    include_governor::Bool        = false
    # Valve saturation treatment for the turbine governor (only read when include_governor).
    governor_limiter::GovernorLimiter = GOV_NO_LIMIT

    # Mixed SG + GFM fleet (FULL_BUS + DQ_4TH). GFM params live in a separate
    # `gfm_dynamic_data.csv` (see TransientConfig.gfm_dynamic_filename), not in
    # gen_dynamic_data*.csv. Phase G0 loads/validates; transient physics is G1+.
    allow_gfm::Bool               = false

    # Do grid-forming converters take part in the Centre of Inertia?
    #
    #   false (default) — COI weights and the `:coi_box` corridors span synchronous
    #                     machines only. A converter moves neither the reference nor
    #                     carries a corridor row: the historical behaviour.
    #   true            — converters join both, weighted by the equivalent inertia of
    #                     their frequency loop, H_eff = Tf / (2·mp), and bounded by the
    #                     same `δ_tol_deg` / `Δω_tol_pu` as the machines.
    #
    # H_eff follows from eliminating P_meas between the measurement filter and the droop
    # law; the pair collapses to a swing-shaped ODE with D_eff = 1/mp. It is a property
    # of the frequency loop alone, so it is unchanged when the current limiter is active
    # — the limiter only moves Pe, which sits on the right-hand side.
    #
    # Governs `bound_style_δ = :coi_box` and `bound_style_Δω = :coi_box` together, and
    # nothing else: `:swing_propagated` needs a per-machine H *and* D, and `:highest_H`
    # still ranks its reference by machine inertia. Inert unless `allow_gfm = true` and
    # a COI-referenced corridor is switched on (warned in `validate_dyn_config!`).
    #
    # Weights are typically lopsided — H_eff is 0.5 s on the shipped `9bus_gfm` against
    # machine inertias of 3–24 s — so expect the new converter corridor rows, not a large
    # shift in the reference.
    gfm_in_coi::Bool              = false

    # REGFM_A1 Fig. 3 P and Q overload limiters.
    #
    #   false (default) — the P–f droop and the Q–V voltage error are exactly the
    #                     loops shipped without these branches. No new states.
    #   true            — both one-sided PI pairs are stamped: Pmax/Pmin added to
    #                     per-unit Δω, Qmax/Qmin added to the voltage error. Eight
    #                     states per converter-step (raw and clipped, each branch).
    #
    # The radial current limiter is unchanged either way. `DGFM` ratings are the
    # comparator inputs; `DGEN` still owns the dispatch box. The integrators follow
    # `gfm_integrator` (same resolver as the Q–V PI) and advance from the clipped
    # previous state. Inert unless `allow_gfm = true` (warned in `validate_dyn_config!`).
    gfm_pq_overload::Bool         = false

    # Opt-in for a converter-only fleet (every active generator listed in
    # `gfm_dynamic_data.csv`). Default `false` keeps the historical hard error, because a
    # case that lists every unit in the GFM file would otherwise *silently* become an
    # all-converter run whose "COI" is a defined weighting rather than a physical centre
    # of inertia. Setting it true is a statement that this was intended.
    #
    # The supported combination is `gfm_in_coi = true` with `bound_style_δ` in
    # (`:coi_box`, `:ref_gen`, `:pairwise`). `:swing_propagated` and `:highest_H` stay
    # rejected: the first needs a per-machine H *and* D, the second ranks its reference
    # by machine inertia, and a converter fleet supplies neither. AVR and governor are
    # rejected too — both attach to `sg_gens` only, so they would be silent no-ops.
    # All of this is enforced by `validate_all_gfm_fleet!` before the warm-start solve.
    allow_all_gfm::Bool           = false

    # First ODE step of each fault/post-fault window for FULL_BUS swing (δ/Δω),
    # EMF, AVR, and governor. `:trapezoidal` (default) matches the pinned SG
    # AVR+TG FULL_BUS path; `:backward_euler` matches the reference implementation
    # BE-at-t=1 (optional for investigation). Remainder steps stay trapezoidal.
    ode_first_step::Symbol        = :trapezoidal

    # Time discretisation of the GFM measurement filters (P/Q/V) and the Q–V PI
    # integrator. Scoped to those two families only: GFM δ and every SG family keep
    # following `ode_first_step` above.
    #   :backward_euler        — BE at every step (DEFAULT). Matches the reference
    #                            implementation, and is the only scheme that stays
    #                            affordable as the horizon grows (see below).
    #   :follow_ode_first_step — BE at t=1 of each window iff `ode_first_step` is
    #                            `:backward_euler`, trapezoidal otherwise and always
    #                            for t ≥ 2 — i.e. the same rule the SG families use.
    #   :trapezoidal           — trapezoidal at every step.
    #
    # Why BE is the default despite trapezoidal being second-order accurate: the
    # trapezoidal row carries `u_{t-1}` (Pe/Qe, nonlinear in V, θ, δ, Id, Iq) into every
    # filter equation, and the resulting Jacobian/Hessian fill-in compounds with horizon
    # length. Measured on 9bus_gfm, SC at bus 7 cleared in 150 ms, t_step = 0.01 s:
    #
    #   horizon   :backward_euler        :follow_ode_first_step
    #   0.6 s     139 it                 97 it,  0.074 s/it
    #   3.0 s     205 it, 0.198 s/it     161 it, 17.6  s/it   (≈70× total solve time)
    #
    # Accuracy cost of the default is small: worst rotor-angle gap between the two was
    # 2.5° absolute / 1.8° relative to COI over 3 s (≈2 % of swing), objective 1.6e-5
    # relative. Use the trapezoidal variants for accuracy studies on short horizons.
    #
    # Stability note: BE is L-stable and never rings; trapezoidal oscillates
    # step-to-step once `t_step > 2·Tf` (warned in the builder).
    gfm_integrator::Symbol        = :backward_euler

    # ZIP load model for the FULL_BUS KCL path. Active and reactive demand carry
    # INDEPENDENT splits, because TSOs do not generally use the same one for both
    # (e.g. REE recommends constant-current P with constant-admittance Q).
    # Coefficient order is (Z, I, P) in BOTH vectors:
    #   [1] = constant impedance (∝ V_t²), [2] = constant current (∝ V_t·V_nom),
    #   [3] = constant power (∝ V_nom²).  Each vector must sum to 1 on its own.
    # Ignored on KRON_REDUCED, where loads are folded into Y_red as constant admittance.
    # WARNING — the reference implementation stores the vector as (P, I, Z) in the KCL
    # formula (`ZIP[1]*V² + ZIP[2]*V_t*V + ZIP[3]*V_t²`). Its `[0,0,1]` is therefore
    # constant-Z, which maps to package `(1,0,0)` — NOT package `(0,0,1)` (const-P).
    zip_load_p::NTuple{3, Float64} = (1.0, 0.0, 0.0)  # active-demand split (Z, I, P)
    zip_load_q::NTuple{3, Float64} = (1.0, 0.0, 0.0)  # reactive-demand split (Z, I, P)

    # --- rotor-angle corridor --------------------------------------------------
    # Half-widths come from TsSimulationConfig.δ_tol_deg*. USE_PM / DQ_4TH reject
    # :swing_propagated (enforced in validate_dyn_config!); :ref_gen requires an id.
    constrain_δ::Bool             = true
    bound_style_δ::Symbol         = :swing_propagated  # :coi_box | :highest_H | :ref_gen | :pairwise
    # Reference machine for bound_style_δ = :ref_gen (a DGEN_DYN generator id). Must be
    # in service, not tripped by a GL disturbance, and not a GFM unit — checked against
    # the actual system data in validate_δ_reference!, before the warm-start solve.
    δ_ref_gen_id::Union{Nothing, Int} = nothing

    # --- speed corridor ----------------------------------------------------------
    # Off by default. :coi_box bounds Δω_i − Δω_COI, :abs bounds the raw Δω of the swing
    # equation and forms no COI. `Δω_tol_pu` is the symmetric half-width; the two
    # overrides set independent below/above limits as positive magnitudes, mirroring
    # δ_tol_deg / δ_tol_deg_lower / δ_tol_deg_upper.
    # e.g. lower=0.05, upper=0.03 → the signed pair (-0.05, +0.03).
    constrain_Δω::Bool            = false
    bound_style_Δω::Symbol        = :coi_box           # :abs
    Δω_tol_pu::Float64            = 0.5
    Δω_tol_pu_lower::Union{Nothing, Float64} = nothing  # below reference [p.u.]; default → Δω_tol_pu
    Δω_tol_pu_upper::Union{Nothing, Float64} = nothing  # above reference [p.u.]; default → Δω_tol_pu

    # Disturbance specification (SC via contingencies.csv; GL / OB via explicit ids)
    fault::FaultConfig            = FaultConfig()
end

# --- FULL_BUS coupling hints (from the mandatory ACOPF warm start) ------------

"""
    SteadyStateHints

Operating-point snapshot used to set `start=` on FULL_BUS dynamic variables.
Populated by `extract_opf_solved_hints` after the mandatory ACOPF pre-solve.
"""
struct SteadyStateHints
    val_V::Dict{Int, Float64}
    val_θ::Dict{Int, Float64}
    val_Pg::Dict{Int, Float64}
    val_Qg::Dict{Int, Float64}
end

"""Extract `JuMP.value` on steady-state V, θ, P_g, Q_g after the warm-start solve."""
function extract_opf_solved_hints(opf_dict::OrderedDict{Symbol, Any})::SteadyStateHints
    return SteadyStateHints(
        Dict(i => Float64(JuMP.value(v)) for (i, v) in opf_dict[:vars][:V]),
        Dict(i => Float64(JuMP.value(v)) for (i, v) in opf_dict[:vars][:θ]),
        Dict(i => Float64(JuMP.value(v)) for (i, v) in opf_dict[:vars][:P_g]),
        Dict(i => Float64(JuMP.value(v)) for (i, v) in opf_dict[:vars][:Q_g]),
    )
end

# --- validation ---------------------------------------------------------------

"""
Columns required by `dyn` that `DGEN_DYN` does not usably provide — either absent
altogether or present with a non-finite entry on some generator row.

The two failures are one bucket on purpose: `Parse_Gen_Dynamic_DataFrame` pre-fills
every optional column with `NaN` and overwrites only the ones the CSV actually
carries, so a missing column is indistinguishable from a `NaN`-valued one by the
time anything downstream reads it — and both are equally fatal.
"""
function missing_dyn_columns(dyn::DynModelConfig, DGEN_DYN::DataFrame)::Vector{Symbol}
    return Symbol[col for col in required_dyn_column_names(dyn)
                  if !(col in propertynames(DGEN_DYN)) || !all(isfinite, DGEN_DYN[!, col])]
end

"""Return `true` when every optional machine column needed by `dyn` is finite."""
function has_required_dyn_columns(dyn::DynModelConfig, DGEN_DYN::DataFrame)::Bool
    return isempty(missing_dyn_columns(dyn, DGEN_DYN))
end

"""Canonical DataFrame column names required for the chosen dynamic configuration."""
function required_dyn_column_names(dyn::DynModelConfig)::Vector{Symbol}
    # E_int is never read from CSV — internal voltage E is a JuMP variable from init.
    cols = Symbol[:bus, :Xd_tr, :H, :D]
    if dyn.gen_order == DQ_4TH
        append!(cols, [:Xq_tr, :Xd, :Xq, :Td, :Tq, :Ra])
    end
    if dyn.include_avr
        # Ta_exc / Tb_exc are the SEXS lead-lag stage; `0` for both bypasses it, but the
        # columns still have to be there so a blank cell cannot silently mean "bypass".
        append!(cols, [:T_exc, :K_exc, :Ta_exc, :Tb_exc])
    end
    if dyn.include_governor
        append!(cols, [:R, :T1, :T2, :T3])
    end
    return unique(cols)
end

# Why each machine column is required, read in the opposite direction from
# `required_dyn_column_names` above. Kept adjacent so the two never drift.
const _DYN_COLUMN_REASON = Dict{Symbol, String}(
    :bus => "the classical 2nd-order machine", :Xd_tr => "the classical 2nd-order machine",
    :H => "the classical 2nd-order machine",   :D => "the classical 2nd-order machine",
    :Xq_tr => "gen_order=DQ_4TH", :Xd => "gen_order=DQ_4TH", :Xq => "gen_order=DQ_4TH",
    :Td => "gen_order=DQ_4TH",    :Tq => "gen_order=DQ_4TH", :Ra => "gen_order=DQ_4TH",
    :T_exc => "include_avr=true", :K_exc => "include_avr=true",
    :Ta_exc => "include_avr=true", :Tb_exc => "include_avr=true",
    :R => "include_governor=true", :T1 => "include_governor=true",
    :T2 => "include_governor=true", :T3 => "include_governor=true",
)

"""
Distinct reasons `cols` are required, in first-appearance order — e.g.
`["include_governor=true"]`, or `["gen_order=DQ_4TH", "include_avr=true"]` when a
run is short of both families. Used to explain *which knob* demands a column the
machine CSV does not carry.
"""
function dyn_column_reasons(cols::AbstractVector{Symbol})::Vector{String}
    return unique(String[get(_DYN_COLUMN_REASON, col, "the selected dynamic model")
                         for col in cols])
end

# validate_dyn_config!(cfg::RunConfig) lives in engine.jl (after RunConfig is defined).

"""Symmetric speed-deviation tolerance tuple `(−Δω_tol_pu, +Δω_tol_pu)` in p.u."""
function Δω_tol_tuple(Δω_tol_pu::Float64)::Tuple{Float64, Float64}
    return (-Δω_tol_pu, Δω_tol_pu)
end

"""
    Δω_tol_tuple(dyn::DynModelConfig) -> (lower, upper)

Signed speed-deviation tolerance pair in p.u., resolved from `Δω_tol_pu` and the optional
`Δω_tol_pu_lower` / `Δω_tol_pu_upper` overrides. Both overrides are **positive magnitudes**,
so `lower = 0.05`, `upper = 0.03` yields `(-0.05, 0.03)` — the same convention as
`δ_tol_deg_lower` / `δ_tol_deg_upper` in `common_ts_parameters`.
"""
function Δω_tol_tuple(dyn::DynModelConfig)::Tuple{Float64, Float64}
    lo = something(dyn.Δω_tol_pu_lower, dyn.Δω_tol_pu)
    hi = something(dyn.Δω_tol_pu_upper, dyn.Δω_tol_pu)
    lo > 0 || throw(ArgumentError("Δω_tol_pu_lower / Δω_tol_pu must be positive."))
    hi > 0 || throw(ArgumentError("Δω_tol_pu_upper / Δω_tol_pu must be positive."))
    return (-lo, hi)
end

"""
    mechanical_power_for_swing(mech_power_mode, P_g, dyn_model_dict)

JuMP variables used as mechanical power in swing and propagated COI bounds.
`USE_PG` → dispatch `P_g`; `USE_PM` → explicit `P_m` (requires prior init).
"""
function mechanical_power_for_swing(
    mech_power_mode::MechPowerMode,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    dyn_model_dict::OrderedDict{Symbol, Any},
)::OrderedDict{Int, JuMP.VariableRef}
    if mech_power_mode == USE_PG
        return P_g
    else
        haskey(dyn_model_dict, :vars) && haskey(dyn_model_dict[:vars], :P_m) ||
            throw(ArgumentError("mech_power_mode=USE_PM but :P_m is not in dyn_model_dict[:vars]."))
        return dyn_model_dict[:vars][:P_m]
    end
end

"""Store `P_g` and swing mechanical-power refs after pre-fault init."""
function register_mech_power_refs!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    mech_power_mode::MechPowerMode,
)
    dyn_model_dict[:refs] = OrderedDict{Symbol, Any}(
        :P_g    => P_g,
        :P_mech => mechanical_power_for_swing(mech_power_mode, P_g, dyn_model_dict),
    )
    dyn_model_dict[:mech_power_mode] = mech_power_mode
    return dyn_model_dict
end
