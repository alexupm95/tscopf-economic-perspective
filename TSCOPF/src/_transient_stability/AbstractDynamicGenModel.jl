#=
================================================================================
 AbstractDynamicGenModel.jl  —  factory dispatch for TS dynamic builders
================================================================================
 Phase 1.2: `RunConfig.dyn_model` selects which physics/network path to build.
 `dynamic_gen_model` maps `DynModelConfig` → a concrete `AbstractDynamicGenModel`
 subtype; `Build_Dynamic_Model!` dispatches on that type.

 Default path: `ClassicalKronModel` → existing `Make_Dynamic_Model_tsred!` or
 `Make_Dynamic_Model_tsredlinear!`.  Later phases add full-bus and dq types here.
================================================================================
=#

# --- abstract root: every dynamic formulation implements trait queries below ---

abstract type AbstractDynamicGenModel end

# ==================================================================================
# ClassicalKronModel — 2nd-order generator + Kron-reduced Ybus (production path)
# ==================================================================================

"""
    ClassicalKronModel

Snapshot of the `DynModelConfig` fields that the Kron builders actually consume.
`linearize=true` routes to the Taylor-reduced Pe formulation (TSC-DCOPF);
`linearize=false` keeps the nonlinear `sin/cos` Pe path (TSC-ACOPF).
"""
struct ClassicalKronModel <: AbstractDynamicGenModel
    linearize::Bool               # false → tsred (ACOPF); true → tsredlinear (DCOPF)
    mech_power_mode::MechPowerMode # USE_PG: swing uses P_g; USE_PM: explicit P_m
    constrain_δ::Bool             # build the rotor-angle corridor
    bound_style_δ::Symbol         # :swing_propagated | :coi_box
    δ_ref_gen_id::Union{Nothing, Int} # reference machine id (only for :ref_gen)
    constrain_Δω::Bool            # build the speed corridor
    bound_style_Δω::Symbol        # :coi_box
    Δω_tol::Tuple{Float64, Float64} # signed (lower, upper) of that box [p.u.]; may be asymmetric
end

"""Trait: network reduction used by this model (always Kron for this type)."""
network_form(::ClassicalKronModel) = KRON_REDUCED

"""Trait: generator dynamic order (classical swing, no dq states)."""
gen_order(::ClassicalKronModel) = CLASSICAL_2ND

# ==================================================================================
# ClassicalFullBusModel — 2nd-order generator + full Ybus (Phase 3)
# ==================================================================================

"""
    ClassicalFullBusModel

Classical swing on the full sparse admittance matrix (no Kron reduction).
Reference full-network TSC-ACOPF formulation.  TSC-ACOPF only; default
`mech_power_mode=USE_PM`.  The two stability corridors follow the same
conventions as `ClassicalKronModel`.  `ode_first_step` selects the
integration rule for the first row of each window (the Kron path cannot).
"""
struct ClassicalFullBusModel <: AbstractDynamicGenModel
    mech_power_mode::MechPowerMode
    constrain_δ::Bool
    bound_style_δ::Symbol         # :swing_propagated | :coi_box
    δ_ref_gen_id::Union{Nothing, Int}
    constrain_Δω::Bool
    bound_style_Δω::Symbol
    Δω_tol::Tuple{Float64, Float64} # signed (lower, upper) [p.u.]; may be asymmetric
    zip_load_p::NTuple{3, Float64}  # active-demand (Z, I, P) split
    zip_load_q::NTuple{3, Float64}  # reactive-demand (Z, I, P) split; independent of P
    include_governor::Bool        # optional TGOV1 turbine governor (time-varying P_mech)
    governor_limiter::GovernorLimiter # valve saturation treatment (only if include_governor)
    ode_first_step::Symbol        # :trapezoidal | :backward_euler (first row of each window)
end

network_form(::ClassicalFullBusModel) = FULL_BUS
gen_order(::ClassicalFullBusModel) = CLASSICAL_2ND

# ==================================================================================
# DqFullBusModel — 4th-order dq machine + full Ybus (Phase 4)
# ==================================================================================

"""
    DqFullBusModel

Fourth-order dq-axis generator on the full sparse admittance matrix.
Requires `gen_dynamic_data_full.csv` (or equivalent columns).  Optional
first-order AVR (`include_avr`) drives time-varying `E_fd`; optional TGOV1
governor (`include_governor`) drives time-varying `P_mech(t)`.
"""
struct DqFullBusModel <: AbstractDynamicGenModel
    mech_power_mode::MechPowerMode
    constrain_δ::Bool
    bound_style_δ::Symbol
    δ_ref_gen_id::Union{Nothing, Int}
    constrain_Δω::Bool
    bound_style_Δω::Symbol
    Δω_tol::Tuple{Float64, Float64} # signed (lower, upper) [p.u.]; may be asymmetric
    zip_load_p::NTuple{3, Float64}  # active-demand (Z, I, P) split
    zip_load_q::NTuple{3, Float64}  # reactive-demand (Z, I, P) split; independent of P
    dq_speed_dev_in_algebra::Bool
    include_avr::Bool
    include_governor::Bool
    governor_limiter::GovernorLimiter
    ode_first_step::Symbol
    gfm_integrator::Symbol          # GFM filters + Q–V PI only; see DynModelConfig
    gfm_in_coi::Bool                # converters in the COI + `:coi_box` corridors
    gfm_pq_overload::Bool           # REGFM_A1 Fig. 3 P and Q overload PIs
    allow_all_gfm::Bool             # permit a converter-only fleet; see DynModelConfig
end

network_form(::DqFullBusModel) = FULL_BUS
gen_order(::DqFullBusModel) = DQ_4TH

# ==================================================================================
# Factory — DynModelConfig → concrete model (throws for unimplemented combos)
# ==================================================================================

"""
    dynamic_gen_model(dyn::DynModelConfig; linearize) -> AbstractDynamicGenModel

Central selector called from `Build_Dynamic_Model!`.

| `linearize` | Steady-state OPF | Pe constraint in swing |
|-------------|------------------|------------------------|
| `false`     | TSC-ACOPF        | Nonlinear on `Yred`    |
| `true`      | TSC-DCOPF        | Taylor around `δ_ref`  |

Unimplemented `(gen_order, network_form)` pairs throw `ArgumentError` so invalid
configurations fail before any JuMP variables are created.
"""
function dynamic_gen_model(dyn::DynModelConfig; linearize::Bool)::AbstractDynamicGenModel
    if dyn.gen_order == CLASSICAL_2ND && dyn.network_form == KRON_REDUCED
        return ClassicalKronModel(
            linearize,
            dyn.mech_power_mode,
            dyn.constrain_δ,
            dyn.bound_style_δ,
            dyn.δ_ref_gen_id,
            dyn.constrain_Δω,
            dyn.bound_style_Δω,
            Δω_tol_tuple(dyn),
        )
    elseif dyn.gen_order == CLASSICAL_2ND && dyn.network_form == FULL_BUS
        linearize && throw(ArgumentError(
            "FULL_BUS classical dynamics require TSC-ACOPF (linearize=false)."))
        return ClassicalFullBusModel(
            dyn.mech_power_mode,
            dyn.constrain_δ,
            dyn.bound_style_δ,
            dyn.δ_ref_gen_id,
            dyn.constrain_Δω,
            dyn.bound_style_Δω,
            Δω_tol_tuple(dyn),
            dyn.zip_load_p,
            dyn.zip_load_q,
            dyn.include_governor,
            dyn.governor_limiter,
            dyn.ode_first_step,
        )
    elseif dyn.gen_order == DQ_4TH && dyn.network_form == FULL_BUS
        return DqFullBusModel(
            dyn.mech_power_mode,
            dyn.constrain_δ,
            dyn.bound_style_δ,
            dyn.δ_ref_gen_id,
            dyn.constrain_Δω,
            dyn.bound_style_Δω,
            Δω_tol_tuple(dyn),
            dyn.zip_load_p,
            dyn.zip_load_q,
            dyn.dq_speed_dev_in_algebra,
            dyn.include_avr,
            dyn.include_governor,
            dyn.governor_limiter,
            dyn.ode_first_step,
            dyn.gfm_integrator,
            dyn.gfm_in_coi,
            dyn.gfm_pq_overload,
            dyn.allow_all_gfm,
        )
    elseif dyn.gen_order == DQ_4TH
        throw(ArgumentError(
            "DQ_4TH requires FULL_BUS network_form (Kron reduction is inconsistent with dq dynamics)."))
    else
        throw(ArgumentError("Unsupported DynModelConfig: $dyn"))
    end
end

# ==================================================================================
# Metadata — record which factory product was built (results / debugging)
# ==================================================================================

"""
    register_dyn_model_meta!(dyn_model_dict, gen_model)

Augment `dyn_model_dict[:meta]` after assembly so saved results and tests can
read back the model type without re-inspecting `RunConfig`.
"""
function register_dyn_model_meta!(
    dyn_model_dict::OrderedDict{Symbol, Any},
    gen_model::AbstractDynamicGenModel,
)
    haskey(dyn_model_dict, :meta) || (dyn_model_dict[:meta] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:meta][:gen_model_type] = string(nameof(typeof(gen_model)))
    dyn_model_dict[:meta][:gen_order] = string(gen_order(gen_model))
    dyn_model_dict[:meta][:network_form] = string(network_form(gen_model))
    # Stability corridors — the same five fields on every model type, so they are written
    # once here rather than repeated (and, previously, half-forgotten on the Kron branch).
    for key in (:constrain_δ, :bound_style_δ, :δ_ref_gen_id, :constrain_Δω,
                :bound_style_Δω, :Δω_tol)
        dyn_model_dict[:meta][key] = getfield(gen_model, key)
    end

    if gen_model isa ClassicalKronModel
        dyn_model_dict[:meta][:linearize] = gen_model.linearize
    elseif gen_model isa ClassicalFullBusModel
        dyn_model_dict[:meta][:zip_load_p] = collect(gen_model.zip_load_p)
        dyn_model_dict[:meta][:zip_load_q] = collect(gen_model.zip_load_q)
        dyn_model_dict[:meta][:include_governor] = gen_model.include_governor
        dyn_model_dict[:meta][:governor_limiter] = string(gen_model.governor_limiter)
        dyn_model_dict[:meta][:ode_first_step] = gen_model.ode_first_step
    elseif gen_model isa DqFullBusModel
        dyn_model_dict[:meta][:zip_load_p] = collect(gen_model.zip_load_p)
        dyn_model_dict[:meta][:zip_load_q] = collect(gen_model.zip_load_q)
        dyn_model_dict[:meta][:dq_speed_dev_in_algebra] = gen_model.dq_speed_dev_in_algebra
        dyn_model_dict[:meta][:include_avr] = gen_model.include_avr
        dyn_model_dict[:meta][:include_governor] = gen_model.include_governor
        dyn_model_dict[:meta][:governor_limiter] = string(gen_model.governor_limiter)
        dyn_model_dict[:meta][:ode_first_step] = gen_model.ode_first_step
        dyn_model_dict[:meta][:gfm_integrator] = gen_model.gfm_integrator
        dyn_model_dict[:meta][:gfm_in_coi] = gen_model.gfm_in_coi
        dyn_model_dict[:meta][:gfm_pq_overload] = gen_model.gfm_pq_overload
        dyn_model_dict[:meta][:allow_all_gfm] = gen_model.allow_all_gfm
    end
    return dyn_model_dict
end
