#=
================================================================================
 TsConfig.jl  —  transient-stability simulation + builder configuration
================================================================================
 Replaces `Define_Inputs_4TS_Optim`, `Time_Parameters_4_TS_*`, and
 `Common_Parameters_4_TS` from `config_file.jl`.
================================================================================
=#

# --- simulation timing & stability limits -------------------------------------

"""
    TsSimulationConfig

Transient-stability simulation timing and stability limits
(`RunConfig.transient.simulation`).

Controls the trapezoidal time grid (`t_step`, `t_start_sim`, `t_end_sim`),
disturbance timing (`t_start_fault`, `clearing_time` for SC only),
synchronous frequency `f_syn`, and the rotor-angle stability tolerance
`δ_tol_deg` (symmetric default, in degrees).  Optional `δ_tol_deg_lower` /
`δ_tol_deg_upper` set independent below/above COI limits (degrees, positive).

`δ_pair_tol_deg` is a separate, symmetric limit read only by
`bound_style_δ = :pairwise`, which bounds |δ_i − δ_j| for every unordered pair
of synchronised units instead of measuring each unit against a common reference.

See the parameter reference, section 3.1.
"""
Base.@kwdef struct TsSimulationConfig
    δ_tol_deg::Float64 = 90.0         # symmetric |δ − δ_COI| limit [deg] when overrides unset
    δ_tol_deg_lower::Union{Nothing, Float64} = nothing  # below COI [deg]; default → δ_tol_deg
    δ_tol_deg_upper::Union{Nothing, Float64} = nothing  # above COI [deg]; default → δ_tol_deg
    # Pairwise separation limit, used only by `bound_style_δ = :pairwise`.  Symmetric by
    # construction: the corridor bounds |δ_i − δ_j|, which has no "above"/"below" reading —
    # the pair (i,j) and the pair (j,i) are the same constraint seen from the two sides, and
    # both sides are already built.  Kept apart from `δ_tol_deg` because the two numbers live
    # on different scales: a COI corridor is tens of degrees, a loss-of-synchronism
    # separation limit is of the order of 180.
    δ_pair_tol_deg::Float64 = 180.0   # symmetric |δ_i − δ_j| limit [deg], :pairwise only
    t_start_sim::Float64 = 0.0
    t_end_sim::Float64 = 5.0
    t_step::Float64 = 0.01
    t_start_fault::Float64 = 0.01
    clearing_time::Float64 = 0.3      # SC only [s]
    f_syn::Float64 = 50.0             # synchronous frequency [Hz]
    Δω_0::Float64 = 0.0                 # initial speed deviation [p.u.]
end

function copy_ts_simulation(sim::TsSimulationConfig)::TsSimulationConfig
    return TsSimulationConfig(; (f => getfield(sim, f) for f in fieldnames(TsSimulationConfig))...)
end

"""Rotor-angle tolerance tuple and synchronous parameters for TS builders."""
function common_ts_parameters(sim::TsSimulationConfig)::Tuple{Tuple{Float64, Float64}, Float64, Float64, Float64}
    lo_deg = something(sim.δ_tol_deg_lower, sim.δ_tol_deg)
    hi_deg = something(sim.δ_tol_deg_upper, sim.δ_tol_deg)
    lo_deg > 0 || throw(ArgumentError("δ_tol_deg_lower / δ_tol_deg must be positive."))
    hi_deg > 0 || throw(ArgumentError("δ_tol_deg_upper / δ_tol_deg must be positive."))
    δ_tol_tuple = (deg2rad(-lo_deg), deg2rad(hi_deg))
    ω_syn = 2π * sim.f_syn
    return δ_tol_tuple, sim.f_syn, ω_syn, sim.Δω_0
end

"""
    δ_pair_tol_tuple(sim) -> Tuple{Float64, Float64}

Signed radian band for the pairwise separation corridor, `(-Δ, +Δ)`.

Deliberately *not* part of `common_ts_parameters`: that function's return tuple is
destructured positionally in all four builders, so widening it would touch every builder for
a value only one corridor style reads.  The builders put this straight into `meta` instead.
"""
function δ_pair_tol_tuple(sim::TsSimulationConfig)::Tuple{Float64, Float64}
    sim.δ_pair_tol_deg > 0 || throw(ArgumentError("δ_pair_tol_deg must be positive."))
    return (deg2rad(-sim.δ_pair_tol_deg), deg2rad(sim.δ_pair_tol_deg))
end

"""SC disturbance time windows from `TsSimulationConfig`."""
function time_windows_sc(sim::TsSimulationConfig)::Tuple{
    Float64, Float64, Float64, Float64, Float64, Float64,
    Vector{Float64}, Vector{Float64}, Vector{Float64},
}
    t_clear_fault = sim.clearing_time + sim.t_start_fault
    t_window_fault = collect(sim.t_start_fault:sim.t_step:t_clear_fault)
    t_window_postf = collect(t_clear_fault + sim.t_step:sim.t_step:sim.t_end_sim)
    t_window_total = vcat(t_window_fault, t_window_postf)
    sim.t_end_sim <= t_clear_fault &&
        throw(ArgumentError("t_end_sim must exceed t_clear_fault (SC timing)."))
    return (sim.t_start_sim, sim.t_end_sim, sim.t_step, sim.t_start_fault,
            sim.clearing_time, t_clear_fault, t_window_fault, t_window_postf, t_window_total)
end

"""GL disturbance time windows from `TsSimulationConfig`."""
function time_windows_gld(sim::TsSimulationConfig)::Tuple{
    Float64, Float64, Float64, Float64, Vector{Float64}, Vector{Float64},
}
    t_window_fault = collect(sim.t_start_fault:sim.t_step:sim.t_end_sim)
    sim.t_end_sim <= sim.t_start_fault &&
        throw(ArgumentError("t_end_sim must exceed t_start_fault (GLD timing)."))
    return (sim.t_start_sim, sim.t_end_sim, sim.t_step, sim.t_start_fault,
            t_window_fault, t_window_fault)
end

# --- TS model builder toggles -------------------------------------------------

"""
    TsBuilderConfig

Toggles which transient-stability inequality families and explicit variable
bounds the model builder constructs (`RunConfig.transient.builder`).

Pre-fault init equalities that couple the OPF operating point to the initial
dynamic state (`P_g`/`Q_g`/`P_m` links) are **always built by the path** — they
are physics, not user toggles:

- `eq_const_P_init` — always (AC and linearized Kron)
- `eq_const_Q_init` — AC classical / DQ only (omitted on TSC-DCOPF / Kron linear)
- `eq_const_Pm_init` — only when `mech_power_mode = USE_PM` (omitted under `USE_PG`)

Fault- and post-fault-window equalities (swing, Pe, COI, nodal balance) are
likewise always built.

COI-referenced stability inequalities and optional explicit variable bounds
control what limits the solution space.

See the parameter reference, section 3.2.
"""
Base.@kwdef struct TsBuilderConfig
    # Variable bounds as explicit inequalities (mostly off in production)
    bound_E::Bool = false
    bound_δ::Bool = false
    bound_P_m::Bool = false
    bound_δ_tf::Bool = false
    bound_Δω_tf::Bool = false
    bound_Pe_tf::Bool = false
    bound_Qe_tf::Bool = false
    bound_δCOI_tf::Bool = false
    bound_δ_tpf::Bool = false
    bound_Δω_tpf::Bool = false
    bound_Pe_tpf::Bool = false
    bound_Qe_tpf::Bool = false
    bound_δCOI_tpf::Bool = false

    # DQ pre-fault algebraic (DQ_4TH only; default off)
    bound_Ed::Bool = false
    bound_Eq::Bool = false
    bound_Id::Bool = false
    bound_Iq::Bool = false

    # AVR / governor pre-fault set-points (default off)
    bound_V_ref::Bool = false
    bound_P_ref::Bool = false

    # DQ machine states — fault window (default off)
    bound_Ed_tf::Bool = false
    bound_Eq_tf::Bool = false
    bound_Id_tf::Bool = false
    bound_Iq_tf::Bool = false
    bound_Te_tf::Bool = false

    # DQ machine states — post-fault window (default off)
    bound_Ed_tpf::Bool = false
    bound_Eq_tpf::Bool = false
    bound_Id_tpf::Bool = false
    bound_Iq_tpf::Bool = false
    bound_Te_tpf::Bool = false

    # AVR research boxes on E_fd_unlim (default off; saturation physics unchanged)
    bound_E_fd_unlim_tf::Bool = false
    bound_E_fd_unlim_tpf::Bool = false

    # Governor research boxes on valve/mechanical trajectories (default off)
    bound_Pv_tf::Bool = false
    bound_Pv_tpf::Bool = false
    bound_Pm_tf::Bool = false
    bound_Pm_tpf::Bool = false

    # GFM boxes (read only when allow_gfm). Default **true**, unlike the SG
    # families above: these are guard-rails the nonconvex current limiter needs
    # to converge, so the shipped behaviour keeps them on. Switch one off to see
    # what it was holding — expect a harder solve, not a different optimum.
    # Values come from the `gfm_*` fields of `TsBoundLimitsConfig`.
    bound_gfm_δ::Bool = true
    bound_gfm_P_meas_tf::Bool = true
    bound_gfm_Q_meas_tf::Bool = true
    bound_gfm_V_meas_tf::Bool = true
    bound_gfm_E_int_raw_tf::Bool = true
    bound_gfm_E_int_tf::Bool = true
    bound_gfm_E_droop_raw_tf::Bool = true
    bound_gfm_E_droop_tf::Bool = true
    bound_gfm_Id_tf::Bool = true
    bound_gfm_Iq_tf::Bool = true
    bound_gfm_P_meas_tpf::Bool = true
    bound_gfm_Q_meas_tpf::Bool = true
    bound_gfm_V_meas_tpf::Bool = true
    bound_gfm_E_int_raw_tpf::Bool = true
    bound_gfm_E_int_tpf::Bool = true
    bound_gfm_E_droop_raw_tpf::Bool = true
    bound_gfm_E_droop_tpf::Bool = true
    bound_gfm_Id_tpf::Bool = true
    bound_gfm_Iq_tpf::Bool = true
    # P/Q overload raw and clipped states. Read only when `gfm_pq_overload` is on.
    # One switch per window covers every branch: the four raw columns share a box,
    # and the four clipped columns share the other.
    bound_gfm_pq_raw_tf::Bool = true
    bound_gfm_pq_raw_tpf::Bool = true
    bound_gfm_pq_clip_tf::Bool = true
    bound_gfm_pq_clip_tpf::Bool = true

    # COI-referenced stability bounds
    ineq_δ_COI_tf_lower::Bool = true
    ineq_δ_COI_tf_upper::Bool = true
    ineq_δ_COI_tpf_lower::Bool = true
    ineq_δ_COI_tpf_upper::Bool = true

    # Transient bound limit values (physical + placeholder); see `TsBoundLimitsConfig`
    limits::TsBoundLimitsConfig = TsBoundLimitsConfig()

    # How simple box limits are attached (see `BoundEncoding`). δ-COI stability
    # bounds and expression targets (e.g. Kron Qe) always use CONSTRAINT form.
    bound_encoding::BoundEncoding = CONSTRAINT
end

function copy_ts_builder(b::TsBuilderConfig)::TsBuilderConfig
    return TsBuilderConfig(; (f => getfield(b, f) for f in fieldnames(TsBuilderConfig))...)
end

"""
    build_ts_input_param(builder::TsBuilderConfig) -> OrderedDict

Legacy container consumed by TS model builders (mirrors `Define_Inputs_4TS_Optim`).
"""
function build_ts_input_param(builder::TsBuilderConfig)::OrderedDict{Symbol, Any}
    ts_input_param = OrderedDict{Symbol, Any}()
    ts_input_param[:var_names] = OrderedDict(
        :E => "E", :δ => "δ", :P_m => "P_m",
        :Ed => "Ed", :Eq => "Eq", :Id => "Id", :Iq => "Iq",
        :V_ref => "V_ref", :P_ref => "P_ref",
        :δ_tf => "δ_tf", :Δω_tf => "Δω_tf", :Pe_tf => "Pe_tf", :Qe_tf => "Qe_tf",
        :δCOI_tf => "δCOI_tf",
        :Ed_tf => "Ed_tf", :Eq_tf => "Eq_tf", :Id_tf => "Id_tf", :Iq_tf => "Iq_tf", :Te_tf => "Te_tf",
        :E_fd_unlim_tf => "E_fd_unlim_tf", :Pv_tf => "Pv_tf", :Pm_tf => "Pm_tf",
        :δ_tpf => "δ_tpf", :Δω_tpf => "Δω_tpf", :Pe_tpf => "Pe_tpf", :Qe_tpf => "Qe_tpf",
        :δCOI_tpf => "δCOI_tpf",
        :Ed_tpf => "Ed_tpf", :Eq_tpf => "Eq_tpf", :Id_tpf => "Id_tpf", :Iq_tpf => "Iq_tpf", :Te_tpf => "Te_tpf",
        :E_fd_unlim_tpf => "E_fd_unlim_tpf", :Pv_tpf => "Pv_tpf", :Pm_tpf => "Pm_tpf",
    )
    ts_input_param[:var_bounds] = OrderedDict(
        :E => builder.bound_E, :δ => builder.bound_δ, :P_m => builder.bound_P_m,
        :Ed => builder.bound_Ed, :Eq => builder.bound_Eq,
        :Id => builder.bound_Id, :Iq => builder.bound_Iq,
        :V_ref => builder.bound_V_ref, :P_ref => builder.bound_P_ref,
        :δ_tf => builder.bound_δ_tf, :Δω_tf => builder.bound_Δω_tf,
        :Pe_tf => builder.bound_Pe_tf, :Qe_tf => builder.bound_Qe_tf,
        :δCOI_tf => builder.bound_δCOI_tf,
        :Ed_tf => builder.bound_Ed_tf, :Eq_tf => builder.bound_Eq_tf,
        :Id_tf => builder.bound_Id_tf, :Iq_tf => builder.bound_Iq_tf, :Te_tf => builder.bound_Te_tf,
        :E_fd_unlim_tf => builder.bound_E_fd_unlim_tf,
        :Pv_tf => builder.bound_Pv_tf, :Pm_tf => builder.bound_Pm_tf,
        :δ_tpf => builder.bound_δ_tpf, :Δω_tpf => builder.bound_Δω_tpf,
        :Pe_tpf => builder.bound_Pe_tpf, :Qe_tpf => builder.bound_Qe_tpf,
        :δCOI_tpf => builder.bound_δCOI_tpf,
        :Ed_tpf => builder.bound_Ed_tpf, :Eq_tpf => builder.bound_Eq_tpf,
        :Id_tpf => builder.bound_Id_tpf, :Iq_tpf => builder.bound_Iq_tpf, :Te_tpf => builder.bound_Te_tpf,
        :E_fd_unlim_tpf => builder.bound_E_fd_unlim_tpf,
        :Pv_tpf => builder.bound_Pv_tpf, :Pm_tpf => builder.bound_Pm_tpf,
        # GFM boxes (default true; only read when allow_gfm)
        :gfm_δ => builder.bound_gfm_δ,
        :gfm_P_meas_tf => builder.bound_gfm_P_meas_tf,
        :gfm_Q_meas_tf => builder.bound_gfm_Q_meas_tf,
        :gfm_V_meas_tf => builder.bound_gfm_V_meas_tf,
        :gfm_E_int_raw_tf => builder.bound_gfm_E_int_raw_tf,
        :gfm_E_int_tf => builder.bound_gfm_E_int_tf,
        :gfm_E_droop_raw_tf => builder.bound_gfm_E_droop_raw_tf,
        :gfm_E_droop_tf => builder.bound_gfm_E_droop_tf,
        :gfm_Id_tf => builder.bound_gfm_Id_tf,
        :gfm_Iq_tf => builder.bound_gfm_Iq_tf,
        :gfm_P_meas_tpf => builder.bound_gfm_P_meas_tpf,
        :gfm_Q_meas_tpf => builder.bound_gfm_Q_meas_tpf,
        :gfm_V_meas_tpf => builder.bound_gfm_V_meas_tpf,
        :gfm_E_int_raw_tpf => builder.bound_gfm_E_int_raw_tpf,
        :gfm_E_int_tpf => builder.bound_gfm_E_int_tpf,
        :gfm_E_droop_raw_tpf => builder.bound_gfm_E_droop_raw_tpf,
        :gfm_E_droop_tpf => builder.bound_gfm_E_droop_tpf,
        :gfm_Id_tpf => builder.bound_gfm_Id_tpf,
        :gfm_Iq_tpf => builder.bound_gfm_Iq_tpf,
        :gfm_pq_raw_tf => builder.bound_gfm_pq_raw_tf,
        :gfm_pq_raw_tpf => builder.bound_gfm_pq_raw_tpf,
        :gfm_pq_clip_tf => builder.bound_gfm_pq_clip_tf,
        :gfm_pq_clip_tpf => builder.bound_gfm_pq_clip_tpf,
    )
    ts_input_param[:ineq_cons] = OrderedDict(
        :ineq_const_δ_COI_tf_lower => builder.ineq_δ_COI_tf_lower,
        :ineq_const_δ_COI_tf_upper => builder.ineq_δ_COI_tf_upper,
        :ineq_const_δ_COI_tpf_lower => builder.ineq_δ_COI_tpf_lower,
        :ineq_const_δ_COI_tpf_upper => builder.ineq_δ_COI_tpf_upper,
    )
    ts_input_param[:limits] = builder.limits
    ts_input_param[:bound_encoding] = builder.bound_encoding
    return ts_input_param
end

# --- transient bundle (avenue 2) ----------------------------------------------

"""
    TransientConfig

Transient-stability analysis bundle (`RunConfig.transient`, avenue 2),
required when `RunConfig.trans_stab = true`.

Nests `simulation::TsSimulationConfig` (timing and δ_tol),
`builder::TsBuilderConfig` (which TS constraints are built), and
`dyn_model::DynModelConfig` (physics, network form, and disturbance).
`gen_dynamic_filename` names the SG machine-parameter CSV; when
`dyn_model.allow_gfm=true`, `gfm_dynamic_filename` names the separate GFM
parameter CSV (default `gfm_dynamic_data.csv`).

See the parameter reference, section 3, and the configuration map, section 3.
"""
Base.@kwdef struct TransientConfig
    simulation::TsSimulationConfig = TsSimulationConfig()
    builder::TsBuilderConfig = TsBuilderConfig()
    dyn_model::DynModelConfig = DynModelConfig()
    gen_dynamic_filename::String = "gen_dynamic_data.csv"
    gfm_dynamic_filename::String = "gfm_dynamic_data.csv"
    # Further short-circuit contingencies (rows of `contingencies.csv`) enforced in the
    # same optimization, on the same pre-fault operating point as `dyn_model.fault`.
    # Each adds its own transient block; results go to `Transient_Stability_c<id>/`.
    extra_contingency_ids::Vector{Int} = Int[]
end

"""Defaults matching the legacy integrated TSC-ACOPF run (contingency 2, δ_tol 90°)."""
function default_transient_config()::TransientConfig
    return TransientConfig()
end

function copy_transient_config(tc::TransientConfig)::TransientConfig
    dm = tc.dyn_model
    return TransientConfig(
        simulation=copy_ts_simulation(tc.simulation),
        builder=copy_ts_builder(tc.builder),
        dyn_model=DynModelConfig(;
            (f => f == :fault ? copy_fault_config(dm.fault) : getfield(dm, f)
             for f in fieldnames(DynModelConfig))...),
        gen_dynamic_filename=tc.gen_dynamic_filename,
        gfm_dynamic_filename=tc.gfm_dynamic_filename,
        extra_contingency_ids=copy(tc.extra_contingency_ids),
    )
end

"""Copy `tc` overriding selected nested fields (used by δ_tol sweeps)."""
function reconfigure_transient(tc::TransientConfig; kwargs...)
    base = (; (f => getfield(tc, f) for f in fieldnames(TransientConfig))...)
    return TransientConfig(; base..., kwargs...)
end

function validate_transient_config!(tc::TransientConfig)
    sim = tc.simulation
    sim.δ_tol_deg > 0 ||
        throw(ArgumentError("δ_tol_deg must be positive."))
    if sim.δ_tol_deg_lower !== nothing
        sim.δ_tol_deg_lower > 0 ||
            throw(ArgumentError("δ_tol_deg_lower must be positive when set."))
    end
    if sim.δ_tol_deg_upper !== nothing
        sim.δ_tol_deg_upper > 0 ||
            throw(ArgumentError("δ_tol_deg_upper must be positive when set."))
    end
    # Checked unconditionally, not only under `bound_style_δ = :pairwise`: the value is
    # written into every builder's meta, and a run that flips the style later should not
    # discover the knob is unusable only then.
    sim.δ_pair_tol_deg > 0 ||
        throw(ArgumentError("δ_pair_tol_deg must be positive."))
    tc.simulation.t_step > 0 ||
        throw(ArgumentError("t_step must be positive."))
    validate_ts_bound_limits!(tc.builder.limits)
    return nothing
end
