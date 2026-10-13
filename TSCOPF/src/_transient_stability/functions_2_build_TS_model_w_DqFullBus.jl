#=
================================================================================
 functions_2_build_TS_model_w_DqFullBus.jl — TSC-ACOPF FULL_BUS dq 4th-order
================================================================================
 Milestone 1: 4th-order dq machine on the full sparse Ybus (no Kron reduction).
 Optional AVR (time-varying E_fd) and optional TGOV1 governor (time-varying P_mech).

 Architecture mirrors `functions_2_build_TS_model_w_FullBus.jl`:
 - Same SC/GL fault windows, COI bounds, swing equations, and ZIP nodal KCL shell
   (independent `zip_load_p` / `zip_load_q` splits, as on the classical FULL_BUS path,
    refined per bus by an optional `load_data.csv` through `DLOAD`)
 - Machine-specific pieces live in `functions_4_TS_dq_{helpers,variables,eqconst}.jl`
 - Generator injection enters KCL via dq currents (Id, Iq), not classical Pe = E·V·sin(δ−θ)/Xd′
 - `dq_speed_dev_in_algebra` toggles (1+Δω) in Pe and stator Vd/Vq (RMS-style); swing unchanged
================================================================================
=#

"""
Append the FULL_BUS dq transient-stability sub-model onto an existing ACOPF `model`.

See `Make_Dynamic_Model_fullbus!` for the SC/GL window layout. DQ-specific requirements:
`mech_power_mode=USE_PM`, `bound_style_δ=:coi_box`, and full machine columns in
`gen_dynamic_data` (use `gen_dynamic_data_full.csv`).
"""
function Make_Dynamic_Model_dqfullbus!(
    model::JuMP.Model,
    opf_dict::OrderedDict{Symbol, Any},
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    DCIR::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nBUS::Int64,
    nGEN::Int64,
    nCIR::Int64,
    ts_fault_details::OrderedDict{Symbol, Any};
    simulation::TsSimulationConfig=TsSimulationConfig(),
    ts_builder::TsBuilderConfig=TsBuilderConfig(),
    mech_power_mode::MechPowerMode=USE_PM,
    constrain_δ::Bool=true,
    bound_style_δ::Symbol=:coi_box,
    δ_ref_gen_id::Union{Nothing, Int}=nothing,
    constrain_Δω::Bool=false,
    bound_style_Δω::Symbol=:coi_box,
    Δω_tol::Tuple{Float64, Float64}=(-0.5, 0.5),
    zip_load_p::NTuple{3, Float64}=(1.0, 0.0, 0.0),
    zip_load_q::NTuple{3, Float64}=(1.0, 0.0, 0.0),
    dq_speed_dev_in_algebra::Bool=true,
    include_avr::Bool=false,
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
    ode_first_step::Symbol=:trapezoidal,
    gfm_integrator::Symbol=:backward_euler,
    gfm_in_coi::Bool=false,
    gfm_pq_overload::Bool=false,
    allow_all_gfm::Bool=false,
    steady_state_hints::SteadyStateHints,
    DGFM::Union{Nothing, DataFrame}=nothing,
    DLOAD::Union{Nothing, DataFrame}=nothing,
)
    mech_power_mode == USE_PM || throw(ArgumentError(
        "DQ FULL_BUS path requires mech_power_mode=USE_PM."))
    ode_first_step ∈ (:trapezoidal, :backward_euler) || throw(ArgumentError(
        "ode_first_step must be :trapezoidal or :backward_euler (got $ode_first_step)."))
    gfm_integrator ∈ (:follow_ode_first_step, :backward_euler, :trapezoidal) ||
        throw(ArgumentError("gfm_integrator must be :follow_ode_first_step, " *
                            ":backward_euler or :trapezoidal (got $gfm_integrator)."))
    dyn_cfg = DynModelConfig(
        gen_order=DQ_4TH, include_avr=include_avr, include_governor=include_governor,
        ode_first_step=ode_first_step, gfm_integrator=gfm_integrator)
    has_required_dyn_columns(dyn_cfg, DGEN_DYN) || throw(ArgumentError(
        "DQ_4TH requires finite columns in gen_dynamic_data: " *
        join(string.(required_dyn_column_names(dyn_cfg)), ", ") * "."))

    ts_input_param = build_ts_input_param(ts_builder)
    # Independent (Z, I, P) splits for the active and reactive nodal balances, one row
    # per bus (`zip_by_bus`: system-wide tuple, or the demand-weighted mean of the loads
    # at that bus when a `load_data.csv` table is present).
    ZIP_P, ZIP_Q = zip_by_bus(DBUS, DLOAD, zip_load_p, zip_load_q; nBUS=nBUS)
    dyn_model_dict = OrderedDict{Symbol, Any}()
    dyn_model_dict[:vars] = OrderedDict{Symbol, Any}()
    dyn_model_dict[:eq_const] = OrderedDict{Symbol, Any}()
    dyn_model_dict[:ineq_const] = OrderedDict{Symbol, Any}()
    dyn_model_dict[:meta] = OrderedDict{Symbol, Any}(
        :mech_power_mode => mech_power_mode,
        :constrain_δ => constrain_δ,
        :bound_style_δ => bound_style_δ,
        :δ_ref_gen_id => δ_ref_gen_id,
        # Read only by `bound_style_δ = :pairwise`; carried for every style so the
        # details export can report the band without knowing which one is active.
        :δ_pair_tol => δ_pair_tol_tuple(simulation),
        :constrain_Δω => constrain_Δω,
        :bound_style_Δω => bound_style_Δω,
        :Δω_tol => Δω_tol,
        :network_form => "FULL_BUS",
        :gen_order => "DQ_4TH",
        :load_table_rows => DLOAD === nothing ? 0 : nrow(DLOAD),
        :coupling_init_source => "acopf_warmstart",   # the only FULL_BUS coupling path
        :dq_speed_dev_in_algebra => dq_speed_dev_in_algebra,
        :include_avr => include_avr,
        :include_governor => include_governor,
        :governor_limiter => string(governor_limiter),
        :ode_first_step => ode_first_step,
        :gfm_integrator => gfm_integrator,
        # Read by `_coi_gens` during assembly, so it must be in this literal:
        # `register_dyn_model_meta!` only runs once the model is already built.
        :gfm_in_coi => gfm_in_coi,
        # Read by `_attach_gfm_window!` while the model is built, so it has to be
        # in this literal: `register_dyn_model_meta!` only runs once the model is
        # already built. Same reason for the SG-presence guards below.
        :gfm_pq_overload => gfm_pq_overload,
        :allow_all_gfm => allow_all_gfm,
        # System base, read by `_gfm_pq_eps` to put the overload-clip ε back on the
        # converter's own base (the DGFM gains are already system base by then).
        :base_MVA => base_MVA,
        :ineq_cons => ts_input_param[:ineq_cons],
        :var_bounds => ts_input_param[:var_bounds],
        :bound_encoding => bound_encoding_from_param(ts_input_param),
    )

    val_V = steady_state_hints.val_V
    val_θ = steady_state_hints.val_θ
    val_Pg = steady_state_hints.val_Pg
    val_Qg = steady_state_hints.val_Qg

    DBUS_mod = deepcopy(DBUS)
    DGEN_mod = deepcopy(DGEN)
    # Row g = generator id g, so the positional `DGEN_DYN.<col>[gen]` reads below stay
    # correct when a converter sits mid-fleet (see `dgen_dyn_by_id`).
    DGEN_DYN_mod = dgen_dyn_by_id(DGEN_DYN, DGEN)
    DCIR_mod = deepcopy(DCIR)

    δ_tol, f_syn, ω_syn, Δω_0 = common_ts_parameters(simulation)
    active_gen = findall(x -> x == 1, DGEN_mod.g_status)
    bus_gen_circ_on = bus_gen_circ_dict_ON
    register_gfm_meta!(dyn_model_dict, DGFM)

    dyn_model_dict[:meta][:fault_type] = ts_fault_details[:fault_type]
    if ts_fault_details[:fault_type] == "GL"
        dyn_model_dict[:meta][:gl_element] = ts_fault_details[:gl][:element_2_disconnect]
    elseif ts_fault_details[:fault_type] == "OB"
        dyn_model_dict[:meta][:ob_branch_ids] = copy(ts_fault_details[:ob][:branch_id])
    end

    if ts_fault_details[:fault_type] == "SC"
        t_start_sim, t_end_sim, t_step, t_start_fault, clearing_time, t_clear_fault,
            t_window_fault, t_window_postf, t_window_total = time_windows_sc(simulation)

        Ybus_pref = Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA)
        sc_bus = ts_fault_details[:sc][:bus]
        bus_fault_id = sc_bus[:bus_id]
        # y_f and the stamp come from FaultConfig. The default adds the inductive
        # solid short; the resistive option overwrites with the real FAULT_BUS_SHUNT.
        Ybus_fault = Calculate_Ybus_fullbus_dynamics(Ybus_pref;
            bus_fault=bus_fault_id,
            y_fault=sc_bus[:admittance_pu],
            stamp=sc_bus[:stamp])
        # Seed the fault-window V and θ at the shorted bus from its current balance
        # under the stamped shunt. The angle follows y_f (near the pre-fault phasor
        # for the default inductive short). See fault_bus_voltage_seed.
        val_V_fault, val_θ_fault = fault_bus_voltage_seed(Ybus_fault, val_V, val_θ, bus_fault_id)
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
            "Ybus_fault", Ybus_fault)

        model, dyn_model_dict = Define_Initial_Condition_4_dq!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode,
            include_avr=include_avr,
            DGFM=DGFM)

        model, dyn_model_dict = Define_Fault_Dynamic_Model_dq!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V_fault, val_θ_fault, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            dq_speed_dev_in_algebra=dq_speed_dev_in_algebra,
            include_avr=include_avr,
            include_governor=include_governor,
            governor_limiter=governor_limiter,
            DGFM=DGFM)

        if ts_fault_details[:sc][:bus][:disconnect_branch]
            DCIR_mod.l_status[ts_fault_details[:sc][:bus][:branch_id_2_disconnect]] .= 0
            Ybus_postf = Calculate_Ybus_fullbus_dynamics(
                Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
            Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
                "Ybus_postfault", Ybus_postf)

            last_var_Pe_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pe_tf])
            last_var_δ_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:δ_tf])
            last_var_Δω_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Δω_tf])

            model, dyn_model_dict = Define_PostFault_Dynamic_Model_dq!(
                model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
                active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_postf, Ybus_postf,
                opf_dict[:vars][:V], opf_dict[:vars][:θ], last_var_Pe_tf, last_var_δ_tf,
                last_var_Δω_tf, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
                bus_ids=collect(Int, DBUS_mod.bus),
                dq_speed_dev_in_algebra=dq_speed_dev_in_algebra,
                include_avr=include_avr,
                include_governor=include_governor,
                governor_limiter=governor_limiter,
                DGFM=DGFM)
        else
            throw(ArgumentError("DQ FULL_BUS SC requires branch trip after clearing."))
        end

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}(
            :t_start_sim => t_start_sim, :t_end_sim => t_end_sim, :t_step => t_step,
            :t_start_fault => t_start_fault, :clearing_time => clearing_time,
            :t_clear_fault => t_clear_fault, :t_window_fault => t_window_fault,
            :t_window_postf => t_window_postf, :t_window_total => t_window_total,
        )

    elseif ts_fault_details[:fault_type] == "GL"
        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        if ts_fault_details[:gl][:element_2_disconnect] == "gen"
            gen_ids = ts_fault_details[:gl][:gen][:gen_id]
            for gen_id in gen_ids
                DGEN_mod.g_status[gen_id] = 0
            end
            ts_fault_details[:gl][:gen][:bus_id] = [Int(DGEN_mod.bus[g]) for g in gen_ids]
            active_gen = findall(x -> x == 1, DGEN_mod.g_status)
            _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)
        elseif ts_fault_details[:gl][:element_2_disconnect] == "load"
            bus_ids = ts_fault_details[:gl][:load][:bus_id]
            apply_gl_load_scaling!(
                DBUS_mod, bus_ids, ts_fault_details[:gl][:load][:percent_power])
            _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)
        else
            throw(ArgumentError("Unknown GL element: $(ts_fault_details[:gl][:element_2_disconnect])"))
        end

        Ybus_fault = Calculate_Ybus_fullbus_dynamics(
            Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
            "Ybus_fault", Ybus_fault)

        model, dyn_model_dict = Define_Initial_Condition_4_dq!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode,
            include_avr=include_avr,
            DGFM=DGFM)

        model, dyn_model_dict = Define_Fault_Dynamic_Model_dq!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            dq_speed_dev_in_algebra=dq_speed_dev_in_algebra,
            include_avr=include_avr,
            include_governor=include_governor,
            governor_limiter=governor_limiter,
            DGFM=DGFM)

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}(
            :t_start_sim => t_start_sim, :t_end_sim => t_end_sim, :t_step => t_step,
            :t_start_fault => t_start_fault, :t_window_fault => t_window_fault,
            :t_window_total => t_window_total,
        )
    elseif ts_fault_details[:fault_type] == "OB"
        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        branch_ids = ts_fault_details[:ob][:branch_id]
        DCIR_mod.l_status[branch_ids] .= 0
        _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)

        Ybus_fault = Calculate_Ybus_fullbus_dynamics(
            Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
            "Ybus_fault", Ybus_fault)

        model, dyn_model_dict = Define_Initial_Condition_4_dq!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode,
            include_avr=include_avr,
            DGFM=DGFM)

        model, dyn_model_dict = Define_Fault_Dynamic_Model_dq!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            dq_speed_dev_in_algebra=dq_speed_dev_in_algebra,
            include_avr=include_avr,
            include_governor=include_governor,
            governor_limiter=governor_limiter,
            DGFM=DGFM)

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}(
            :t_start_sim => t_start_sim, :t_end_sim => t_end_sim, :t_step => t_step,
            :t_start_fault => t_start_fault, :t_window_fault => t_window_fault,
            :t_window_total => t_window_total,
        )
    else
        throw(ArgumentError("Unknown fault_type: $(ts_fault_details[:fault_type])"))
    end

    dyn_model_dict[:active_gen] = active_gen
    dyn_model_dict[:mech_power_mode] = mech_power_mode
    dyn_parameters_dict[:common] = OrderedDict{Symbol, Any}(
        :δ_tol => δ_tol, :δ_pair_tol => δ_pair_tol_tuple(simulation),
        :f_syn => f_syn, :ω_syn => ω_syn, :Δω_0 => Δω_0,
        :constrain_Δω => constrain_Δω,
    )
    if constrain_Δω
        dyn_parameters_dict[:common][:Δω_tol] = Δω_tol
    end

    build_dual_registry!(dyn_model_dict)
    return model, dyn_model_dict, dyn_parameters_dict
end

"""
Pre-fault dq initial conditions at t = 0.

`E_fd` and `δ` reuse the classical FULL_BUS bound builders (`var_tsred_*`, `var_kron_*`)
because `lims[:E]` and `lims[:δ]` are per-generator vectors, not scalars.
`Ed/Eq/Id/Iq` are free algebraic states pinned by `eq_const_dq_init_steady_state!`.

When `DGFM` is present, SG machines use this path on `sg_gens` only; GFM units are
attached via `Attach_GFM_init!` (Phase G1).
"""
function Define_Initial_Condition_4_dq!(
    model::JuMP.Model,
    ts_input_param::OrderedDict{Symbol, Any},
    dyn_model_dict::OrderedDict{Symbol, Any},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64},
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef},
    val_V::Dict,
    val_θ::Dict,
    val_Pg::Dict,
    val_Qg::Dict;
    mech_power_mode::MechPowerMode=USE_PM,
    include_avr::Bool=false,
    DGFM::Union{Nothing, DataFrame}=nothing,
)
    sg_gens = Vector{Int64}(sg_active_gens(active_gen, DGFM))
    gfm_gens = Vector{Int64}(gfm_active_gens(active_gen, DGFM))
    # Converter-only fleets are opt-in; see `resolve_sg_gfm_gens!` and
    # `validate_all_gfm_fleet!`, which screens the unsupported corridor/control combos
    # before any model is built.
    if isempty(sg_gens) && !get(dyn_model_dict[:meta], :allow_all_gfm, false)
        throw(ArgumentError(
            "DQ FULL_BUS init requires at least one active SG (sg_gens empty; all-GFM " *
            "fleets unsupported); set DynModelConfig.allow_all_gfm=true to permit a " *
            "converter-only fleet."))
    end
    dyn_model_dict[:meta][:sg_gens] = sg_gens
    dyn_model_dict[:meta][:gfm_gens] = gfm_gens

    register_ts_bound_limit_specs!(
        ts_input_param, dyn_model_dict, ts_input_param[:limits],
        Vector{Int}(sg_gens), DGEN, base_MVA;
        DGFM = DGFM, gfm_gens = gfm_gens)
    lims = ts_input_param[:var_limit_specs]
    E_min, E_max = lims[:E]
    δ_min, δ_max = lims[:δ]

    # Field voltage E_fd: constant in milestone 1 (AVR would make it a state later).
    begin
        enc = bound_encoding_from_meta(dyn_model_dict[:meta])
        result = var_tsred_gen_voltage_magnitude!(
            model, sg_gens, "E_fd";
            bounded=ts_input_param[:var_bounds][:E],
            min_lim=E_min, max_lim=E_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_E_fd_lower, export_key_upper=:ineq_const_E_fd_upper)
        E_fd = store_ts_scalar_var_bounds!(dyn_model_dict, :E_fd, result, :ineq_const_E_fd_lower, :ineq_const_E_fd_upper)
    end

    begin
        result = var_kron_gen_rotor_angle!(
            model, sg_gens, ts_input_param[:var_names][:δ];
            bounded=ts_input_param[:var_bounds][:δ], min_lim=δ_min, max_lim=δ_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_δ_lower, export_key_upper=:ineq_const_δ_upper)
        δ = store_ts_scalar_var_bounds!(dyn_model_dict, :δ, result, :ineq_const_δ_lower, :ineq_const_δ_upper)
    end

    (; Ed, Eq, Id, Iq) = var_dq_prefault_algebraic!(model, sg_gens)
    dyn_model_dict[:vars][:Ed] = Ed
    dyn_model_dict[:vars][:Eq] = Eq
    dyn_model_dict[:vars][:Id] = Id
    dyn_model_dict[:vars][:Iq] = Iq

    P_m = nothing
    if mech_power_mode == USE_PM
        Pm_min, Pm_max = lims[:P_m]
        result = var_kron_gen_mech_power!(
            model, sg_gens, ts_input_param[:var_names][:P_m];
            bounded=ts_input_param[:var_bounds][:P_m], min_lim=Pm_min, max_lim=Pm_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_P_m_lower, export_key_upper=:ineq_const_P_m_upper)
        P_m = store_ts_scalar_var_bounds!(dyn_model_dict, :P_m, result, :ineq_const_P_m_lower, :ineq_const_P_m_upper)
    end

    _set_dq_init_warm_starts!(E_fd, δ, Ed, Eq, Id, Iq, P_m, sg_gens, DGEN, DGEN_DYN,
        val_V, val_θ, val_Pg, val_Qg)

    (; eq_Ed, eq_Eq, eq_Vd, eq_Vq, eq_P, eq_Q, Vd_init, Vq_init) =
        eq_const_dq_init_steady_state!(
            model, sg_gens, DGEN, DGEN_DYN, V, θ, P_g, Q_g, E_fd, δ, Ed, Eq, Id, Iq)
    dyn_model_dict[:eq_const][:eq_const_Ed_init] = eq_Ed
    dyn_model_dict[:eq_const][:eq_const_Eq_init] = eq_Eq
    dyn_model_dict[:eq_const][:eq_const_Vd_init] = eq_Vd
    dyn_model_dict[:eq_const][:eq_const_Vq_init] = eq_Vq
    dyn_model_dict[:eq_const][:eq_const_P_init] = eq_P
    dyn_model_dict[:eq_const][:eq_const_Q_init] = eq_Q
    haskey(dyn_model_dict, :expressions) ||
        (dyn_model_dict[:expressions] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:expressions][:Vd_init] = Vd_init
    dyn_model_dict[:expressions][:Vq_init] = Vq_init

    if mech_power_mode == USE_PM
        dyn_model_dict[:eq_const][:eq_const_Pm_init] = eq_const_kron_initial_mechanical_power!(
            model, dyn_model_dict[:vars][:P_m], P_g, sg_gens)
    end

    register_mech_power_refs!(dyn_model_dict, P_g, mech_power_mode)

    if include_avr
        Attach_Avr_init!(model, dyn_model_dict, sg_gens, DGEN, DGEN_DYN, V, val_V)
    end

    if !isempty(gfm_gens)
        DGFM === nothing && throw(ArgumentError("gfm_gens non-empty but DGFM is nothing."))
        Attach_GFM_init!(
            model, dyn_model_dict, gfm_gens, DGEN, DGFM,
            V, θ, P_g, Q_g, val_V, val_θ, val_Pg, val_Qg)
    end

    attach_dq_prefault_algebraic_bounds!(model, dyn_model_dict)
    attach_avr_prefault_bounds!(model, dyn_model_dict)

    return model, dyn_model_dict
end

function Define_Fault_Dynamic_Model_dq!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    DBUS::DataFrame,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64},
    nBUS::Int64,
    base_MVA::Float64,
    ZIP_P::AbstractMatrix{Float64},
    ZIP_Q::AbstractMatrix{Float64},
    Δt::Float64,
    time_window::Vector{Float64},
    Ybus::SparseMatrixCSC,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Δω_0::Float64,
    ω_syn::Float64,
    δ_tol::Tuple{Float64, Float64},
    val_V::Dict,
    val_θ::Dict,
    val_Pg::Dict,
    val_Qg::Dict;
    bus_ids::Vector{Int},
    dq_speed_dev_in_algebra::Bool=true,
    include_avr::Bool=false,
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
    DGFM::Union{Nothing, DataFrame}=nothing,
)
    sg_gens, gfm_gens = resolve_sg_gfm_gens!(dyn_model_dict, active_gen, DGFM)

    E_fd = dyn_model_dict[:vars][:E_fd]
    δ_0 = dyn_model_dict[:vars][:δ]
    Ed_0 = dyn_model_dict[:vars][:Ed]
    Eq_0 = dyn_model_dict[:vars][:Eq]
    Id_0 = dyn_model_dict[:vars][:Id]
    Iq_0 = dyn_model_dict[:vars][:Iq]
    P_mech = dyn_model_dict[:refs][:P_mech]

    V_tf = var_fullbus_bus_voltage_magnitude_time!(model, bus_ids, time_window, val_V; suffix="tf")
    θ_tf = var_fullbus_bus_angle_time!(model, bus_ids, time_window, val_θ; suffix="tf")
    Pe_tf = var_fullbus_gen_Pe_time!(model, active_gen, time_window, val_Pg; suffix="tf")
    Qe_tf = var_fullbus_gen_Qe_time!(model, active_gen, time_window, val_Qg; suffix="tf")
    δ_tf = var_fullbus_gen_rotor_angle_time!(model, active_gen, time_window, δ_0; suffix="tf")
    Δω_tf = var_fullbus_gen_speed_dev_time!(model, active_gen, time_window; suffix="tf")
    dq_tf = var_dq_gen_state_time!(
        model, active_gen, time_window, δ_0, Ed_0, Eq_0, Id_0, Iq_0, val_Pg;
        suffix="tf", emf_gens=sg_gens)
    Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf = dq_tf.Ed, dq_tf.Eq, dq_tf.Id, dq_tf.Iq, dq_tf.Te

    dyn_model_dict[:vars][:V_tf] = V_tf
    dyn_model_dict[:vars][:θ_tf] = θ_tf
    V_min = dyn_model_dict[:meta][:var_limit_specs][:V_bus_min][1]
    attach_fullbus_V_lower_bounds!(model, dyn_model_dict, V_tf, V_min, :ineq_const_V_tf_lower)
    dyn_model_dict[:vars][:Pe_tf] = Pe_tf
    dyn_model_dict[:vars][:Qe_tf] = Qe_tf
    dyn_model_dict[:vars][:δ_tf] = δ_tf
    dyn_model_dict[:vars][:Δω_tf] = Δω_tf
    dyn_model_dict[:vars][:Ed_tf] = Ed_tf
    dyn_model_dict[:vars][:Eq_tf] = Eq_tf
    dyn_model_dict[:vars][:Id_tf] = Id_tf
    dyn_model_dict[:vars][:Iq_tf] = Iq_tf
    dyn_model_dict[:vars][:Te_tf] = Te_tf

    # The corridors receive the full active set and narrow it themselves: the
    # machine-referenced δ styles and the absolute Δω style span the converters too, while
    # every inertia-weighted quantity is rebuilt over the SGs alone (`_sg_ids`).
    δCOI_tf = _attach_δCOI!(model, dyn_model_dict, δ_tf, active_gen, DGEN_DYN,
        time_window, :tf)
    attach_fault_tf_var_bounds!(model, dyn_model_dict)
    _add_δ_bounds_fault!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tf, δCOI_tf,
        Δω_tf, Pe_tf, time_window, δ_tol, δ_0, Δω_0, ω_syn, Δt)

    # Optional speed corridor (reference set by `bound_style_Δω`).
    _add_Δω_bounds!(model, dyn_model_dict, :tf, active_gen, DGEN_DYN, Δω_tf, time_window)

    alg_tf = eq_const_dq_machine_algebra!(
        model, sg_gens, DGEN, DGEN_DYN, Pe_tf, Qe_tf, δ_tf, Δω_tf,
        Ed_tf, Eq_tf, Id_tf, Iq_tf, Te_tf, V_tf, θ_tf, time_window;
        speed_in_algebra=dq_speed_dev_in_algebra)
    dyn_model_dict[:eq_const][:eq_const_Pe_tf] = alg_tf.eq_Pe
    dyn_model_dict[:eq_const][:eq_const_Qe_tf] = alg_tf.eq_Qe
    dyn_model_dict[:eq_const][:eq_const_Vd_tf] = alg_tf.eq_Vd
    dyn_model_dict[:eq_const][:eq_const_Vq_tf] = alg_tf.eq_Vq
    dyn_model_dict[:eq_const][:eq_const_Te_tf] = alg_tf.eq_Te

    # Nodal admittance terms from Ybus; dq currents provide the generator injection.
    terms_Pb, terms_Qb = _build_fullbus_nodal_injection_terms!(
        model, nBUS, bus_gen_circ_dict, V_tf, θ_tf, time_window, Ybus)
    haskey(dyn_model_dict, :expressions) ||
        (dyn_model_dict[:expressions] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:expressions][:Vd_tf] = alg_tf.Vd_out
    dyn_model_dict[:expressions][:Vq_tf] = alg_tf.Vq_out
    dyn_model_dict[:expressions][:P_inj_tf] = terms_Pb
    dyn_model_dict[:expressions][:Q_inj_tf] = terms_Qb
    dyn_model_dict[:eq_const][:eq_const_Pbalance_tf] = eq_const_dq_Pbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, δ_tf, Id_tf, Iq_tf,
        terms_Pb, V_tf, θ_tf, time_window, ZIP_P)
    dyn_model_dict[:eq_const][:eq_const_Qbalance_tf] = eq_const_dq_Qbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, δ_tf, Id_tf, Iq_tf,
        terms_Qb, V_tf, θ_tf, time_window, ZIP_Q)

    ode_fs = dyn_model_dict[:meta][:ode_first_step]
    dyn_model_dict[:eq_const][:eq_const_δ_tf] = eq_const_fullbus_δ_swing_fault!(
        model, sg_gens, δ_tf, Δω_tf, δ_0, Δω_0, time_window, ω_syn, Δt;
        ode_first_step=ode_fs)
    if include_governor
        P_m0 = dyn_model_dict[:vars][:P_m]
        Pm_tf = Attach_Governor_fault!(
            model, dyn_model_dict, sg_gens, DGEN, DGEN_DYN, base_MVA,
            Δω_tf, time_window, Δt, governor_limiter)
        dyn_model_dict[:eq_const][:eq_const_Δω_tf] = eq_const_fullbus_Δω_swing_fault!(
            model, sg_gens, DGEN_DYN, Pm_tf, P_m0, Pe_tf, Δω_tf, Δω_0, P_g,
            time_window, Δt; ode_first_step=ode_fs)
    else
        dyn_model_dict[:eq_const][:eq_const_Δω_tf] = eq_const_fullbus_Δω_swing_fault!(
            model, sg_gens, DGEN_DYN, P_mech, Pe_tf, Δω_tf, Δω_0, P_g, time_window, Δt;
            ode_first_step=ode_fs)
    end

    E_fd_tf = nothing
    if include_avr
        E_fd_tf = Attach_Avr_fault!(
            model, dyn_model_dict, sg_gens, DGEN, DGEN_DYN, V, V_tf, time_window, Δt)
    end

    emf_kw = if include_avr
        (E_fd_time=E_fd_tf, E_fd_prev0=E_fd,
         ode_first_step=dyn_model_dict[:meta][:ode_first_step])
    else
        (ode_first_step=dyn_model_dict[:meta][:ode_first_step],)
    end
    dyn_model_dict[:eq_const][:eq_const_Ed_tf], dyn_model_dict[:eq_const][:eq_const_Eq_tf] =
        eq_const_dq_emf_dynamics!(
            model, sg_gens, DGEN_DYN, E_fd, Ed_0, Eq_0, Id_0, Iq_0,
            Ed_tf, Eq_tf, Id_tf, Iq_tf, time_window, Δt; emf_kw...)

    if !isempty(gfm_gens)
        DGFM === nothing && throw(ArgumentError("gfm_gens non-empty but DGFM is nothing."))
        Attach_GFM_fault!(
            model, dyn_model_dict, gfm_gens, DGEN, DGFM, Δt, ω_syn, time_window)
    end

    return model, dyn_model_dict
end

function Define_PostFault_Dynamic_Model_dq!(
    model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    DBUS::DataFrame,
    bus_gen_circ_dict::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64},
    nBUS::Int64,
    base_MVA::Float64,
    ZIP_P::AbstractMatrix{Float64},
    ZIP_Q::AbstractMatrix{Float64},
    Δt::Float64,
    time_window::Vector{Float64},
    Ybus::SparseMatrixCSC,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    Pe_ant::OrderedDict{Int64, JuMP.VariableRef},
    δ_ant::OrderedDict{Int64, JuMP.VariableRef},
    Δω_ant::OrderedDict{Int64, JuMP.VariableRef},
    ω_syn::Float64,
    δ_tol::Tuple{Float64, Float64},
    val_V::Dict,
    val_θ::Dict,
    val_Pg::Dict,
    val_Qg::Dict;
    bus_ids::Vector{Int},
    dq_speed_dev_in_algebra::Bool=true,
    include_avr::Bool=false,
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
    DGFM::Union{Nothing, DataFrame}=nothing,
)
    sg_gens, gfm_gens = resolve_sg_gfm_gens!(dyn_model_dict, active_gen, DGFM)

    E_fd = dyn_model_dict[:vars][:E_fd]
    P_mech = dyn_model_dict[:refs][:P_mech]
    Pe_tf = dyn_model_dict[:vars][:Pe_tf]
    Ed_tf = dyn_model_dict[:vars][:Ed_tf]
    Eq_tf = dyn_model_dict[:vars][:Eq_tf]
    Id_tf = dyn_model_dict[:vars][:Id_tf]
    Iq_tf = dyn_model_dict[:vars][:Iq_tf]

    Ed_ant = OrderedDict(g => last(inner).second for (g, inner) in Ed_tf)
    Eq_ant = OrderedDict(g => last(inner).second for (g, inner) in Eq_tf)
    Id_ant = OrderedDict(g => last(inner).second for (g, inner) in Id_tf)
    Iq_ant = OrderedDict(g => last(inner).second for (g, inner) in Iq_tf)

    V_tpf = var_fullbus_bus_voltage_magnitude_time!(model, bus_ids, time_window, val_V; suffix="tpf")
    θ_tpf = var_fullbus_bus_angle_time!(model, bus_ids, time_window, val_θ; suffix="tpf")
    Pe_tpf = var_fullbus_gen_Pe_time!(model, active_gen, time_window, val_Pg; suffix="tpf")
    Qe_tpf = var_fullbus_gen_Qe_time!(model, active_gen, time_window, val_Qg; suffix="tpf")
    δ_tpf = var_fullbus_gen_rotor_angle_time!(
        model, active_gen, time_window, dyn_model_dict[:vars][:δ]; suffix="tpf")
    Δω_tpf = var_fullbus_gen_speed_dev_time!(model, active_gen, time_window; suffix="tpf")
    dq_tpf = var_dq_gen_state_time!(
        model, active_gen, time_window, dyn_model_dict[:vars][:δ],
        dyn_model_dict[:vars][:Ed], dyn_model_dict[:vars][:Eq],
        dyn_model_dict[:vars][:Id], dyn_model_dict[:vars][:Iq], val_Pg;
        suffix="tpf", emf_gens=sg_gens)
    Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf = dq_tpf.Ed, dq_tpf.Eq, dq_tpf.Id, dq_tpf.Iq, dq_tpf.Te

    dyn_model_dict[:vars][:V_tpf] = V_tpf
    dyn_model_dict[:vars][:θ_tpf] = θ_tpf
    V_min = dyn_model_dict[:meta][:var_limit_specs][:V_bus_min][1]
    attach_fullbus_V_lower_bounds!(model, dyn_model_dict, V_tpf, V_min, :ineq_const_V_tpf_lower)
    dyn_model_dict[:vars][:Pe_tpf] = Pe_tpf
    dyn_model_dict[:vars][:Qe_tpf] = Qe_tpf
    dyn_model_dict[:vars][:δ_tpf] = δ_tpf
    dyn_model_dict[:vars][:Δω_tpf] = Δω_tpf
    dyn_model_dict[:vars][:Ed_tpf] = Ed_tpf
    dyn_model_dict[:vars][:Eq_tpf] = Eq_tpf
    dyn_model_dict[:vars][:Id_tpf] = Id_tpf
    dyn_model_dict[:vars][:Iq_tpf] = Iq_tpf
    dyn_model_dict[:vars][:Te_tpf] = Te_tpf

    # As in the fault window: full active set in, narrowing decided by the corridor style.
    δCOI_tpf = _attach_δCOI!(model, dyn_model_dict, δ_tpf, active_gen, DGEN_DYN,
        time_window, :tpf)
    attach_postfault_tpf_var_bounds!(model, dyn_model_dict)
    _add_δ_bounds_postf!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tpf, δCOI_tpf,
        Δω_tpf, Pe_tpf, time_window, δ_tol, δ_ant, Δω_ant, Pe_ant, ω_syn, Δt)

    _add_Δω_bounds!(model, dyn_model_dict, :tpf, active_gen, DGEN_DYN, Δω_tpf, time_window)

    alg_tpf = eq_const_dq_machine_algebra!(
        model, sg_gens, DGEN, DGEN_DYN, Pe_tpf, Qe_tpf, δ_tpf, Δω_tpf,
        Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, Te_tpf, V_tpf, θ_tpf, time_window;
        speed_in_algebra=dq_speed_dev_in_algebra)
    dyn_model_dict[:eq_const][:eq_const_Pe_tpf] = alg_tpf.eq_Pe
    dyn_model_dict[:eq_const][:eq_const_Qe_tpf] = alg_tpf.eq_Qe
    dyn_model_dict[:eq_const][:eq_const_Vd_tpf] = alg_tpf.eq_Vd
    dyn_model_dict[:eq_const][:eq_const_Vq_tpf] = alg_tpf.eq_Vq
    dyn_model_dict[:eq_const][:eq_const_Te_tpf] = alg_tpf.eq_Te
    dyn_model_dict[:expressions][:Vd_tpf] = alg_tpf.Vd_out
    dyn_model_dict[:expressions][:Vq_tpf] = alg_tpf.Vq_out

    terms_Pb, terms_Qb = _build_fullbus_nodal_injection_terms!(
        model, nBUS, bus_gen_circ_dict, V_tpf, θ_tpf, time_window, Ybus)
    dyn_model_dict[:expressions][:P_inj_tpf] = terms_Pb
    dyn_model_dict[:expressions][:Q_inj_tpf] = terms_Qb
    dyn_model_dict[:eq_const][:eq_const_Pbalance_tpf] = eq_const_dq_Pbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, δ_tpf, Id_tpf, Iq_tpf,
        terms_Pb, V_tpf, θ_tpf, time_window, ZIP_P)
    dyn_model_dict[:eq_const][:eq_const_Qbalance_tpf] = eq_const_dq_Qbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, δ_tpf, Id_tpf, Iq_tpf,
        terms_Qb, V_tpf, θ_tpf, time_window, ZIP_Q)

    ode_fs = dyn_model_dict[:meta][:ode_first_step]
    dyn_model_dict[:eq_const][:eq_const_δ_tpf] = eq_const_fullbus_δ_swing_postf!(
        model, sg_gens, δ_tpf, Δω_tpf, δ_ant, Δω_ant, time_window, ω_syn, Δt;
        ode_first_step=ode_fs)
    if include_governor
        Pm_tf_last = OrderedDict(
            g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pm_tf])
        Pm_tpf = Attach_Governor_postf!(
            model, dyn_model_dict, sg_gens, DGEN, DGEN_DYN, base_MVA,
            Δω_tpf, time_window, Δt, governor_limiter)
        dyn_model_dict[:eq_const][:eq_const_Δω_tpf] = eq_const_fullbus_Δω_swing_postf!(
            model, sg_gens, DGEN_DYN, Pm_tpf, Pm_tf_last, Pe_tpf, Δω_tpf, Δω_ant, Pe_ant,
            Pe_tf, time_window, Δt; ode_first_step=ode_fs)
    else
        dyn_model_dict[:eq_const][:eq_const_Δω_tpf] = eq_const_fullbus_Δω_swing_postf!(
            model, sg_gens, DGEN_DYN, P_mech, Pe_tpf, Δω_tpf, Δω_ant, Pe_ant,
            Pe_tf, time_window, Δt; ode_first_step=ode_fs)
    end

    E_fd_tpf = nothing
    if include_avr
        E_fd_tpf = Attach_Avr_postfault!(
            model, dyn_model_dict, sg_gens, DGEN, DGEN_DYN,
            dyn_model_dict[:vars][:V_tf], V_tpf, time_window, Δt)
        E_fd_last = OrderedDict(
            g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:E_fd_tf])
        emf_kw = (E_fd_time=E_fd_tpf, E_fd_prev0=E_fd_last,
                  ode_first_step=dyn_model_dict[:meta][:ode_first_step])
    else
        emf_kw = (ode_first_step=dyn_model_dict[:meta][:ode_first_step],)
    end
    dyn_model_dict[:eq_const][:eq_const_Ed_tpf], dyn_model_dict[:eq_const][:eq_const_Eq_tpf] =
        eq_const_dq_emf_dynamics!(
            model, sg_gens, DGEN_DYN, E_fd, Ed_ant, Eq_ant, Id_ant, Iq_ant,
            Ed_tpf, Eq_tpf, Id_tpf, Iq_tpf, time_window, Δt; emf_kw...)

    if !isempty(gfm_gens)
        DGFM === nothing && throw(ArgumentError("gfm_gens non-empty but DGFM is nothing."))
        Attach_GFM_postf!(
            model, dyn_model_dict, gfm_gens, DGEN, DGFM, Δt, ω_syn, time_window)
    end

    return model, dyn_model_dict
end
