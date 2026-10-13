#=
================================================================================
 functions_2_build_TS_model_w_FullBus.jl  —  TSC-ACOPF full-network classical 2nd
================================================================================
 Phase 3: classical generator on the full sparse Ybus (no Kron reduction).
 Physics aligned with the reference full-network TSC-ACOPF formulation.

 Architecture mirrors the Kron path:
 - `build_ts_input_param(ts_builder)` + explicit `@constraint` bounds when `var_bounds[k]=true`
 - One builder per variable / constraint family (`var_fullbus_*`, `eq_const_fullbus_*`)
 - `bound_style_δ` dispatches δ-COI stability; `constrain_Δω` gates Δω-COI bounds
 - Mandatory ACOPF warm start (`SteadyStateHints`); model TXT export before `optimize!`
================================================================================
=#

# ==================================================================================
# Top-level orchestrator — SC bus fault + GL gen/load + OB open-branch
# ==================================================================================

"""
Append the FULL_BUS classical transient-stability sub-model onto an existing ACOPF `model`.

Three fault families are handled, both producing per-time-step JuMP variables/constraints
that are stored (by family) in the returned `dyn_model_dict`:

- `"SC"` short-circuit bus fault: pre-fault initial conditions → fault-on window (faulted
  Ybus) → post-fault window (faulted branch tripped). A branch trip after clearing is
  mandatory here.
- `"GL"` generation/load loss: pre-fault initial conditions → a single perturbed window
  after the generator(s) or load are disconnected/scaled. There is no post-fault window.
- `"OB"` open branch: same single-window timeline as GL, but the disturbance is opening
  one or more in-service lines/transformers (`DCIR.l_status → 0`) with no short-circuit.

Key choices:
- `mech_power_mode` must be `USE_PM` (explicit mechanical-power variable P_m).
- `bound_style_δ` selects the δ-COI stability bound flavour; `constrain_Δω` optionally
  adds speed-deviation COI bounds.
- `zip_load_p` / `zip_load_q` = independent (Z, I, P) load splits for the active and the
  reactive nodal balance (a TSO may model P as constant current and Q as constant admittance).
- `DLOAD` = optional per-load table (`load_data.csv`). When given, each bus uses the
  demand-weighted mean of the shares of the loads sitting on it instead of the system-wide
  split; see `zip_by_bus`.
- `steady_state_hints` carries the solved ACOPF point (V, θ, P_g, Q_g) used for warm starts.

Returns `(model, dyn_model_dict, dyn_parameters_dict)`.
"""
function Make_Dynamic_Model_fullbus!(
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
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
    ode_first_step::Symbol=:trapezoidal,
    steady_state_hints::SteadyStateHints,
    DLOAD::Union{Nothing, DataFrame}=nothing,
)
    mech_power_mode == USE_PM || throw(ArgumentError(
        "FULL_BUS classical path currently requires mech_power_mode=USE_PM."))
    include_governor && mech_power_mode != USE_PM && throw(ArgumentError(
        "include_governor requires mech_power_mode=USE_PM."))
    ode_first_step ∈ (:trapezoidal, :backward_euler) || throw(ArgumentError(
        "ode_first_step must be :trapezoidal or :backward_euler (got $ode_first_step)."))

    # `ts_input_param` holds the per-family toggles (var names, which bounds/eq-constraints
    # to emit) decoded from the builder config. ZIP_P / ZIP_Q are `nBUS × 3` (Z, I, P)
    # matrices — one row per bus — for the active and reactive balances respectively.
    # Without `DLOAD` every row is the system-wide split; with it, a bus carrying several
    # loads gets their demand-weighted mean (exact, see `zip_by_bus`).
    ts_input_param = build_ts_input_param(ts_builder)
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
        :load_table_rows => DLOAD === nothing ? 0 : nrow(DLOAD),
        :coupling_init_source => "acopf_warmstart",   # the only FULL_BUS coupling path
        :ode_first_step => ode_first_step,
        :ineq_cons => ts_input_param[:ineq_cons],
        :var_bounds => ts_input_param[:var_bounds],
        :bound_encoding => bound_encoding_from_param(ts_input_param),
    )

    # Solved ACOPF operating point — used as warm starts and as the steady-state reference
    # voltage V[i] in the ZIP load scaling of the nodal balances.
    val_V  = steady_state_hints.val_V
    val_θ  = steady_state_hints.val_θ
    val_Pg = steady_state_hints.val_Pg
    val_Qg = steady_state_hints.val_Qg

    # Work on private copies: GL faults mutate generator status / bus loads, OB/SC
    # trip branches (l_status). The caller's original DataFrames must stay intact.
    DBUS_mod = deepcopy(DBUS)
    DGEN_mod = deepcopy(DGEN)
    DGEN_DYN_mod = deepcopy(DGEN_DYN)
    DCIR_mod = deepcopy(DCIR)

    # δ_tol = (lower, upper) rotor-angle COI tolerance; ω_syn synchronous speed [rad/s];
    # Δω_0 initial speed deviation (0 at equilibrium).
    δ_tol, f_syn, ω_syn, Δω_0 = common_ts_parameters(simulation)

    active_gen = findall(x -> x == 1, DGEN_mod.g_status)
    bus_gen_circ_on = bus_gen_circ_dict_ON

    dyn_model_dict[:meta][:fault_type] = ts_fault_details[:fault_type]
    if ts_fault_details[:fault_type] == "GL"
        dyn_model_dict[:meta][:gl_element] = ts_fault_details[:gl][:element_2_disconnect]
    elseif ts_fault_details[:fault_type] == "OB"
        dyn_model_dict[:meta][:ob_branch_ids] = copy(ts_fault_details[:ob][:branch_id])
    end

    # ------------------------------------------------------------------ SC: bus fault
    if ts_fault_details[:fault_type] == "SC"
        if ts_fault_details[:sc][:fault_location] != "bus"
            throw(ArgumentError("FULL_BUS SC supports bus faults only."))
        end

        # Split the simulation horizon into fault-on and post-fault sub-windows.
        t_start_sim, t_end_sim, t_step, t_start_fault, clearing_time, t_clear_fault,
            t_window_fault, t_window_postf, t_window_total = time_windows_sc(simulation)

        # Fault-on admittance: pre-fault sparse Ybus with a large shunt at the faulted bus
        # (a bus short-circuit clamps that bus voltage). Network only — loads stay in the
        # ZIP balance, never folded into Ybus.
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

        # Pre-fault block: internal EMF E, rotor angle δ, mechanical power P_m + their links
        # to the ACOPF terminal variables (initial-condition equality constraints).
        model, dyn_model_dict = Define_Initial_Condition_4_fullbus!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode)

        # Fault-on window: per-step bus/gen states + swing dynamics on the faulted Ybus,
        # initialized from the pre-fault equilibrium (δ_0, Δω_0).
        model, dyn_model_dict = Define_Fault_Dynamic_Model_fullbus!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V_fault, val_θ_fault, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            include_governor=include_governor, governor_limiter=governor_limiter)

        if ts_fault_details[:sc][:bus][:disconnect_branch]
            # Clear the fault by tripping the faulted branch, then rebuild the (network-only)
            # post-fault Ybus on the new topology.
            DCIR_mod.l_status[ts_fault_details[:sc][:bus][:branch_id_2_disconnect]] .= 0
            Ybus_postf = Calculate_Ybus_fullbus_dynamics(
                Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
            Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
                "Ybus_postfault", Ybus_postf)

            # The post-fault window starts from the last fault-on step, so grab the final
            # Pe/δ/Δω of each generator as its initial condition.
            last_var_Pe_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pe_tf])
            last_var_δ_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:δ_tf])
            last_var_Δω_tf = OrderedDict(
                g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Δω_tf])

            model, dyn_model_dict = Define_PostFault_Dynamic_Model_fullbus!(
                model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
                active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_postf, Ybus_postf,
                opf_dict[:vars][:V], opf_dict[:vars][:θ], last_var_Pe_tf, last_var_δ_tf,
                last_var_Δω_tf, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
                bus_ids=collect(Int, DBUS_mod.bus),
                include_governor=include_governor, governor_limiter=governor_limiter)
        else
            # A self-clearing SC (same topology before/after) would leave the post-fault
            # network identical to the pre-fault one; not supported on this path.
            throw(ArgumentError("FULL_BUS SC requires branch trip after clearing (disconnect_branch=true)."))
        end

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}(
            :t_start_sim => t_start_sim, :t_end_sim => t_end_sim, :t_step => t_step,
            :t_start_fault => t_start_fault, :clearing_time => clearing_time,
            :t_clear_fault => t_clear_fault, :t_window_fault => t_window_fault,
            :t_window_postf => t_window_postf, :t_window_total => t_window_total,
        )

    # --------------------------------------------------- GL: generation / load loss
    elseif ts_fault_details[:fault_type] == "GL"
        # GL has no clearing/branch-trip stage, so only one perturbed window is needed.
        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        if ts_fault_details[:gl][:element_2_disconnect] == "gen"
            # Disconnect the tripped generator(s): zero their status, recompute the active
            # set, and rebuild the bus/gen/circuit topology map.
            gen_ids = ts_fault_details[:gl][:gen][:gen_id]
            for gen_id in gen_ids
                DGEN_mod.g_status[gen_id] = 0
            end
            ts_fault_details[:gl][:gen][:bus_id] = [Int(DGEN_mod.bus[g]) for g in gen_ids]
            active_gen = findall(x -> x == 1, DGEN_mod.g_status)
            _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)

        elseif ts_fault_details[:gl][:element_2_disconnect] == "load"
            # Scale the load at the affected bus(es) by `percent_power` (full drop = 0%).
            bus_ids = ts_fault_details[:gl][:load][:bus_id]
            apply_gl_load_scaling!(
                DBUS_mod, bus_ids, ts_fault_details[:gl][:load][:percent_power])
            _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)
        else
            throw(ArgumentError("Unknown GL element: $(ts_fault_details[:gl][:element_2_disconnect])"))
        end

        # GL keeps the network intact; the perturbation lives in the load/gen data, so the
        # "fault" Ybus is just the (network-only) Ybus of the post-event system.
        Ybus_fault = Calculate_Ybus_fullbus_dynamics(
            Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
            "Ybus_fault", Ybus_fault)

        model, dyn_model_dict = Define_Initial_Condition_4_fullbus!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode)

        model, dyn_model_dict = Define_Fault_Dynamic_Model_fullbus!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            include_governor=include_governor, governor_limiter=governor_limiter)

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}(
            :t_start_sim => t_start_sim, :t_end_sim => t_end_sim, :t_step => t_step,
            :t_start_fault => t_start_fault, :t_window_fault => t_window_fault,
            :t_window_total => t_window_total,
        )

    # --------------------------------------------------- OB: open branch (no short-circuit)
    elseif ts_fault_details[:fault_type] == "OB"
        # Same single-window timeline as GL; topology change only (no fault-on shunt).
        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        branch_ids = ts_fault_details[:ob][:branch_id]
        ts_fault_details[:ob][:from_bus] = [Int(DCIR_mod.from_bus[b]) for b in branch_ids]
        ts_fault_details[:ob][:to_bus]   = [Int(DCIR_mod.to_bus[b]) for b in branch_ids]
        DCIR_mod.l_status[branch_ids] .= 0
        _, bus_gen_circ_on = Organize_Bus_Gen_Circ(DBUS_mod, DGEN_mod, DCIR_mod)

        Ybus_fault = Calculate_Ybus_fullbus_dynamics(
            Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA))
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices],
            "Ybus_fault", Ybus_fault)

        model, dyn_model_dict = Define_Initial_Condition_4_fullbus!(
            model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen,
            base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], val_V, val_θ, val_Pg, val_Qg;
            mech_power_mode=mech_power_mode)

        model, dyn_model_dict = Define_Fault_Dynamic_Model_fullbus!(
            model, dyn_model_dict, DBUS_mod, bus_gen_circ_on, DGEN_mod, DGEN_DYN_mod,
            active_gen, nBUS, base_MVA, ZIP_P, ZIP_Q, t_step, t_window_fault, Ybus_fault,
            opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            Δω_0, ω_syn, δ_tol, val_V, val_θ, val_Pg, val_Qg;
            bus_ids=collect(Int, DBUS_mod.bus),
            include_governor=include_governor, governor_limiter=governor_limiter)

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

    # Index every stored constraint family so duals can be exported by name downstream.
    build_dual_registry!(dyn_model_dict)
    return model, dyn_model_dict, dyn_parameters_dict
end

# ===================================================================================
# Pre-fault: E, δ, P_m — mirrors `Define_Initial_Condition_4_tsred!`
# ===================================================================================

"""
Create the pre-fault (t=0) generator state and tie it to the ACOPF operating point.

Adds the classical-model decision variables — internal EMF magnitude `E`, rotor angle `δ`,
and (under `USE_PM`) mechanical power `P_m` — optionally with explicit bound constraints,
then warm-starts them from the solved ACOPF point. The equality constraints link these to
the ACOPF terminal variables (V, θ, P_g, Q_g): active/reactive power injection at the
internal node and P_m = P_g at equilibrium. Everything is stored in `dyn_model_dict`.
"""
function Define_Initial_Condition_4_fullbus!(
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
)
    register_ts_bound_limit_specs!(
        ts_input_param, dyn_model_dict, ts_input_param[:limits],
        active_gen, DGEN, base_MVA)
    lims = ts_input_param[:var_limit_specs]
    E_min, E_max = lims[:E]
    δ_min, δ_max = lims[:δ]

    #-------------------------------------------------------------------------
    #                                Variables
    #-------------------------------------------------------------------------
    # Internal EMF magnitude E, bounded when var_bounds[:E] is on.
    begin
        enc = bound_encoding_from_meta(dyn_model_dict[:meta])
        result = var_tsred_gen_voltage_magnitude!(
            model, active_gen, ts_input_param[:var_names][:E];
            bounded=ts_input_param[:var_bounds][:E],
            min_lim=E_min, max_lim=E_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_E_lower, export_key_upper=:ineq_const_E_upper)
        E = store_ts_scalar_var_bounds!(dyn_model_dict, :E, result, :ineq_const_E_lower, :ineq_const_E_upper)
    end

    begin
        result = var_kron_gen_rotor_angle!(
            model, active_gen, ts_input_param[:var_names][:δ];
            bounded=ts_input_param[:var_bounds][:δ], min_lim=δ_min, max_lim=δ_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_δ_lower, export_key_upper=:ineq_const_δ_upper)
        δ = store_ts_scalar_var_bounds!(dyn_model_dict, :δ, result, :ineq_const_δ_lower, :ineq_const_δ_upper)
    end

    P_m = nothing
    if mech_power_mode == USE_PM
        Pm_min, Pm_max = lims[:P_m]
        result = var_kron_gen_mech_power!(
            model, active_gen, ts_input_param[:var_names][:P_m];
            bounded=ts_input_param[:var_bounds][:P_m], min_lim=Pm_min, max_lim=Pm_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_P_m_lower, export_key_upper=:ineq_const_P_m_upper)
        P_m = store_ts_scalar_var_bounds!(dyn_model_dict, :P_m, result, :ineq_const_P_m_lower, :ineq_const_P_m_upper)
    end

    _set_fullbus_init_warm_starts!(E, δ, P_m, active_gen, DGEN, DGEN_DYN,
        val_V, val_θ, val_Pg, val_Qg)

    #-------------------------------------------------------------------------
    #                          Equality Constraints (always-on physics)
    #-------------------------------------------------------------------------
    # Couple the internal EMF (E, δ) to the ACOPF terminal injection so the dynamic state
    # starts exactly on the steady-state operating point: P_g, Q_g must equal the classical
    # active/reactive injection of E∠δ behind Xd_tr. Always built.
    dyn_model_dict[:eq_const][:eq_const_P_init] = eq_const_tsred_initial_active_power!(
        model, V, θ, P_g, E, δ, DGEN, DGEN_DYN, active_gen)
    dyn_model_dict[:eq_const][:eq_const_Q_init] = eq_const_tsred_initial_reactive_power!(
        model, V, θ, Q_g, E, δ, DGEN, DGEN_DYN, active_gen)
    # Pin mechanical power to the dispatched P_g (equilibrium: P_m = P_e = P_g).
    if mech_power_mode == USE_PM
        dyn_model_dict[:eq_const][:eq_const_Pm_init] = eq_const_kron_initial_mechanical_power!(
            model, dyn_model_dict[:vars][:P_m], P_g, active_gen)
    end

    # Record which JuMP refs play the "mechanical power" role (P_m here) for the swing eqs.
    register_mech_power_refs!(dyn_model_dict, P_g, mech_power_mode)
    return model, dyn_model_dict
end

# ===================================================================================
# Fault-on window — thin orchestrator calling per-family builders
# ===================================================================================

"""
Build the fault-on time window: per-step bus phasors, generator electrical/rotor states,
and the swing dynamics on the faulted `Ybus`.

Pipeline (each family stored in `dyn_model_dict`):
1. Time-indexed variables: bus V_tf/θ_tf, generator Pe_tf/Qe_tf/δ_tf/Δω_tf (the "tf"
   suffix tags the fault-on window).
2. Center-of-inertia angle δCOI_tf and the δ-COI stability bounds (`bound_style_δ`).
3. Optional Δω-COI variable + bounds when `constrain_Δω` is set.
4. Network: classical Pe/Qe at each internal node, the nodal injection expressions from
   Ybus, and the ZIP active/reactive power-balance equalities.
5. Swing equations (trapezoidal): δ integrates ω·Δω, Δω integrates (P_m − P_e)/2H − D·Δω.

The first step (t=1) is seeded from the pre-fault equilibrium (δ_0, Δω_0, P_g).
"""
function Define_Fault_Dynamic_Model_fullbus!(
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
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
)
    # Pre-fault state carried in from the initial-condition block: E and δ_0 are the t=0
    # generator state, P_mech is the mechanical-power ref used by the swing equation.
    E = dyn_model_dict[:vars][:E]
    δ_0 = dyn_model_dict[:vars][:δ]
    P_mech = dyn_model_dict[:refs][:P_mech]

    # ---- Time-indexed decision variables for the fault-on window ("tf" suffix) ----
    V_tf = var_fullbus_bus_voltage_magnitude_time!(
        model, bus_ids, time_window, val_V; suffix="tf")
    θ_tf = var_fullbus_bus_angle_time!(
        model, bus_ids, time_window, val_θ; suffix="tf")
    Pe_tf = var_fullbus_gen_Pe_time!(
        model, active_gen, time_window, val_Pg; suffix="tf")
    Qe_tf = var_fullbus_gen_Qe_time!(
        model, active_gen, time_window, val_Qg; suffix="tf")
    δ_tf = var_fullbus_gen_rotor_angle_time!(
        model, active_gen, time_window, δ_0; suffix="tf")
    Δω_tf = var_fullbus_gen_speed_dev_time!(
        model, active_gen, time_window; suffix="tf")

    dyn_model_dict[:vars][:V_tf] = V_tf
    dyn_model_dict[:vars][:θ_tf] = θ_tf
    V_min = dyn_model_dict[:meta][:var_limit_specs][:V_bus_min][1]
    attach_fullbus_V_lower_bounds!(model, dyn_model_dict, V_tf, V_min, :ineq_const_V_tf_lower)
    dyn_model_dict[:vars][:Pe_tf] = Pe_tf
    dyn_model_dict[:vars][:Qe_tf] = Qe_tf
    dyn_model_dict[:vars][:δ_tf] = δ_tf
    dyn_model_dict[:vars][:Δω_tf] = Δω_tf

    # Center-of-inertia (COI) rotor angle: inertia-weighted mean of δ. Stability is judged
    # on each machine's deviation from this COI reference, not on absolute angles.
    δCOI_tf = _attach_δCOI!(model, dyn_model_dict, δ_tf, active_gen, DGEN_DYN,
        time_window, :tf)

    attach_fault_tf_var_bounds!(model, dyn_model_dict)

    # δ-COI stability bounds (flavour set by `bound_style_δ`).
    _add_δ_bounds_fault!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tf, δCOI_tf,
        Δω_tf, Pe_tf, time_window, δ_tol, δ_0, Δω_0, ω_syn, Δt)

    # Optional speed-deviation COI bounds (|Δω_g − Δω_COI| ≤ tol).
    # Optional speed corridor (reference set by `bound_style_Δω`).
    _add_Δω_bounds!(model, dyn_model_dict, :tf, active_gen, DGEN_DYN, Δω_tf, time_window)

    # Classical machine electrical injections (Pe, Qe) as functions of E, δ and terminal V, θ.
    dyn_model_dict[:eq_const][:eq_const_Pe_tf] = eq_const_fullbus_gen_Pe!(
        model, active_gen, DGEN, DGEN_DYN, E, Pe_tf, δ_tf, V_tf, θ_tf, time_window)
    dyn_model_dict[:eq_const][:eq_const_Qe_tf] = eq_const_fullbus_gen_Qe!(
        model, active_gen, DGEN, DGEN_DYN, E, Qe_tf, δ_tf, V_tf, θ_tf, time_window)

    # Nodal power injections from the network (Ybus) at every bus and time step, kept as
    # named expressions so the balance constraints and dual export can reuse them.
    terms_Pb, terms_Qb = _build_fullbus_nodal_injection_terms!(
        model, nBUS, bus_gen_circ_dict, V_tf, θ_tf, time_window, Ybus)
    haskey(dyn_model_dict, :expressions) ||
        (dyn_model_dict[:expressions] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:expressions][:P_inj_tf] = terms_Pb
    dyn_model_dict[:expressions][:Q_inj_tf] = terms_Qb
    # KCL at each bus: network injection = generator injection − ZIP load demand.
    dyn_model_dict[:eq_const][:eq_const_Pbalance_tf] = eq_const_fullbus_Pbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, Pe_tf, terms_Pb,
        V_tf, time_window, ZIP_P)
    dyn_model_dict[:eq_const][:eq_const_Qbalance_tf] = eq_const_fullbus_Qbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, Qe_tf, terms_Qb,
        V_tf, time_window, ZIP_Q)

    # Swing equations, seeded from the pre-fault equilibrium at t=1
    # (trapezoidal or BE via ode_first_step).
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)
    dyn_model_dict[:eq_const][:eq_const_δ_tf] = eq_const_fullbus_δ_swing_fault!(
        model, active_gen, δ_tf, Δω_tf, δ_0, Δω_0, time_window, ω_syn, Δt;
        ode_first_step=ode_fs)

    # Mechanical power in the Δω swing: constant P_m (default) or the governor's
    # time-varying P_mech(t). The governor consumes Δω_tf and returns Pm_tf; the t=1
    # step is anchored to the pre-fault equilibrium P_m.
    if include_governor
        P_m0 = dyn_model_dict[:vars][:P_m]
        Pm_tf = Attach_Governor_fault!(
            model, dyn_model_dict, active_gen, DGEN, DGEN_DYN, base_MVA,
            Δω_tf, time_window, Δt, governor_limiter)
        dyn_model_dict[:eq_const][:eq_const_Δω_tf] = eq_const_fullbus_Δω_swing_fault!(
            model, active_gen, DGEN_DYN, Pm_tf, P_m0, Pe_tf, Δω_tf, Δω_0, P_g,
            time_window, Δt; ode_first_step=ode_fs)
    else
        dyn_model_dict[:eq_const][:eq_const_Δω_tf] = eq_const_fullbus_Δω_swing_fault!(
            model, active_gen, DGEN_DYN, P_mech, Pe_tf, Δω_tf, Δω_0, P_g,
            time_window, Δt; ode_first_step=ode_fs)
    end

    return model, dyn_model_dict
end

# ===================================================================================
# Post-fault window
# ===================================================================================

"""
Build the post-fault time window ("tpf" suffix) on the cleared topology `Ybus`.

Structurally identical to the fault-on window, with two differences:
- The window's initial condition is the *last fault-on step*, passed in as `δ_ant`,
  `Δω_ant`, `Pe_ant` (and `Pe_tf` for the Δω trapezoid's first-step electrical term).
- The network admittance is the post-fault (branch-tripped) Ybus.

The internal EMF `E` and mechanical power `P_mech` are unchanged (constant across the
classical transient), so they are reused from `dyn_model_dict`.
"""
function Define_PostFault_Dynamic_Model_fullbus!(
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
    include_governor::Bool=false,
    governor_limiter::GovernorLimiter=GOV_NO_LIMIT,
)
    # Reused (constant) state from earlier blocks: internal EMF E, mechanical power P_mech,
    # and the full fault-on Pe trajectory (its last point links the two windows' swing eqs).
    E = dyn_model_dict[:vars][:E]
    P_mech = dyn_model_dict[:refs][:P_mech]
    Pe_tf = dyn_model_dict[:vars][:Pe_tf]

    V_tpf = var_fullbus_bus_voltage_magnitude_time!(
        model, bus_ids, time_window, val_V; suffix="tpf")
    θ_tpf = var_fullbus_bus_angle_time!(
        model, bus_ids, time_window, val_θ; suffix="tpf")
    Pe_tpf = var_fullbus_gen_Pe_time!(
        model, active_gen, time_window, val_Pg; suffix="tpf")
    Qe_tpf = var_fullbus_gen_Qe_time!(
        model, active_gen, time_window, val_Qg; suffix="tpf")
    δ_tpf = var_fullbus_gen_rotor_angle_time!(
        model, active_gen, time_window, dyn_model_dict[:vars][:δ]; suffix="tpf")
    Δω_tpf = var_fullbus_gen_speed_dev_time!(
        model, active_gen, time_window; suffix="tpf")

    dyn_model_dict[:vars][:V_tpf] = V_tpf
    dyn_model_dict[:vars][:θ_tpf] = θ_tpf
    V_min = dyn_model_dict[:meta][:var_limit_specs][:V_bus_min][1]
    attach_fullbus_V_lower_bounds!(model, dyn_model_dict, V_tpf, V_min, :ineq_const_V_tpf_lower)
    dyn_model_dict[:vars][:Pe_tpf] = Pe_tpf
    dyn_model_dict[:vars][:Qe_tpf] = Qe_tpf
    dyn_model_dict[:vars][:δ_tpf] = δ_tpf
    dyn_model_dict[:vars][:Δω_tpf] = Δω_tpf

    δCOI_tpf = _attach_δCOI!(model, dyn_model_dict, δ_tpf, active_gen, DGEN_DYN,
        time_window, :tpf)

    attach_postfault_tpf_var_bounds!(model, dyn_model_dict)

    _add_δ_bounds_postf!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tpf, δCOI_tpf,
        Δω_tpf, Pe_tpf, time_window, δ_tol, δ_ant, Δω_ant, Pe_ant, ω_syn, Δt)

    _add_Δω_bounds!(model, dyn_model_dict, :tpf, active_gen, DGEN_DYN, Δω_tpf, time_window)

    dyn_model_dict[:eq_const][:eq_const_Pe_tpf] = eq_const_fullbus_gen_Pe!(
        model, active_gen, DGEN, DGEN_DYN, E, Pe_tpf, δ_tpf, V_tpf, θ_tpf, time_window)
    dyn_model_dict[:eq_const][:eq_const_Qe_tpf] = eq_const_fullbus_gen_Qe!(
        model, active_gen, DGEN, DGEN_DYN, E, Qe_tpf, δ_tpf, V_tpf, θ_tpf, time_window)

    terms_Pb, terms_Qb = _build_fullbus_nodal_injection_terms!(
        model, nBUS, bus_gen_circ_dict, V_tpf, θ_tpf, time_window, Ybus)
    dyn_model_dict[:expressions][:P_inj_tpf] = terms_Pb
    dyn_model_dict[:expressions][:Q_inj_tpf] = terms_Qb
    dyn_model_dict[:eq_const][:eq_const_Pbalance_tpf] = eq_const_fullbus_Pbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, Pe_tpf, terms_Pb,
        V_tpf, time_window, ZIP_P)
    dyn_model_dict[:eq_const][:eq_const_Qbalance_tpf] = eq_const_fullbus_Qbalance!(
        model, DBUS, nBUS, bus_gen_circ_dict, base_MVA, V, Qe_tpf, terms_Qb,
        V_tpf, time_window, ZIP_Q)

    # Swing equations seeded from the last fault-on step (δ_ant, Δω_ant, Pe_ant).
    ode_fs = get(dyn_model_dict[:meta], :ode_first_step, :trapezoidal)
    dyn_model_dict[:eq_const][:eq_const_δ_tpf] = eq_const_fullbus_δ_swing_postf!(
        model, active_gen, δ_tpf, Δω_tpf, δ_ant, Δω_ant, time_window, ω_syn, Δt;
        ode_first_step=ode_fs)

    # Governor post-fault window continues the fault-on governor states; the t=1 step is
    # anchored to the last fault-on mechanical power. Non-governor runs keep constant P_m.
    if include_governor
        Pm_tf_last = OrderedDict(
            g => last(inner).second for (g, inner) in dyn_model_dict[:vars][:Pm_tf])
        Pm_tpf = Attach_Governor_postf!(
            model, dyn_model_dict, active_gen, DGEN, DGEN_DYN, base_MVA,
            Δω_tpf, time_window, Δt, governor_limiter)
        dyn_model_dict[:eq_const][:eq_const_Δω_tpf] = eq_const_fullbus_Δω_swing_postf!(
            model, active_gen, DGEN_DYN, Pm_tpf, Pm_tf_last, Pe_tpf, Δω_tpf, Δω_ant, Pe_ant,
            Pe_tf, time_window, Δt; ode_first_step=ode_fs)
    else
        dyn_model_dict[:eq_const][:eq_const_Δω_tpf] = eq_const_fullbus_Δω_swing_postf!(
            model, active_gen, DGEN_DYN, P_mech, Pe_tpf, Δω_tpf, Δω_ant, Pe_ant,
            Pe_tf, time_window, Δt; ode_first_step=ode_fs)
    end

    return model, dyn_model_dict
end
