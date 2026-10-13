#=
================================================================================
 functions_2_build_TS_model_w_Kron.jl  —  TSC-ACOPF (nonlinear Pe on Yred)
================================================================================
 Phase 2: shared var/COI/swing-bound builders live in
 `functions_2_build_TS_model_w_Kron_shared.jl`.  This file keeps TSC-ACOPF-only
 pieces: internal voltage E, P/Q init, nonlinear Pe, Qe expressions (diagnostic),
 and nonlinear Δω swing.
================================================================================
=#

""" tsred = Transient Stability using the Kron-reduced network."""

# ==================================================================================
# Function to create the dynamic model
# ==================================================================================
function Make_Dynamic_Model_tsred!(
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
    mech_power_mode::MechPowerMode=USE_PG,
    constrain_δ::Bool=true,
    bound_style_δ::Symbol=:swing_propagated,
    δ_ref_gen_id::Union{Nothing, Int}=nothing,
    constrain_Δω::Bool=false,
    bound_style_Δω::Symbol=:coi_box,
    Δω_tol::Tuple{Float64, Float64}=(-0.5, 0.5),
    )

    ts_input_param = build_ts_input_param(ts_builder)
    dyn_model_dict = OrderedDict{Symbol, Any}()       # Define a Dicitionary to Save Variables, Constraints and other relevant data
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
        :fault_type => ts_fault_details[:fault_type],
        :ineq_cons => ts_input_param[:ineq_cons],
        :var_bounds => ts_input_param[:var_bounds],
        :bound_encoding => bound_encoding_from_param(ts_input_param),
    )


    DBUS_mod     = deepcopy(DBUS)
    DGEN_mod     = deepcopy(DGEN)
    DGEN_DYN_mod = deepcopy(DGEN_DYN)
    DCIR_mod     = deepcopy(DCIR)

    δ_tol, f_syn, ω_syn, Δω_0 = common_ts_parameters(simulation)

    if ts_fault_details[:fault_type] == "SC"
        active_gen      = findall(x -> x == 1, DGEN_mod.g_status)     # Indices of active generators
        active_branches = findall(x -> x == 1, DCIR_mod.l_status)     # Indices of active branches
        nGEN_active = length(active_gen) # Count how many generators are active

        t_start_sim, t_end_sim, t_step, t_start_fault, clearing_time, t_clear_fault,
            t_window_fault, t_window_postf, t_window_total = time_windows_sc(simulation)

        if ts_fault_details[:sc][:fault_location] == "bus"
            if !haskey(ts_fault_details[:sc][:bus], :bus_id) throw(ArgumentError("You must define a bus ID to apply a fault.")) end
            begin
                # Calculate the new admittance matrix for the fault period and its reduced version 
                sc_bus = ts_fault_details[:sc][:bus]
                # Fault shunt from FaultConfig. The default adds the inductive solid short.
                Ybus_fault = Calculate_Ybus_fault_SC_Kron(
                    Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA),
                    DBUS_mod, DGEN_mod, DGEN_DYN_mod, nBUS, active_gen, base_MVA,
                    sc_bus[:bus_id];
                    y_fault=sc_bus[:admittance_pu], stamp=sc_bus[:stamp])
                Yred_fault = Reduce_Matrix(Ybus_fault, nBUS)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_fault", Ybus_fault)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_fault", Yred_fault)

                # =================================================================================
                # Link the steady-state operating point with transient stability decision variables
                # =================================================================================
                model, dyn_model_dict = Define_Initial_Condition_4_tsred!(model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen, base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g], opf_dict[:vars][:Q_g]; mech_power_mode=mech_power_mode)

                # =================================================================================
                # Model the Fault
                # =================================================================================
                model, dyn_model_dict = Define_Fault_Dynamic_Model_tsred!(model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step, t_window_fault, Yred_fault, dyn_model_dict[:vars][:δ], Δω_0, ω_syn, δ_tol)
            end

            if ts_fault_details[:sc][:bus][:disconnect_branch] == true # A branch will be disconnected
                # =================================================================================
                # Model the Post-Fault
                # =================================================================================
                begin
                    DCIR_mod.l_status[ts_fault_details[:sc][:bus][:branch_id_2_disconnect]] .= 0

                    # Calculate the new admittance matrix for the post-fault period and its reduced version 
                    Ybus_postf = Calculate_Ybus_postf_ClearFault_Kron(DBUS_mod, DCIR_mod, DGEN_mod, DGEN_DYN_mod, nBUS, nCIR, active_gen, base_MVA)
                    Yred_postf = Reduce_Matrix(Ybus_postf, nBUS)
                    Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_postfault", Ybus_postf)
                    Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_postfault", Yred_postf)

                    # Extract the last variable for each generator
                    last_var_Pe_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:Pe_tf])
                    last_var_δ_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:δ_tf])
                    last_var_Δω_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:Δω_tf])
                    
                    model, dyn_model_dict = Define_PostFault_Dynamic_Model_tsred!(model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step, t_window_postf, Yred_postf, last_var_Pe_tf, last_var_δ_tf, last_var_Δω_tf, ω_syn, δ_tol)
                end

            else # The fault will extinguish by itself
                # =================================================================================
                # Model the Post-Fault
                # =================================================================================
                begin
                    # Extract the last variable for each generator
                    last_var_Pe_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:Pe_tf])
                    last_var_δ_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:δ_tf])
                    last_var_Δω_tf = OrderedDict(gen_id => last(inner_dict).second for (gen_id, inner_dict) in dyn_model_dict[:vars][:Δω_tf])

                    # Calculate the new admittance matrix for the post-fault period and its reduced version 
                    Ybus_postf = Calculate_Ybus_postf_ClearFault_Kron(DBUS_mod, DCIR_mod, DGEN_mod, DGEN_DYN_mod, nBUS, nCIR, active_gen, base_MVA)
                    Yred_postf = Reduce_Matrix(Ybus_postf, nBUS)
                    Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_postfault", Ybus_postf)
                    Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_postfault", Yred_postf)
                    
                    model, dyn_model_dict = Define_PostFault_Dynamic_Model_tsred!(model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step, t_window_postf, Yred_postf, last_var_Pe_tf, last_var_δ_tf, last_var_Δω_tf, ω_syn, δ_tol)
                end
            end


        # Branch SC faults are out of scope (bus faults only). See roadmap.
        elseif ts_fault_details[:sc][:fault_location] == "branch"
            throw(ArgumentError("Branch SC faults are not implemented (bus faults only)."))
        end

    elseif ts_fault_details[:fault_type] == "GL"

        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        if ts_fault_details[:gl][:element_2_disconnect] == "gen"
            begin
                DGEN_mod.g_status[ts_fault_details[:gl][:gen][:gen_id]] .= 0 # Modify the status of the generator
                active_gen = findall(x -> x == 1, DGEN_mod.g_status)         # Indices of active generators

                # Calculate the new admittance matrix for the fault period and its reduced version 
                Ybus_fault = Calculate_Ybus_fault_LG_Kron(Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA), DBUS_mod, DGEN_mod, DGEN_DYN_mod, nBUS, active_gen, base_MVA)  # Admittance matrix augmented for the fault
                Yred_fault = Reduce_Matrix(Ybus_fault, nBUS)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_fault", Ybus_fault)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_fault", Yred_fault)

                # =================================================================================
                # Link the steady-state operating point with transient stability decision variables
                # =================================================================================
                model, dyn_model_dict = Define_Initial_Condition_4_tsred!(model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen, base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g], opf_dict[:vars][:Q_g]; mech_power_mode=mech_power_mode)

                # =================================================================================
                # Model the Fault
                # =================================================================================
                model, dyn_model_dict = Define_Fault_Dynamic_Model_tsred!(model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step, t_window_fault, Yred_fault, dyn_model_dict[:vars][:δ], Δω_0, ω_syn, δ_tol)
            end
            
        elseif ts_fault_details[:gl][:element_2_disconnect] == "load"
            begin
                bus_ids = ts_fault_details[:gl][:load][:bus_id]
                apply_gl_load_scaling!(
                    DBUS_mod, bus_ids, ts_fault_details[:gl][:load][:percent_power])

                active_gen = findall(x -> x == 1, DGEN_mod.g_status)         # Indices of active generators

                # Calculate the new admittance matrix for the fault period and its reduced version 
                Ybus_fault = Calculate_Ybus_fault_LG_Kron(Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA), DBUS_mod, DGEN_mod, DGEN_DYN_mod, nBUS, active_gen, base_MVA)  # Admittance matrix augmented for the fault
                Yred_fault = Reduce_Matrix(Ybus_fault, nBUS)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_fault", Ybus_fault)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_fault", Yred_fault)

                # =================================================================================
                # Link the steady-state operating point with transient stability decision variables
                # =================================================================================
                model, dyn_model_dict = Define_Initial_Condition_4_tsred!(model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod, active_gen, base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g], opf_dict[:vars][:Q_g]; mech_power_mode=mech_power_mode)

                # =================================================================================
                # Model the Fault
                # =================================================================================
                model, dyn_model_dict = Define_Fault_Dynamic_Model_tsred!(model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step, t_window_fault, Yred_fault, dyn_model_dict[:vars][:δ], Δω_0, ω_syn, δ_tol)

            end
            
        end

    elseif ts_fault_details[:fault_type] == "OB"

        t_start_sim, t_end_sim, t_step, t_start_fault, t_window_fault, t_window_total =
            time_windows_gld(simulation)

        begin
            branch_ids = ts_fault_details[:ob][:branch_id]
            DCIR_mod.l_status[branch_ids] .= 0

            active_gen = findall(x -> x == 1, DGEN_mod.g_status)

            # OB mirrors the GL single-window path, but the perturbation is purely
            # topological: open the requested branches, rebuild Ybus/Yred, and
            # simulate one post-event window with no prior short-circuit stage.
            Ybus_fault = Calculate_Ybus_fault_LG_Kron(
                Calculate_Ybus_sparse(DBUS_mod, DCIR_mod, nBUS, nCIR, base_MVA),
                DBUS_mod, DGEN_mod, DGEN_DYN_mod, nBUS, active_gen, base_MVA)
            Yred_fault = Reduce_Matrix(Ybus_fault, nBUS)
            Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_fault", Ybus_fault)
            Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus_red_fault", Yred_fault)

            model, dyn_model_dict = Define_Initial_Condition_4_tsred!(
                model, ts_input_param, dyn_model_dict, DGEN_mod, DGEN_DYN_mod,
                active_gen, base_MVA, opf_dict[:vars][:V], opf_dict[:vars][:θ],
                opf_dict[:vars][:P_g], opf_dict[:vars][:Q_g];
                mech_power_mode=mech_power_mode)

            model, dyn_model_dict = Define_Fault_Dynamic_Model_tsred!(
                model, dyn_model_dict, active_gen, DGEN_DYN_mod, t_step,
                t_window_fault, Yred_fault, dyn_model_dict[:vars][:δ], Δω_0,
                ω_syn, δ_tol)
        end
    else
        throw(ArgumentError("Unknown fault_type: $(ts_fault_details[:fault_type])"))
    end

    # Save time parameters and common parameters in a dictionary
    begin
        dyn_model_dict[:active_gen] = active_gen

        dyn_parameters_dict = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time] = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:time][:t_start_sim] = t_start_sim
        dyn_parameters_dict[:time][:t_end_sim] = t_end_sim
        dyn_parameters_dict[:time][:t_step] = t_step
        dyn_parameters_dict[:time][:t_start_fault] = t_start_fault
        if @isdefined(clearing_time) dyn_parameters_dict[:time][:clearing_time] = clearing_time end
        if @isdefined(t_clear_fault) dyn_parameters_dict[:time][:t_clear_fault] = t_clear_fault end
        dyn_parameters_dict[:time][:t_window_fault] = t_window_fault
        if @isdefined(t_window_postf) dyn_parameters_dict[:time][:t_window_postf] = t_window_postf end
        dyn_parameters_dict[:time][:t_window_total] = t_window_total

        dyn_parameters_dict[:common] = OrderedDict{Symbol, Any}()
        dyn_parameters_dict[:common][:δ_tol] = δ_tol
        dyn_parameters_dict[:common][:δ_pair_tol] = δ_pair_tol_tuple(simulation)
        dyn_parameters_dict[:common][:f_syn] = f_syn
        dyn_parameters_dict[:common][:ω_syn] = ω_syn
        dyn_parameters_dict[:common][:Δω_0] = Δω_0
        dyn_parameters_dict[:common][:constrain_Δω] = constrain_Δω
        if constrain_Δω
            dyn_parameters_dict[:common][:Δω_tol] = Δω_tol
        end
    end

    # Record which constraint families exist so the dual save layer can export
    # them without hard-coded haskey checks (see DynDualRegistry.jl, Phase 1.3).
    build_dual_registry!(dyn_model_dict)

    return model, dyn_model_dict, dyn_parameters_dict
 
end

# ===================================================================================
#                 DEFINE INITIAL VARIABLES IN THE PRE-FAULT
# ===================================================================================
# Function used to define E, δ and Pm variables of each active generator
function Define_Initial_Condition_4_tsred!(model::JuMP.Model,
    ts_input_param::OrderedDict{Symbol, Any},
    dyn_model_dict::OrderedDict{Symbol, Any},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64},
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef};
    mech_power_mode::MechPowerMode=USE_PG
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
    # Initial Internal Voltage Magnitude
    begin
        enc = bound_encoding_from_meta(dyn_model_dict[:meta])
        result = var_tsred_gen_voltage_magnitude!(model, active_gen, ts_input_param[:var_names][:E]; bounded = ts_input_param[:var_bounds][:E], min_lim = E_min, max_lim = E_max, encoding=enc, meta=dyn_model_dict[:meta], export_key_lower=:ineq_const_E_lower, export_key_upper=:ineq_const_E_upper)
        E = store_ts_scalar_var_bounds!(dyn_model_dict, :E, result, :ineq_const_E_lower, :ineq_const_E_upper)
    end

    # Initial Rotor Angle
    begin
        result = var_kron_gen_rotor_angle!(model, active_gen, ts_input_param[:var_names][:δ]; bounded = ts_input_param[:var_bounds][:δ], min_lim = δ_min, max_lim = δ_max, encoding=enc, meta=dyn_model_dict[:meta], export_key_lower=:ineq_const_δ_lower, export_key_upper=:ineq_const_δ_upper)
        δ = store_ts_scalar_var_bounds!(dyn_model_dict, :δ, result, :ineq_const_δ_lower, :ineq_const_δ_upper)
    end

    # Initial Mechanical Power (only when mech_power_mode == USE_PM)
    if mech_power_mode == USE_PM
        Pm_min, Pm_max = lims[:P_m]
        result = var_kron_gen_mech_power!(model, active_gen, ts_input_param[:var_names][:P_m];
            bounded=ts_input_param[:var_bounds][:P_m], min_lim=Pm_min, max_lim=Pm_max,
            encoding=enc, meta=dyn_model_dict[:meta],
            export_key_lower=:ineq_const_P_m_lower, export_key_upper=:ineq_const_P_m_upper)
        P_m = store_ts_scalar_var_bounds!(dyn_model_dict, :P_m, result, :ineq_const_P_m_lower, :ineq_const_P_m_upper)
    end

    #-------------------------------------------------------------------------
    #                          Equality Constraints (always-on physics)
    #-------------------------------------------------------------------------
    # Couple OPF terminal injections to the classical EMF behind X'd; always
    # built. Pin P_m = P_g only when mechanical power is an independent variable.
    dyn_model_dict[:eq_const][:eq_const_P_init] = eq_const_tsred_initial_active_power!(
        model, V, θ, P_g, E, δ, DGEN, DGEN_DYN, active_gen)
    dyn_model_dict[:eq_const][:eq_const_Q_init] = eq_const_tsred_initial_reactive_power!(
        model, V, θ, Q_g, E, δ, DGEN, DGEN_DYN, active_gen)
    if mech_power_mode == USE_PM
        dyn_model_dict[:eq_const][:eq_const_Pm_init] = eq_const_kron_initial_mechanical_power!(
            model, dyn_model_dict[:vars][:P_m], P_g, active_gen)
    end

    register_mech_power_refs!(dyn_model_dict, P_g, mech_power_mode)

    return model, dyn_model_dict
    
end

# ===================================================================================
#                                  FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the fault period
function Define_Fault_Dynamic_Model_tsred!(model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Δt::Float64,
    time_window::Vector{Float64},
    Yred::Matrix,
    δ_0::OrderedDict{Int64, JuMP.VariableRef},
    Δω_0::Float64,
    ω_syn::Float64,
    δ_tol::Tuple{Float64, Float64}
	)

    P_g    = dyn_model_dict[:refs][:P_g]
    P_mech = dyn_model_dict[:refs][:P_mech]

    δ_tf  = var_kron_gen_time_generic!(model, active_gen, "δ_tf",  time_window)
    Δω_tf = var_kron_gen_time_generic!(model, active_gen, "Δω_tf", time_window)
    Pe_tf = var_kron_gen_time_generic!(model, active_gen, "Pe_tf",  time_window)

    dyn_model_dict[:vars][:δ_tf]    = δ_tf
    dyn_model_dict[:vars][:Δω_tf]   = Δω_tf
    dyn_model_dict[:vars][:Pe_tf]   = Pe_tf

    # eq_const_δCOI_tf = eq_const_kron_COI_generic!(model, δ_tf, δCOI_tf, active_gen, DGEN_DYN, time_window)
    # eq_const_Pe_tf   = eq_const_tsred_Pe_generic!(model, E, Pe_tf, δ_tf, active_gen, t_window, Yred)
    # eq_const_δ_tf    = eq_const_kron_δ_swingeq_generic!(model, active_gen, δ_tf, Δω_tf, δ_0, Δω_0, t_window, ω_syn, Δt)
    # eq_const_Δω_tf   = eq_const_tsred_Δω_swingeq_generic!(model, active_gen, DGEN_DYN, Pg, Pm, Pe_tf, Δω_tf, Δω_0, t_window, Δt)

    δCOI_tf = _attach_δCOI!(model, dyn_model_dict, δ_tf, active_gen, DGEN_DYN,
        time_window, :tf)
    dyn_model_dict[:eq_const][:eq_const_Pe_tf]   = eq_const_tsred_Pe_generic!(model, dyn_model_dict[:vars][:E], Pe_tf, δ_tf, active_gen, time_window, Yred)
    dyn_model_dict[:eq_const][:eq_const_δ_tf]    = eq_const_kron_δ_swingeq_generic!(model, active_gen, δ_tf, Δω_tf, δ_0, Δω_0, time_window, ω_syn, Δt)
    dyn_model_dict[:eq_const][:eq_const_Δω_tf]   = eq_const_tsred_Δω_swingeq_generic!(model, active_gen, DGEN_DYN, P_g, P_mech, Pe_tf, Δω_tf, Δω_0, time_window, Δt)

    haskey(dyn_model_dict, :expressions) ||
        (dyn_model_dict[:expressions] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:expressions][:Qe_tf] = expr_tsred_Qe_generic!(
        dyn_model_dict[:vars][:E], δ_tf, active_gen, time_window, Yred)

    attach_fault_tf_var_bounds!(model, dyn_model_dict)

    # δ-COI stability bounds (flavour set by `bound_style_δ`).
    _add_δ_bounds_fault!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tf, δCOI_tf,
        Δω_tf, Pe_tf, time_window, δ_tol, δ_0, Δω_0, ω_syn, Δt)

    # Optional speed corridor (reference set by `bound_style_Δω`).
    _add_Δω_bounds!(model, dyn_model_dict, :tf, active_gen, DGEN_DYN, Δω_tf, time_window)

    return model, dyn_model_dict

end

# ===================================================================================
#                               POST-FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the post-fault period
function Define_PostFault_Dynamic_Model_tsred!(model::JuMP.Model,
    dyn_model_dict::OrderedDict{Symbol, Any},
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Δt::Float64,
    time_window::Vector{Float64},
    Yred::Matrix,
    Pe_ant::OrderedDict{Int64, JuMP.VariableRef},
    δ_ant::OrderedDict{Int64, JuMP.VariableRef},
    Δω_ant::OrderedDict{Int64, JuMP.VariableRef},
    ω_syn::Float64,
    δ_tol::Tuple{Float64, Float64}
	)

    P_mech = dyn_model_dict[:refs][:P_mech]

    δ_tpf  = var_kron_gen_time_generic!(model, active_gen, "δ_tpf",  time_window)
    Δω_tpf = var_kron_gen_time_generic!(model, active_gen, "Δω_tpf", time_window)
    Pe_tpf = var_kron_gen_time_generic!(model, active_gen, "Pe_tpf",  time_window)

    dyn_model_dict[:vars][:δ_tpf]    = δ_tpf
    dyn_model_dict[:vars][:Δω_tpf]   = Δω_tpf
    dyn_model_dict[:vars][:Pe_tpf]   = Pe_tpf

    δCOI_tpf = _attach_δCOI!(model, dyn_model_dict, δ_tpf, active_gen, DGEN_DYN,
        time_window, :tpf)
    dyn_model_dict[:eq_const][:eq_const_Pe_tpf]   = eq_const_tsred_Pe_generic!(model, dyn_model_dict[:vars][:E], Pe_tpf, δ_tpf, active_gen, time_window, Yred)
    dyn_model_dict[:eq_const][:eq_const_δ_tpf]    = eq_const_kron_δ_swingeq_generic!(model, active_gen, δ_tpf, Δω_tpf, δ_ant, Δω_ant, time_window, ω_syn, Δt)
    dyn_model_dict[:eq_const][:eq_const_Δω_tpf]   = eq_const_tsred_Δω_swingeq_generic!(model, active_gen, DGEN_DYN, Pe_ant, P_mech, Pe_tpf, Δω_tpf, Δω_ant, time_window, Δt)

    haskey(dyn_model_dict, :expressions) ||
        (dyn_model_dict[:expressions] = OrderedDict{Symbol, Any}())
    dyn_model_dict[:expressions][:Qe_tpf] = expr_tsred_Qe_generic!(
        dyn_model_dict[:vars][:E], δ_tpf, active_gen, time_window, Yred)

    attach_postfault_tpf_var_bounds!(model, dyn_model_dict)

    _add_δ_bounds_postf!(
        model, dyn_model_dict, active_gen, DGEN_DYN, P_mech, δ_tpf, δCOI_tpf,
        Δω_tpf, Pe_tpf, time_window, δ_tol, δ_ant, Δω_ant, Pe_ant, ω_syn, Δt)

    _add_Δω_bounds!(model, dyn_model_dict, :tpf, active_gen, DGEN_DYN, Δω_tpf, time_window)

    return model, dyn_model_dict

end

