# Function to build the AC OPF model including additional decision variables to represent the power flows in the branches
function Make_ACOPF_Model!(
    model::Model, 
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64,
    opf_input_param::OrderedDict{Symbol, Any}
    )

    """
    This function builds a general AC OPF Model with the objective function, variables and constraints
    """
    active_gen               = findall(x -> x == 1, DGEN.g_status)     # Indices of active generators
    active_branches          = findall(x -> x == 1, DCIR.l_status)     # Indices of active branches
    pair_info, pair_circ_map = get_angle_limits(DCIR, active_branches) # Limits for angular differences between adjacent buses
    limits = opf_input_param[:limits]                                  # DispatchLimitsConfig (θ box, clamps)

    # Rated branches: ±l_cap_1/base_MVA; unrated (l_cap_1 == 0): ±Inf → skipped at attach
    lower_bounds_circ, upper_bounds_circ = branch_flow_limit_vectors(DCIR, base_MVA)

    # Defining a dictionary to save opf related data
    opf_dict = init_opf_dict!()
    enc = bound_encoding_from_param(opf_input_param)

    #------------------------------------------
    # Check if there is at least one swing bus
    #------------------------------------------
    SW = findall(x -> x == 3, DBUS.type)
    if isempty(SW) 
        throw(ArgumentError("You must define one bus as the SLACK BUS (type 3).")) 
    elseif length(SW) > 1 
        throw(ArgumentError("This code still does not support more than one SLACK BUS (type 3).")) 
    end

    # ========================================================================
    #              DEFINE THE PRIMAL VARIABLES OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                                Buses
    #-------------------------------------------------------------------------
    # Voltage Magnitude
    begin
        result = var_voltage_magnitude!(model, DBUS.bus, opf_input_param[:var_names][:V]; bounded = opf_input_param[:var_bounds][:V], min_lim = DBUS.v_min, max_lim = DBUS.v_max, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_volt_mag_lower, export_key_upper=:ineq_const_volt_mag_upper)
        V = store_scalar_var_bounds!(opf_dict, :V, result, :ineq_const_volt_mag_lower, :ineq_const_volt_mag_upper)
    end

    # Voltage Angle
    begin
        result = var_voltage_angle!(model, DBUS.bus, opf_input_param[:var_names][:θ]; bounded = opf_input_param[:var_bounds][:θ], min_lim=limits.θ_min_rad, max_lim=limits.θ_max_rad, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_volt_ang_lower, export_key_upper=:ineq_const_volt_ang_upper)
        θ = store_scalar_var_bounds!(opf_dict, :θ, result, :ineq_const_volt_ang_lower, :ineq_const_volt_ang_upper)
    end

    #-------------------------------------------------------------------------
    #                            Generators
    #-------------------------------------------------------------------------
    # Active power generated
    begin
        result = var_gen_power_active!(model, DGEN.id[active_gen], opf_input_param[:var_names][:P_g]; bounded = opf_input_param[:var_bounds][:P_g], min_lim=DGEN.pg_min[active_gen] ./ base_MVA, max_lim=DGEN.pg_max[active_gen] ./ base_MVA, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_pg_lower, export_key_upper=:ineq_const_pg_upper)
        P_g = store_scalar_var_bounds!(opf_dict, :P_g, result, :ineq_const_pg_lower, :ineq_const_pg_upper)
    end

    # Reactive power generated
    begin
        result = var_gen_power_reactive!(model, DGEN.id[active_gen], opf_input_param[:var_names][:Q_g]; bounded = opf_input_param[:var_bounds][:Q_g], min_lim=DGEN.qg_min[active_gen] ./ base_MVA, max_lim=DGEN.qg_max[active_gen] ./ base_MVA, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_qg_lower, export_key_upper=:ineq_const_qg_upper)
        Q_g = store_scalar_var_bounds!(opf_dict, :Q_g, result, :ineq_const_qg_lower, :ineq_const_qg_upper)
    end
  
    #-------------------------------------------------------------------------
    #                                Branches
    #-------------------------------------------------------------------------
    # Active power flow from bus "i" to bus "k"
    begin
        result = var_powerflow_ik_active!(model, DCIR.id[active_branches], DCIR.from_bus[active_branches], DCIR.to_bus[active_branches], opf_input_param[:var_names][:P_ik]; bounded = opf_input_param[:var_bounds][:P_ik], min_lim = lower_bounds_circ[active_branches], max_lim = upper_bounds_circ[active_branches], encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_pik_lower, export_key_upper=:ineq_const_pik_upper)
        P_ik = store_scalar_var_bounds!(opf_dict, :P_ik, result, :ineq_const_pik_lower, :ineq_const_pik_upper)
    end

    # Reactive power flow from bus "i" to bus "k"
    begin
        result = var_powerflow_ik_reactive!(model, DCIR.id[active_branches], DCIR.from_bus[active_branches], DCIR.to_bus[active_branches], opf_input_param[:var_names][:Q_ik]; bounded = opf_input_param[:var_bounds][:Q_ik], min_lim = lower_bounds_circ[active_branches], max_lim = upper_bounds_circ[active_branches], encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_qik_lower, export_key_upper=:ineq_const_qik_upper)
        Q_ik = store_scalar_var_bounds!(opf_dict, :Q_ik, result, :ineq_const_qik_lower, :ineq_const_qik_upper)
    end

    # Active power flow from bus "k" to bus "i"
    begin
        result = var_powerflow_ki_active!(model, DCIR.id[active_branches], DCIR.from_bus[active_branches], DCIR.to_bus[active_branches], opf_input_param[:var_names][:P_ki]; bounded = opf_input_param[:var_bounds][:P_ki], min_lim=lower_bounds_circ[active_branches], max_lim=upper_bounds_circ[active_branches], encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_pki_lower, export_key_upper=:ineq_const_pki_upper)
        P_ki = store_scalar_var_bounds!(opf_dict, :P_ki, result, :ineq_const_pki_lower, :ineq_const_pki_upper)
    end

    # Reactive power flow from bus "k" to bus "i"
    begin
        result = var_powerflow_ki_reactive!(model, DCIR.id[active_branches], DCIR.from_bus[active_branches], DCIR.to_bus[active_branches], opf_input_param[:var_names][:Q_ki]; bounded = opf_input_param[:var_bounds][:Q_ki], min_lim=lower_bounds_circ[active_branches], max_lim=upper_bounds_circ[active_branches], encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_qki_lower, export_key_upper=:ineq_const_qki_upper)
        Q_ki = store_scalar_var_bounds!(opf_dict, :Q_ki, result, :ineq_const_qki_lower, :ineq_const_qki_upper)
    end
    
    # ========================================================================
    #            DEFINE THE OBJECTIVE FUNCTION OF THE MODEL
    # ========================================================================
    if opf_input_param[:obj_function][:type] == "quadratic"
        # Minimize total fuel cost
        obj_function, obj_function_MVA = obj_minimize_fuel_cost!(model, P_g, DGEN.id[active_gen], DGEN.g_cost_2[active_gen], DGEN.g_cost_1[active_gen], DGEN.g_cost_0[active_gen], base_MVA)

    elseif opf_input_param[:obj_function][:type] == "linear"
        # Minimize total marginal cost
        obj_function, obj_function_MVA = obj_minimize_marginal_cost!(model, P_g, DGEN.id[active_gen], DGEN.g_cost_1[active_gen], base_MVA)

    end
    
    # ========================================================================
    #            DEFINE THE CONSTRAINTS OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                  Equality and Inequality Constraints 
    #-------------------------------------------------------------------------

    # Mandatory equalities (branch-flow formulation)
    eq_const_angle_sw = const_angle_slack_bus!(model, θ, SW[1])

    # ****************************************
    # Constraint thermal limits of generators
    # Reproducing the capability curve
    # ****************************************
    if opf_input_param[:ineq_cons][:sg_upper]
        ineq_const_sg_upper = const_gen_power_apparent!(model, P_g, Q_g, DGEN.id[active_gen], DGEN.pg_max[active_gen] ./ base_MVA)
    end
    
    # Optional: promote the bus active demand to a pinned variable, so the pin's dual is
    # ∂L/∂P_d(k) directly (see `const_pin_load_active!`). `nothing` keeps the demand an
    # inlined Float64 and the model byte-identical to before the flag existed.
    P_d, eq_const_pd_pin = if get(opf_input_param, :pd_as_variable, false)
        const_pin_load_active!(model, bus_gen_circ_dict_ON, base_MVA)
    else
        nothing, nothing
    end

    eq_const_p_balance = ac_const_active_power_balance!(model, V, P_g, P_ik, P_ki, bus_gen_circ_dict_ON, DCIR, base_MVA; P_d = P_d)
    eq_const_q_balance = ac_const_reactive_power_balance!(model, V, Q_g, Q_ik, Q_ki, bus_gen_circ_dict_ON, DCIR, base_MVA)
    eq_const_p_ik = ac_const_powerflow_ik_active!(model, V, θ, P_ik, DCIR, active_branches)
    eq_const_q_ik = ac_const_powerflow_ik_reactive!(model, V, θ, Q_ik, DCIR, active_branches)
    eq_const_p_ki = ac_const_powerflow_ki_active!(model, V, θ, P_ki, DCIR, active_branches)
    eq_const_q_ki = ac_const_powerflow_ki_reactive!(model, V, θ, Q_ki, DCIR, active_branches)
    
    # **********************************************
    # Constraints for Thermal limits at the branches
    # **********************************************
    if opf_input_param[:ineq_cons][:sbranch_upper]
        ineq_const_s_ik, ineq_const_s_ki = ac_const_branch_thermal_limit!(model, P_ik, Q_ik, P_ki, Q_ki, upper_bounds_circ) # Dicts to save inequality constraints of apparent power flow
    end

    # *************************************************
    # Constraints for angle differences at the branches
    # *************************************************
    if opf_input_param[:ineq_cons][:ang_diff_branch]
        ineq_const_ang_diff_lower, ineq_const_ang_diff_upper = ac_const_angle_differences!(model, θ, pair_info, pair_circ_map; clamp_deg=limits.ang_diff_clamp_ac_deg) # Vector to save inequality constraints of angle difference between adjacent buses
    end
    
    # ========================================================================
    #           SAVING THE REST OF CONSTRAINTS OF THE MODEL
    # ========================================================================
    # Saving the constraints in the DICTIONARY
    begin
        # First, save the constraints that must be defined anyway
        opf_dict[:eq_const][:eq_const_angle_sw] = OrderedDict(1 => eq_const_angle_sw)
        opf_dict[:eq_const][:eq_const_p_balance] = eq_const_p_balance
        if P_d !== nothing
            opf_dict[:vars][:P_d] = P_d
            opf_dict[:eq_const][:eq_const_pd_pin] = eq_const_pd_pin
        end
        opf_dict[:eq_const][:eq_const_q_balance] = eq_const_q_balance
        opf_dict[:eq_const][:eq_const_p_ik] = eq_const_p_ik
        opf_dict[:eq_const][:eq_const_q_ik] = eq_const_q_ik
        opf_dict[:eq_const][:eq_const_p_ki] = eq_const_p_ki
        opf_dict[:eq_const][:eq_const_q_ki] = eq_const_q_ki

        # Then, verify whether the constraints are defined and save them
        if @isdefined(ineq_const_sg_upper)
            opf_dict[:ineq_const][:ineq_const_sg_upper] = ineq_const_sg_upper
        end

        if @isdefined(ineq_const_s_ik)
            opf_dict[:ineq_const][:ineq_const_s_ik] = ineq_const_s_ik
        end

        if @isdefined(ineq_const_s_ki)
            opf_dict[:ineq_const][:ineq_const_s_ki] = ineq_const_s_ki
        end

        if @isdefined(ineq_const_ang_diff_lower)
            opf_dict[:ineq_const][:ineq_const_ang_diff_lower] = ineq_const_ang_diff_lower
        end

        if @isdefined(ineq_const_ang_diff_upper)
            opf_dict[:ineq_const][:ineq_const_ang_diff_upper] = ineq_const_ang_diff_upper
        end
    end

    #-------------------------------------------------------------------------------------
    #                         SAVE MODEL SUMMARY AND DETAILS
    #-------------------------------------------------------------------------------------
    begin
        println("--------------------------------------------------------------------------------------------------------------------------------------")
        Export_OPF_Model(model, path_names, obj_function_MVA, opf_dict; model_label="AC-OPF")
        println("--------------------------------------------------------------------------------------------------------------------------------------")
    end

    return model, obj_function, obj_function_MVA, opf_dict
end

# Function to build the AC OPF model using the bus admittance matrix in the power balance constraints
function Make_ACOPF_Model_w_Ybus!(
    model::Model, 
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64,
    opf_input_param::OrderedDict{Symbol, Any};
    save_matrices::Bool=true,
)

    """
    This function builds a general AC OPF Model with the objective function, variables and constraints
    """
    #--------------------------------------
    # Calculating the bus admittance matrix
    #--------------------------------------
    Ybus = Calculate_Ybus_sparse(DBUS, DCIR, nBUS, nCIR, base_MVA) # Admittance matrix for the steady-state condition
    if save_matrices
        Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus", Ybus)
    end

    active_gen               = findall(x -> x == 1, DGEN.g_status)     # Indices of active generators
    active_branches          = findall(x -> x == 1, DCIR.l_status)     # Indices of active branches
    pair_info, pair_circ_map = get_angle_limits(DCIR, active_branches) # Limits for angular differences between adjacent buses
    limits = opf_input_param[:limits]                                  # DispatchLimitsConfig (θ box, clamps)

    # Rated branches: ±l_cap_1/base_MVA; unrated (l_cap_1 == 0): ±Inf → skipped at attach
    lower_bounds_circ, upper_bounds_circ = branch_flow_limit_vectors(DCIR, base_MVA)

    # Defining a dictionary to save opf related data
    opf_dict = init_opf_dict!(:Ybus => Ybus)
    enc = bound_encoding_from_param(opf_input_param)

    #------------------------------------------
    # Check if there is at least one swing bus
    #------------------------------------------
    SW = findall(x -> x == 3, DBUS.type)
    if isempty(SW) 
        throw(ArgumentError("You must define one bus as the SLACK BUS (type 3).")) 
    elseif length(SW) > 1 
        throw(ArgumentError("This code still does not support more than one SLACK BUS (type 3).")) 
    end

    # ========================================================================
    #              DEFINE THE PRIMAL VARIABLES OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                                Buses
    #-------------------------------------------------------------------------
    # Voltage Magnitude
    begin
        result = var_voltage_magnitude!(model, DBUS.bus, opf_input_param[:var_names][:V]; bounded = opf_input_param[:var_bounds][:V], min_lim = DBUS.v_min, max_lim = DBUS.v_max, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_volt_mag_lower, export_key_upper=:ineq_const_volt_mag_upper)
        V = store_scalar_var_bounds!(opf_dict, :V, result, :ineq_const_volt_mag_lower, :ineq_const_volt_mag_upper)
    end

    # Voltage Angle
    begin
        result = var_voltage_angle!(model, DBUS.bus, opf_input_param[:var_names][:θ]; bounded = opf_input_param[:var_bounds][:θ], min_lim=limits.θ_min_rad, max_lim=limits.θ_max_rad, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_volt_ang_lower, export_key_upper=:ineq_const_volt_ang_upper)
        θ = store_scalar_var_bounds!(opf_dict, :θ, result, :ineq_const_volt_ang_lower, :ineq_const_volt_ang_upper)
    end

    #-------------------------------------------------------------------------
    #                            Generators
    #-------------------------------------------------------------------------
    # Active power generated
    begin
        result = var_gen_power_active!(model, DGEN.id[active_gen], opf_input_param[:var_names][:P_g]; bounded = opf_input_param[:var_bounds][:P_g], min_lim=DGEN.pg_min[active_gen] ./ base_MVA, max_lim=DGEN.pg_max[active_gen] ./ base_MVA, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_pg_lower, export_key_upper=:ineq_const_pg_upper)
        P_g = store_scalar_var_bounds!(opf_dict, :P_g, result, :ineq_const_pg_lower, :ineq_const_pg_upper)
    end

    # Reactive power generated
    begin
        result = var_gen_power_reactive!(model, DGEN.id[active_gen], opf_input_param[:var_names][:Q_g]; bounded = opf_input_param[:var_bounds][:Q_g], min_lim=DGEN.qg_min[active_gen] ./ base_MVA, max_lim=DGEN.qg_max[active_gen] ./ base_MVA, encoding=enc, meta=opf_dict[:meta], export_key_lower=:ineq_const_qg_lower, export_key_upper=:ineq_const_qg_upper)
        Q_g = store_scalar_var_bounds!(opf_dict, :Q_g, result, :ineq_const_qg_lower, :ineq_const_qg_upper)
    end
    
    # ========================================================================
    #            DEFINE THE OBJECTIVE FUNCTION OF THE MODEL
    # ========================================================================
    if opf_input_param[:obj_function][:type] == "quadratic"
        # Minimize total fuel cost
        obj_function, obj_function_MVA = obj_minimize_fuel_cost!(model, P_g, DGEN.id[active_gen], DGEN.g_cost_2[active_gen], DGEN.g_cost_1[active_gen], DGEN.g_cost_0[active_gen], base_MVA)

    elseif opf_input_param[:obj_function][:type] == "linear"
        # Minimize total marginal cost
        obj_function, obj_function_MVA = obj_minimize_marginal_cost!(model, P_g, DGEN.id[active_gen], DGEN.g_cost_1[active_gen], base_MVA)

    end
    
    # ========================================================================
    #            DEFINE THE CONSTRAINTS OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                  Equality and Inequality Constraints 
    #-------------------------------------------------------------------------

    # Mandatory equalities (Ybus formulation — no explicit branch flow eq)
    eq_const_angle_sw = const_angle_slack_bus!(model, θ, SW[1])

    # ****************************************
    # Constraint thermal limits of generators
    # Reproducing the capability curve
    # ****************************************
    if opf_input_param[:ineq_cons][:sg_upper]
        ineq_const_sg_upper = const_gen_power_apparent!(model, P_g, Q_g, DGEN.id[active_gen], DGEN.pg_max[active_gen] ./ base_MVA)
    end 
    
    # Optional: promote the bus active demand to a pinned variable, so the pin's dual is
    # ∂L/∂P_d(k) directly (see `const_pin_load_active!`). `nothing` keeps the demand an
    # inlined Float64 and the model byte-identical to before the flag existed.
    P_d, eq_const_pd_pin = if get(opf_input_param, :pd_as_variable, false)
        const_pin_load_active!(model, bus_gen_circ_dict_ON, base_MVA)
    else
        nothing, nothing
    end

    eq_const_p_balance = ac_const_active_power_balance_Ybus!(model, Ybus, V, θ, P_g, bus_gen_circ_dict_ON, DCIR, base_MVA; P_d = P_d)
    eq_const_q_balance = ac_const_reactive_power_balance_Ybus!(model, Ybus, V, θ, Q_g, bus_gen_circ_dict_ON, DCIR, base_MVA)
    
    # **********************************************
    # Constraints for Thermal limits at the branches
    # **********************************************
    if opf_input_param[:ineq_cons][:sbranch_upper]
        ineq_const_s_ik, ineq_const_s_ki = ac_const_branch_thermal_limit_Ybus!(model, Ybus, V, θ, DCIR, upper_bounds_circ) # Dicts to save inequality constraints of apparent power flow
    end

    # *************************************************
    # Constraints for angle differences at the branches
    # *************************************************
    if opf_input_param[:ineq_cons][:ang_diff_branch]
        ineq_const_ang_diff_lower, ineq_const_ang_diff_upper = ac_const_angle_differences!(model, θ, pair_info, pair_circ_map; clamp_deg=limits.ang_diff_clamp_ac_deg) # Vector to save inequality constraints of angle difference between adjacent buses
    end
    
    # ========================================================================
    #           SAVING THE REST OF CONSTRAINTS OF THE MODEL
    # ========================================================================
    # Saving the constraints in the DICTIONARY
    begin
        # First, save the constraints that must be defined anyway
        opf_dict[:eq_const][:eq_const_angle_sw] = OrderedDict(1 => eq_const_angle_sw)
        opf_dict[:eq_const][:eq_const_p_balance] = eq_const_p_balance
        if P_d !== nothing
            opf_dict[:vars][:P_d] = P_d
            opf_dict[:eq_const][:eq_const_pd_pin] = eq_const_pd_pin
        end
        opf_dict[:eq_const][:eq_const_q_balance] = eq_const_q_balance

        # Then, verify whether the constraints are defined and save them
        if @isdefined(ineq_const_sg_upper)
            opf_dict[:ineq_const][:ineq_const_sg_upper] = ineq_const_sg_upper
        end

        if @isdefined(ineq_const_s_ik)
            opf_dict[:ineq_const][:ineq_const_s_ik] = ineq_const_s_ik
        end

        if @isdefined(ineq_const_s_ki)
            opf_dict[:ineq_const][:ineq_const_s_ki] = ineq_const_s_ki
        end

        if @isdefined(ineq_const_ang_diff_lower)
            opf_dict[:ineq_const][:ineq_const_ang_diff_lower] = ineq_const_ang_diff_lower
        end

        if @isdefined(ineq_const_ang_diff_upper)
            opf_dict[:ineq_const][:ineq_const_ang_diff_upper] = ineq_const_ang_diff_upper
        end
    end

    #-------------------------------------------------------------------------------------
    #                         SAVE MODEL SUMMARY AND DETAILS
    #-------------------------------------------------------------------------------------
    begin
        println("--------------------------------------------------------------------------------------------------------------------------------------")
        Export_OPF_Model(model, path_names, obj_function_MVA, opf_dict; model_label="AC-OPF (Ybus)")
        println("--------------------------------------------------------------------------------------------------------------------------------------")
    end

    return model, obj_function, obj_function_MVA, opf_dict

end