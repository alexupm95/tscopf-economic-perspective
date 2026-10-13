# Function to build the Economic Dispatch model
function Make_ED_Model!(
    model::Model, 
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    opf_input_param::OrderedDict{Symbol, Any}
    )

    """
    Build the Economic Dispatch (ED) model: minimise generation cost subject to
    active-power balance and generator output limits (no network constraints).
    """
    active_gen = findall(x -> x == 1, DGEN.g_status)     # Indices of active generators

    # Defining a dictionary to save opf related data
    opf_dict = init_opf_dict!()
    enc = bound_encoding_from_param(opf_input_param)

    # ========================================================================
    #              DEFINE THE PRIMAL VARIABLES OF THE MODEL
    # ========================================================================

    #-------------------------------------------------------------------------
    #                            Generators
    #-------------------------------------------------------------------------
    # Active power generated
    begin
        result = var_gen_power_active!(
            model, DGEN.id[active_gen], opf_input_param[:var_names][:P_g];
            bounded = opf_input_param[:var_bounds][:P_g],
            min_lim = DGEN.pg_min[active_gen] ./ base_MVA,
            max_lim = DGEN.pg_max[active_gen] ./ base_MVA,
            encoding = enc,
            meta = opf_dict[:meta],
            export_key_lower = :ineq_const_pg_lower,
            export_key_upper = :ineq_const_pg_upper,
        )
        P_g = store_scalar_var_bounds!(
            opf_dict, :P_g, result,
            :ineq_const_pg_lower, :ineq_const_pg_upper,
        )
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
    
    # Mandatory equality (ED)
    eq_const_p_balance = ed_const_power_balance!(model, P_g, bus_gen_circ_dict_ON, active_gen, base_MVA)
    
    # ========================================================================
    #           SAVING THE REST OF CONSTRAINTS OF THE MODEL
    # ========================================================================
    # Saving the constraints in the DICTIONARY
    begin
        # First, save the constraints that must be defined anyway
        opf_dict[:eq_const][:eq_const_p_balance] = eq_const_p_balance
    end

    #-------------------------------------------------------------------------------------
    #                         SAVE MODEL SUMMARY AND DETAILS
    #-------------------------------------------------------------------------------------
    begin
        println("--------------------------------------------------------------------------------------------------------------------------------------")
        Export_OPF_Model(model, path_names, obj_function_MVA, opf_dict; model_label="ED")
        println("--------------------------------------------------------------------------------------------------------------------------------------")
    end

    return model, obj_function, obj_function_MVA, opf_dict
end