# Function to decide which model to use
function Build_SteadyState_Model!(
    path_names::OrderedDict{Symbol, String},
    model::Model,
    dispatch::DispatchConfig,
    DBUS::DataFrame,
    DGEN::DataFrame,
    DCIR::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nBUS::Int64,
    nGEN::Int64,
    nCIR::Int64;
    save_matrices::Bool=true,
)
    opf_input_param = build_opf_input_param(dispatch)
    type_model = dispatch.type_model
    use_matrix = dispatch.use_matrix

    if type_model == "ACOPF"
        if use_matrix
            # Build the optimization model
            model, obj_function, obj_function_MVA, opf_dict = Make_ACOPF_Model_w_Ybus!(
                model, path_names, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON,
                base_MVA, nBUS, nGEN, nCIR, opf_input_param; save_matrices=save_matrices)

        else
            # Branch formulation: Ybus is still useful for inspection when requested.
            if save_matrices
                Ybus = Calculate_Ybus_sparse(DBUS, DCIR, nBUS, nCIR, base_MVA)
                Save_Admittance_Matrix_XLSX(path_names[:pf_main], path_names[:pf_bus_matrices], "Ybus", Ybus)
            end
            # Build the optimization model
            model, obj_function, obj_function_MVA, opf_dict = Make_ACOPF_Model!(model, path_names, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR,opf_input_param)

        end

    elseif type_model == "DCOPF"
        if use_matrix
            # Build the optimization model
            model, obj_function, obj_function_MVA, opf_dict = Make_DCOPF_Model_w_Bbus!(
                model, path_names, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON,
                base_MVA, nBUS, nGEN, nCIR, opf_input_param; save_matrices=save_matrices)

        else
            # Build the optimization model
            model, obj_function, obj_function_MVA, opf_dict = Make_DCOPF_Model!(model, path_names, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR,opf_input_param)
            
        end

    elseif type_model == "ED"
            # Build the optimization model
            model, obj_function, obj_function_MVA, opf_dict = Make_ED_Model!(model, path_names, DBUS, DGEN, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, opf_input_param)
    elseif type_model == "UC"
            model, obj_function, obj_function_MVA, opf_dict = Make_UC_Model!(
                model, path_names, DBUS, DGEN, bus_gen_circ_dict_ON,
                base_MVA, nBUS, nGEN, opf_input_param)
    end
    

    return model, obj_function, obj_function_MVA, opf_dict
    
end
