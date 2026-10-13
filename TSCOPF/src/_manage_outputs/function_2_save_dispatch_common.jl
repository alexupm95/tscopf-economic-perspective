# ===================================================================================
#                  PRINT THE REPORTS IN TXT AND CSV
# ===================================================================================

# Function to save the solution of the Optimal Dispatch
function Save_Solution_Optimal_Dispatch(
    path_names::OrderedDict{Symbol, String},
    model::Model,
    type_model::String,
    use_matrix::Bool,
    obj_function::JuMP.AbstractJuMPScalar,
    opf_dict::OrderedDict{Symbol, Any},
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict;
    _save_duals::Bool=true
    )

    if type_model == "ACOPF"
        if use_matrix 
            RBUS, RGEN, RCIR = Save_Solution_ACOPF_Model_w_Ybus(path_names, obj_function, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping)

            if _save_duals Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; model=model, model_label="AC-OPF") end

        else
            RBUS, RGEN, RCIR = Save_Solution_ACOPF_Model(path_names, obj_function, opf_dict[:vars][:V], opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:Q_g], opf_dict[:vars][:P_ik], opf_dict[:vars][:Q_ik], opf_dict[:vars][:P_ki], opf_dict[:vars][:Q_ki], 
            bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping)

            if _save_duals Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; model=model, model_label="AC-OPF") end

        end
    elseif type_model == "DCOPF"
        if use_matrix 
            RBUS, RGEN, RCIR = Save_Solution_DCOPF_Model_w_Bbus(path_names, opf_dict[:Bbus], obj_function, opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping)

            if _save_duals Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; model=model, model_label="DC-OPF") end

        else
            RBUS, RGEN, RCIR = Save_Solution_DCOPF_Model(path_names, obj_function, opf_dict[:vars][:θ], opf_dict[:vars][:P_g],
            opf_dict[:vars][:P_ik], opf_dict[:vars][:P_ki], bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping)

            if _save_duals Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; model=model, model_label="DC-OPF") end

        end
        
    elseif type_model == "ED"
        RBUS, RGEN, RCIR = Save_Solution_ED_Model(path_names, obj_function, opf_dict[:vars][:P_g],
        bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping, reverse_bus_mapping)

        if _save_duals Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; model=model, model_label="ED") end

    elseif type_model == "UC"
        RBUS, RGEN, RCIR = Save_Solution_UC_Model(
            path_names, model, obj_function, opf_dict,
            bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR,
            bus_mapping, reverse_bus_mapping; _save_duals = _save_duals)

    end

    return RBUS, RGEN, RCIR
    
end