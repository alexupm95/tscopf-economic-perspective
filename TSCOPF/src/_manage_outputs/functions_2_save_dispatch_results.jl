# ===========================================================
#                       AC-OPF
# ===========================================================
# Function to print and save the solution of the ACOPF model
function Save_Solution_ACOPF_Model(
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    P_ik::OrderedDict{Int, VariableRef}, 
    Q_ik::OrderedDict{Int, VariableRef}, 
    P_ki::OrderedDict{Int, VariableRef}, 
    Q_ki::OrderedDict{Int, VariableRef},
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict
    )

    obj_value = JuMP.value.(obj_function)

    println("===================================================")
    println("Objective Function: € "*string(round(obj_value, digits=2))*"")
    println("===================================================")


    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g
    Q_g_optim = [JuMP.value(v) for (i, v) in Q_g]   # Get the results of the optimization process -> Variable Q_g
    S_g_optim = abs.(P_g_optim .+ 1im .* Q_g_optim) # Calculate the output complex power of the generator
    V_optim   = [JuMP.value(v) for (i, v) in V]     # Get the results of the optimization process -> Variable V
    θ_optim   = [JuMP.value(v) for (i, v) in θ]     # Get the results of the optimization process -> Variable θ

    Pik_optim = [JuMP.value(v) for (i, v) in P_ik] # Get the active power flow from bus i to bus k
    Qik_optim = [JuMP.value(v) for (i, v) in Q_ik] # Get the reactive power flow from bus i to bus k
    Pki_optim = [JuMP.value(v) for (i, v) in P_ki] # Get the active power flow from bus k to bus i
    Qki_optim = [JuMP.value(v) for (i, v) in Q_ki] # Get the reactive power flow from bus k to bus i

    Sik_optim = abs.(Pik_optim .+ 1im .* Qik_optim) # Apparent power flow from bus i to bus k
    Ski_optim = abs.(Pki_optim .+ 1im .* Qki_optim) # Apparent power flow from bus k to bus i

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()
    Q_g_dict = Dict{Int, Float64}()
    S_g_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
            Q_g_dict[id] = Float64(Q_g_optim[aux_count])
            S_g_dict[id] = Float64(S_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
            Q_g_all[i] = Q_g_dict[id]  # Use optimized value
            S_g_all[i] = S_g_dict[id]
        end
    end

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power
    for i in 1:nGEN
        if Q_g_all[i] >= 0.0 # Check if the generator is providing reactive power to the system
            gen_loading_q[i] = DGEN.g_status[i] * (Q_g_all[i] / (DGEN.qg_max[i] / base_MVA))
        else # Or if the generator is consuming reactive power from the system
            gen_loading_q[i] = -DGEN.g_status[i] * (abs(Q_g_all[i]) / (abs(DGEN.qg_min[i]) / base_MVA))
        end
    end

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    Q_g_bus = zeros(Float64, nBUS) # Vector of reactive power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
            Q_g_bus[i] = sum(Q_g_all[indices_bus_gen])
        end
    end

    # =============================================================================
    #                                 Branches
    # =============================================================================
    # Initialize Vector for all Branches
    P_ik_all = zeros(Float64, nCIR)
    Q_ik_all = zeros(Float64, nCIR)
    S_ik_all = zeros(Float64, nCIR)
    P_ki_all = zeros(Float64, nCIR)
    Q_ki_all = zeros(Float64, nCIR)
    S_ki_all = zeros(Float64, nCIR)

    P_ik_dict = Dict{Int, Float64}()
    Q_ik_dict = Dict{Int, Float64}()
    S_ik_dict = Dict{Int, Float64}()
    P_ki_dict = Dict{Int, Float64}()
    Q_ki_dict = Dict{Int, Float64}()
    S_ki_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DCIR.id)
        if DCIR.l_status[i] == 1
            aux_count += 1
            P_ik_dict[id] = Float64(Pik_optim[aux_count])
            Q_ik_dict[id] = Float64(Qik_optim[aux_count])
            S_ik_dict[id] = Float64(Sik_optim[aux_count])
            P_ki_dict[id] = Float64(Pki_optim[aux_count])
            Q_ki_dict[id] = Float64(Qki_optim[aux_count])
            S_ki_dict[id] = Float64(Ski_optim[aux_count])
        end
    end

    for (i, id) in enumerate(DCIR.id)
        if haskey(P_ik_dict, id)
            P_ik_all[i] = P_ik_dict[id]  # Use optimized value
            Q_ik_all[i] = Q_ik_dict[id]  # Use optimized value
            S_ik_all[i] = S_ik_dict[id]
            P_ki_all[i] = P_ki_dict[id]  # Use optimized value
            Q_ki_all[i] = Q_ki_dict[id]  # Use optimized value
            S_ki_all[i] = S_ki_dict[id]
        end
    end

    Plosses = P_ik_all + P_ki_all  # Active power losses in the branches
    Qlosses = Q_ik_all + Q_ki_all  # Reactive power losses in the branches

    circ_loading = [DCIR.l_status[lin] * abs(max(S_ik_all[lin], S_ki_all[lin])) / (DCIR.l_cap_1[lin] / base_MVA) for lin in eachindex(DCIR.id)]

    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             # Bus identifies  
        v    = V_optim,                                             # Voltage magnitude                            [p.u.]     
        θ    = rad2deg.(θ_optim),                 # Voltage angle                                [deg]  
        p    = (P_g_bus .* base_MVA) .- DBUS.p_d, # Net Active power                             [MW] 
        q    = (Q_g_bus .* base_MVA) .- DBUS.q_d, # Net Reactive power                           [MVAr] 
        p_g  = P_g_bus .* base_MVA,               # Active power generated                       [MW]  
        q_g  = Q_g_bus .* base_MVA,               # Reactive power generated                     [MVAr]  
        p_d  = DBUS.p_d,                          # Active power generated                       [MW]  
        q_d  = DBUS.q_d,                          # Reactive power demanded by load              [MVAr]  
        p_sh = DBUS.g_sh .* (V_optim.^2),         # Active power demanded by shunt conductance   [MW]   
        q_sh = DBUS.b_sh .* (V_optim.^2)          # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        id      = DCIR.id,                              # Circuit identifier
        from_bus  = from_bus,                               # From bus identifier
        to_bus    = to_bus,                                 # To bus identifier
        p_ik      = P_ik_all .* base_MVA, # Circuit active power flow from i to k   [MW]
        q_ik      = Q_ik_all .* base_MVA, # Circuit reactive power flow from i to k [MVAr]
        s_ik      = S_ik_all .* base_MVA, # Circuit apparent power flow from i to k [MVA]
        p_ki      = P_ki_all .* base_MVA, # Circuit active power flow from k to i   [MW]
        q_ki      = Q_ki_all .* base_MVA, # Circuit reactive power flow from k to i [MVAr]
        s_ki      = S_ki_all .* base_MVA, # Circuit apparent power flow from k to i [MVA]
        p_losses  = Plosses .* base_MVA,  # Losses of active power                  [MW]
        q_losses  = Qlosses .* base_MVA,  # Losses of reactive power                [MVAr]
        s_cap     = DCIR.l_cap_1,                           # Circuit maximum power capacity          [MVA]
        loading   = circ_loading                            # Circuit loading
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id        = DGEN.id,                                 # Generator ID
        bus       = gen_bus,                                 # Bus in which the generator is connected
        p_g       = (P_g_all .* base_MVA), # Active power generated           [MW]
        q_g       = (Q_g_all .* base_MVA), # Reactive power generated         [MVAr]
        s_g       = (S_g_all .* base_MVA), # Apparent power generated         [MVA]
        loading_p = gen_loading_p, # Generator active power loading
        loading_q = gen_loading_q  # Generator reactive power loading
    )

    Save_ResultsTXT_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a TXT file
    Save_ResultsCSV_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a CSV file
    Save_ResultsXLSX_File(path_names, RBUS, RGEN, RCIR, Float64(obj_value))# Save the results in a XLSX file


    return RBUS, RGEN, RCIR

end

# Function to print and save the solution of the ACOPF modelusnig the bus admittance matrix
function Save_Solution_ACOPF_Model_w_Ybus(
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict
    )

    obj_value = JuMP.value.(obj_function)

    println("===================================================")
    println("Objective Function: € "*string(round(obj_value, digits=2))*"")
    println("===================================================")

    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g
    Q_g_optim = [JuMP.value(v) for (i, v) in Q_g]   # Get the results of the optimization process -> Variable Q_g
    S_g_optim = abs.(P_g_optim .+ 1im .* Q_g_optim) # Calculate the output complex power of the generator
    V_optim   = [JuMP.value(v) for (i, v) in V]     # Get the results of the optimization process -> Variable V
    θ_optim   = [JuMP.value(v) for (i, v) in θ]     # Get the results of the optimization process -> Variable θ

    P_ik, Q_ik, S_ik, P_ki, Q_ki, S_ki, Plosses, Qlosses, circ_loading = Calculate_AC_Power_Flow(DCIR, nCIR, V_optim, θ_optim, base_MVA)
    # P_ik, Q_ik, S_ik, P_ki, Q_ki, S_ki, Plosses, Qlosses, circ_loading = Calculate_AC_Power_Flow(Ybus, DCIR, nCIR, V_optim, θ_optim, base_MVA)

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()
    Q_g_dict = Dict{Int, Float64}()
    S_g_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
            Q_g_dict[id] = Float64(Q_g_optim[aux_count])
            S_g_dict[id] = Float64(S_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
            Q_g_all[i] = Q_g_dict[id]  # Use optimized value
            S_g_all[i] = S_g_dict[id]
        end
    end

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power
    for i in 1:nGEN
        if Q_g_all[i] >= 0.0 # Check if the generator is providing reactive power to the system
            gen_loading_q[i] = DGEN.g_status[i] * (Q_g_all[i] / (DGEN.qg_max[i] / base_MVA))
        else # Or if the generator is consuming reactive power from the system
            gen_loading_q[i] = -DGEN.g_status[i] * (abs(Q_g_all[i]) / (abs(DGEN.qg_min[i]) / base_MVA))
        end
    end

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    Q_g_bus = zeros(Float64, nBUS) # Vector of reactive power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
            Q_g_bus[i] = sum(Q_g_all[indices_bus_gen])
        end
    end

    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             # Bus identifies  
        v    = V_optim,                                             # Voltage magnitude                            [p.u.]     
        θ    = rad2deg.(θ_optim),                 # Voltage angle                                [deg]  
        p    = (P_g_bus .* base_MVA) .- DBUS.p_d, # Net Active power                             [MW] 
        q    = (Q_g_bus .* base_MVA) .- DBUS.q_d, # Net Reactive power                           [MVAr] 
        p_g  = P_g_bus .* base_MVA,               # Active power generated                       [MW]  
        q_g  = Q_g_bus .* base_MVA,               # Reactive power generated                     [MVAr]  
        p_d  = DBUS.p_d,                          # Active power generated                       [MW]  
        q_d  = DBUS.q_d,                          # Reactive power demanded by load              [MVAr]  
        p_sh = DBUS.g_sh .* (V_optim.^2),         # Active power demanded by shunt conductance   [MW]   
        q_sh = DBUS.b_sh .* (V_optim.^2)          # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        id      = DCIR.id,                            # Circuit identifier
        from_bus  = from_bus,                           # From bus identifier
        to_bus    = to_bus,                             # To bus identifier
        p_ik      = P_ik .* base_MVA, # Circuit active power flow from i to k   [MW]
        q_ik      = Q_ik .* base_MVA, # Circuit reactive power flow from i to k [MVAr]
        s_ik      = S_ik .* base_MVA, # Circuit apparent power flow from i to k [MVA]
        p_ki      = P_ki .* base_MVA, # Circuit active power flow from k to i   [MW]
        q_ki      = Q_ki .* base_MVA, # Circuit reactive power flow from k to i [MVAr]
        s_ki      = S_ki .* base_MVA, # Circuit apparent power flow from k to i [MVA]
        p_losses  = Plosses .* base_MVA,  # Losses of active power                  [MW]
        q_losses  = Qlosses .* base_MVA,  # Losses of reactive power                [MVAr]
        s_cap     = DCIR.l_cap_1,                           # Circuit maximum power capacity          [MVA]
        loading   = circ_loading                            # Circuit loading
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id        = DGEN.id,                                 # Generator ID
        bus       = gen_bus,                                 # Bus in which the generator is connected
        p_g       = (P_g_all .* base_MVA), # Active power generated           [MW]
        q_g       = (Q_g_all .* base_MVA), # Reactive power generated         [MVAr]
        s_g       = (S_g_all .* base_MVA), # Apparent power generated         [MVA]
        loading_p = gen_loading_p, # Generator active power loading
        loading_q = gen_loading_q  # Generator reactive power loading
    )

    Save_ResultsTXT_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a TXT file
    Save_ResultsCSV_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a CSV file
    Save_ResultsXLSX_File(path_names, RBUS, RGEN, RCIR, Float64(obj_value))# Save the results in a XLSX file


    return RBUS, RGEN, RCIR
    
end

# ===========================================================
#                       DC-OPF
# ===========================================================
# Function to print and save the solution of the DCOPF model
function Save_Solution_DCOPF_Model(
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    P_ik::OrderedDict{Int, VariableRef}, 
    P_ki::OrderedDict{Int, VariableRef}, 
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict
    )

    obj_value = JuMP.value.(obj_function)

    println("===================================================")
    println("Objective Function: € "*string(round(obj_value, digits=2))*"")
    println("===================================================")


    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g
    θ_optim   = [JuMP.value(v) for (i, v) in θ]     # Get the results of the optimization process -> Variable θ
    Pik_optim = [JuMP.value(v) for (i, v) in P_ik] # Get the active power flow from bus i to bus k
    Pki_optim = [JuMP.value(v) for (i, v) in P_ki] # Get the active power flow from bus k to bus i

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()

    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
        end
    end
    S_g_all = P_g_all

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    Q_g_bus = zeros(Float64, nBUS) # Vector of reactive power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
        end
    end

    # =============================================================================
    #                                 Branches
    # =============================================================================
    # Initialize Vector for all Branches
    P_ik_all = zeros(Float64, nCIR)
    Q_ik_all = zeros(Float64, nCIR)
    S_ik_all = zeros(Float64, nCIR)
    P_ki_all = zeros(Float64, nCIR)
    Q_ki_all = zeros(Float64, nCIR)
    S_ki_all = zeros(Float64, nCIR)

    P_ik_dict = Dict{Int, Float64}()
    P_ki_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DCIR.id)
        if DCIR.l_status[i] == 1
            aux_count += 1
            P_ik_dict[id] = Float64(Pik_optim[aux_count])
            P_ki_dict[id] = Float64(Pki_optim[aux_count])
        end
    end

    for (i, id) in enumerate(DCIR.id)
        if haskey(P_ik_dict, id)
            P_ik_all[i] = P_ik_dict[id]  # Use optimized value
            P_ki_all[i] = P_ki_dict[id]  # Use optimized value
        end
    end
    S_ik_all = P_ik_all
    S_ki_all = P_ki_all

    Plosses = zeros(Float64, nCIR)  # Active power losses in the branches
    Qlosses = zeros(Float64, nCIR)  # Reactive power losses in the branches

    circ_loading = [DCIR.l_status[lin] * S_ik_all[lin] / (DCIR.l_cap_1[lin] / base_MVA) for lin in eachindex(DCIR.id)]

    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             # Bus identifies  
        v    = ones(Float64, nBUS),                                 # Voltage magnitude                            [p.u.]     
        θ    = rad2deg.(θ_optim),                 # Voltage angle                                [deg]  
        p    = (P_g_bus .* base_MVA) .- DBUS.p_d, # Net Active power                             [MW] 
        q    = (Q_g_bus .* base_MVA) .- DBUS.q_d, # Net Reactive power                           [MVAr] 
        p_g  = P_g_bus .* base_MVA,               # Active power generated                       [MW]  
        q_g  = Q_g_bus .* base_MVA,               # Reactive power generated                     [MVAr]  
        p_d  = DBUS.p_d,                          # Active power generated                       [MW]  
        q_d  = DBUS.q_d,                          # Reactive power demanded by load              [MVAr]  
        p_sh = zeros(Float64, nBUS),                                # Active power demanded by shunt conductance   [MW]   
        q_sh = zeros(Float64, nBUS)                                 # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        id      = DCIR.id,                              # Circuit identifier
        from_bus  = from_bus,                               # From bus identifier
        to_bus    = to_bus,                                 # To bus identifier
        p_ik      = P_ik_all .* base_MVA, # Circuit active power flow from i to k   [MW]
        q_ik      = Q_ik_all .* base_MVA, # Circuit reactive power flow from i to k [MVAr]
        s_ik      = S_ik_all .* base_MVA, # Circuit apparent power flow from i to k [MVA]
        p_ki      = P_ki_all .* base_MVA, # Circuit active power flow from k to i   [MW]
        q_ki      = Q_ki_all .* base_MVA, # Circuit reactive power flow from k to i [MVAr]
        s_ki      = S_ki_all .* base_MVA, # Circuit apparent power flow from k to i [MVA]
        p_losses  = Plosses .* base_MVA,  # Losses of active power                  [MW]
        q_losses  = Qlosses .* base_MVA,  # Losses of reactive power                [MVAr]
        s_cap     = DCIR.l_cap_1,                           # Circuit maximum power capacity          [MVA]
        loading   = circ_loading                            # Circuit loading
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id        = DGEN.id,                                 # Generator ID
        bus       = gen_bus,                                 # Bus in which the generator is connected
        p_g       = (P_g_all .* base_MVA), # Active power generated           [MW]
        q_g       = (Q_g_all .* base_MVA), # Reactive power generated         [MVAr]
        s_g       = (S_g_all .* base_MVA), # Apparent power generated         [MVA]
        loading_p = gen_loading_p, # Generator active power loading
        loading_q = gen_loading_q  # Generator reactive power loading
    )

    Save_ResultsTXT_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a TXT file
    Save_ResultsCSV_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a CSV file
    Save_ResultsXLSX_File(path_names, RBUS, RGEN, RCIR, Float64(obj_value))# Save the results in a XLSX file


    return RBUS, RGEN, RCIR

end

# Function to print and save the solution of the DCOPF modelusnig the bus susceptance matrix
function Save_Solution_DCOPF_Model_w_Bbus(
    path_names::OrderedDict{Symbol, String},
    Bbus::SparseMatrixCSC,
    obj_function::JuMP.AbstractJuMPScalar,
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict
    )

    obj_value = JuMP.value.(obj_function)

    println("===================================================")
    println("Objective Function: € "*string(round(obj_value, digits=2))*"")
    println("===================================================")

    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g
    θ_optim   = [JuMP.value(v) for (i, v) in θ]     # Get the results of the optimization process -> Variable θ

    # Recompute branch flows from the same Bbus that built the model, so the report honours
    # the selected susceptance convention (SIMPLE vs POWERMODELS) instead of assuming 1/x.
    P_ik, Q_ik, S_ik, P_ki, Q_ki, S_ki, Plosses, Qlosses, circ_loading = Calculate_DC_Power_Flow(Bbus, DCIR, nCIR, θ_optim, base_MVA)

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
        end
    end
    S_g_all = P_g_all

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    Q_g_bus = zeros(Float64, nBUS) # Vector of reactive power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
        end
    end

    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             # Bus identifies  
        v    = ones(Float64, nBUS),                                             # Voltage magnitude                            [p.u.]     
        θ    = rad2deg.(θ_optim),                 # Voltage angle                                [deg]  
        p    = (P_g_bus .* base_MVA) .- DBUS.p_d, # Net Active power                             [MW] 
        q    = (Q_g_bus .* base_MVA) .- DBUS.q_d, # Net Reactive power                           [MVAr] 
        p_g  = P_g_bus .* base_MVA,               # Active power generated                       [MW]  
        q_g  = Q_g_bus .* base_MVA,               # Reactive power generated                     [MVAr]  
        p_d  = DBUS.p_d,                          # Active power generated                       [MW]  
        q_d  = DBUS.q_d,                          # Reactive power demanded by load              [MVAr]  
        p_sh = zeros(Float64, nBUS),                                # Active power demanded by shunt conductance   [MW]   
        q_sh = zeros(Float64, nBUS)                                 # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        id      = DCIR.id,                            # Circuit identifier
        from_bus  = from_bus,                           # From bus identifier
        to_bus    = to_bus,                             # To bus identifier
        p_ik      = P_ik .* base_MVA, # Circuit active power flow from i to k   [MW]
        q_ik      = Q_ik .* base_MVA, # Circuit reactive power flow from i to k [MVAr]
        s_ik      = S_ik .* base_MVA, # Circuit apparent power flow from i to k [MVA]
        p_ki      = P_ki .* base_MVA, # Circuit active power flow from k to i   [MW]
        q_ki      = Q_ki .* base_MVA, # Circuit reactive power flow from k to i [MVAr]
        s_ki      = S_ki .* base_MVA, # Circuit apparent power flow from k to i [MVA]
        p_losses  = Plosses .* base_MVA,  # Losses of active power                  [MW]
        q_losses  = Qlosses .* base_MVA,  # Losses of reactive power                [MVAr]
        s_cap     = DCIR.l_cap_1,                           # Circuit maximum power capacity          [MVA]
        loading   = circ_loading                            # Circuit loading
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id        = DGEN.id,                                 # Generator ID
        bus       = gen_bus,                                 # Bus in which the generator is connected
        p_g       = (P_g_all .* base_MVA), # Active power generated           [MW]
        q_g       = (Q_g_all .* base_MVA), # Reactive power generated         [MVAr]
        s_g       = (S_g_all .* base_MVA), # Apparent power generated         [MVA]
        loading_p = gen_loading_p, # Generator active power loading
        loading_q = gen_loading_q  # Generator reactive power loading
    )

    Save_ResultsTXT_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a TXT file
    Save_ResultsCSV_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a CSV file
    Save_ResultsXLSX_File(path_names, RBUS, RGEN, RCIR, Float64(obj_value))# Save the results in a XLSX file


    return RBUS, RGEN, RCIR
    
end

# ===========================================================
#                   Economic Dispatch
# ===========================================================
# Function to print and save the solution of the DCOPF model
function Save_Solution_ED_Model(
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    P_g::OrderedDict{Int, VariableRef},  
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict
    )

    obj_value = JuMP.value.(obj_function)

    println("===================================================")
    println("Objective Function: € "*string(round(obj_value, digits=2))*"")
    println("===================================================")


    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()

    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
        end
    end
    S_g_all = P_g_all

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
        end
    end


    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # Struct to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                               # Bus identifies  
        v    = zeros(Float64, nBUS),                                  # Voltage magnitude                            [p.u.]     
        θ    = zeros(Float64, nBUS),                                  # Voltage angle                                [deg]  
        p    = ((P_g_bus .* base_MVA) .- DBUS.p_d), # Net Active power                             [MW] 
        q    = zeros(Float64, nBUS),                                  # Net Reactive power                           [MVAr] 
        p_g  = (P_g_bus .* base_MVA),               # Active power generated                       [MW]  
        q_g  = zeros(Float64, nBUS),                                  # Reactive power generated                     [MVAr]  
        p_d  = DBUS.p_d,                            # Active power generated                       [MW]  
        q_d  = zeros(Float64, nBUS),                                  # Reactive power demanded by load              [MVAr]  
        p_sh = zeros(Float64, nBUS),                                  # Active power demanded by shunt conductance   [MW]   
        q_sh = zeros(Float64, nBUS)                                   # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # Struct to save the results related to the circuits
    RCIR = DataFrame(
        id        = DCIR.id,              # Circuit identifier
        from_bus  = from_bus,             # From bus identifier
        to_bus    = to_bus,               # To bus identifier
        p_ik      = zeros(Float64, nCIR), # Circuit active power flow from i to k   [MW]
        q_ik      = zeros(Float64, nCIR), # Circuit reactive power flow from i to k [MVAr]
        s_ik      = zeros(Float64, nCIR), # Circuit apparent power flow from i to k [MVA]
        p_ki      = zeros(Float64, nCIR), # Circuit active power flow from k to i   [MW]
        q_ki      = zeros(Float64, nCIR), # Circuit reactive power flow from k to i [MVAr]
        s_ki      = zeros(Float64, nCIR), # Circuit apparent power flow from k to i [MVA]
        p_losses  = zeros(Float64, nCIR), # Losses of active power                  [MW]
        q_losses  = zeros(Float64, nCIR), # Losses of reactive power                [MVAr]
        s_cap     = zeros(Float64, nCIR), # Circuit maximum power capacity          [MVA]
        loading   = zeros(Float64, nCIR)  # Circuit loading
    ) 

    # Struct to save the results related to the generators
    RGEN = DataFrame(
        id        = DGEN.id,                                 # Generator ID
        bus       = gen_bus,                                 # Bus in which the generator is connected
        p_g       = (P_g_all .* base_MVA), # Active power generated           [MW]
        q_g       = (Q_g_all .* base_MVA), # Reactive power generated         [MVAr]
        s_g       = (S_g_all .* base_MVA), # Apparent power generated         [MVA]
        loading_p = gen_loading_p, # Generator active power loading
        loading_q = gen_loading_q  # Generator reactive power loading
    )

    Save_ResultsTXT_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a TXT file
    Save_ResultsCSV_Files(path_names, RBUS, RGEN, RCIR, Float64(obj_value)) # Save the results in a CSV file
    Save_ResultsXLSX_File(path_names, RBUS, RGEN, RCIR, Float64(obj_value))# Save the results in a XLSX file


    return RBUS, RGEN, RCIR

end

# ===========================================================
#                   Unit Commitment
# ===========================================================

# Re-write generator CSV/XLSX with u_commit (ED saver has no binary variables).
function _save_uc_generator_exports!(
    path_names::OrderedDict{Symbol, String},
    RGEN::DataFrame,
    RBUS::DataFrame,
    RCIR::DataFrame,
    objective::Float64,
)
    df_generators = DataFrame(
        ID = RGEN.id, BUS = RGEN.bus, u_commit = RGEN.u_commit,
        P_MW = RGEN.p_g, Q_MVAr = RGEN.q_g, S_MVA = RGEN.s_g,
        Loading_P = RGEN.loading_p, Loading_Q = RGEN.loading_q,
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "generators_report.csv"), df_generators; delim = ';')

    file_path = joinpath(path_names[:pf_dispatch], "OPF_Dispatch_Results.xlsx")
    df_buses = DataFrame(
        BUS = RBUS.bus, V_pu = RBUS.v, Theta_deg = RBUS.θ,
        P_MW = RBUS.p, Q_MVAr = RBUS.q, PG_MW = RBUS.p_g, QG_MVAr = RBUS.q_g,
        PD_MW = RBUS.p_d, QD_MVAr = RBUS.q_d, Psh_MW = RBUS.p_sh, Qsh_MVAr = RBUS.q_sh,
    )
    df_circuits = DataFrame(
        ID_CIRC = RCIR.id, FROM_BUS = RCIR.from_bus, TO_BUS = RCIR.to_bus,
        Pik_MW = RCIR.p_ik, Qik_MVAr = RCIR.q_ik, Sik_MVA = RCIR.s_ik,
        Pki_MW = RCIR.p_ki, Qki_MVAr = RCIR.q_ki, Ski_MVA = RCIR.s_ki,
        Cap_MVA = RCIR.s_cap, Loading = RCIR.loading,
        Ploss_MW = RCIR.p_losses, Qloss_MVAr = RCIR.q_losses,
    )
    df_optimization = DataFrame(Metric = ["Total Cost (Euros)"], Value = [objective])
    df_uc = DataFrame(Gen_ID = RGEN.id, BUS = RGEN.bus, u_commit = RGEN.u_commit, P_MW = RGEN.p_g)
    XLSX.writetable(file_path,
        "Buses" => df_buses, "Generators" => df_generators, "Circuits" => df_circuits,
        "Optimization" => df_optimization, "UC_Commitment" => df_uc,
        overwrite = true,
    )
end

function _save_uc_primal_snapshot!(
    path_names::OrderedDict{Symbol, String},
    opf_dict::OrderedDict{Symbol, Any},
    objective::Float64,
)
    u = opf_dict[:vars][:u_commit]
    P_g = opf_dict[:vars][:P_g]
    open(joinpath(path_names[:pf_dispatch], "uc_primal_solution.txt"), "w") do io
        println(io, "=== UC MILP primal solution ===")
        println(io, "objective (Euros): ", objective)
        println(io)
        for g in sort(collect(keys(u)))
            println(io, "gen $g: u = ", round(Int, JuMP.value(u[g])),
                ", P_g = ", JuMP.value(P_g[g]), " pu")
        end
    end
end

# Primal: MILP dispatch + commitment tables. Duals: restricted LP with u* fixed (not MILP KKT).
function Save_Solution_UC_Model(
    path_names::OrderedDict{Symbol, String},
    model::Model,
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
    _save_duals::Bool = true,
)
    mkpath(path_names[:pf_dispatch])
    mkpath(path_names[:pf_dispatch_CSV])
    mkpath(path_names[:pf_dispatch_CSV_duals])

    P_g = opf_dict[:vars][:P_g]
    u = opf_dict[:vars][:u_commit]
    obj_value = JuMP.objective_value(model)

    # ED-style TXT/CSV/XLSX for P_g (duals skipped here; UC path handles them below).
    RBUS, RGEN, RCIR = Save_Solution_ED_Model(
        path_names, obj_function, P_g,
        bus_gen_circ_dict, DBUS, DGEN, DCIR, base_MVA, nBUS, nGEN, nCIR,
        bus_mapping, reverse_bus_mapping)

    @assert all(haskey(u, g) for g in RGEN.id) "UC: u_commit is missing generator(s) present in RGEN.id; cannot align commitment with report rows."
    u_vals = [round(Int, JuMP.value(u[g])) for g in RGEN.id]
    RGEN.u_commit = u_vals
    _save_uc_generator_exports!(path_names, RGEN, RBUS, RCIR, obj_value)
    _save_uc_primal_snapshot!(path_names, opf_dict, obj_value)

    CSV.write(
        joinpath(path_names[:pf_dispatch_CSV], "uc_commitment.csv"),
        DataFrame(Gen_ID = RGEN.id, u_commit = u_vals);
        delim = ';',
    )

    if _save_duals
        u_star = uc_commitment_values(u)
        lp_model, lp_dict = Solve_UC_Restricted_Pricing!(
            DGEN, bus_gen_circ_dict, base_MVA, u_star; silent_solver = true)
        Export_UC_Restricted_LP_Model!(path_names, lp_model, lp_dict, u_star)
        Save_Duals_UC_Model(path_names, lp_dict, base_MVA, u_star)
    end

    return RBUS, RGEN, RCIR
end

# ===========================================================
#                       Common
# ===========================================================
# Save the reports of the power flow in TXT files
function Save_ResultsTXT_Files(
    path_names::OrderedDict{Symbol, String},
    RBUS::DataFrame,
    RGEN::DataFrame,
    RCIR::DataFrame,
    objective::Float64
    )

    # Summation of power generated and demmanded
    sum_Pg = 0.0
    sum_Qg = 0.0
    sum_Pd = 0.0
    sum_Qd = 0.0
    sum_Psh = 0.0
    sum_Qsh = 0.0
    sum_losses_P = 0.0
    sum_losses_Q = 0.0
    sum_gen_Pg = 0.0
    sum_gen_Qg = 0.0
    sum_gen_Sg = 0.0

    io = open(joinpath(path_names[:pf_dispatch], "buses_report.txt"), "w")
    @printf(io, "BUSES REPORT\n")
    @printf(io, "============================================================================================================================= \n")
    @printf(io, "   BUS    V (pu)     Θ (º)     P (MW)     Q (MVAr)      PG (MW)    QG (MVAr)    PD (MW)   QD (MVAr)     Psh (MW)   Qsh (MVAr) \n")
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------- \n")
    for i in eachindex(RBUS.bus)
        @printf(io, " %4d    %6.4f    %6.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n", RBUS.bus[i], RBUS.v[i], RBUS.θ[i], RBUS.p[i], RBUS.q[i], RBUS.p_g[i], RBUS.q_g[i], RBUS.p_d[i], RBUS.q_d[i], RBUS.p_sh[i], RBUS.q_sh[i])
        sum_Pg += RBUS.p_g[i]
        sum_Qg += RBUS.q_g[i]
        sum_Pd += RBUS.p_d[i]
        sum_Qd += RBUS.q_d[i]
        sum_Psh += RBUS.p_sh[i]
        sum_Qsh += RBUS.q_sh[i]
    end
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                               %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f\n", sum_Pg, sum_Qg, sum_Pd , sum_Qd, sum_Psh , sum_Qsh)
    @printf(io, "============================================================================================================================= \n")
    @printf(io, "\n")
    close(io)  

    io = open(joinpath(path_names[:pf_dispatch], "generators_report.txt"), "w")
    @printf(io, "GENERATORS REPORT\n")
    @printf(io, "========================================================================== \n")
    @printf(io, "   ID     BUS     P_g (MW)  Q_g (MVAr)   S_g (MVA)   Loading_P   Loading_Q \n")
    @printf(io, "-------------------------------------------------------------------------- \n")
    for i in eachindex(RGEN.id)
        @printf(io, " %4d   %4d    %8.2f    %8.2f    %8.2f    %8.4f    %8.4f   \n", RGEN.id[i], RGEN.bus[i], RGEN.p_g[i], RGEN.q_g[i], RGEN.s_g[i], RGEN.loading_p[i], RGEN.loading_q[i])
        sum_gen_Pg += RGEN.p_g[i]
        sum_gen_Qg += RGEN.q_g[i]
        sum_gen_Sg += RGEN.s_g[i]
    end
    @printf(io, "-------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:         %8.4f   %8.4f   %8.4f\n", round(sum_gen_Pg, digits=3), round(sum_gen_Qg, digits=3), round(sum_gen_Sg, digits=3))
    @printf(io, "========================================================================== \n")
    close(io)
    
    io = open(joinpath(path_names[:pf_dispatch], "circuits_report.txt"), "w")
    @printf(io, "CIRCUITS REPORT\n")
    @printf(io, "=============================================================================================================================================== \n")
    @printf(io, "   ID     FROM    TO      Pik (MW)  Qik (MVAr)   Sik (MVA)   Pki (MW)  Qki (MVAr)    Ski (MVA)  Cap (MVA)    Loading    Ploss (MW)  Qloss (MVAr)\n")
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    for i in eachindex(RCIR.id)
        @printf(io, " %4d   %4d   %4d    %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f   %8.4f     %8.4f     %8.3f     %8.3f\n", RCIR.id[i], RCIR.from_bus[i], RCIR.to_bus[i], RCIR.p_ik[i], RCIR.q_ik[i], RCIR.s_ik[i], RCIR.p_ki[i], RCIR.q_ki[i], RCIR.s_ki[i], RCIR.s_cap[i], RCIR.loading[i], RCIR.p_losses[i], RCIR.q_losses[i])
        sum_losses_P += RCIR.p_losses[i]
        sum_losses_Q += RCIR.q_losses[i]
    end
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                                                                                                  %8.3f     %8.3f\n", sum_losses_P , sum_losses_Q)
    @printf(io, "=============================================================================================================================================== \n")
    close(io)
    
    io = open(joinpath(path_names[:pf_dispatch], "optimization_report.txt"), "w")
    @printf(io, "OBJECTIVE\n")
    @printf(io, "================================ \n")
    @printf(io, "Total cost: (Euros) %8.2f \n", round(objective, digits = 2))
    @printf(io, "================================ \n")
    close(io)

    println("OPF results successfully saved as TXT files in: ", path_names[:pf_dispatch])


end

# Save the reports of the power flow in CSV files
function Save_ResultsCSV_Files(
    path_names::OrderedDict{Symbol, String},
    RBUS::DataFrame,
    RGEN::DataFrame,
    RCIR::DataFrame,
    objective::Float64
    )

    # Save Bus Report as CSV
    df_buses = DataFrame(
        BUS = RBUS.bus,
        V_pu = RBUS.v,
        Theta_deg = RBUS.θ,
        P_MW = RBUS.p,
        Q_MVAr = RBUS.q,
        PG_MW = RBUS.p_g,
        QG_MVAr = RBUS.q_g,
        PD_MW = RBUS.p_d,
        QD_MVAr = RBUS.q_d,
        Psh_MW =  RBUS.p_sh,
        Qsh_MVAr =  RBUS.q_sh
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "buses_report.csv"), df_buses; delim=';', header=true)

    # Save Generators Report as CSV
    df_generators = DataFrame(
        ID = RGEN.id,
        BUS = RGEN.bus,
        P_MW = RGEN.p_g,
        Q_MVAr = RGEN.q_g,
        S_MVA = RGEN.s_g,
        Loading_P = RGEN.loading_p,
        Loading_Q = RGEN.loading_q
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "generators_report.csv"), df_generators; delim=';', header=true)

    # Save Circuit Report as CSV
    df_circuits = DataFrame(
        ID_CIRC = RCIR.id,
        FROM_BUS = RCIR.from_bus,
        TO_BUS = RCIR.to_bus,
        Pik_MW = RCIR.p_ik,
        Qik_MVAr = RCIR.q_ik,
        Sik_MVA = RCIR.s_ik,
        Pki_MW = RCIR.p_ki,
        Qki_MVAr = RCIR.q_ki,
        Ski_MVA = RCIR.s_ki,
        Cap_MVA = RCIR.s_cap,
        Loading = RCIR.loading,
        Ploss_MW = RCIR.p_losses,
        Qloss_MVAr = RCIR.q_losses
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "circuits_report.csv"), df_circuits; delim=';', header=true)

    # Save Optimization Report as CSV
    df_optimization = DataFrame(
        Metric = ["Total Cost (Euros)"],
        Value = [objective]
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "optimization_report.csv"), df_optimization; delim=';', header=true)

    println("OPF results successfully saved as CSV files in: ", path_names[:pf_dispatch_CSV])

end

# Save the reports of the power flow in XLSX files
function Save_ResultsXLSX_File(
    path_names::OrderedDict{Symbol, String},
    RBUS::DataFrame,
    RGEN::DataFrame,
    RCIR::DataFrame,
    objective::Float64
    )

    # Define the output path (saved in the Dispatch folder)
    file_path = joinpath(path_names[:pf_dispatch], "OPF_Dispatch_Results.xlsx")
  
    # 1. Prepare Bus Report DataFrame
    df_buses = DataFrame(
        BUS       = RBUS.bus,
        V_pu      = RBUS.v,
        Theta_deg = RBUS.θ,
        P_MW      = RBUS.p,
        Q_MVAr    = RBUS.q,
        PG_MW     = RBUS.p_g,
        QG_MVAr   = RBUS.q_g,
        PD_MW     = RBUS.p_d,
        QD_MVAr   = RBUS.q_d,
        Psh_MW    = RBUS.p_sh,
        Qsh_MVAr  = RBUS.q_sh
    )

    # 2. Prepare Generators Report DataFrame
    df_generators = DataFrame(
        ID        = RGEN.id,
        BUS       = RGEN.bus,
        P_MW      = RGEN.p_g,
        Q_MVAr    = RGEN.q_g,
        S_MVA     = RGEN.s_g,
        Loading_P = RGEN.loading_p,
        Loading_Q = RGEN.loading_q
    )

    # 3. Prepare Circuit Report DataFrame
    df_circuits = DataFrame(
        ID_CIRC    = RCIR.id,
        FROM_BUS   = RCIR.from_bus,
        TO_BUS     = RCIR.to_bus,
        Pik_MW     = RCIR.p_ik,
        Qik_MVAr   = RCIR.q_ik,
        Sik_MVA    = RCIR.s_ik,
        Pki_MW     = RCIR.p_ki,
        Qki_MVAr   = RCIR.q_ki,
        Ski_MVA    = RCIR.s_ki,
        Cap_MVA    = RCIR.s_cap,
        Loading    = RCIR.loading,
        Ploss_MW   = RCIR.p_losses,
        Qloss_MVAr = RCIR.q_losses
    )

    # 4. Prepare Optimization Summary DataFrame
    df_optimization = DataFrame(
        Metric = ["Total Cost (Euros)"],
        Value  = [objective]
    )

    # Save all DataFrames as separate sheets in one file
    XLSX.writetable(file_path, 
        "Buses"        => df_buses, 
        "Generators"   => df_generators, 
        "Circuits"     => df_circuits, 
        "Optimization" => df_optimization,
        overwrite = true
    )

    println("OPF results successfully saved in: $file_path")

end