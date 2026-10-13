# =============================================================================
# -------------------------- Equality Constraints -----------------------------

# This file contains several functions to define different equality constraints
# =============================================================================

"""
    assert_bus_not_islanded(bus, circ_ids)

Throw if `circ_ids` is empty (no incident line/transformer at `bus`).
Shared by ACOPF/DCOPF balance builders and FULL_BUS nodal injection.
"""
function assert_bus_not_islanded(bus, circ_ids)
    if isempty(circ_ids)
        throw(ArgumentError(
            "The bus $bus is islanded, i.e., there is no line or transformer connected to it."))
    end
    return nothing
end

# Function to force the angle of the slack bus to be zero at the steady state
function const_angle_slack_bus!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    idx_slack_bus::Int
    )

    eq_const_angle_sw = JuMP.@constraint(model, θ[idx_slack_bus] == 0.0)   # Set the constraint -> Angle == 0

    return eq_const_angle_sw
end

"""
    const_pin_load_active!(model, bus_gen_circ_dict_ON, base_MVA) -> (P_d, eq_const_pd_pin)

Promote the bus active demand from an inlined `Float64` to a JuMP variable held at its
data value by an explicit equality, one per bus:

    P_d[k]  free variable,      eq_const_pd_pin[k] :  P_d[k] - p_d⁰[k] == 0   : σ_k

The feasible set is unchanged (the pin fixes the variable), so the optimum, the objective
and every other dual are identical to the inlined build. What is gained is a real Jacobian
column for the demand.

WHY THIS EXISTS. The demand-side price is

    π_k = ∂f/∂P_d(k) - Σ_r λ_r ∂h_r/∂P_d(k) - Σ_s μ_s ∂g_s/∂P_d(k)

and stationarity of the Lagrangian in the new variable `P_d[k]` is exactly that sum minus
σ_k, so

    σ_k = π_k.

The solver therefore accumulates the parameter derivative over EVERY row that happens to
contain the demand, whatever those rows are. No family list, no census, nothing to go stale
when a builder changes — which is the failure mode a hand-written table of "rows containing
P_d" invites.

The alternative is to recover the same column by rebuilding the model at `p_d ± ε` and
differencing the constraint residuals at the frozen solution. That route remains the only
one available where the demand is absorbed into a numeric matrix before JuMP ever sees it
(the Kron-reduced `Y_red` constant-impedance shunt), but where the demand survives as an
inlined coefficient this pin is exact, costs one solve instead of three per bus, and needs
no `ε`.

Deliberately an explicit `@constraint` rather than `JuMP.fix()`: a fixed-variable dual is a
variable-bound dual, and the CONSTRAINT bound encoding exists precisely so `JuMP.dual` stays
clean and uniform across families.

Every bus is pinned, zero-demand buses included — that is what makes π_k well defined at a
bus with no load.
"""
function const_pin_load_active!(
    model::Model,
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}},
    base_MVA::Float64
    )

    P_d             = OrderedDict{Int, JuMP.VariableRef}()
    eq_const_pd_pin = OrderedDict{Int, JuMP.ConstraintRef}()

    for bus in keys(bus_gen_circ_dict_ON)
        p_d0 = bus_gen_circ_dict_ON[bus][:pd_tot] / base_MVA  # Load active power in the bus (p.u.)

        # `start` at the pinned value: the variable is determined, so this hands the solver a
        # feasible point for the new column instead of the default zero.
        P_d[bus] = JuMP.@variable(model, base_name = "P_d[$bus]", start = p_d0)
        eq_const_pd_pin[bus] = JuMP.@constraint(model, P_d[bus] - p_d0 == 0.0)
    end

    return P_d, eq_const_pd_pin
end

# (ACOPF) Function to create the active power balance constraint using the variables Pik and Pki
function ac_const_active_power_balance!(
    model::Model, 
    V::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    P_ik::OrderedDict{Int, JuMP.VariableRef}, 
    P_ki::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64;
    P_d::Union{Nothing, OrderedDict{Int, JuMP.VariableRef}} = nothing
    )

    # Create expressions for the inflow/outflow of power in all buses according to the model variables
    p_flow_terms_dict = OrderedDict{Int, Vector{JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary of active power inflow/outflow

    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        terms_p = JuMP.GenericAffExpr{Float64, JuMP.VariableRef}[]  # List of expressions

        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus i
        assert_bus_not_islanded(bus, indices_circ_connected)

        for id_branch in bus_gen_circ_dict_ON[bus][:circ]

            from_bus = DCIR.from_bus[id_branch]
            to_bus = DCIR.to_bus[id_branch]

            if from_bus == bus
                push!(terms_p, P_ik[id_branch])  # Power flows out of bus
            elseif to_bus == bus
                push!(terms_p, P_ki[id_branch])  # Power flows into bus
            end
        end
        p_flow_terms_dict[bus] = terms_p
    end
   
    # ************************************
    # Constraints for Active Power Balance
    # ************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        indices_bus_gen        = bus_gen_circ_dict_ON[bus][:gen_ids] # Generators at bus
        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus
        terms_p                = p_flow_terms_dict[bus]              # Active power flow terms for bus

        assert_bus_not_islanded(bus, indices_circ_connected)

        # Get generator variables from the dictionary for those ON at this bus
        Pg_terms = [P_g[gen] for gen in indices_bus_gen if haskey(P_g, gen)]

        Pg_sum = isempty(Pg_terms) ? 0.0 : sum(Pg_terms) # Sum P_g terms; zero otherwise

        # Demand: the data value, or the pinned variable when `const_pin_load_active!` ran.
        # Both branches write the SAME row — the pin holds the variable at `p_d⁰` — so the
        # solution is invariant and only the Jacobian gains a column.
        p_d = P_d === nothing ? bus_gen_circ_dict_ON[bus][:pd_tot] / base_MVA : P_d[bus]
        g_sh = bus_gen_circ_dict_ON[bus][:gsh_tot] / base_MVA # Shunt conductance in the bus

        if g_sh == 0.0 # Check if the bus has shunt elements connected to it 
            eq_const_p_balance[bus] = JuMP.@constraint(model, Pg_sum - p_d - sum(terms_p) == 0.0)

        else
            eq_const_p_balance[bus] = JuMP.@constraint(model, Pg_sum - p_d - g_sh * V[bus]^2 - sum(terms_p) == 0.0)

        end
    end

    return eq_const_p_balance
end

# (ACOPF) Function to create the Equality Constraint reactive power balance constraint using the variables Qik and Qki
function ac_const_reactive_power_balance!(
    model::Model, 
    V::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef}, 
    Q_ik::OrderedDict{Int, JuMP.VariableRef}, 
    Q_ki::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64
    )

    # Create expressions for the inflow/outflow of power in all buses according to the model variables
    q_flow_terms_dict = OrderedDict{Int, Vector{JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary of reactive power inflow/outflow

    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        terms_q = JuMP.GenericAffExpr{Float64, JuMP.VariableRef}[]  # List of expressions

        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus i
        assert_bus_not_islanded(bus, indices_circ_connected)

        for id_branch in bus_gen_circ_dict_ON[bus][:circ]

            from_bus = DCIR.from_bus[id_branch]
            to_bus = DCIR.to_bus[id_branch]

            if from_bus == bus
                push!(terms_q, Q_ik[id_branch])  # Power flows out of bus
            elseif to_bus == bus
                push!(terms_q, Q_ki[id_branch])  # Power flows into bus
            end
        end
        q_flow_terms_dict[bus] = terms_q
    end
   
    # ************************************
    # Constraints for Rective Power Balance
    # ************************************
    eq_const_q_balance = OrderedDict{Int, JuMP.ConstraintRef}()


    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        indices_bus_gen        = bus_gen_circ_dict_ON[bus][:gen_ids] # Generators at bus
        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus
        terms_q                = q_flow_terms_dict[bus]              # Reactive power flow terms for bus

        assert_bus_not_islanded(bus, indices_circ_connected)

        # Get generator variables from the dictionary for those ON at this bus
        Qg_terms = [Q_g[gen] for gen in indices_bus_gen if haskey(Q_g, gen)]

        Qg_sum = isempty(Qg_terms) ? 0.0 : sum(Qg_terms) # Sum Q_g terms; zero otherwise

        q_d = bus_gen_circ_dict_ON[bus][:qd_tot] / base_MVA   # Load reactive power in the bus
        b_sh = bus_gen_circ_dict_ON[bus][:bsh_tot] / base_MVA # Shunt susceptance in the bus

        if b_sh == 0.0 # Check if the bus has shunt elements connected to it 
            eq_const_q_balance[bus] = JuMP.@constraint(model, Qg_sum - q_d - sum(terms_q) == 0.0)

        else
            eq_const_q_balance[bus] = JuMP.@constraint(model, Qg_sum - q_d + b_sh * V[bus]^2 - sum(terms_q) == 0.0)

        end
    end

    return eq_const_q_balance
end

# (ACOPF) Function to create the active power balance constraint using the bus admittance matrix
function ac_const_active_power_balance_Ybus!(
    model::Model, 
    Ybus::SparseMatrixCSC,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64;
    P_d::Union{Nothing, OrderedDict{Int, JuMP.VariableRef}} = nothing
    )

    # ******************************************************************************************************
    # Create expressions for the net inflow/outflow of power in all buses according to the admittance matrix
    # ******************************************************************************************************
    # Net active injection at bus i:  P_i = V_i · Σ_k V_k (G_ik cos θ_ik + B_ik sin θ_ik).
    # We build it in ONE pass over the nonzeros of the sparse Ybus (O(nnz)), instead
    # of probing Ybus[i,k] over every (i,k) pair inside an outer per-bus loop that
    # rebuilt and overwrote the whole dictionary nBUS times (was O(nBUS³) and defeated
    # the sparsity). Ybus is column-major: iterating column k and reading row i = the
    # term that lands in bus i's sum, so accumulating by row reconstructs each row dot
    # product exactly. The k = i (diagonal) term is included (θ_ii = 0 ⇒ V_i·G_ii).
    buses = keys(bus_gen_circ_dict_ON)

    # Fail fast on an islanded bus (no incident branch) before assembling anything.
    for bus in buses
        assert_bus_not_islanded(bus, bus_gen_circ_dict_ON[bus][:circ])
    end

    terms_p_net = OrderedDict{Int, Vector{JuMP.NonlinearExpr}}(i => JuMP.NonlinearExpr[] for i in buses)
    rows = SparseArrays.rowvals(Ybus)
    vals = SparseArrays.nonzeros(Ybus)
    for k in buses                                   # column index k of Ybus
        for p in SparseArrays.nzrange(Ybus, k)
            i   = rows[p]                            # row index i (nonzero Ybus[i,k])
            Yik = vals[p]
            Gik = real(Yik); Bik = imag(Yik)
            angik = θ[i] - θ[k]                      # angular difference between bus i and k
            push!(terms_p_net[i], V[k] * (Gik * cos(angik) + Bik * sin(angik)))
        end
    end

    p_net_dict = OrderedDict{Int, JuMP.NonlinearExpr}() # Net active power at the bus
    for i in buses
        p_net_dict[i] = V[i] * sum(terms_p_net[i])   # diagonal guarantees a nonempty sum
    end

    # ************************************
    # Constraints for Active Power Balance
    # ************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in buses
        indices_bus_gen = bus_gen_circ_dict_ON[i][:gen_ids] # Generators at bus i

        # Get generator variables from the dictionary for those ON at this bus
        Pg_terms = [P_g[g] for g in indices_bus_gen if haskey(P_g, g)]

        Pg_sum = isempty(Pg_terms) ? 0.0 : sum(Pg_terms) # Sum P_g terms; zero otherwise

        # Same substitution as the branch-flow form: data value, or the pinned variable.
        p_d = P_d === nothing ? bus_gen_circ_dict_ON[i][:pd_tot] / base_MVA : P_d[i]

        eq_const_p_balance[i] = JuMP.@constraint(model, Pg_sum - p_d - p_net_dict[i]== 0.0)

    end

    return eq_const_p_balance
end

# (ACOPF) Function to create the Equality Constraint reactive power balance constraint using the bus admittance matrix
function ac_const_reactive_power_balance_Ybus!(
    model::Model, 
    Ybus::SparseMatrixCSC,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64
    )

    # ******************************************************************************************************
    # Create expressions for the net inflow/outflow of power in all buses according to the admittance matrix
    # ******************************************************************************************************
    # Net reactive injection at bus i:  Q_i = V_i · Σ_k V_k (G_ik sin θ_ik − B_ik cos θ_ik).
    # Single O(nnz) pass over the sparse Ybus, accumulating by row (column k, row i),
    # mirroring the active-power builder. (The previous version also allocated an
    # unused `p_net_dict` and re-built everything nBUS times — both removed.)
    buses = keys(bus_gen_circ_dict_ON)

    # Fail fast on an islanded bus (no incident branch) before assembling anything.
    for bus in buses
        assert_bus_not_islanded(bus, bus_gen_circ_dict_ON[bus][:circ])
    end

    terms_q_net = OrderedDict{Int, Vector{JuMP.NonlinearExpr}}(i => JuMP.NonlinearExpr[] for i in buses)
    rows = SparseArrays.rowvals(Ybus)
    vals = SparseArrays.nonzeros(Ybus)
    for k in buses                                   # column index k of Ybus
        for p in SparseArrays.nzrange(Ybus, k)
            i   = rows[p]                            # row index i (nonzero Ybus[i,k])
            Yik = vals[p]
            Gik = real(Yik); Bik = imag(Yik)
            angik = θ[i] - θ[k]                      # angular difference between bus i and k
            push!(terms_q_net[i], V[k] * (Gik * sin(angik) - Bik * cos(angik)))
        end
    end

    q_net_dict = OrderedDict{Int, JuMP.NonlinearExpr}() # Net reactive power at the bus
    for i in buses
        q_net_dict[i] = V[i] * sum(terms_q_net[i])   # diagonal guarantees a nonempty sum
    end

    # ************************************
    # Constraints for Reactive Power Balance
    # ************************************
    eq_const_q_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in buses
        indices_bus_gen = bus_gen_circ_dict_ON[i][:gen_ids] # Generators at bus i

        # Get generator variables from the dictionary for those ON at this bus
        Qg_terms = [Q_g[g] for g in indices_bus_gen if haskey(Q_g, g)]

        Qg_sum = isempty(Qg_terms) ? 0.0 : sum(Qg_terms) # Sum Q_g terms; zero otherwise

        q_d = bus_gen_circ_dict_ON[i][:qd_tot] / base_MVA   # Load reactive power in the bus

        eq_const_q_balance[i] = JuMP.@constraint(model, Qg_sum - q_d - q_net_dict[i] == 0.0)

    end

    return eq_const_q_balance
end

# (ACOPF) Function to create the Equality Constraint active power flow from bus i to bus k
function ac_const_powerflow_ik_active!(
    model::Model, 
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_ik::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    )
    eq_const_p_ik = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        yik     = 1 / (DCIR.l_res[branch] + 1im*DCIR.l_reac[branch]) # Series admittance
        bik_sh  = DCIR.l_sh_susp[branch] / 2                      # Shunt suscpetance
        g       = real(yik)                                    # Branch series conductance
        b       = imag(yik)                                    # Branch series susceptance

        t_tap   = DCIR.t_tap[branch]                              # Transformer tap ratio (tap:1)
        t_shift = deg2rad(DCIR.t_shift[branch])                   # Transformer shift angle

        t_r     = t_tap * cos(t_shift)                         # Real number for transformers
        t_i     = t_tap * sin(t_shift)                         # Imaginary number for transformers

        angik   = θ[i] - θ[k]                                  # Angular difference between bus i and k

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************

        # Line flow from i to k
        # Same model used in POWERMODELS
        eq_const_p_ik[branch] = JuMP.@constraint(model, P_ik[branch] == g*((1/t_tap) * V[i])^2 + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*sin(angik)) )  
    end

    return eq_const_p_ik
end

# (ACOPF) Function to create the Equality Constraint reactive power flow from bus i to bus k
function ac_const_powerflow_ik_reactive!(
    model::Model, 
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    Q_ik::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    )

    eq_const_q_ik = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        yik     = 1 / (DCIR.l_res[branch] + 1im*DCIR.l_reac[branch]) # Series admittance
        bik_sh  = DCIR.l_sh_susp[branch] / 2                      # Shunt suscpetance
        g       = real(yik)                                    # Branch series conductance
        b       = imag(yik)                                    # Branch series susceptance

        t_tap   = DCIR.t_tap[branch]                              # Transformer tap ratio (tap:1)
        t_shift = deg2rad(DCIR.t_shift[branch])                   # Transformer shift angle

        t_r     = t_tap * cos(t_shift)                         # Real number for transformers
        t_i     = t_tap * sin(t_shift)                         # Imaginary number for transformers

        angik   = θ[i] - θ[k]                                  # Angular difference between bus i and k

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************

        # Line flow from i to k
        # Same model used in POWERMODELS
        eq_const_q_ik[branch] = JuMP.@constraint(model, Q_ik[branch] == -(b + bik_sh)*((1/t_tap) * V[i])^2 - (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*sin(angik)) )
    end

    return eq_const_q_ik
end

# (ACOPF) Function to create the Equality Constraint active power flow from bus k to bus i
function ac_const_powerflow_ki_active!(
    model::Model,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_ki::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    )

    eq_const_p_ki = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        yik     = 1 / (DCIR.l_res[branch] + 1im*DCIR.l_reac[branch]) # Series admittance
        bik_sh  = DCIR.l_sh_susp[branch] / 2                         # Shunt suscpetance
        g       = real(yik)                                          # Branch series conductance
        b       = imag(yik)                                          # Branch series susceptance

        t_tap   = DCIR.t_tap[branch]                              # Transformer tap ratio (tap:1)
        t_shift = deg2rad(DCIR.t_shift[branch])                   # Transformer shift angle

        t_r     = t_tap * cos(t_shift)                         # Real number for transformers
        t_i     = t_tap * sin(t_shift)                         # Imaginary number for transformers

        angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************
        # Line flow from k to i
        # Same model used in POWERMODELS
        eq_const_p_ki[branch] = JuMP.@constraint(model, P_ki[branch] == g*(V[k]^2) + (-g*t_r-b*t_i)/t_tap^2*(V[i] * V[k] * cos(angki)) + (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * sin(angki)) )
           
    end

    return eq_const_p_ki
end

# (ACOPF) Function to create the Equality Constraint reactive power flow from bus k to bus i
function ac_const_powerflow_ki_reactive!(
    model::Model,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    Q_ki::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    )

    eq_const_q_ki = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        yik     = 1 / (DCIR.l_res[branch] + 1im*DCIR.l_reac[branch]) # Series admittance
        bik_sh  = DCIR.l_sh_susp[branch] / 2                      # Shunt suscpetance
        g       = real(yik)                                    # Branch series conductance
        b       = imag(yik)                                    # Branch series susceptance

        t_tap   = DCIR.t_tap[branch]                              # Transformer tap ratio (tap:1)
        t_shift = deg2rad(DCIR.t_shift[branch])                   # Transformer shift angle

        t_r     = t_tap * cos(t_shift)                         # Real number for transformers
        t_i     = t_tap * sin(t_shift)                         # Imaginary number for transformers

        angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************

        # Line flow from k to i
        # Same model used in POWERMODELS
        eq_const_q_ki[branch] = JuMP.@constraint(model, Q_ki[branch] == -(b + bik_sh)*(V[k]^2)  - (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * cos(angki)) + (-g*t_r-b*t_i)/t_tap^2*( V[i] * V[k] * sin(angki)) )
           
    end

    return eq_const_q_ki
end


# ========================================================================================
#                  Exclusive functions for the DC-OPF model
# ========================================================================================


# (DCOPF) Function to create the active power balance constraint using the variables Pik and Pki
function dc_const_active_power_balance!(
    model::Model, 
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    P_ik::OrderedDict{Int, JuMP.VariableRef}, 
    P_ki::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64
    )

    # Create expressions for the inflow/outflow of power in all buses according to the model variables
    p_flow_terms_dict = OrderedDict{Int, Vector{JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary of active power inflow/outflow

    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        terms_p = JuMP.GenericAffExpr{Float64, JuMP.VariableRef}[]  # List of expressions

        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus i
        assert_bus_not_islanded(bus, indices_circ_connected)

        for id_branch in bus_gen_circ_dict_ON[bus][:circ]

            from_bus = DCIR.from_bus[id_branch]
            to_bus = DCIR.to_bus[id_branch]

            if from_bus == bus
                push!(terms_p, P_ik[id_branch])  # Power flows out of bus
            elseif to_bus == bus
                push!(terms_p, P_ki[id_branch])  # Power flows into bus
            end
        end
        p_flow_terms_dict[bus] = terms_p
    end
   
    # ************************************
    # Constraints for Active Power Balance
    # ************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for bus in keys(bus_gen_circ_dict_ON) # Loop in all buses
        indices_bus_gen        = bus_gen_circ_dict_ON[bus][:gen_ids] # Generators at bus
        indices_circ_connected = bus_gen_circ_dict_ON[bus][:circ]    # Circuits connected to bus
        terms_p                = p_flow_terms_dict[bus]              # Active power flow terms for bus

        assert_bus_not_islanded(bus, indices_circ_connected)

        # Get generator variables from the dictionary for those ON at this bus
        Pg_terms = [P_g[gen] for gen in indices_bus_gen if haskey(P_g, gen)]

        Pg_sum = isempty(Pg_terms) ? 0.0 : sum(Pg_terms) # Sum P_g terms; zero otherwise

        p_d = bus_gen_circ_dict_ON[bus][:pd_tot] / base_MVA   # Load active power in the bus

        eq_const_p_balance[bus] = JuMP.@constraint(model, Pg_sum - p_d - sum(terms_p) == 0.0)

    end

    return eq_const_p_balance
end

# (DCOPF) Function to create the active power balance constraint using the bus susceptance matrix
function dc_const_active_power_balance_Bbus!(
    model::Model, 
    Bbus::SparseMatrixCSC,
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    DCIR::DataFrame,
    base_MVA::Float64
    )

    # ******************************************************************************************************
    # Create expressions for the net inflow/outflow of power in all buses according to the susceptance matrix
    # ******************************************************************************************************
    # Net active injection at bus i (DC):  P_i = Σ_k B_ik (θ_i − θ_k).
    # Single O(nnz) pass over the sparse Bbus, accumulating by row (column k, row i),
    # mirroring the AC builders. (The previous nested form was O(nBUS³) and probed
    # Bbus[i,k] over every pair, discarding the sparsity.) The k = i term is 0 (θ_ii = 0).
    buses = keys(bus_gen_circ_dict_ON)

    # Fail fast on an islanded bus (no incident branch) before assembling anything.
    for bus in buses
        assert_bus_not_islanded(bus, bus_gen_circ_dict_ON[bus][:circ])
    end

    terms_p_net = OrderedDict{Int, Vector{JuMP.AffExpr}}(i => JuMP.AffExpr[] for i in buses)
    rows = SparseArrays.rowvals(Bbus)
    vals = SparseArrays.nonzeros(Bbus)
    for k in buses                                   # column index k of Bbus
        for p in SparseArrays.nzrange(Bbus, k)
            i   = rows[p]                            # row index i (nonzero Bbus[i,k])
            Bik = vals[p]
            angik = θ[i] - θ[k]                      # angular difference between bus i and k
            push!(terms_p_net[i], Bik * angik)
        end
    end

    p_net_dict = OrderedDict{Int, JuMP.AffExpr}() # Net active power at the bus
    for i in buses
        p_net_dict[i] = sum(terms_p_net[i])          # diagonal guarantees a nonempty sum
    end

    # ************************************
    # Constraints for Active Power Balance
    # ************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    for i in buses
        indices_bus_gen = bus_gen_circ_dict_ON[i][:gen_ids] # Generators at bus i

        # Get generator variables from the dictionary for those ON at this bus
        Pg_terms = [P_g[g] for g in indices_bus_gen if haskey(P_g, g)]

        Pg_sum = isempty(Pg_terms) ? 0.0 : sum(Pg_terms) # Sum P_g terms; zero otherwise

        p_d = bus_gen_circ_dict_ON[i][:pd_tot] / base_MVA   # Load active power in the bus

        eq_const_p_balance[i] = JuMP.@constraint(model, Pg_sum - p_d - p_net_dict[i]== 0.0)

    end

    return eq_const_p_balance
end

# (DCOPF) Function to create the Equality Constraint active power flow from bus i to bus k
function dc_const_powerflow_ik_active!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_ik::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    susc_model::SusceptanceModel,
    )
    eq_const_p_ik = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        b       = dc_branch_susceptance(DCIR.l_res[branch], DCIR.l_reac[branch], susc_model) # Branch series susceptance

        angik   = θ[i] - θ[k]                                  # Angular difference between bus i and k

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************

        # Line flow from i to k
        eq_const_p_ik[branch] = JuMP.@constraint(model, P_ik[branch] == b * (angik) )
    end

    return eq_const_p_ik
end

# (DCOPF) Function to create the Equality Constraint active power flow from bus k to bus i
function dc_const_powerflow_ki_active!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_ki::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    active_branches::Vector{Int64},
    susc_model::SusceptanceModel,
    )
    eq_const_p_ki = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop in all branches
    for branch in active_branches
        i       = DCIR.from_bus[branch]                           # Bus i (from)
        k       = DCIR.to_bus[branch]                             # Bus k (to)

        b       = dc_branch_susceptance(DCIR.l_res[branch], DCIR.l_reac[branch], susc_model) # Branch series susceptance

        angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i

        # ***************************************************
        # Equality Constraints for power flow in the branches
        # ***************************************************

        # Line flow from k to i
        eq_const_p_ki[branch] = JuMP.@constraint(model, P_ki[branch] == b * (angki) )
    end

    return eq_const_p_ki
end

# ========================================================================================
#                Exclusive functions for the Economic Dispatch model
# ========================================================================================


# (ED) Function to create the active power balance constraint using the variables Pik and Pki
function ed_const_power_balance!(
    model::Model, 
    P_g::OrderedDict{Int, JuMP.VariableRef},  
    bus_gen_circ_dict_ON::OrderedDict{Int, Dict{Symbol, Any}}, 
    active_gen_ids::Vector{Int64},
    base_MVA::Float64
    )

    # ************************************
    # Constraints for Active Power Balance
    # ************************************
    eq_const_p_balance = OrderedDict{Int, JuMP.ConstraintRef}()

    p_d_total = sum(bus_gen_circ_dict_ON[bus][:pd_tot] for bus in keys(bus_gen_circ_dict_ON)) / base_MVA

    eq_const_p_balance[1] = JuMP.@constraint(model, sum(P_g[gen] for gen in active_gen_ids) - p_d_total == 0.0)

    return eq_const_p_balance
end
