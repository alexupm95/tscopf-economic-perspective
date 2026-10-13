# ===============================================================================
# -------------------------- Inequality Constraints -----------------------------

# This file contains several functions to define different inequality constraints
# ===============================================================================

# Function to create the Inequality Constraint maximum apparent power generated
function const_gen_power_apparent!(
    model::Model,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef},
    gen_id::AbstractVector{Int},
    max_lim::AbstractVector{Float64}
    )

    upper_cons = OrderedDict{Int, JuMP.ConstraintRef}()
    for (idx, gen) in enumerate(gen_id)
        # Logic: (P_g[gen])² + (Q_g[gen])² <= (max_lim[idx])²
        upper_cons[gen] = @constraint(model, (P_g[gen])^2 + (Q_g[gen])^2 - (max_lim[idx])^2 ≤ 0.0 )
    end

    return upper_cons
end

# (ACOPF) Function to create the Inequality Constraint maximum capacity of the branch
function ac_const_branch_thermal_limit!(
    model::Model,
    P_ik::OrderedDict{Int, JuMP.VariableRef},
    Q_ik::OrderedDict{Int, JuMP.VariableRef},
    P_ki::OrderedDict{Int, JuMP.VariableRef},
    Q_ki::OrderedDict{Int, JuMP.VariableRef},
    limit::AbstractVector{Float64}
    )
    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    for branch in keys(P_ik)
        isfinite(limit[branch]) || continue
        ineq_const_s_ik[branch] = JuMP.@constraint(model, P_ik[branch]^2 + Q_ik[branch]^2 - limit[branch]^2 ≤ 0.0 )
        ineq_const_s_ki[branch] = JuMP.@constraint(model, P_ki[branch]^2 + Q_ki[branch]^2 - limit[branch]^2 ≤ 0.0 )
    end

    return ineq_const_s_ik, ineq_const_s_ki
end

# (ACOPF) Function to create the Inequality Constraint maximum capacity of the branch using the Admittance Matrix
function ac_const_branch_thermal_limit_Ybus!(
    model::Model,
    Ybus::SparseMatrixCSC,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    limit::AbstractVector{Float64}
    )
    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    # Loop in all branches
    for branch in DCIR.id
        if DCIR.l_status[branch] == 1 # Check if the branch is ON

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
            angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i

            # ***************************************************
            # Equality Constraints for power flow in the branches
            # ***************************************************

            # Branch flow from i to k
            # Same model used in POWERMODELS
            P_ik_expr = g*((1/t_tap) * V[i])^2 + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*sin(angik))
            Q_ik_expr = -(b + bik_sh)*((1/t_tap) * V[i])^2 - (-b*t_r-g*t_i)/t_tap^2*(V[i]*V[k]*cos(angik)) + (-g*t_r+b*t_i)/t_tap^2*(V[i]*V[k]*sin(angik))


            # Branch flow from k to i
            # Same model used in POWERMODELS
            P_ki_expr = g*(V[k]^2) + (-g*t_r-b*t_i)/t_tap^2*(V[i] * V[k] * cos(angki)) + (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * sin(angki))
            Q_ki_expr = -(b + bik_sh)*(V[k]^2)  - (-b*t_r+g*t_i)/t_tap^2*( V[i] * V[k] * cos(angki)) + (-g*t_r-b*t_i)/t_tap^2*( V[i] * V[k] * sin(angki))

            # ****************************************************************
            # Inequality Constraint for capacity/thermal limits of power flows
            # ****************************************************************
            isfinite(limit[branch]) || continue
            ineq_const_s_ik[branch] = JuMP.@constraint(model, P_ik_expr^2 + Q_ik_expr^2 - limit[branch]^2 ≤ 0.0 )
            ineq_const_s_ki[branch] = JuMP.@constraint(model, P_ki_expr^2 + Q_ki_expr^2 - limit[branch]^2 ≤ 0.0 )
         
        end
    end

    return ineq_const_s_ik, ineq_const_s_ki
end

# (ACOPF) Function to create the Inequality Constraint maximum angular deviation between adjacent buses.
# `clamp_deg` (DispatchLimitsConfig.ang_diff_clamp_ac_deg) tightens any input limit beyond ±clamp_deg.
function ac_const_angle_differences!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    sorted_pair_info::Vector,
    pair_circ_map::Vector;
    clamp_deg::Float64=60.0,
    )
    return _const_angle_differences!(model, θ, sorted_pair_info, pair_circ_map, clamp_deg)
end

# ========================================================================================
#                  Exclusive functions for the DC-OPF model
# ========================================================================================

# (DCOPF) Function to create the Inequality Constraint maximum capacity of the branch
function dc_const_branch_thermal_limit!(
    model::Model,
    P_ik::OrderedDict{Int, JuMP.VariableRef},
    P_ki::OrderedDict{Int, JuMP.VariableRef},
    limit::AbstractVector{Float64}
    )
    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    for branch in keys(P_ik)
        isfinite(limit[branch]) || continue
        ineq_const_s_ik[branch] = JuMP.@constraint(model, P_ik[branch] - limit[branch] ≤ 0.0 )
        ineq_const_s_ki[branch] = JuMP.@constraint(model, P_ki[branch] - limit[branch] ≤ 0.0 )
    end

    return ineq_const_s_ik, ineq_const_s_ki
end

# (DCOPF) Function to create the Inequality Constraint maximum capacity of the branch using the Susceptance Matrix
function dc_const_branch_thermal_limit_Bbus!(
    model::Model,
    Bbus::SparseMatrixCSC,
    θ::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    limit::AbstractVector{Float64}
    )
    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    # Loop in all branches
    for branch in DCIR.id
        if DCIR.l_status[branch] == 1 # Check if the branch is ON

            i  = DCIR.from_bus[branch]                              # Bus i (from)
            k  = DCIR.to_bus[branch]                                # Bus k (to)

            # Read the off-diagonal susceptance straight from Bbus so the thermal limit
            # uses exactly the same convention (SIMPLE 1/x or POWERMODELS x/(r²+x²)) that
            # built Bbus — Bbus[i,k] = +b for a branch i→k under either susceptance model.
            b  = Bbus[i, k]

            angik   = θ[i] - θ[k]                                  # Angular difference between bus i and k
            angki   = θ[k] - θ[i]                                  # Angular difference between bus k and i

            # ***************************************************
            # Equality Constraints for power flow in the branches
            # ***************************************************

            # Branch flow from i to k
            P_ik_expr = b * (angik)

            # Branch flow from k to i
            P_ki_expr = b * (angki)

            # ****************************************************************
            # Inequality Constraint for capacity/thermal limits of power flows
            # ****************************************************************
            isfinite(limit[branch]) || continue
            ineq_const_s_ik[branch] = JuMP.@constraint(model, P_ik_expr - limit[branch] ≤ 0.0 )
            ineq_const_s_ki[branch] = JuMP.@constraint(model, P_ki_expr - limit[branch] ≤ 0.0 )
         
        end
    end

    return ineq_const_s_ik, ineq_const_s_ki
end

# (DCOPF) Function to create the Inequality Constraint maximum capacity of the branch
function dc_const_branch_thermal_limit_Bbus!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    DCIR::DataFrame,
    limit::AbstractVector{Float64}
    )
    ineq_const_s_ik = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus i to bus k
    ineq_const_s_ki = OrderedDict{Int, JuMP.ConstraintRef}() # Dict to save inequality constraints of apparent power flow from bus k to bus i

    # Loop in all branches
    for branch in DCIR.id
        if DCIR.l_status[branch] == 1 # Check if the branch is ON

            i     = DCIR.from_bus[branch]      # Bus i (from)
            k     = DCIR.to_bus[branch]        # Bus k (to)

            b     = 1 / (DCIR.l_reac[branch])  # Series admittance

            angik = θ[i] - θ[k]                # Angular difference between bus i and k
            angki = θ[k] - θ[i]                # Angular difference between bus k and i

            # ***************************************************
            # Equality Constraints for power flow in the branches
            # ***************************************************

            # Branch flow from i to k
            P_ik_expr = b * (angik)

            # Branch flow from k to i
            P_ki_expr = b * (angki)

            # ****************************************************************
            # Inequality Constraint for capacity/thermal limits of power flows
            # ****************************************************************
            isfinite(limit[branch]) || continue
            ineq_const_s_ik[branch] = JuMP.@constraint(model, P_ik_expr - limit[branch] ≤ 0.0 )
            ineq_const_s_ki[branch] = JuMP.@constraint(model, P_ki_expr - limit[branch] ≤ 0.0 )
         
        end
    end

    return ineq_const_s_ik, ineq_const_s_ki
end

# (DCOPF) Function to create the Inequality Constraint maximum angular deviation between adjacent buses.
# `clamp_deg` (DispatchLimitsConfig.ang_diff_clamp_dc_deg) tightens any input limit beyond ±clamp_deg.
function dc_const_angle_differences!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    sorted_pair_info::Vector,
    pair_circ_map::Vector;
    clamp_deg::Float64=30.0,
    )
    return _const_angle_differences!(model, θ, sorted_pair_info, pair_circ_map, clamp_deg)
end

"""Clamp `DCIR.ang_min`/`ang_max` to ±`clamp_deg` (safety clamp; see `DispatchLimitsConfig`)."""
function clamp_angle_diff_limits(min_ang::Float64, max_ang::Float64, clamp_deg::Float64,
                                 bus_from::Int, bus_to::Int)
    if min_ang < -clamp_deg
        println("Correcting angle constraints between adjacent buses ($bus_from, $bus_to): setting ang_min to -$(clamp_deg)°.")
        min_ang = -clamp_deg
    end
    if max_ang > clamp_deg
        println("Correcting angle constraints between adjacent buses ($bus_from, $bus_to): setting ang_max to +$(clamp_deg)°.")
        max_ang = clamp_deg
    end
    return min_ang, max_ang
end

# Shared body of the AC / DC angle-difference constraints (they differ only in the clamp value).
function _const_angle_differences!(
    model::Model,
    θ::OrderedDict{Int, JuMP.VariableRef},
    sorted_pair_info::Vector,
    pair_circ_map::Vector,
    clamp_deg::Float64,
    )

    ineq_const_ang_diff_lower  = OrderedDict{Int, JuMP.ConstraintRef}() # Vector to save inequality constraints of angle difference between adjacent buses
    ineq_const_ang_diff_upper  = OrderedDict{Int, JuMP.ConstraintRef}() # Vector to save inequality constraints of angle difference between adjacent buses

    pair_circ_map_dict = OrderedDict(pair_circ_map)

    for (pair, pair_data) in sorted_pair_info
        bus_from, bus_to = pair
        ang_ik = θ[bus_from] - θ[bus_to]

        circ_id = pair_circ_map_dict[pair]

        min_ang, max_ang = clamp_angle_diff_limits(
            pair_data.min_ang, pair_data.max_ang, clamp_deg, bus_from, bus_to)

        ineq_const_ang_diff_lower[circ_id] = JuMP.@constraint(model, deg2rad(min_ang) - ang_ik ≤ 0.0)
        ineq_const_ang_diff_upper[circ_id] = JuMP.@constraint(model, ang_ik - deg2rad(max_ang) ≤ 0.0)
    end

    return ineq_const_ang_diff_lower, ineq_const_ang_diff_upper
end

