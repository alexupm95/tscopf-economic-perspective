# ============================================================================
# ------------------------------- Objective ----------------------------------

# This file contains several functions to define different objective functions
# ============================================================================

# Function to create the objective function minimize fuel cost (it can be a quadratic function)
function obj_minimize_fuel_cost!(
    model::Model, 
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    gen_id::AbstractVector{Int},
    a::AbstractVector{Float64},
    b::AbstractVector{Float64},
    c::AbstractVector{Float64},
    base_MVA::Float64,
)::Tuple{JuMP.AbstractJuMPScalar, JuMP.AbstractJuMPScalar}
    # obj_function_MVA is a QuadExpr when any quadratic cost coefficient a[idx] ≠ 0
    # (e.g. case39) and an AffExpr otherwise (e.g. the 9-bus cases). The return type
    # must therefore stay AbstractJuMPScalar — annotating it as GenericAffExpr forces a
    # convert(AffExpr, QuadExpr) that throws on any truly quadratic cost.
    # Ipopt minimizes the scaled objective (P_g in p.u.); reported cost uses obj_function_MVA [€].
    objective_scaling = 1.0 / base_MVA

    total_cost_MVA = 0.0
    for (idx, gen) in enumerate(gen_id)
        if haskey(P_g, gen)
            P_MW = P_g[gen] * base_MVA
            c[idx] != 0.0 && (total_cost_MVA += c[idx])
            b[idx] != 0.0 && (total_cost_MVA += b[idx] * P_MW)
            a[idx] != 0.0 && (total_cost_MVA += a[idx] * P_MW^2)
        end
    end

    obj_function = JuMP.@objective(model, Min, objective_scaling * total_cost_MVA)
    obj_function_MVA = JuMP.@expression(model, total_cost_MVA)

    return obj_function, obj_function_MVA

end

# Function to create the objective function minimize marginal cost (it is a linear function)
function obj_minimize_marginal_cost!(
    model::Model, 
    P_g::OrderedDict{Int, JuMP.VariableRef}, 
    gen_id::AbstractVector{Int},
    marginal_cost::AbstractVector{Float64},
    base_MVA::Float64,
)::Tuple{JuMP.AbstractJuMPScalar, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}
    
    objective_scaling = 1.0 / base_MVA # Scaling factor that multiplies the objective function

    total_cost_MVA = 0.0
    for (idx, gen) in enumerate(gen_id)
        if haskey(P_g, gen)
            if marginal_cost[idx] != 0.0
                total_cost_MVA += marginal_cost[idx] * P_g[gen] * base_MVA
            end
        end
    end

    obj_function = JuMP.@objective(model, Min, objective_scaling * total_cost_MVA) # This is needed to set the objective function w.r.t. per-unit values
    # obj_function = JuMP.@objective(model, Min, total_cost_MVA)
    obj_function_MVA = JuMP.@expression(model, total_cost_MVA) # This expression is only used to recover the true cost of dispatch

    return obj_function, obj_function_MVA    

end