# ==============================================================================
#                    UNIT COMMITMENT (UC) — decision variables
# ==============================================================================
# Variable builders for the UC model, split out from functions_make_uc_model.jl
# so new UC variables (multi-period u[g,t] / P[g,t], reserves, startup flags, …)
# have a single home. Each returns an OrderedDict{Int, VariableRef} keyed by
# generator id, matching the dispatch-container convention.

"""Continuous active-power dispatch P_g (pu) for each active generator."""
function var_uc_dispatch!(model::Model, gen_ids::AbstractVector{Int})
    P_g = OrderedDict{Int, JuMP.VariableRef}()
    for gen in gen_ids
        P_g[gen] = JuMP.@variable(model, base_name = "P_g[$gen]")
    end
    return P_g
end

"""Binary commitment variables u_g ∈ {0,1} for each active generator."""
function var_uc_commitment!(model::Model, gen_ids::AbstractVector{Int})
    u = OrderedDict{Int, JuMP.VariableRef}()
    for gen in gen_ids
        # base_name keeps JuMP model human-readable in model_summary.txt exports.
        u[gen] = JuMP.@variable(model, binary = true, base_name = "u[$gen]")
    end
    return u
end
