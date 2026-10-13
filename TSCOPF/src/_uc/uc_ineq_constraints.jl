# ==============================================================================
#                 UNIT COMMITMENT (UC) — inequality constraints
# ==============================================================================
# Inequality-constraint builders for the UC model. Bounds are written in explicit
# ≤ form (ED convention) so JuMP.dual() on the restricted-pricing LP exposes the
# scarcity rents cleanly. New UC limits (ramp up/down, min up/down time, reserve
# contribution, …) belong here.

"""
Coupled on/off generation limits for the MILP (bilinear in u):
  P_min u_g - P_g ≤ 0   and   P_g - P_max u_g ≤ 0.

u_g = 0 forces P_g = 0; u_g = 1 reduces to ED bounds. Returns `(lower, upper)`
as OrderedDict{Int, ConstraintRef} keyed by generator id.
"""
function ineq_const_uc_pg_limits!(
    model::Model,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    u::OrderedDict{Int, JuMP.VariableRef},
    gen_ids::AbstractVector{Int},
    pg_min_pu::AbstractVector{Float64},
    pg_max_pu::AbstractVector{Float64},
)
    ineq_lower = OrderedDict{Int, JuMP.ConstraintRef}()
    ineq_upper = OrderedDict{Int, JuMP.ConstraintRef}()
    for (idx, gen) in enumerate(gen_ids)
        ineq_lower[gen] = JuMP.@constraint(model, pg_min_pu[idx] * u[gen] - P_g[gen] ≤ 0.0)
        ineq_upper[gen] = JuMP.@constraint(model, P_g[gen] - pg_max_pu[idx] * u[gen] ≤ 0.0)
    end
    return ineq_lower, ineq_upper
end

"""
Generation limits for the restricted-pricing LP with commitment fixed at u*:
  on  (u*_g = 1):  P_min ≤ P_g ≤ P_max
  off (u*_g = 0):  0 ≤ P_g ≤ 0

Linear (no binaries) so JuMP.dual() yields the SMP and scarcity rents. Returns
`(lower, upper)` keyed by generator id.
"""
function ineq_const_uc_pg_limits_fixed!(
    model::Model,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    gen_ids::AbstractVector{Int},
    pg_min_pu::AbstractVector{Float64},
    pg_max_pu::AbstractVector{Float64},
    u_star::OrderedDict{Int, Float64},
)
    ineq_lower = OrderedDict{Int, JuMP.ConstraintRef}()
    ineq_upper = OrderedDict{Int, JuMP.ConstraintRef}()
    for (idx, gen) in enumerate(gen_ids)
        if u_star[gen] > 0.5
            ineq_lower[gen] = JuMP.@constraint(model, pg_min_pu[idx] - P_g[gen] ≤ 0.0)
            ineq_upper[gen] = JuMP.@constraint(model, P_g[gen] - pg_max_pu[idx] ≤ 0.0)
        else
            ineq_lower[gen] = JuMP.@constraint(model, -P_g[gen] ≤ 0.0)
            ineq_upper[gen] = JuMP.@constraint(model, P_g[gen] ≤ 0.0)
        end
    end
    return ineq_lower, ineq_upper
end
