# ==============================================================================
#                  UNIT COMMITMENT (UC) — equality constraints
# ==============================================================================
# Equality-constraint builders for the UC model. The single-period UC reuses the
# ED system power balance; isolating it here gives future inter-temporal
# equalities (energy balance across periods, storage state-of-charge, …) a home
# without editing the orchestrator.

"""
System active-power balance Σ_g P_g = demand (per bus), keyed by bus id.

Delegates to the shared ED builder so single-period UC and ED stay numerically
identical. The seam exists so multi-period / storage balances can replace this
without touching `Make_UC_Model!`.
"""
function eq_const_uc_power_balance!(
    model::Model,
    P_g::OrderedDict{Int, JuMP.VariableRef},
    bus_gen_circ_dict_ON::OrderedDict,
    gen_ids::AbstractVector{Int},
    base_MVA::Float64,
)
    # ponytail: thin delegation today; the wrapper is the extension point, add
    # multi-period balance here when n_periods > 1.
    return ed_const_power_balance!(model, P_g, bus_gen_circ_dict_ON, gen_ids, base_MVA)
end
