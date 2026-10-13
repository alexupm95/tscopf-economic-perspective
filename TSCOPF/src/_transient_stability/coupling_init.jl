#=
================================================================================
 coupling_init.jl — steady-state / dynamics coupling initialisation
================================================================================
 FULL_BUS TSC-ACOPF seeds every coupling variable from a solved ACOPF warm start
 (`extract_opf_solved_hints`). That pre-solve is mandatory — see
 `resolve_fullbus_coupling_hints!` — so there is no alternative recipe here.
================================================================================
=#

"""Copy coupling hints onto steady-state OPF JuMP variables (`start=` only)."""
function apply_steady_state_hints_to_opf!(
    opf_dict::OrderedDict{Symbol, Any},
    hints::SteadyStateHints,
)
    for (bus, v) in opf_dict[:vars][:V]
        haskey(hints.val_V, bus) && JuMP.set_start_value(v, hints.val_V[bus])
    end
    for (bus, v) in opf_dict[:vars][:θ]
        haskey(hints.val_θ, bus) && JuMP.set_start_value(v, hints.val_θ[bus])
    end
    for (gen, v) in opf_dict[:vars][:P_g]
        haskey(hints.val_Pg, gen) && JuMP.set_start_value(v, hints.val_Pg[gen])
    end
    for (gen, v) in opf_dict[:vars][:Q_g]
        haskey(hints.val_Qg, gen) && JuMP.set_start_value(v, hints.val_Qg[gen])
    end
    return nothing
end
