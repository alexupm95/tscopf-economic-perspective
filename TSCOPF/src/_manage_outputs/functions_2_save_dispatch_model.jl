# ===================================================================================
#                   PRINT THE OPTIMIZATION MODEL IN TXT FILE
# ===================================================================================
# One generic exporter (`Export_OPF_Model`) replaces the former
# Export_{ACOPF, ACOPF_w_Ybus, DCOPF, DCOPF_w_Bbus, ED}_Model trio-plus, which
# were ~340-line copies differing only in which variables / constraint families
# their model contained. The variable order and the constraint print order +
# banner titles now live once, in `MODEL_VAR_ORDER` and `MODEL_EXPORT_SPECS`;
# the writer skips anything absent from `opf_dict` (same `haskey` filtering as the
# dual exporter in DispatchDualRegistry.jl), so all model types share one function.

# Canonical variable print order; only those present in opf_dict[:vars] are written.
const MODEL_VAR_ORDER = Symbol[:V, :θ, :P_g, :Q_g, :P_ik, :Q_ik, :P_ki, :Q_ki]

# Constraint print order + banner titles. `:section` rows are unconditional
# dividers (no constraint lookup); the rest are skipped when the key is absent.
const MODEL_EXPORT_SPECS = Tuple{Symbol, Symbol, String}[
    (:eq_const, :eq_const_angle_sw, "Equality Constraint Angle Swing"),
    (:eq_const, :eq_const_p_balance, "Equality Constraints Active Power Balance for Buses"),
    (:eq_const, :eq_const_q_balance, "Equality Constraints Reactive Power Balance for Buses"),
    (:eq_const, :eq_const_p_ik, "Equality Constraints Active Power Flow from Line i to k"),
    (:eq_const, :eq_const_q_ik, "Equality Constraints Reactive Power Flow from Line i to k"),
    (:eq_const, :eq_const_p_ki, "Equality Constraints Active Power Flow from Line k to i"),
    (:eq_const, :eq_const_q_ki, "Equality Constraints Reactive Power Flow from Line k to i"),

    (:ineq_const, :ineq_const_sg_upper, "Inequality Constraints Capability Curve Generators"),
    (:ineq_const, :ineq_const_s_ik, "Inequality Constraints Capacity Power Flow from Line i to k"),
    (:ineq_const, :ineq_const_s_ki, "Inequality Constraints Capacity Power Flow from Line k to i"),
    (:ineq_const, :ineq_const_ang_diff_lower, "Inequality Constraints Voltage Angle Differences between Buses - Lower Bound"),
    (:ineq_const, :ineq_const_ang_diff_upper, "Inequality Constraints Voltage Angle Differences between Buses - Upper Bound"),
]

const MODEL_BOUND_SECTION_TITLE =
    "Inequality Constraints Inferior Limits Decision Variables"

# Box limits: explicit ≤-constraints (CONSTRAINT) or JuMP variable bounds (VARIABLE).
const MODEL_BOUND_SPECS = Tuple{Symbol, Symbol, String}[
    (:ineq_const, :ineq_const_volt_mag_lower, "Voltage Magnitudes - Lower Bound"),
    (:ineq_const, :ineq_const_volt_mag_upper, "Voltage Magnitudes - Upper Bound"),
    (:ineq_const, :ineq_const_volt_ang_lower, "Voltage Angles - Lower Bound"),
    (:ineq_const, :ineq_const_volt_ang_upper, "Voltage Angles - Upper Bound"),
    (:ineq_const, :ineq_const_pg_lower, "Active Power Generated - Lower Bound"),
    (:ineq_const, :ineq_const_pg_upper, "Active Power Generated - Upper Bound"),
    (:ineq_const, :ineq_const_qg_lower, "Reactive Power Generated - Lower Bound"),
    (:ineq_const, :ineq_const_qg_upper, "Reactive Power Generated - Upper Bound"),
    (:ineq_const, :ineq_const_pik_lower, "Active Power Flow ik - Lower Bound"),
    (:ineq_const, :ineq_const_pik_upper, "Active Power Flow ik - Upper Bound"),
    (:ineq_const, :ineq_const_qik_lower, "Reactive Power Flow ik - Lower Bound"),
    (:ineq_const, :ineq_const_qik_upper, "Reactive Power Flow ik - Upper Bound"),
    (:ineq_const, :ineq_const_pki_lower, "Active Power Flow ki - Lower Bound"),
    (:ineq_const, :ineq_const_pki_upper, "Active Power Flow ki - Upper Bound"),
    (:ineq_const, :ineq_const_qki_lower, "Reactive Power Flow ki - Lower Bound"),
    (:ineq_const, :ineq_const_qki_upper, "Reactive Power Flow ki - Upper Bound"),
]

"""Print a `===` / title / `===` banner sized to the title."""
function _print_section_banner(io::IO, title::AbstractString)
    bar = "=" ^ (length(title) + 1)
    println(io, bar)
    println(io, title, " ")
    println(io, bar)
end

"""`true` when box limits exist as constraints or as JuMP variable bounds."""
function _model_bound_present(opf_dict::OrderedDict{Symbol, Any}, key::Symbol)::Bool
    if haskey(opf_dict, :ineq_const) && haskey(opf_dict[:ineq_const], key)
        return true
    end
    return bound_manifest_entry_present(opf_dict, key)
end

"""One line for a JuMP variable bound, using the same ≤-form as constraint export."""
function _format_variable_bound_line(var::JuMP.VariableRef, side::Symbol)::String
    if side == :lower
        JuMP.has_lower_bound(var) ||
            throw(ArgumentError("expected lower bound on $var"))
        lb = JuMP.lower_bound(var)
        return string("-", var, " <= ", -lb)
    elseif side == :upper
        JuMP.has_upper_bound(var) ||
            throw(ArgumentError("expected upper bound on $var"))
        return string(var, " <= ", JuMP.upper_bound(var))
    end
    throw(ArgumentError("side must be :lower or :upper, got $side"))
end

"""Write one box-limit family from constraints or from the bound manifest."""
function _print_bound_entry(
    io::IO,
    opf_dict::OrderedDict{Symbol, Any},
    entry::Tuple{Symbol, Symbol, String},
)
    _, key, title = entry
    _model_bound_present(opf_dict, key) || return
    _print_section_banner(io, title)
    if haskey(opf_dict, :ineq_const) && haskey(opf_dict[:ineq_const], key)
        for (i, info) in opf_dict[:ineq_const][key]
            println(io, "$i: ", info)
        end
    else
        manifest_entry = find_bound_manifest_entry(opf_dict, key)
        manifest_entry === nothing && return
        for (i, var) in manifest_entry.vars
            println(io, "$i: ", _format_variable_bound_line(var, manifest_entry.side))
        end
    end
    println(io, "\n")
end

"""Write one constraint family (banner + entries). Section dividers are handled separately."""
function _print_model_entry(io::IO, opf_dict::OrderedDict{Symbol, Any}, entry::Tuple{Symbol, Symbol, String})
    source, key, title = entry
    source == :section && return
    (haskey(opf_dict, source) && haskey(opf_dict[source], key)) || return
    _print_section_banner(io, title)
    for (i, info) in opf_dict[source][key]
        println(io, "$i: ", info)
    end
    println(io, "\n")
end

function _print_bound_block(io::IO, opf_dict::OrderedDict{Symbol, Any})
    any(_model_bound_present(opf_dict, key) for (_, key, _) in MODEL_BOUND_SPECS) || return
    _print_section_banner(io, MODEL_BOUND_SECTION_TITLE)
    for entry in MODEL_BOUND_SPECS
        _print_bound_entry(io, opf_dict, entry)
    end
end

"""
    Export_OPF_Model(model, path_names, obj_function, opf_dict; model_label="OPF")

Dump the solved/built model to `model_summary.txt` (JuMP `show`) and a structured
`model_details.txt` (objective, variables, constraint families) under
`path_names[:pf_dispatch]`. Variable and constraint sets are filtered from
`MODEL_VAR_ORDER` / `MODEL_EXPORT_SPECS` by what exists in `opf_dict`, so AC-OPF
(with or without Ybus), DC-OPF (with or without Bbus) and ED all use this one
function. `model_label` only tags the console message.
"""
function Export_OPF_Model(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    opf_dict::OrderedDict{Symbol, Any};
    model_label::String = "OPF",
    )

    open(joinpath(path_names[:pf_dispatch], "model_summary.txt"), "w") do io
        show(io, model)
    end

    open(joinpath(path_names[:pf_dispatch], "model_details.txt"), "w") do io
        # --- objective ---
        println(io, "=========")
        println(io, "Objective ")
        println(io, "=========")
        println(io, obj_function)
        println(io, "\n")

        # --- variables (canonical order, present only) ---
        println(io, "=========")
        println(io, "Variables")
        println(io, "=========")
        for key in MODEL_VAR_ORDER
            haskey(opf_dict[:vars], key) || continue
            for (j, info) in opf_dict[:vars][key]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # --- constraints (registry-ordered, present only) ---
        for entry in MODEL_EXPORT_SPECS
            _print_model_entry(io, opf_dict, entry)
        end
        _print_bound_block(io, opf_dict)
    end

    println("$model_label Model successfully saved as TXT file in: ", path_names[:pf_dispatch])
end
