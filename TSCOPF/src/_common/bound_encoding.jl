#=
================================================================================
 bound_encoding.jl — how simple box limits are attached to decision variables
================================================================================
 CONSTRAINT (default): explicit ≤-form inequalities (dual-friendly ConstraintRef export).
 VARIABLE: JuMP variable lower_bound / upper_bound at creation (fewer rows for the solver).
=#

"""How simple box limits are encoded in the JuMP model."""
@enum BoundEncoding CONSTRAINT VARIABLE

"""Read bound encoding from a builder param dict; default CONSTRAINT."""
function bound_encoding_from_param(param::OrderedDict{Symbol, Any})::BoundEncoding
    return get(param, :bound_encoding, CONSTRAINT)
end

function bound_encoding_from_meta(meta::OrderedDict{Symbol, Any})::BoundEncoding
    return get(meta, :bound_encoding, CONSTRAINT)
end
