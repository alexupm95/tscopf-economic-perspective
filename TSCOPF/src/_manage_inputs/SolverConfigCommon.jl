"""Apply user-supplied raw MOI optimizer attributes after typed options."""
function _apply_solver_raw_options!(model::JuMP.Model, raw_options::Dict{String, Any})::JuMP.Model
    for (name, value) in raw_options
        JuMP.set_optimizer_attribute(model, name, value)
    end
    return model
end
