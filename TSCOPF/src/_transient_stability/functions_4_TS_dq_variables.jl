#=
================================================================================
 functions_4_TS_dq_variables.jl — DQ 4th-order machine variables
================================================================================
 Pre-fault: E_fd and δ are created in the orchestrator (bounded, per-gen limits).
 Here: Ed/Eq/Id/Iq at t=0 and time-indexed Ed/Eq/Id/Iq/Te during fault/post-fault windows.
================================================================================
=#

"""
Pre-fault subtransient emfs and dq currents.

Unbounded algebraic states; steady-state is enforced by `eq_const_dq_init_steady_state!`.
`E_fd` and `δ` are declared separately via the classical FULL_BUS var builders.

Returns a `NamedTuple` `(; Ed, Eq, Id, Iq)`. The dq builders return named fields rather than
bare tuples so a caller cannot silently bind the wrong state to the wrong name when a sibling
builder carries one field more (see `var_dq_gen_state_time!`, which also returns `Te`).
"""
function var_dq_prefault_algebraic!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
)
    Ed = OrderedDict{Int, JuMP.VariableRef}()
    Eq = OrderedDict{Int, JuMP.VariableRef}()
    Id = OrderedDict{Int, JuMP.VariableRef}()
    Iq = OrderedDict{Int, JuMP.VariableRef}()
    for gen in active_gen
        Ed[gen] = JuMP.@variable(model, base_name="Ed_p[$gen]")
        Eq[gen] = JuMP.@variable(model, base_name="Eq_p[$gen]")
        Id[gen] = JuMP.@variable(model, base_name="Id[$gen]")
        Iq[gen] = JuMP.@variable(model, base_name="Iq[$gen]")
    end
    return (; Ed, Eq, Id, Iq)
end

"""Per-generator, per-step dq machine states.

`Id`/`Iq` are created for every generator in `active_gen` (SG and GFM need
terminal currents for KCL). `Ed`/`Eq`/`Te` are created only for `emf_gens`
(default: all of `active_gen`). The reference implementation never declares time-window Ed/Eq/Te for
GFMs — leaving them free breaks Ipopt L-BFGS.

Returns a `NamedTuple` `(; Ed, Eq, Id, Iq, Te)`.
"""
function var_dq_gen_state_time!(
    model::JuMP.Model,
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    δ_ref::OrderedDict{Int, JuMP.VariableRef},
    Ed_ref::OrderedDict{Int, JuMP.VariableRef},
    Eq_ref::OrderedDict{Int, JuMP.VariableRef},
    Id_ref::OrderedDict{Int, JuMP.VariableRef},
    Iq_ref::OrderedDict{Int, JuMP.VariableRef},
    val_Pg::Dict;
    suffix::String,
    emf_gens::AbstractVector{Int} = active_gen,
)
    Ed = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Eq = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Id = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Iq = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    Te = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}()
    emf_set = Set{Int}(Int(g) for g in emf_gens)
    for gen in active_gen
        Id[gen] = OrderedDict{Int, JuMP.VariableRef}()
        Iq[gen] = OrderedDict{Int, JuMP.VariableRef}()
        Id0 = _opf_scalar_hint(Id_ref[gen], 0.0)
        Iq0 = _opf_scalar_hint(Iq_ref[gen], 0.0)
        make_emf = gen in emf_set
        if make_emf
            Ed[gen] = OrderedDict{Int, JuMP.VariableRef}()
            Eq[gen] = OrderedDict{Int, JuMP.VariableRef}()
            Te[gen] = OrderedDict{Int, JuMP.VariableRef}()
            Ed0 = _opf_scalar_hint(Ed_ref[gen], 0.0)
            Eq0 = _opf_scalar_hint(Eq_ref[gen], 0.0)
            pg0 = get(val_Pg, gen, 0.0)
        end
        for t in eachindex(time_window)
            Id[gen][t] = JuMP.@variable(model, base_name="Id_$suffix[$gen,$t]", start=Id0)
            Iq[gen][t] = JuMP.@variable(model, base_name="Iq_$suffix[$gen,$t]", start=Iq0)
            if make_emf
                Ed[gen][t] = JuMP.@variable(model, base_name="Ed_$suffix[$gen,$t]", start=Ed0)
                Eq[gen][t] = JuMP.@variable(model, base_name="Eq_$suffix[$gen,$t]", start=Eq0)
                Te[gen][t] = JuMP.@variable(model, base_name="Te_$suffix[$gen,$t]", start=pg0)
            end
        end
    end
    return (; Ed, Eq, Id, Iq, Te)
end
