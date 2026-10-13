#=
================================================================================
 functions_4_TS_kron_eqconst.jl — Kron TS equality builders
================================================================================
=#

# ===================================================================================
#                           Equality Constraints (shared)
# ===================================================================================

function eq_const_kron_initial_mechanical_power!(model::JuMP.Model,
    P_m::OrderedDict,
    P_g::OrderedDict,
    active_gen::Vector{Int64}
    )

    eq_const_Pm_init = OrderedDict{Int64, JuMP.ConstraintRef}()

    for gen in active_gen
        eq_const_Pm_init[gen] = JuMP.@constraint(model, P_m[gen] - P_g[gen] == 0.0)
    end

    return eq_const_Pm_init
end

# Function to create equality constraints to determine the COI equations
#
# `weights` is the inertia weight per unit, built by `_coi_weights`. It used to be read
# here as `DGEN_DYN.H[gen]`, which is *row-position* indexing and therefore only correct
# while generator ids happen to equal row numbers. Taking the weights as an argument keeps
# that lookup in one place and lets grid-forming converters — absent from `DGEN_DYN`
# entirely — carry their equivalent inertia `Tf/(2·mp)` when `gfm_in_coi` is on.
function eq_const_kron_COI_generic!(model::JuMP.Model,
    _var::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}},
    _COI::OrderedDict{Int, JuMP.VariableRef},
    active_gen::Vector{Int64},
    weights::OrderedDict{Int, Float64},
    time_window::Vector{Float64}
    )

    # Synchronized-set inertia (surviving generators only). For GL gen-trip,
    # active_gen excludes disconnected machines; SC runs keep the full set.
    H_total = sum(weights[gen] for gen in active_gen)
    H_total > 0.0 || throw(ArgumentError(
        "COI inertia denominator is $(H_total); the synchronised set carries no inertia."))

    # Terms that will be summed to calculate the dynamics of the COI
    terms_COI_tf_dict  = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary

    # Loop to create the variables
    for gen in active_gen
        terms_COI_tf_dict[gen]  = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()

        for t in eachindex(time_window)
            terms_COI_tf_dict[gen][t]  = weights[gen] * _var[gen][t]
        end
    end

    # Gives the expressions to calculate the COI variable
    expr_COI_per_time  = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_COI_tf_dict))  for t in eachindex(time_window))

    eq_const_COI_tf  = OrderedDict{Int, JuMP.ConstraintRef}()
    for t in eachindex(time_window)
        eq_const_COI_tf[t]  = JuMP.@constraint(model, _COI[t] - (expr_COI_per_time[t]  / H_total)  == 0.0)

    end

    return eq_const_COI_tf
end

function eq_const_kron_δ_swingeq_generic!(model::JuMP.Model,
    active_gen::Vector{Int64},
    δ::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    δ_0::OrderedDict{Int64, JuMP.VariableRef},
    Δω_0::Float64,
    time_window::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_δ[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        for t in eachindex(time_window)
            if t == 1
                eq_const_δ[gen][t] = JuMP.@constraint(model, δ[gen][t] - δ_0[gen] - (ω_syn * (Δt/2) * (Δω[gen][t] + Δω_0)) == 0.0)
            else
                eq_const_δ[gen][t] = JuMP.@constraint(model, δ[gen][t] - δ[gen][t-1] - (ω_syn * (Δt/2) * (Δω[gen][t] + Δω[gen][t-1])) == 0.0)
            end
        end
    end

    return eq_const_δ

end

# Post-Fault
function eq_const_kron_δ_swingeq_generic!(model::JuMP.Model,
    active_gen::Vector{Int64},
    δ::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    δ_0::OrderedDict{Int64, JuMP.VariableRef},
    Δω_0::OrderedDict{Int64, JuMP.VariableRef},
    time_window::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_δ[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        for t in eachindex(time_window)
            if t == 1
                eq_const_δ[gen][t] = JuMP.@constraint(model, δ[gen][t] - δ_0[gen] - (ω_syn * (Δt/2) * (Δω[gen][t] + Δω_0[gen])) == 0.0)
            else
                eq_const_δ[gen][t] = JuMP.@constraint(model, δ[gen][t] - δ[gen][t-1] - (ω_syn * (Δt/2) * (Δω[gen][t] + Δω[gen][t-1])) == 0.0)
            end
        end
    end

    return eq_const_δ

end
# ===================================================================================
#                           Equality Constraints (TSC-ACOPF specific)
# ===================================================================================

# Function to create the equality constraint initial active power (base for dynamic simulation)
function eq_const_tsred_initial_active_power!(model::JuMP.Model,
    V::OrderedDict,
    θ::OrderedDict,
    P_g::OrderedDict,
    E::OrderedDict,
    δ::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64}
    )

    eq_const_P_init = OrderedDict{Int64, JuMP.ConstraintRef}()

    for gen in active_gen
        bus = DGEN.bus[gen]
        Xd_prime = DGEN_DYN.Xd_tr[gen] # Transient reactance of the generator

        # Constraint for active power
        eq_const_P_init[gen] = JuMP.@constraint(model, P_g[gen] - ((E[gen] * V[bus] * sin(δ[gen] - θ[bus])) / Xd_prime) == 0.0 )
    end

    return eq_const_P_init
    
end

# Function to create the equality constraint initial reactive power (base for dynamic simulation)
function eq_const_tsred_initial_reactive_power!(model::JuMP.Model,
    V::OrderedDict,
    θ::OrderedDict,
    Q_g::OrderedDict,
    E::OrderedDict,
    δ::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64}
    )

    eq_const_Q_init = OrderedDict{Int64, JuMP.ConstraintRef}()

    for gen in active_gen
        bus = DGEN.bus[gen]
        Xd_prime = DGEN_DYN.Xd_tr[gen] # Transient reactance of the generator

        # Constraint for active power
        eq_const_Q_init[gen] = JuMP.@constraint(model, Q_g[gen] - ((E[gen] * V[bus] * cos(δ[gen] - θ[bus])) / Xd_prime) + ((V[bus]^2) / Xd_prime) == 0.0 )
    end

    return eq_const_Q_init
    
end

# Function to create the equality constraint initial mechanical power (base for dynamic simulation)

# Function to create equality constraints for the electrical power as a function of δ
function eq_const_tsred_Pe_generic!(model::JuMP.Model,
    E::OrderedDict{Int64, JuMP.VariableRef},
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    δ::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    Yred::Matrix
    )

    nGEN_active = length(active_gen) # Number of active generators
    eq_const_Pe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    terms_Pe_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() # Dictionary of Terms that will be summed to calculate Pe

    for (i, id_gen1) in enumerate(active_gen)

        eq_const_Pe[id_gen1] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints

        terms_Pe_dict[id_gen1] = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of terms to form a nonlinear expression

        for t in eachindex(time_window) # Loop over the time
            terms_to_sum = JuMP.NonlinearExpr[] # List of expressions

            for (j, id_gen2) in enumerate(active_gen) #in eachindex(active_gen) # Loop over active generators
                G = real(Yred[i,j]) # Conductance of the reduced admittance matrix
                B = imag(Yred[i,j]) # Susceptance of the reduced admittance matrix

                expression_Pe = E[id_gen2] * (G*cos(δ[id_gen1][t] - δ[id_gen2][t]) + B*sin(δ[id_gen1][t] - δ[id_gen2][t]))
                push!(terms_to_sum, expression_Pe)
            end
            terms_Pe_dict[id_gen1][t] = E[id_gen1] * sum(terms_to_sum)
            eq_const_Pe[id_gen1][t] = JuMP.@constraint(model, Pe[id_gen1][t] - terms_Pe_dict[id_gen1][t] == 0.0) # Equality constraint for the electrical power
        end
    end

    return eq_const_Pe

end

# ==================================================================================
# Kron-reduced Qe(δ) — post-solve diagnostics as JuMP expressions (not constraints)
# ==================================================================================
# Qe does not enter the classical swing equations or TSC objective on Yred, so we
# avoid equality constraints (and dual_Qe) and evaluate trajectories after solve.

"""
    expr_tsred_Qe_generic!(E, δ, active_gen, time_window, Yred)

Generator reactive power on the Kron-reduced network:

``Q_{e,g}(t) = E_g \\sum_{h} E_h \\bigl(G_{gh}\\sin(\\delta_g-\\delta_h)
- B_{gh}\\cos(\\delta_g-\\delta_h)\\bigr)``,

with `(G,B)` from `real/imag(Yred)` in the same generator ordering as `Pe`.
"""
function expr_tsred_Qe_generic!(
    E::OrderedDict{Int64, JuMP.VariableRef},
    δ::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    Yred::Matrix,
)
    Qe_expr = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}()

    for (i, id_gen1) in enumerate(active_gen)
        Qe_expr[id_gen1] = OrderedDict{Int, JuMP.NonlinearExpr}()

        for t in eachindex(time_window)
            terms_to_sum = JuMP.NonlinearExpr[]

            for (j, id_gen2) in enumerate(active_gen)
                G = real(Yred[i, j])
                B = imag(Yred[i, j])
                push!(terms_to_sum,
                    E[id_gen2] * (G * sin(δ[id_gen1][t] - δ[id_gen2][t]) -
                                  B * cos(δ[id_gen1][t] - δ[id_gen2][t])))
            end
            Qe_expr[id_gen1][t] = E[id_gen1] * sum(terms_to_sum)
        end
    end

    return Qe_expr
end

# Function to define the equality constraints of rotor according to the swing equation
# Fault

# Function to define the equality constraints of speed deviation according to the swing equation
# Fault
function eq_const_tsred_Δω_swingeq_generic!(model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pe_0::OrderedDict{Int64, JuMP.VariableRef},
    Pm::OrderedDict{Int64, JuMP.VariableRef},
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω_0::Float64,
    time_window::Vector{Float64},
    Δt::Float64
    )

    eq_const_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        H = DGEN_DYN.H[gen]
        D = DGEN_DYN.D[gen]

        for t in eachindex(time_window)
            if t == 1
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_0 - (Δt / (4*H))*(2*Pm[gen] - Pe[gen][t] - Pe_0[gen])  == 0.0)

            else
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω[gen][t-1] - (Δt / (4*H))*(2*Pm[gen] - Pe[gen][t] - Pe[gen][t-1])  == 0.0)
            end
        end
    end

    return eq_const_Δω

end

# Post-Fault
function eq_const_tsred_Δω_swingeq_generic!(model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pe_0::OrderedDict{Int64, JuMP.VariableRef},
    Pm::OrderedDict{Int64, JuMP.VariableRef},
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω_0::OrderedDict{Int64, JuMP.VariableRef},
    time_window::Vector{Float64},
    Δt::Float64
    )

    eq_const_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        H = DGEN_DYN.H[gen]
        D = DGEN_DYN.D[gen]

        for t in eachindex(time_window)
            if t == 1
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_0[gen] - (Δt / (4*H))*(2*Pm[gen] - Pe[gen][t] - Pe_0[gen])  == 0.0)

            else
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω[gen][t-1] - (Δt / (4*H))*(2*Pm[gen] - Pe[gen][t] - Pe[gen][t-1])  == 0.0)
            end
        end
    end

    return eq_const_Δω

end
# ===================================================================================
#                     Equality Constraints (TSC-DCOPF specific)
# ===================================================================================
# Linearised P–δ init, Taylor Pe, and modified Δω swing (Pg-based first time step).

# Pre-fault active power: P_g − (δ − θ)/X'd = 0 (E = V = 1 p.u., small-angle).
function eq_const_tsredlinear_initial_active_power!(model::JuMP.Model,
    θ::OrderedDict,
    P_g::OrderedDict,
    δ::OrderedDict,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    active_gen::Vector{Int64}
    )

    eq_const_P_init = OrderedDict{Int64, JuMP.ConstraintRef}()

    for gen in active_gen
        bus = DGEN.bus[gen]
        Xd_prime = DGEN_DYN.Xd_tr[gen] # Transient reactance of the generator

        # Constraint for active power
        eq_const_P_init[gen] = JuMP.@constraint(model, P_g[gen] - ((1.0 * 1.0 * (δ[gen] - θ[bus])) / Xd_prime) == 0.0 )
    end

    return eq_const_P_init
    
end

# Function to create the equality constraint initial mechanical power (base for dynamic simulation)
function eq_const_tsredlinear_Δω_swingeq_generic_modified!(model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pg::OrderedDict{Int64, JuMP.VariableRef},
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω_0::Float64,
    time_window::Vector{Float64},
    Δt::Float64
    )

    eq_const_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        H = DGEN_DYN.H[gen]
        D = DGEN_DYN.D[gen]

        for t in eachindex(time_window)
            if t == 1
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_0 - (Δt / (4*H))*(Pg[gen] - Pe[gen][t])  == 0.0)

            else
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω[gen][t-1] - (Δt / (4*H))*(2*Pg[gen] - Pe[gen][t] - Pe[gen][t-1])  == 0.0)
            end
        end
    end

    return eq_const_Δω

end

# Post-Fault
function eq_const_tsredlinear_Δω_swingeq_generic_modified!(model::JuMP.Model,
    active_gen::Vector{Int64},
    DGEN_DYN::DataFrame,
    Pe_0::OrderedDict{Int64, JuMP.VariableRef},
    Pg::OrderedDict{Int64, JuMP.VariableRef},
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    Δω_0::OrderedDict{Int64, JuMP.VariableRef},
    time_window::Vector{Float64},
    Δt::Float64
    )

    eq_const_Δω = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for gen in active_gen
        eq_const_Δω[gen] = OrderedDict{Int, JuMP.ConstraintRef}()

        H = DGEN_DYN.H[gen]
        D = DGEN_DYN.D[gen]

        for t in eachindex(time_window)
            if t == 1
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω_0[gen] - (Δt / (4*H))*(2*Pg[gen] - Pe[gen][t] - Pe_0[gen])  == 0.0)

            else
                eq_const_Δω[gen][t] = JuMP.@constraint(model, (1.0 + ((D*Δt) / (4*H))) * Δω[gen][t] - (1.0 - ((D*Δt) / (4*H))) * Δω[gen][t-1] - (Δt / (4*H))*(2*Pg[gen] - Pe[gen][t] - Pe[gen][t-1])  == 0.0)
            end
        end
    end

    return eq_const_Δω

end

# ===================================================================================
#                         Inequality Constraints
# ===================================================================================
# Function to define the transient stability constraint
# Bounds the angular motion of each generator with respect to the Center of Inertia
function eq_const_tsredlinear_Pe_taylor!(
    model::JuMP.Model,
    Pe::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    δ::OrderedDict{Int64, OrderedDict{Int64, JuMP.VariableRef}},
    active_gen::Vector{Int64},
    time_window::Vector{Float64},
    Yred::Matrix,
    δ_ref::OrderedDict{Int64, Float64}    # Reference angles from OPF [rad]
    )

    eq_const_Pe = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    for (i, id_gen1) in enumerate(active_gen)

        eq_const_Pe[id_gen1] = OrderedDict{Int, JuMP.ConstraintRef}()

        # ------------------------------------------------------------------
        # Pre-compute Taylor coefficients for generator i vs all j.
        # These are pure scalars — computed once outside the time loop.
        # For each pair (i,j):
        #   C_const[j] = G·cos(Δδ⁰) + B·sin(Δδ⁰)         [constant term]
        #              + (G·sin(Δδ⁰) - B·cos(Δδ⁰)) · Δδ⁰  [offset absorbed]
        #   C_coeff[j] = -G·sin(Δδ⁰) + B·cos(Δδ⁰)         [coefficient of (δᵢ-δⱼ)]
        #
        # So Pe_i[t] = Σⱼ [ C_const[j] + C_coeff[j]·(δᵢ[t] - δⱼ[t]) ]
        # ------------------------------------------------------------------
        C_const = zeros(Float64, length(active_gen))
        C_coeff = zeros(Float64, length(active_gen))

        for (j, id_gen2) in enumerate(active_gen)
            G   = real(Yred[i, j])
            B   = imag(Yred[i, j])
            Δδ⁰ = δ_ref[id_gen1] - δ_ref[id_gen2]   # Reference angle difference [rad]

            cosΔδ⁰  = cos(Δδ⁰)
            sinΔδ⁰  = sin(Δδ⁰)

            # Coefficient on the linear term (δᵢ - δⱼ)
            C_coeff[j] = -G * sinΔδ⁰ + B * cosΔδ⁰

            # Constant term: Pe at reference point minus the linear correction offset
            # i.e.  (G·cos⁰ + B·sin⁰) - C_coeff·Δδ⁰
            C_const[j] = (G * cosΔδ⁰ + B * sinΔδ⁰) - C_coeff[j] * Δδ⁰
        end

        # ------------------------------------------------------------------
        # Build one affine Pe constraint per time step.
        # All time steps share the same Taylor coefficients since the
        # linearization point δ_ref is fixed (from OPF, not time-varying).
        # ------------------------------------------------------------------
        for t in eachindex(time_window)

            Pe_expr = JuMP.AffExpr(0.0)

            for (j, id_gen2) in enumerate(active_gen)
                # Constant contribution from pair (i,j)
                JuMP.add_to_expression!(Pe_expr, C_const[j])

                # Linear contribution: +C_coeff · δᵢ[t]
                JuMP.add_to_expression!(Pe_expr,  C_coeff[j], δ[id_gen1][t])

                # Linear contribution: -C_coeff · δⱼ[t]
                JuMP.add_to_expression!(Pe_expr, -C_coeff[j], δ[id_gen2][t])
            end

            # Pe[i][t] = linearized expression → single equality per (i,t)
            eq_const_Pe[id_gen1][t] = JuMP.@constraint(model,
                Pe[id_gen1][t] - Pe_expr == 0.0)
        end
    end

    return eq_const_Pe

end
