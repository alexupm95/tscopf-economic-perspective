#=
================================================================================
 functions_4_TS_dq_helpers.jl — DQ 4th-order machine helpers
================================================================================
 ACOPF warm-start projection (complex V, I → dq steady state) and the optional
 (1+Δω) scale used in stator algebra. Swing equations are unchanged.
================================================================================
=#

"""
    dgen_dyn_by_id(DGEN_DYN, DGEN) -> DataFrame

Copy of the SG dynamic table re-laid so that **row `g` holds generator id `g`**, for
`g ∈ 1:nrow(DGEN)`.

Every builder reads machine data positionally, `DGEN_DYN.<col>[gen]`, with `gen` a
generator id. That is only right when the table is id-aligned, which the SG-only paths
guarantee (`validate_dgen_dyn_row_ids!`). Under GFM the table holds SG rows only
(`validate_dgen_dyn_subset_ids!` accepts any id subset), so a converter at an interior
id shifts every later machine up one row: with converter id 2 of 10, G3 read G4's
reactances and inertia, and G10 threw a `BoundsError`. Aligning once here fixes all
of those readers at the one place they share.

Converter ids get a placeholder row: `id` and `bus` from `DGEN`, `NaN` in every other
column. Nothing should read a converter's machine data; a `NaN` that reaches the model
fails loudly at the solver instead of passing a neighbour's numbers.

An already-aligned table (ids `1:n`, or no `id` column) comes back as a plain copy, and
so does an empty one: a converter-only fleet has no machine to misalign, and padding it
with `NaN` rows would make it look like a fleet of synchronous machines.
"""
function dgen_dyn_by_id(DGEN_DYN::DataFrame, DGEN::DataFrame)::DataFrame
    nrow(DGEN_DYN) == 0 && return deepcopy(DGEN_DYN)
    n = nrow(DGEN)
    ids = :id in propertynames(DGEN_DYN) ? Int.(DGEN_DYN.id) : collect(1:nrow(DGEN_DYN))
    ids == 1:n && return deepcopy(DGEN_DYN)

    row_of = Dict(g => i for (i, g) in enumerate(ids))
    out = DataFrame()
    for col in propertynames(DGEN_DYN)
        src = DGEN_DYN[!, col]
        out[!, col] = if col === :id
            collect(1:n)
        elseif col === :bus
            Int64[haskey(row_of, g) ? src[row_of[g]] : DGEN.bus[g] for g in 1:n]
        else
            Float64[haskey(row_of, g) ? src[row_of[g]] : NaN for g in 1:n]
        end
    end
    return out
end

"""
    dq_machine_warmstart(v, θ, P_g, Q_g, Xd_tr, Xq_tr, Xd, Xq, Ra) -> (E_fd, δ, Ed, Eq, Id, Iq)

Map a solved ACOPF terminal injection to the dq steady-state.

Matches the reference `Define_Dyn_Var_EδPm!`: angle from `Eq_c = V + I*(Ra + j*Xq)`
(synchronous `Xq`), then Ed/Eq from transient reactances.
"""
function dq_machine_warmstart(
    v::Real,
    θ::Real,
    P_g::Real,
    Q_g::Real,
    Xd_tr::Real,
    Xq_tr::Real,
    Xd::Real,
    Xq::Real,
    Ra::Real,
)::NTuple{6, Float64}
    V_c = Float64(v) * exp(1im * Float64(θ))
    I_c = (Float64(P_g) - 1im * Float64(Q_g)) / conj(V_c)
    # Reference: angle from phasor behind Ra + j Xq (synchronous), not Xq_tr
    Eq_c = V_c + I_c * (Ra + 1im * Float64(Xq))
    δ = angle(Eq_c)
    I_mag = abs(I_c)
    I_ang = angle(I_c)
    Id = I_mag * sin(δ - I_ang)
    Iq = I_mag * cos(δ - I_ang)
    Vd = Float64(v) * sin(δ - θ)
    Vq = Float64(v) * cos(δ - θ)
    Ed = Vd + Ra * Id - Xq_tr * Iq
    Eq = Vq + Ra * Iq + Xd_tr * Id
    E_fd = Eq + (Xd - Xd_tr) * Id
    return (Float64(E_fd), Float64(δ), Float64(Ed), Float64(Eq), Float64(Id), Float64(Iq))
end

"""
    dq_rotor_scale(Δω, speed_in_algebra) -> JuMP expression

When `speed_in_algebra=true` (default, RMS reference style), stator equations use
`(1+Δω)` on Ed/Eq and `Pe = (1+Δω)·Te`. When false, speed is neglected in those
algebraic links only; Δω still enters the swing equation.
"""
dq_rotor_scale(Δω::JuMP.VariableRef, speed_in_algebra::Bool) =
    speed_in_algebra ? (1.0 + Δω) : 1.0

"""Copy solved ACOPF (V, θ, P_g, Q_g) into JuMP starts for all dq pre-fault variables."""
function _set_dq_init_warm_starts!(
    E_fd::OrderedDict{Int, JuMP.VariableRef},
    δ::OrderedDict{Int, JuMP.VariableRef},
    Ed::OrderedDict{Int, JuMP.VariableRef},
    Eq::OrderedDict{Int, JuMP.VariableRef},
    Id::OrderedDict{Int, JuMP.VariableRef},
    Iq::OrderedDict{Int, JuMP.VariableRef},
    P_m::Union{Nothing, OrderedDict{Int, JuMP.VariableRef}},
    active_gen::Vector{Int64},
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    val_V::Dict,
    val_θ::Dict,
    val_Pg::Dict,
    val_Qg::Dict,
)
    for gen in active_gen
        bus = DGEN.bus[gen]
        start_E_fd, start_δ, start_Ed, start_Eq, start_Id, start_Iq = dq_machine_warmstart(
            val_V[bus], val_θ[bus], val_Pg[gen], val_Qg[gen],
            DGEN_DYN.Xd_tr[gen], DGEN_DYN.Xq_tr[gen], DGEN_DYN.Xd[gen],
            DGEN_DYN.Xq[gen], DGEN_DYN.Ra[gen])
        JuMP.set_start_value(E_fd[gen], start_E_fd)
        JuMP.set_start_value(δ[gen], start_δ)
        JuMP.set_start_value(Ed[gen], start_Ed)
        JuMP.set_start_value(Eq[gen], start_Eq)
        JuMP.set_start_value(Id[gen], start_Id)
        JuMP.set_start_value(Iq[gen], start_Iq)
        if P_m !== nothing
            JuMP.set_start_value(P_m[gen], val_Pg[gen])
        end
    end
    return nothing
end
