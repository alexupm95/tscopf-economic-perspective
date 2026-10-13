#=
================================================================================
 functions_4_TS_fullbus_helpers.jl — FullBus TS helper routines
================================================================================
=#

"""Read `JuMP.start_value` when set, else `fallback` (used for δ_tf trajectory seeds)."""
function _opf_scalar_hint(v::JuMP.VariableRef, fallback::Float64)::Float64
    sv = JuMP.start_value(v)
    return isnothing(sv) ? fallback : Float64(sv)
end

# The FULL_BUS dynamic-period admittance is built by Calculate_Ybus_fullbus_dynamics in
# _common/functions_4_admittance_matrices.jl (network only + optional fault shunt). Loads
# are NOT folded into the Ybus here — the nodal power balance models them via the ZIP
# model (eq_const_fullbus_Pbalance!), so folding them in would double-count.

"""Apply solved ACOPF values as JuMP start hints on pre-fault E, δ, P_m.

For each active generator, the converged ACOPF terminal phasor (V∠θ, P_g + jQ_g) is
mapped to the classical model's internal EMF E∠δ behind the transient reactance Xd_tr
(`classical_E_delta_warmstart`). Seeding E and δ at this consistent operating point keeps
the nonlinear solve close to the steady state; P_m is seeded at P_g (mechanical balances
electrical injection at t=0). These are only `start` hints, not constraints.
"""
function _set_fullbus_init_warm_starts!(
    E::OrderedDict{Int, JuMP.VariableRef},
    δ::OrderedDict{Int, JuMP.VariableRef},
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
        # Internal EMF phasor E∠δ from the terminal phasor behind Xd_tr (classical model).
        start_E, start_δ = classical_E_delta_warmstart(
            val_V[bus], val_θ[bus], val_Pg[gen], val_Qg[gen], DGEN_DYN.Xd_tr[gen])
        JuMP.set_start_value(E[gen], start_E)
        JuMP.set_start_value(δ[gen], start_δ)
        if P_m !== nothing
            JuMP.set_start_value(P_m[gen], val_Pg[gen])  # P_m ≈ P_g at the equilibrium
        end
    end
    return nothing
end
"""
Build the standard AC power-flow nodal injection expressions on the sparse `Ybus`.

For each bus i and time t, returns the real/reactive power injected into the network:

    P_i = V_i · Σ_j V_j (G_ij·cos θ_ij + B_ij·sin θ_ij)
    Q_i = V_i · Σ_j V_j (G_ij·sin θ_ij − B_ij·cos θ_ij)

with θ_ij = θ_i − θ_j, G = real(Y), B = imag(Y). Only the structural non-zeros of row i
(`row_Y.nzind/nzval`) are summed, which is why the FULL_BUS path keeps the network sparse
rather than Kron-reducing it. Islanded buses (no incident branch) are rejected up front.
Returned as named expressions reused by both the balance constraints and dual export.
"""
function _build_fullbus_nodal_injection_terms!(
    model::JuMP.Model,
    nBUS::Int64,
    bus_gen_circ_dict::OrderedDict,
    V_t::OrderedDict,
    θ_t::OrderedDict,
    time_window::Vector{Float64},
    Ybus::SparseMatrixCSC,
)
    terms_Pb = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}()
    terms_Qb = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}()
    for i in 1:nBUS
        assert_bus_not_islanded(i, bus_gen_circ_dict[i][:circ])
        terms_Pb[i] = OrderedDict{Int, JuMP.NonlinearExpr}()
        terms_Qb[i] = OrderedDict{Int, JuMP.NonlinearExpr}()
        row_Y = Ybus[i, :]  # only this row's non-zeros contribute to bus i's injection
        for t in eachindex(time_window)
            terms_Pb[i][t] = JuMP.@expression(model,
                V_t[i][t] * sum(
                    V_t[j][t] * (real(y) * cos(θ_t[i][t] - θ_t[j][t]) +
                                 imag(y) * sin(θ_t[i][t] - θ_t[j][t]))
                    for (j, y) in zip(row_Y.nzind, row_Y.nzval)))
            terms_Qb[i][t] = JuMP.@expression(model,
                V_t[i][t] * sum(
                    V_t[j][t] * (real(y) * sin(θ_t[i][t] - θ_t[j][t]) -
                                 imag(y) * cos(θ_t[i][t] - θ_t[j][t]))
                    for (j, y) in zip(row_Y.nzind, row_Y.nzval)))
        end
    end
    return terms_Pb, terms_Qb
end

"""
    zip_by_bus(DBUS, DLOAD, zip_load_p, zip_load_q; nBUS) -> (Matrix, Matrix)

Per-bus ZIP rows for the FULL_BUS nodal balances: `nBUS × 3` matrices in (Z, I, P) order,
one row per bus, for the active and the reactive demand.

Without a load table (`DLOAD === nothing`) every row is the system-wide tuple and the
balances behave exactly as before. With one, a bus carrying several loads gets the
**demand-weighted mean** of their shares, which is not an approximation: every load at
bus *i* sees the same trajectory voltage `V_t[i][t]` and the same nominal `V[i]`, so

    Σ_l p_d,l·(a_l·V_t² + b_l·V_t·V + c_l·V²) ≡ p_d,bus·(ā·V_t² + b̄·V_t·V + c̄·V²)

with `ā = Σ_l (p_d,l / p_d,bus)·a_l`. The individual loads therefore need no rows of their
own in the model; their separate prices are recovered afterwards by contracting each
load's own polynomial with the bus's balance dual.

Weights are `p_d` for the active matrix and `q_d` for the reactive one — a load may well
be a different share of its bus's MW than of its MVAr. A bus whose rows cancel to ≈0
falls back to the system-wide tuple; that balance row is dropped by the zero-load branch
of the builders anyway.
"""
function zip_by_bus(
    DBUS::DataFrame,
    DLOAD::Union{Nothing, DataFrame},
    zip_load_p::NTuple{3, Float64},
    zip_load_q::NTuple{3, Float64};
    nBUS::Int = nrow(DBUS),
)::Tuple{Matrix{Float64}, Matrix{Float64}}
    ZIP_P = repeat(collect(Float64, zip_load_p)', nBUS)
    ZIP_Q = repeat(collect(Float64, zip_load_q)', nBUS)
    (DLOAD === nothing || nrow(DLOAD) == 0) && return (ZIP_P, ZIP_Q)

    # Accumulate Σ_l p_d,l·(Z, I, P)_l and Σ_l p_d,l per bus, then divide.
    acc_p = zeros(Float64, nBUS, 3); w_p = zeros(Float64, nBUS)
    acc_q = zeros(Float64, nBUS, 3); w_q = zeros(Float64, nBUS)
    for l in 1:nrow(DLOAD)
        i = Int(DLOAD.bus[l])
        1 <= i <= nBUS || throw(ArgumentError(
            "load_data id $(Int(DLOAD.id[l])) sits at bus index $i, outside 1:$nBUS."))
        p, q = DLOAD.p_d[l], DLOAD.q_d[l]
        acc_p[i, :] .+= p .* (DLOAD.aZ_p[l], DLOAD.aI_p[l], DLOAD.aP_p[l])
        acc_q[i, :] .+= q .* (DLOAD.aZ_q[l], DLOAD.aI_q[l], DLOAD.aP_q[l])
        w_p[i] += p
        w_q[i] += q
    end
    for i in 1:nBUS
        abs(w_p[i]) > 1e-12 && (ZIP_P[i, :] .= acc_p[i, :] ./ w_p[i])
        abs(w_q[i]) > 1e-12 && (ZIP_Q[i, :] .= acc_q[i, :] ./ w_q[i])
    end
    return (ZIP_P, ZIP_Q)
end
